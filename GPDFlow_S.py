import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import normflows as nf


class DataTransform(nn.Module):
    """
    Elementwise transform:
      If gamma_j = 0:
        forward_transform(y_j): x_j = sigma_j * y_j
        inverse_transform(x_j): y_j = x_j / sigma_j
      Otherwise (gamma_j != 0):
        forward_transform(y_j): x_j = sigma_j * [exp(gamma_j * y_j) - 1] / gamma_j
        inverse_transform(x_j): y_j = (1 / gamma_j)* log(1 + gamma_j * x_j / sigma_j)
    """

    def __init__(self, dim, device, fix_margin):
        super().__init__()
        self.dim = dim
        # Store log_sigma so sigma = exp(log_sigma) > 0
        if fix_margin:
            self.log_sigma = torch.zeros(dim, device=device)
            self.theta = torch.zeros(dim, device=device)
        else:
            self.log_sigma = nn.Parameter(torch.zeros(dim, device=device))
            self.theta = nn.Parameter(torch.zeros(dim, device=device))

    def get_sigma(self):
        sigma = torch.exp(self.log_sigma)
        sigma = torch.clamp(sigma, min=1e-6, max=1e6)
        return sigma

    def get_gamma(self):
        gamma = torch.tanh(self.theta)
        return gamma

    def forward_transform(self, y):
        sigma = self.get_sigma().unsqueeze(0)
        gamma = self.get_gamma().unsqueeze(0)

        sigma_b = sigma.expand_as(y)
        gamma_b = gamma.expand_as(y)

        # Reformulate sigma*expm1(gamma*y)/gamma as sigma*y*expm1(t)/t where t=gamma*y.
        # This avoids dividing by gamma, whose 1/gamma^2 backward gradient caused NaN
        # when |gamma| was near eps and quantile_lambda was large.
        t = gamma_b * y
        # Clamp so expm1(t) cannot overflow to inf (and poison the backward pass
        # with NaN); mirrors the t-clamp in inverse_transform below. Only an
        # upper bound is needed: for very negative t, expm1(t) -> -1 and
        # expm1_over_t -> 0, both already finite and gradient-safe. Kept well
        # below expm1's own overflow point (~88 in float32) so that even the
        # worst-case compounding of sigma's outer clamp (1e6) with y at its
        # _r_to_v-clamped extreme (~2.2e4) can't produce an x large enough to
        # overflow when later squared inside torch.cdist in the MMD kernel.
        #
        # clamp_mask marks samples where t is actually being capped. For those,
        # the *unclamped* y must not be used as the outer multiplier below: doing
        # so previously made x keep growing linearly (and unboundedly) in the raw
        # y with a fixed huge slope (expm1(15)/15), instead of saturating -- and
        # made dx/dsigma explode while dx/dgamma collapsed to exactly 0 (since
        # expm1_over_t no longer depends on gamma once t is pinned). Under MMD
        # training, y_fake can legitimately reach |y| ~ exp(10) via _r_to_v, so
        # this was a real, frequently-triggered bug. Detaching x in the clamped
        # branch makes it a fixed numerical-safety extrapolation (matching every
        # other clamp in this file: zero gradient beyond the safe region) rather
        # than an inconsistent, gradient-exploding mix of clamped-t/unclamped-y.
        clamp_mask = t > 15.0
        t = torch.where(clamp_mask, torch.full_like(t, 15.0), t)
        eps = 1e-3
        is_small = t.abs() < eps
        safe_t = torch.where(is_small, torch.ones_like(t), t)
        expm1_over_t = torch.where(
            is_small,
            1.0 + t / 2.0 + (t ** 2) / 6.0 + (t ** 3) / 24.0,
            torch.expm1(t) / safe_t,
        )
        x = sigma_b * y * expm1_over_t
        return torch.where(clamp_mask, x.detach(), x)

    def inverse_transform(self, x):
        sigma = self.get_sigma().unsqueeze(0)
        gamma = self.get_gamma().unsqueeze(0)

        sigma_b = sigma.expand_as(x)
        gamma_b = gamma.expand_as(x)

        # Reformulate log1p(gamma*x/sigma)/gamma as (x/sigma)*log1p(t)/t where t=gamma*x/sigma.
        # This avoids dividing by gamma, whose 1/gamma^2 backward gradient caused NaN
        # when |gamma| was near eps and quantile_lambda was large.
        w = x / sigma_b
        t = gamma_b * w
        t = torch.clamp(t, min=-1.0 + 1e-7)
        eps = 1e-3
        is_small = t.abs() < eps
        safe_t = torch.where(is_small, torch.ones_like(t), t)
        log1p_over_t = torch.where(
            is_small,
            1.0 - t / 2.0 + (t ** 2) / 3.0 - (t ** 3) / 4.0,
            torch.log1p(t) / safe_t,
        )
        return w * log1p_over_t

    def forward(self, z, reverse=False):
        """
        For normflows compatibility:
          - If reverse=False: we do y->x
          - If reverse=True: we do x->y
        """
        return self.inverse_transform(z) if reverse else self.forward_transform(z)


# ---------------------------------------------------------------------------
# Conditional normalizing flow infrastructure
# ---------------------------------------------------------------------------

class ConditionalMaskedAffineFlow(nn.Module):
    """
    RealNVP coupling layer conditioned on a face-index embedding.

    The scale and shift networks receive torch.cat([z_masked, cond], dim=-1)
    as input, where cond is the (batch, embed_dim) face embedding looked up
    externally and passed in at call time.
    """

    def __init__(self, b, s_net, t_net):
        """
        Args:
            b:     binary mask tensor of shape (latent_size,)
            s_net: scale MLP with input size latent_size + embed_dim
            t_net: shift MLP with input size latent_size + embed_dim
        """
        super().__init__()
        self.register_buffer('b', b.view(1, *b.size()))
        self.s_net = s_net
        self.t_net = t_net

    def forward(self, z, cond):
        """Generative direction (base → data). Returns (z_out, log_det)."""
        z_m = self.b * z
        inp = torch.cat([z_m, cond], dim=-1)
        s = self.s_net(inp)
        nan = torch.tensor(float('nan'), dtype=z.dtype, device=z.device)
        s = torch.where(torch.isfinite(s), s, nan)
        t = self.t_net(inp)
        t = torch.where(torch.isfinite(t), t, nan)
        z_out = z_m + (1 - self.b) * (z * s.exp() + t)
        log_det = ((1 - self.b) * s).sum(dim=1)
        return z_out, log_det

    def inverse(self, x, cond):
        """Normalizing direction (data → base). Returns (x_out, log_det)."""
        x_m = self.b * x
        inp = torch.cat([x_m, cond], dim=-1)
        s = self.s_net(inp)
        nan = torch.tensor(float('nan'), dtype=x.dtype, device=x.device)
        s = torch.where(torch.isfinite(s), s, nan)
        t = self.t_net(inp)
        t = torch.where(torch.isfinite(t), t, nan)
        x_out = x_m + (1 - self.b) * (x - t) * (-s).exp()
        log_det = -((1 - self.b) * s).sum(dim=1)
        return x_out, log_det


class ActNormConditional(nn.Module):
    """
    ActNorm wrapper that accepts (z, cond) to match the conditional flow API.
    The cond argument is ignored; normalization is unconditional.
    """

    def __init__(self, dim):
        super().__init__()
        self.act_norm = nf.flows.ActNorm(dim)

    def forward(self, z, cond):
        return self.act_norm.forward(z)

    def inverse(self, x, cond):
        return self.act_norm.inverse(x)


class ConditionalNormalizingFlow(nn.Module):
    """
    Normalizing flow on R^{d-1} shared across all d faces, conditioned on the
    face index k via a learned embedding injected into every coupling layer.

    A single set of network weights serves all faces; face identity is
    communicated solely through the embedding lookup.
    """

    def __init__(self, d, embed_dim, flows, base):
        """
        Args:
            d:         original data dimension; face indices range over 0..d-1
            embed_dim: embedding dimension for the face index
            flows:     list of ConditionalMaskedAffineFlow / ActNormConditional
            base:      DiagGaussian(d-1) base distribution
        """
        super().__init__()
        self.d = d
        self.embedding = nn.Embedding(d, embed_dim)
        self.flows = nn.ModuleList(flows)
        self.base = base

    def log_prob(self, x, k):
        """
        Evaluate log density of x in R^{d-1} given face index k.

        Args:
            x: (batch, d-1) tensor
            k: (batch,) long tensor of face indices
        Returns:
            (batch,) log densities
        """
        cond = self.embedding(k)                              # (batch, embed_dim)
        z = x
        log_det_acc = torch.zeros(x.shape[0], device=x.device, dtype=x.dtype)
        # Inverting the composition sample() builds (flows[0].forward, then
        # flows[1].forward, ..., flows[-1].forward) requires undoing the
        # layers in the opposite order: flows[-1].inverse first, ...,
        # flows[0].inverse last.
        for flow in reversed(self.flows):
            z, ld = flow.inverse(z, cond)
            log_det_acc = log_det_acc + ld
        return self.base.log_prob(z) + log_det_acc

    def sample(self, n, k):
        """
        Draw samples from the distribution conditioned on face index k.

        Args:
            n: number of samples
            k: (n,) long tensor of face indices
        Returns:
            (n, d-1) samples in R^{d-1}  (before the h mapping to the negative orthant)
        """
        cond = self.embedding(k)                              # (n, embed_dim)
        z = torch.randn(n, self.d - 1, device=k.device)
        for flow in self.flows:
            z, _ = flow.forward(z, cond)
        return z


# ---------------------------------------------------------------------------
# Censoring function for the likelihood-free (MMD) objective
# ---------------------------------------------------------------------------

def censoring_function(z, c=0.0, tau=0.1):
    """
    Smooth (softplus-based) censoring at threshold c:
        C_{c,tau}(z) = c + tau * log(1 + exp((z - c) / tau))
    A soft approximation of max(z, c) that stays differentiable everywhere;
    tau controls how sharp the approximation is (tau -> 0 recovers max(z, c)).
    """
    return c + tau * F.softplus((z - c) / tau)


def soft_indicator(z, c=0.0, tau=0.1):
    """
    Smooth (sigmoid-based) approximation of the exceedance indicator 1(z > c):
        D_{c,tau}(z) = sigmoid((z - c) / tau)
    Shares the threshold c and smoothing scale tau with censoring_function;
    tau -> 0 recovers the hard indicator. Stays differentiable everywhere,
    which the likelihood-free (MMD) path requires since gradients must flow
    through it into the flow network and log_pi_raw.
    """
    return torch.sigmoid((z - c) / tau)


# ---------------------------------------------------------------------------
# Main S-representation mGPD model
# ---------------------------------------------------------------------------

class S_mGPD_NF(nn.Module):
    """
    S-representation of the multivariate GPD with normalizing flows.

    The density of standardized observations Z is:
        log p_Z(z) = log p_S(z - max(z)) - max(z)

    where S = Z - max(Z) lives on the union of d faces
        S = union_k {s : s_k = 0, s_j < 0 for j != k}.

    For face k, S_{-k} is modeled as
        S_{-k} = -exp(R_k),    R_k = flow(U_k, k),    U_k ~ N(0, I_{d-1}).

    The face indicator K ~ Categorical(pi) with learnable pi stored as
    log_pi_raw (softmax-normalized before use).
    """

    def __init__(
        self,
        dim,
        flow,
        device,
        penalty_lambda,
        fix_margin,
        min_exceedances=1,
        marg_lambda=0.0,
        # Likelihood-free (MMD, Eq. 11) training controls
        use_likelihood=True,
        mmd_lambda=0.0,
        mmd_samples_per_face=8,
        mmd_bandwidths=(0.25, 0.5, 1.0, 2.0, 4.0),
        mmd_epsilon=1e-6,
        mmd_censor_c=0.0,
        mmd_censor_tau=0.1,
        mmd_indicator_alpha=1.0,
        mmd_indicator_rho=1.0,
    ):
        super().__init__()
        self.dim = dim
        self.device = device
        self.penalty_lambda = penalty_lambda
        self.min_exceedances = min_exceedances
        self.marg_lambda = marg_lambda

        # Likelihood-free (MMD) training controls
        self.use_likelihood = use_likelihood
        self.mmd_lambda = mmd_lambda
        self.mmd_samples_per_face = mmd_samples_per_face
        self.mmd_bandwidths = mmd_bandwidths
        self.mmd_epsilon = mmd_epsilon
        self.mmd_censor_c = mmd_censor_c
        self.mmd_censor_tau = mmd_censor_tau
        self.mmd_indicator_alpha = mmd_indicator_alpha
        self.mmd_indicator_rho = mmd_indicator_rho

        self.data_transform = DataTransform(dim, device, fix_margin)
        self.flow_model = flow

        # Log-unnormalized face probabilities; pi = softmax(log_pi_raw)
        self.log_pi_raw = nn.Parameter(torch.zeros(dim, device=device))

    def get_sigma(self):
        return self.data_transform.get_sigma()

    def get_gamma(self):
        return self.data_transform.get_gamma()

    # ------------------------------------------------------------------
    # Core density (S representation)
    # ------------------------------------------------------------------

    def log_prob_S_mGPD_std(self, y):
        """
        Log density of standardized mGPD observations under the S representation.

        For observation y with k = argmax(y) and e = max(y):
            s_{-k} = y_{-k} - e  (all < 0)
            r      = log(-s_{-k})          [inverse of h(r) = -exp(r)]
            log p_Z(y) = log pi_k
                         + log p_flow(r, k)
                         + sum_j log(1 / (-s_{-k,j}))   [log |Jac of h^{-1}|]
                         - e

        Args:
            y: (batch, d) standardized observations with max(y) > 0
        Returns:
            (batch,) log densities
        """
        batch, d = y.shape
        e, k = y.max(dim=1)                                    # (batch,) each

        s = y - e.unsqueeze(1)                                 # s_k = 0, s_{-k} < 0

        # Extract s_{-k}: remove the k-th column for each row → (batch, d-1)
        all_idx = torch.arange(d, device=y.device).unsqueeze(0).expand(batch, -1)
        not_k = all_idx != k.unsqueeze(1)                      # (batch, d) bool
        s_neg_k = s.masked_select(not_k).view(batch, d - 1)
        s_neg_k = s_neg_k.clamp(max=-1e-7)                     # guard boundary ties

        # Inverse of h(r) = -exp(r):  r = log(-v)
        r = torch.log(-s_neg_k)                                # (batch, d-1)

        # Log |Jac of h^{-1}|: d/dv[log(-v)] = 1/v → log|1/v| = -log(-v)
        log_jac = -torch.log(-s_neg_k).sum(dim=1)              # (batch,)

        log_pi = torch.log_softmax(self.log_pi_raw, dim=0)     # (d,) numerically stable
        log_p_r = self.flow_model.log_prob(r, k)               # (batch,)

        return log_pi[k] + log_p_r + log_jac - e

    def log_prob(self, x_data, very_negative=-1e30):
        sigma = self.get_sigma().unsqueeze(0)
        gamma = self.get_gamma().unsqueeze(0)
        inside_raw = sigma + gamma * x_data
        valid = (inside_raw > 0).all(dim=1)

        log_prob_x = torch.full(
            (x_data.size(0),),
            float(very_negative),
            device=x_data.device,
            dtype=x_data.dtype,
        )

        if valid.any():
            x_v = x_data[valid]
            inside_v = inside_raw[valid]
            y_v = self.data_transform.inverse_transform(x_v)
            log_prob_y_v = self.log_prob_S_mGPD_std(y_v)
            log_abs_detJ_v = -torch.sum(torch.log(inside_v), dim=1)
            log_prob_x[valid] = log_prob_y_v + log_abs_detJ_v

        return log_prob_x

    # ------------------------------------------------------------------
    # Auxiliary loss term: conditional marginal GPD negative log-likelihood
    # ------------------------------------------------------------------

    def _exceedance_mask(self, x_data, inside_raw):
        mask = (x_data > 0) & (inside_raw > 0)
        active = mask.sum(dim=0) >= self.min_exceedances
        return mask, active

    def marginal_loss(self, x_data, inside_raw):
        """
        Conditional marginal GPD negative log-likelihood:

            L_marg = (1/d) * sum_j (1/N_j) * sum_{i: X_ij>0}
                     [ g_STD(X_ij; sigma_j, gamma_j) + log(sigma_j + gamma_j*X_ij) ]

        where g_STD is DataTransform.inverse_transform.
        """
        dtype = x_data.dtype
        device = x_data.device
        mask, active = self._exceedance_mask(x_data, inside_raw)

        if not active.any():
            return torch.tensor(0.0, device=device, dtype=dtype)

        y = self.data_transform.inverse_transform(x_data)
        log_term = torch.log(inside_raw.clamp(min=1e-12))
        per_element = y + log_term

        loss_per_margin = torch.zeros(self.dim, device=device, dtype=dtype)
        for j in range(self.dim):
            if not active[j]:
                continue
            loss_per_margin[j] = per_element[:, j][mask[:, j]].mean()

        return loss_per_margin[active].mean()

    # ------------------------------------------------------------------
    # Sampling (S representation)
    # ------------------------------------------------------------------

    def _build_s_from_v(self, v, k):
        """
        Construct the d-dimensional S vector by inserting 0 at position k[i]
        and placing v[i] at all other positions, for each row i.

        Uses masked_scatter to preserve the autograd graph through v.

        Args:
            v: (n, d-1) values for the non-k components (typically < 0)
            k: (n,) long face indices
        Returns:
            s: (n, d) with s[i, k[i]] = 0 and s[i, j] = v[i, ...] for j != k[i]
        """
        n, d = v.shape[0], self.dim
        s = torch.zeros(n, d, device=v.device, dtype=v.dtype)
        all_idx = torch.arange(d, device=v.device).unsqueeze(0).expand(n, -1)
        not_k = all_idx != k.unsqueeze(1)                      # (n, d) bool
        # masked_scatter fills True positions row-by-row from v.flatten()
        s = s.masked_scatter(not_k, v.reshape(-1))
        return s

    def _r_to_v(self, r):
        """
        Map flow output r to v = -exp(r), clamped so exp cannot overflow (or
        produce non-finite gradients) if the flow drifts to large |r|, e.g.
        under unconstrained ActNorm scaling during MMD training.
        """
        # exp(10) ~ 2.2e4, already far beyond mmd_bandwidths/mmd_censor_tau scale
        # (both < 5), so the kernel already treats this as a fully saturated
        # extreme point. Kept modest (rather than e.g. 30) so that y = e + s
        # built from v stays small enough for forward_transform's t-clamp to
        # have real headroom, and so x_fake can't overflow when squared inside
        # torch.cdist in the MMD kernel.
        r_max = 10.0
        r_clamped = torch.clamp(r, max=r_max)
        return -torch.exp(r_clamped)

    # ------------------------------------------------------------------
    # Likelihood-free (MMD, extends Eq. 11 with an indicator-augmented kernel) objective
    # ------------------------------------------------------------------

    def sample_with_grad(self, samples_per_face):
        """
        Differentiable sample of Z_theta with the discrete face variable K
        marginalized out exactly (rather than stochastically drawn).

        For every face k = 0..d-1, draws `samples_per_face` reparameterized
        continuous samples from the flow conditioned on k (gradients flow
        into the flow weights and the face embedding). Each face's block of
        samples is assigned weight pi_k / samples_per_face, so the returned
        (points, weights) pair represents the exact mixture
        sum_k pi_k * Q_k as a weighted point cloud, with pi_k = softmax(log_pi_raw)[k]
        differentiable w.r.t. log_pi_raw.

        Returns:
            y_all: (dim*samples_per_face, dim) standardized samples Z_theta
            w_all: (dim*samples_per_face,) mixture weights summing to 1
        """
        n = self.dim * samples_per_face
        k_all = torch.arange(self.dim, device=self.device).repeat_interleave(samples_per_face)

        r_all = self.flow_model.sample(n, k_all)
        v_all = self._r_to_v(r_all)
        s_all = self._build_s_from_v(v_all, k_all)

        e_all = torch.empty(n, 1, device=self.device).exponential_(1.0)
        y_all = e_all + s_all

        pi = torch.softmax(self.log_pi_raw, dim=0)
        w_all = (pi.unsqueeze(1).expand(-1, samples_per_face) / samples_per_face).reshape(-1)

        return y_all, w_all

    def _multiscale_rbf_kernel(self, a, b):
        """
        Sum of RBF kernels exp(-||a_i - b_j||^2 / (2*h^2)) over self.mmd_bandwidths.

        Args:
            a: (n, d), b: (m, d)
        Returns:
            (n, m) kernel matrix
        """
        sq_dist = torch.cdist(a, b) ** 2
        K = sum(torch.exp(-sq_dist / (2.0 * h ** 2)) for h in self.mmd_bandwidths)
        return K

    def _weighted_off_diag_mean(self, K, w):
        """
        Weighted, off-diagonal (i != j) estimate of E[k(X, X')] for X, X'
        drawn independently from a distribution represented by the weighted
        point cloud (K's rows/cols, w). Reduces to the usual unbiased
        U-statistic mean when w is uniform (w_i = 1/n).
        """
        s = w @ K @ w
        diag_term = (torch.diag(K) * w ** 2).sum()
        denom = (1.0 - (w ** 2).sum()).clamp_min(self.mmd_epsilon)
        return (s - diag_term) / denom

    def _indicator_kernel(self, d_a, d_b, rho):
        """
        Dimension-normalized Hamming-style exponential kernel on the soft
        exceedance-indicator vectors:
            k_D(d,d') = exp(-rho * H_norm(d,d')),  H_norm(d,d') = (1/dim) * sum_j |d_j - d'_j|
        |d_j - d'_j| reduces to the Hamming indicator 1(d_j != d'_j) for hard
        0/1 values; normalizing by dim keeps rho's scale comparable across dim.

        Args:
            d_a: (n, dim), d_b: (m, dim)
        Returns:
            (n, m) kernel matrix
        """
        h_norm = (d_a.unsqueeze(1) - d_b.unsqueeze(0)).abs().mean(dim=2)
        return torch.exp(-rho * h_norm)

    def _augmented_kernel(self, z_a, d_a, z_b, d_b):
        """
        K((z,d),(z',d')) = K_Z(z,z'; mmd_bandwidths) + mmd_indicator_alpha * k_D(d,d'; mmd_indicator_rho)

        Additive combination (not gated/multiplicative): z here is a smooth,
        well-defined function of the underlying standardized value everywhere
        (not sentinel-valued below threshold, unlike the older hard-mask
        kernel used for the T representation), so there is no reason to gate
        the value term on indicator agreement. An additive D-block instead
        contributes an independent penalty for co-exceedance-pattern
        mismatches even when censored values already look similar. A sum of
        two PSD kernels stays PSD/characteristic, so this remains a valid
        MMD kernel.
        """
        K_z = self._multiscale_rbf_kernel(z_a, z_b)
        K_d = self._indicator_kernel(d_a, d_b, self.mmd_indicator_rho)
        return K_z + self.mmd_indicator_alpha * K_d

    def mmd_censored_loss(self, x_data, valid):
        """
        Squared MMD between augmented real/generated samples on the
        STANDARDIZED (z) scale, each represented as (C_{c,tau}(Z),
        soft_indicator(Z, c, tau)), so the kernel captures joint
        co-exceedance structure across dimensions, not just similarity of
        the censored continuous values. Matches the paper's literal Eq. 11
        objective (censored continuous values), augmented with the
        exceedance-indicator block described in _augmented_kernel.

        On this scale, exceedance margins are unit-exponential-like by
        construction (Z = E + S, E ~ Exp(1)), i.e. genuinely O(1) -- unlike
        the raw observational (x) scale, whose spread depends on the data's
        arbitrary units/threshold. That's what makes self.mmd_bandwidths'
        fixed default a well-justified multiscale span here, with no
        data-dependent bandwidth heuristic needed.

        The real side is y_real = inverse_transform(x_real), mapping the
        observational data onto the flow's native z-space through
        sigma/gamma. The fake side, y_fake, is already on that z-space scale
        (it's built directly as e_all + s_all in
        sample_with_grad), so no transform is needed there. Only
        the real side stays attached to sigma/gamma in this forward graph
        (no detach), so ordinary autograd carries d(mmd_loss)/d(sigma,
        gamma) through it.

        The discrete face variable is marginalized exactly (see
        sample_with_grad), so gradients flow into log_pi_raw as
        well as the flow weights, through both the censored-value and
        indicator paths of the kernel.
        """
        x_real = x_data[valid]
        y_real = self.data_transform.inverse_transform(x_real)
        yc_real = censoring_function(y_real, self.mmd_censor_c, self.mmd_censor_tau)
        d_real = soft_indicator(y_real, self.mmd_censor_c, self.mmd_censor_tau)

        y_fake, w_fake = self.sample_with_grad(self.mmd_samples_per_face)
        yc_fake = censoring_function(y_fake, self.mmd_censor_c, self.mmd_censor_tau)
        d_fake = soft_indicator(y_fake, self.mmd_censor_c, self.mmd_censor_tau)

        B = yc_real.shape[0]
        w_real = torch.full((B,), 1.0 / B, device=yc_real.device, dtype=yc_real.dtype)

        K_rr = self._augmented_kernel(yc_real, d_real, yc_real, d_real)
        K_ff = self._augmented_kernel(yc_fake, d_fake, yc_fake, d_fake)
        K_rf = self._augmented_kernel(yc_real, d_real, yc_fake, d_fake)

        mmd_sq = (
            self._weighted_off_diag_mean(K_rr, w_real)
            + self._weighted_off_diag_mean(K_ff, w_fake)
            - 2.0 * (w_real @ K_rf @ w_fake)
        )
        return mmd_sq

    # ------------------------------------------------------------------
    # Training objective
    # ------------------------------------------------------------------

    def loss_components(self, x_data, very_negative=-1e30, x_data_aux=None):
        sigma = self.get_sigma().unsqueeze(0)
        gamma = self.get_gamma().unsqueeze(0)

        inside_raw = sigma + gamma * x_data
        valid = (inside_raw > 0).all(dim=1)

        if self.use_likelihood:
            log_prob_x = torch.full(
                (x_data.size(0),),
                float(very_negative),
                device=x_data.device,
                dtype=x_data.dtype,
            )

            if valid.any():
                x_v = x_data[valid]
                inside_v = inside_raw[valid]
                y_v = self.data_transform.inverse_transform(x_v)
                log_prob_y_v = self.log_prob_S_mGPD_std(y_v)
                log_abs_detJ_v = -torch.sum(torch.log(inside_v), dim=1)
                log_prob_x[valid] = log_prob_y_v + log_abs_detJ_v
                nll = -log_prob_x[valid].mean()
            else:
                nll = torch.tensor(0.0, device=x_data.device, dtype=x_data.dtype)
        else:
            # Likelihood-free mode: skip the expensive exact-density computation entirely.
            nll = torch.tensor(0.0, device=x_data.device, dtype=x_data.dtype)

        negative_part = torch.relu(-inside_raw)

        if self.use_likelihood:
            support_loss = (negative_part ** 2).sum(dim=1).mean()
        else:
            # mmd_censored_loss compares samples through a soft censoring at `c`
            # (censoring_function/soft_indicator), which collapses everything below
            # c toward c. gamma_j > 0 violations push the standardized value toward
            # -inf (already collapsed by censoring -> invisible to MMD, no barrier
            # needed). gamma_j < 0 violations push it toward +inf, into the
            # uncensored tail region MMD actually compares -> must be penalized.
            upper_violation = negative_part * (gamma < 0).float()
            support_loss = (upper_violation ** 2).sum(dim=1).mean()

        if x_data_aux is not None:
            x_aux = x_data_aux
            inside_raw_aux = self.get_sigma().unsqueeze(0) + self.get_gamma().unsqueeze(0) * x_aux
        else:
            x_aux = x_data
            inside_raw_aux = inside_raw

        marg_loss = self.marginal_loss(x_aux, inside_raw_aux)

        support_loss_weighted = self.penalty_lambda * support_loss
        marg_loss_weighted    = self.marg_lambda    * marg_loss

        if (not self.use_likelihood) and self.mmd_lambda > 0.0 and valid.any():
            mmd_loss = self.mmd_censored_loss(x_data, valid)
        else:
            mmd_loss = torch.tensor(0.0, device=x_data.device, dtype=x_data.dtype)
        mmd_loss_weighted = self.mmd_lambda * mmd_loss

        if self.use_likelihood:
            total_loss = nll + support_loss_weighted + marg_loss_weighted
        else:
            # Extends Eq. 11: L(eta, theta) = L_margin(eta) + lambda_dep * MMD^2[...],
            # where MMD^2 is now computed on the augmented (censored-value,
            # soft-exceedance-indicator) representation (see mmd_censored_loss /
            # _augmented_kernel), not the literal Eq. 11 censored-value-only kernel.
            total_loss = (
                marg_loss_weighted + mmd_loss_weighted + support_loss_weighted
            )

        return {
            "nll": nll,
            "support_loss": support_loss,
            "support_loss_weighted": support_loss_weighted,
            "marg_loss": marg_loss,
            "marg_loss_weighted": marg_loss_weighted,
            "mmd_loss": mmd_loss,
            "mmd_loss_weighted": mmd_loss_weighted,
            "total_loss": total_loss,
            "num_valid_rows": valid.sum(),
        }

    def forward(self, x_data, very_negative=-1e30, return_components=False, x_data_aux=None):
        comps = self.loss_components(
            x_data, very_negative=very_negative, x_data_aux=x_data_aux
        )
        if return_components:
            return comps
        return comps["total_loss"]

    def sample(self, n_samples=1):
        """
        Sample from the learned distribution in x-space.

        Returns:
            samples_x: (n_samples, d) on the original data scale
            samples_y: (n_samples, d) on the standardized mGPD scale
            samples_s: (n_samples, d) the S component (max is 0, others ≤ 0)
        """
        self.flow_model.eval()
        with torch.no_grad():
            pi = torch.softmax(self.log_pi_raw, dim=0)
            k = torch.multinomial(pi.expand(n_samples, -1), 1).squeeze(1)
            r = self.flow_model.sample(n_samples, k)
            v = self._r_to_v(r)
            samples_s = self._build_s_from_v(v, k)
            samples_E = torch.empty(n_samples, 1, device=self.device).exponential_(1.0)
            samples_y = samples_E + samples_s
            samples_x = self.data_transform.forward_transform(samples_y)
        self.flow_model.train()
        return samples_x, samples_y, samples_s

