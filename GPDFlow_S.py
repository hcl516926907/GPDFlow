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
        gamma = 0.5 * torch.tanh(self.theta)
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
        eps = 1e-3
        is_small = t.abs() < eps
        safe_t = torch.where(is_small, torch.ones_like(t), t)
        expm1_over_t = torch.where(
            is_small,
            1.0 + t / 2.0 + (t ** 2) / 6.0 + (t ** 3) / 24.0,
            torch.expm1(t) / safe_t,
        )
        return sigma_b * y * expm1_over_t

    def inverse_transform(self, x, detach_params=False):
        sigma = self.get_sigma().unsqueeze(0)
        gamma = self.get_gamma().unsqueeze(0)
        if detach_params:
            sigma = sigma.detach()
            gamma = gamma.detach()

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
        for flow in self.flows:
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
        mse_quantile_lambda=0.0,
        quantile_min_exceedances=1,
        # MMD censored-loss hyperparameters
        mmd_lambda=0.0,
        mmd_zeta=-1.0,
        mmd_h_m=1.0,
        mmd_h_x=1.0,
        mmd_alpha=1.0,
        mmd_epsilon=1e-6,
        mmd_update_margins=True,
    ):
        super().__init__()
        self.dim = dim
        self.device = device
        self.penalty_lambda = penalty_lambda
        self.mse_quantile_lambda = mse_quantile_lambda
        self.quantile_min_exceedances = quantile_min_exceedances
        self.mmd_lambda = mmd_lambda
        self.mmd_zeta = mmd_zeta
        self.mmd_h_m = mmd_h_m
        self.mmd_h_x = mmd_h_x
        self.mmd_alpha = mmd_alpha
        self.mmd_epsilon = mmd_epsilon
        self.mmd_update_margins = mmd_update_margins

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
    # Auxiliary loss terms (unchanged from T representation)
    # ------------------------------------------------------------------

    def _exceedance_mask(self, x_data, inside_raw):
        mask = (x_data > 0) & (inside_raw > 0)
        active = mask.sum(dim=0) >= self.quantile_min_exceedances
        return mask, active

    def mse_quantile_loss(self, x_data, inside_raw):
        dtype = x_data.dtype
        device = x_data.device
        mask, active = self._exceedance_mask(x_data, inside_raw)

        if not active.any():
            return torch.tensor(0.0, device=device, dtype=dtype)

        y = self.data_transform.inverse_transform(x_data)

        loss_per_margin = torch.zeros(self.dim, device=device, dtype=dtype)
        for j in range(self.dim):
            if not active[j]:
                continue
            yj_sorted, _ = torch.sort(y[:, j][mask[:, j]])
            n_j = yj_sorted.shape[0]
            p = torch.arange(1, n_j + 1, device=device, dtype=dtype) / (n_j + 1)
            t = -torch.log(1.0 - p)
            loss_per_margin[j] = ((yj_sorted - t) ** 2).mean()

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

    def sample_with_grad(self, n_samples):
        """
        Differentiable sampling on the standardized (y) scale for use in the
        MMD loss.

        Samples K ~ Cat(pi), draws R from the conditional flow for face K,
        maps S_{-K} = -exp(R) to the negative orthant, then returns
        Y = E + S where E ~ Exp(1).  Gradients propagate through the flow
        weights and face embeddings; no gradient through E or K.

        Returns:
            y: (n_samples, d) on the standardized mGPD scale
        """
        pi = torch.softmax(self.log_pi_raw, dim=0)
        k = torch.multinomial(pi.expand(n_samples, -1), 1).squeeze(1)  # (n,)

        r = self.flow_model.sample(n_samples, k)               # (n, d-1) in R^{d-1}
        v = -torch.exp(r)                                      # (n, d-1) in (-inf, 0)

        s = self._build_s_from_v(v, k)                         # (n, d)

        samples_E = torch.empty(n_samples, 1, device=self.device).exponential_(1.0)
        return samples_E + s                                   # (n, d)

    # ------------------------------------------------------------------
    # MMD censored loss (unchanged from T representation)
    # ------------------------------------------------------------------

    @staticmethod
    def _off_diag_mean(K):
        B = K.shape[0]
        return (K.sum() - torch.diag(K).sum()) / (B * (B - 1))

    def _mixed_kernel(self, m1, x1, m2, x2):
        diff_m = (m1.unsqueeze(1) - m2.unsqueeze(0)).abs()
        K_m = torch.exp(-diff_m.sum(dim=2) / self.mmd_h_m)

        joint_mask = m1.unsqueeze(1) * m2.unsqueeze(0)
        S_size = joint_mask.sum(dim=2)

        diff_x = x1.unsqueeze(1) - x2.unsqueeze(0)
        masked_sq_diff = (joint_mask * diff_x.pow(2)).sum(dim=2)

        denom = 2.0 * (self.mmd_h_x ** 2) * (S_size + self.mmd_epsilon)
        K_x = torch.exp(-masked_sq_diff / denom)

        return K_m + self.mmd_alpha * K_m * K_x

    def mmd_censored_loss(self, x_real):
        B = x_real.shape[0]

        y_fake = self.sample_with_grad(B)

        y_real = self.data_transform.inverse_transform(
            x_real, detach_params=not self.mmd_update_margins
        )

        m_real = (y_real > 0).float()
        y_tilde_real = torch.where(
            y_real > 0,
            y_real,
            torch.full_like(y_real, self.mmd_zeta),
        )
        m_fake = (y_fake > 0).float().detach()
        y_tilde_fake = torch.where(
            y_fake > 0,
            y_fake,
            torch.full_like(y_fake, self.mmd_zeta),
        )

        K_rr = self._mixed_kernel(m_real, y_tilde_real, m_real, y_tilde_real)
        K_ff = self._mixed_kernel(m_fake, y_tilde_fake, m_fake, y_tilde_fake)
        K_rf = self._mixed_kernel(m_real, y_tilde_real, m_fake, y_tilde_fake)

        mmd_sq = (
            self._off_diag_mean(K_rr)
            + self._off_diag_mean(K_ff)
            - 2.0 * K_rf.mean()
        )
        return mmd_sq

    # ------------------------------------------------------------------
    # Training objective
    # ------------------------------------------------------------------

    def loss_components(self, x_data, very_negative=-1e30, x_data_quantile=None):
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
            nll = -log_prob_x[valid].mean()
        else:
            nll = torch.tensor(0.0, device=x_data.device, dtype=x_data.dtype)

        negative_part = torch.relu(-inside_raw)
        support_loss = (negative_part ** 2).sum(dim=1).mean()

        if x_data_quantile is not None:
            x_q = x_data_quantile
            inside_raw_q = self.get_sigma().unsqueeze(0) + self.get_gamma().unsqueeze(0) * x_q
        else:
            x_q = x_data
            inside_raw_q = inside_raw

        mse_q_loss = self.mse_quantile_loss(x_q, inside_raw_q)

        support_loss_weighted = self.penalty_lambda       * support_loss
        mse_q_loss_weighted   = self.mse_quantile_lambda  * mse_q_loss

        if self.mmd_lambda > 0.0:
            mmd_loss = self.mmd_censored_loss(x_data)
        else:
            mmd_loss = torch.tensor(0.0, device=x_data.device, dtype=x_data.dtype)
        mmd_loss_weighted = self.mmd_lambda * mmd_loss

        total_loss = nll + support_loss_weighted + mse_q_loss_weighted + mmd_loss_weighted

        return {
            "nll": nll,
            "support_loss": support_loss,
            "support_loss_weighted": support_loss_weighted,
            "mse_quantile_loss": mse_q_loss,
            "mse_quantile_loss_weighted": mse_q_loss_weighted,
            "mmd_loss": mmd_loss,
            "mmd_loss_weighted": mmd_loss_weighted,
            "total_loss": total_loss,
            "num_valid_rows": valid.sum(),
        }

    def forward(self, x_data, very_negative=-1e30, return_components=False, x_data_quantile=None):
        comps = self.loss_components(
            x_data, very_negative=very_negative, x_data_quantile=x_data_quantile
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
            v = -torch.exp(r)
            samples_s = self._build_s_from_v(v, k)
            samples_E = torch.empty(n_samples, 1, device=self.device).exponential_(1.0)
            samples_y = samples_E + samples_s
            samples_x = self.data_transform.forward_transform(samples_y)
        self.flow_model.train()
        return samples_x, samples_y, samples_s


# ---------------------------------------------------------------------------
# Factory function
# ---------------------------------------------------------------------------

def build_S_mGPD_model(dim, device, embed_dim=8, num_layers=16, seed=0):
    """
    Construct an S_mGPD_NF model with the default architecture used in
    the simulation section of the paper.

    Args:
        dim:        data dimension (d >= 2)
        device:     torch device
        embed_dim:  dimension of the face-index embedding (default 8)
        num_layers: number of coupling layers (default 16)
        seed:       random seed for weight initialization

    Returns:
        model: S_mGPD_NF instance on the given device
    """
    latent_size = dim - 1                                       # flow operates on R^{d-1}
    base = nf.distributions.DiagGaussian(latent_size)

    torch.manual_seed(seed)
    b = torch.Tensor([1 if i % 2 == 0 else 0 for i in range(latent_size)])

    flows = []
    for i in range(num_layers):
        mask = b if i % 2 == 0 else (1 - b)
        s_net = nf.nets.MLP(
            [latent_size + embed_dim, 4 * latent_size, latent_size],
            init_zeros=True,
            output_fn='tanh',
        )
        t_net = nf.nets.MLP(
            [latent_size + embed_dim, 4 * latent_size, latent_size],
            init_zeros=True,
            output_fn='tanh',
        )
        flows.append(ConditionalMaskedAffineFlow(mask, s_net, t_net))
        flows.append(ActNormConditional(latent_size))

    flow_model = ConditionalNormalizingFlow(
        d=dim, embed_dim=embed_dim, flows=flows, base=base
    ).to(device)

    model = S_mGPD_NF(
        dim=dim,
        flow=flow_model,
        device=device,
        penalty_lambda=10000,
        fix_margin=False,
    )
    return model
