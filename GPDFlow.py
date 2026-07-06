import math
import torch
import torch.nn as nn
import torch.nn.functional as F


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


class T_mGPD_NF(nn.Module):
    """
    Interpreting:
       x = observed data (size=dim)
       y = standardized data
    We have:

      DataTransform.forward_transform(y)   = x
      DataTransform.inverse_transform(x)   = y

    So in 'forward(x)', to get minus log-likelihood of data x:
       1) We do y = inverse_transform(x)
       2) log_p_flow(y) = flow_model.log_prob(y)
       3) The log-det Jacobian of x->y is the sum of log(d y_i / d x_i).
          We'll compute that from the known formula for inverse_transform.
    """

    def __init__(
        self,
        dim,
        flow,
        device,
        s_min,
        s_max,
        num_integration_points,
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
        # MMD hyperparameters
        self.mmd_lambda = mmd_lambda
        self.mmd_zeta = mmd_zeta
        self.mmd_h_m = mmd_h_m
        self.mmd_h_x = mmd_h_x
        self.mmd_alpha = mmd_alpha
        self.mmd_epsilon = mmd_epsilon
        self.mmd_update_margins = mmd_update_margins

        self.data_transform = DataTransform(dim, device, fix_margin)

        self.flow_model = flow
        self.s_values = torch.linspace(
            s_min, s_max, num_integration_points, device=device
        ).reshape(-1, 1, 1)
        self.s_min = s_min
        self.s_max = s_max
        self.num_integration_points = num_integration_points

    def get_sigma(self):
        return self.data_transform.get_sigma()

    def get_gamma(self):
        return self.data_transform.get_gamma()

    def log_integral_f_T(self, data):
        batch_size = data.shape[0]
        dim = data.shape[1]

        x_expanded = data.unsqueeze(0).expand(self.num_integration_points, -1, -1)
        s_expanded = self.s_values.expand(-1, batch_size, 1)
        x_plus_s = (x_expanded + s_expanded).reshape(-1, dim)

        log_f_T = self.flow_model.log_prob(x_plus_s)
        log_integrand = log_f_T.reshape(self.num_integration_points, batch_size)

        max_vals, _ = torch.max(log_integrand, dim=0, keepdim=True)
        stable_exp = torch.exp(log_integrand - max_vals)

        delta_s = (self.s_max - self.s_min) / (self.num_integration_points - 1)
        sum_exp = torch.trapz(stable_exp, dx=delta_s, dim=0)

        log_integral = max_vals.squeeze(0) + torch.log(sum_exp + 1e-20)
        return log_integral

    def log_prob_T_mGPD_std(self, data):
        log_integral = self.log_integral_f_T(data)
        max_T = torch.max(data, dim=1)[0]
        return log_integral - max_T

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
            log_prob_y_v = self.log_prob_T_mGPD_std(y_v)
            log_abs_detJ_v = -torch.sum(torch.log(inside_v), dim=1)
            log_prob_x[valid] = log_prob_y_v + log_abs_detJ_v

        return log_prob_x

    def _exceedance_mask(self, x_data, inside_raw):
        mask = (x_data > 0 ) & (inside_raw > 0)
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

    def sample_with_grad(self, n_samples):
        """
        Differentiable sampling on the standardized exponential (y) scale.

        Returns y = E + T − max(T) where T comes from the flow and E ~ Exp(1).
        Gradients propagate through the flow weights. No σ/γ dependency — the
        returned tensor is already on the standardized mGPD scale (margins ~
        Exp(1) for exceedances).

        Assumption: self.flow_model follows the normflows NormalizingFlow API
        with attributes .q0 (callable returning (z, log_q)) and .flows
        (iterable of flow layers whose __call__ is the generative direction).
        """
        # Sample noise from base distribution — no model parameters here
        z, _ = self.flow_model.q0(n_samples)
        # Apply each flow layer in the generative direction (same order as sample())
        for flow_layer in self.flow_model.flows:
            z, _ = flow_layer(z)
        samples_T = z  # shape (n_samples, dim) — graph intact, no detach

        # mGPD representation: Z = E + T - max(T)
        samples_T_max = torch.max(samples_T, dim=1, keepdim=True)[0]
        samples_T_1 = samples_T - samples_T_max
        # Exponential noise is independent randomness; no gradient needed through it
        samples_E = torch.empty(n_samples, 1, device=self.device).exponential_(1.0)
        samples_y = samples_E + samples_T_1
        return samples_y  # shape (n_samples, dim) — standardized mGPD scale

    @staticmethod
    def _off_diag_mean(K):
        """
        Mean of the off-diagonal elements of a square matrix K.

        Computes (sum(K) - trace(K)) / (B * (B - 1)), which is the unbiased
        U-statistic estimator used in the MMD² formula.
        """
        B = K.shape[0]
        return (K.sum() - torch.diag(K).sum()) / (B * (B - 1))

    def _mixed_kernel(self, m1, x1, m2, x2):
        """
        Compute the B1 × B2 mixed kernel matrix between two batches of
        augmented censored representations.

        Each augmented representation is split into:
          m  — float exceedance mask  (0.0 or 1.0), shape (B, d)
          x̃  — censored values       (x_j if x_j > 0 else ζ), shape (B, d)

        The kernel is:
          k(A, A') = k_m(m, m') + α · k_m(m, m') · k_x(x̃, x̃'; m, m')

        where:
          k_m(m, m') = exp(-‖m - m'‖₁ / h_m)
          S(m, m')   = {j : m_j = 1 AND m'_j = 1}
          k_x(·)     = exp(-Σ_{j∈S}(x̃_j - x̃'_j)² / (2·h_x²·(|S|+ε)))

        When S(m, m') = ∅, the numerator of k_x is 0, so k_x = exp(0) = 1
        and only the mask kernel differentiates the pair.

        Parameters
        ----------
        m1 : Tensor (B1, d)  — float mask for batch 1
        x1 : Tensor (B1, d)  — censored values for batch 1
        m2 : Tensor (B2, d)  — float mask for batch 2
        x2 : Tensor (B2, d)  — censored values for batch 2

        Returns
        -------
        K  : Tensor (B1, B2) — kernel matrix
        """
        # --- Mask kernel k_m ---
        # L1 distance between masks: (B1, B2, d) → sum → (B1, B2)
        diff_m = (m1.unsqueeze(1) - m2.unsqueeze(0)).abs()
        K_m = torch.exp(-diff_m.sum(dim=2) / self.mmd_h_m)  # (B1, B2)

        # --- Value kernel k_x (only on jointly exceeded coordinates) ---
        # joint_mask[i,j,k] = 1 iff m1[i,k]=1 AND m2[j,k]=1
        joint_mask = m1.unsqueeze(1) * m2.unsqueeze(0)       # (B1, B2, d)
        S_size = joint_mask.sum(dim=2)                        # (B1, B2)

        diff_x = x1.unsqueeze(1) - x2.unsqueeze(0)           # (B1, B2, d)
        masked_sq_diff = (joint_mask * diff_x.pow(2)).sum(dim=2)  # (B1, B2)

        denom = 2.0 * (self.mmd_h_x ** 2) * (S_size + self.mmd_epsilon)
        K_x = torch.exp(-masked_sq_diff / denom)             # (B1, B2)

        return K_m + self.mmd_alpha * K_m * K_x              # (B1, B2)

    def mmd_censored_loss(self, x_real):
        """
        Compute the squared MMD between real and generated augmented censored
        samples on the standardized exponential (y) scale, using the mixed
        mask-and-value kernel.

        Both sides are compared on the y = g_std(x; σ, γ) scale where
        exceedance margins are ~ Exp(1).

        Augmented censored representation A(y):
          m_j   = 1{y_j > 0}          (exceedance mask; equivalent to 1{x_j > 0})
          ỹ_j   = y_j  if y_j > 0
                  ζ    otherwise       (ζ = self.mmd_zeta < 0)

        Generated samples (y_fake) are produced on the y-scale directly via
        sample_with_grad (no σ/γ involved).  Real data (x_real) are mapped to
        y-scale via inverse_transform; whether σ/γ gradients flow through this
        mapping is controlled by self.mmd_update_margins:
          True  → full gradient through σ, γ (inverse_transform uses live params)
          False → σ, γ detached inside inverse_transform; no MMD gradient to margins

        The squared MMD (unbiased U-statistic) is:
          MMD² = (1/B(B-1)) Σ_{i≠i'} k(A_i,A_i')
               + (1/B(B-1)) Σ_{j≠j'} k(Â_j,Â_j')
               - (2/B²)     Σ_{i,j}  k(A_i,Â_j)

        Parameters
        ----------
        x_real : Tensor (B, dim) — observed exceedance data (threshold u = 0)

        Returns
        -------
        mmd_sq : scalar Tensor — the estimated MMD²
        """
        B = x_real.shape[0]

        # --- Generated samples on the exponential (y) scale, gradient intact ---
        y_fake = self.sample_with_grad(B)  # (B, dim)

        # --- Real data mapped to the exponential (y) scale ---
        # detach_params=True when mmd_update_margins=False to block σ/γ gradient
        y_real = self.data_transform.inverse_transform(
            x_real, detach_params=not self.mmd_update_margins
        )  # (B, dim)

        # --- Augmented censored representations on y-scale ---
        # Real
        m_real = (y_real > 0).float()
        y_tilde_real = torch.where(
            y_real > 0,
            y_real,
            torch.full_like(y_real, self.mmd_zeta),
        )
        # Fake (detach mask; gradient flows through y_fake values only)
        m_fake = (y_fake > 0).float().detach()
        y_tilde_fake = torch.where(
            y_fake > 0,
            y_fake,
            torch.full_like(y_fake, self.mmd_zeta),
        )

        # --- Kernel matrices ---
        K_rr = self._mixed_kernel(m_real, y_tilde_real, m_real, y_tilde_real)
        K_ff = self._mixed_kernel(m_fake, y_tilde_fake, m_fake, y_tilde_fake)
        K_rf = self._mixed_kernel(m_real, y_tilde_real, m_fake, y_tilde_fake)

        # --- MMD² (unbiased estimator) ---
        mmd_sq = (
            self._off_diag_mean(K_rr)
            + self._off_diag_mean(K_ff)
            - 2.0 * K_rf.mean()
        )
        return mmd_sq

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
            log_prob_y_v = self.log_prob_T_mGPD_std(y_v)
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

        # --- MMD censored loss (skipped when lambda == 0 to save computation) ---
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
        comps = self.loss_components(x_data, very_negative=very_negative, x_data_quantile=x_data_quantile)

        if return_components:
            return comps

        return comps["total_loss"]

    def sample(self, n_samples=1):
        """
        Sample from the learned distribution in x-space by:
          1) Sample y ~ flow
          2) Convert y->x using forward_transform
        """
        self.flow_model.eval()
        samples_T, _ = self.flow_model.sample(n_samples)
        self.flow_model.train()

        samples_T_max = torch.max(samples_T, axis=1, keepdim=True)[0]
        samples_T_1 = samples_T - samples_T_max

        samples_E = torch.empty(n_samples, device=self.device).exponential_(1.0).unsqueeze(1)
        samples_y = samples_E + samples_T_1

        samples_x = self.data_transform.forward_transform(samples_y)
        return samples_x, samples_y, samples_T
