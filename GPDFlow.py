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
    ):
        super().__init__()
        self.dim = dim
        self.device = device
        self.penalty_lambda = penalty_lambda
        self.mse_quantile_lambda = mse_quantile_lambda
        self.quantile_min_exceedances = quantile_min_exceedances

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

        support_loss_weighted = self.penalty_lambda      * support_loss
        mse_q_loss_weighted   = self.mse_quantile_lambda * mse_q_loss

        total_loss = nll + support_loss_weighted + mse_q_loss_weighted

        return {
            "nll": nll,
            "support_loss": support_loss,
            "support_loss_weighted": support_loss_weighted,
            "mse_quantile_loss": mse_q_loss,
            "mse_quantile_loss_weighted": mse_q_loss_weighted,
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
