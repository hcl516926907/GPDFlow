"""
Parametric mGPD -- T-representation with a fixed parametric generator.

Same mGPD skeleton as GPDFlow (GPDFlow.py::T_mGPD_NF), but the generator T is one
of five parametric families instead of a normalizing flow.  Standardized-mGPD
log-density

    log h(z) = log \\int f_T(z + s 1) ds  -  max(z)                (CLAUDE.md 3.4)

with the 1-D shift integral evaluated by the trapezoidal rule + log-sum-exp
stabiliser, mirroring GPDFlow.py::T_mGPD_NF.log_integral_f_T.  On the
observational scale the marginal Jacobian -sum_j log(sigma_j + gamma_j x_j) is
added, with z = g_std(x; sigma, gamma).

The margin transform and the learnable (sigma, gamma) are reused verbatim from
GPDFlow_S.DataTransform.
"""

import math
import random

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import scipy.stats as stats
from torch.utils.data import DataLoader

from GPDFlow_S import DataTransform

_LOG2PI = math.log(2.0 * math.pi)
GENERATORS = ('gumbel', 'rev_gumbel', 'rev_exp', 'log_gamma', 'gaussian')


def _inv_softplus(y):
    return math.log(math.expm1(y))


# ---------------------------------------------------------------------------
# Main parametric-generator mGPD model (T representation)
# ---------------------------------------------------------------------------

class ParametricMGPD(nn.Module):
    """mGPD with a parametric generator T;  S = T - max(T).

    generator in {'gumbel','rev_gumbel','rev_exp','log_gamma','gaussian'}.
    Independent-component families: per-margin positive scale/shape + per-margin
    location (location[0] pinned to 0 for identifiability -- the integral over the
    common shift s makes the overall generator location unidentified).
    'gaussian': full Cholesky factor L (Sigma = L L^T) + mean (mean[0] pinned).

    Args:
        dim (int): data dimension.
        generator (str): one of GENERATORS.
        device: torch device.
        s_half_width (float): half-width of the fine trapezoid grid laid around
            the per-row integrand mode.
        num_integration_points (int): number of trapezoid grid points.
        penalty_lambda (float): weight of the marginal-support barrier.
        fix_margin (bool): if True, sigma and gamma are held fixed (see
            DataTransform).
    """

    def __init__(self, dim, generator, device, s_half_width=35.0,
                 num_integration_points=1500, penalty_lambda=1e4,
                 fix_margin=False):
        super().__init__()
        assert generator in GENERATORS, generator
        self.dim = int(dim)
        self.device = device
        self.generator = generator
        self.penalty_lambda = float(penalty_lambda)
        self.s_half_width = float(s_half_width)
        self.num_integration_points = int(num_integration_points)

        self.data_transform = DataTransform(self.dim, device, fix_margin)

        # free generator location, length dim-1; component 0 fixed at 0
        self.gen_loc_free = nn.Parameter(torch.zeros(self.dim - 1, device=device))

        if generator in ('gumbel', 'rev_gumbel', 'rev_exp'):
            self.gen_log_scale = nn.Parameter(
                torch.full((self.dim,), _inv_softplus(2.5), device=device))
        elif generator == 'log_gamma':
            self.gen_log_shape = nn.Parameter(
                torch.full((self.dim,), _inv_softplus(1.0), device=device))
        elif generator == 'gaussian':
            self.gen_L_diag_raw = nn.Parameter(
                torch.full((self.dim,), _inv_softplus(2.5), device=device))
            self.gen_L_offdiag = nn.Parameter(
                torch.zeros(self.dim * (self.dim - 1) // 2, device=device))
            self._tril_idx = torch.tril_indices(self.dim, self.dim, offset=-1,
                                                device=device)

        self.to(device)

    # -- accessors -----------------------------------------------------------
    def get_sigma(self):
        return self.data_transform.get_sigma()

    def get_gamma(self):
        return self.data_transform.get_gamma()

    def _loc(self):
        return F.pad(self.gen_loc_free, (1, 0))              # (dim,)  loc[0] = 0

    def _gen_mode(self):
        # per-component mode of T_j, used only to anchor the integration grid
        mu = self._loc()
        if self.generator == 'log_gamma':                   # d/dt: k - e^{t-mu} = 0
            k = F.softplus(self.gen_log_shape).clamp(0.1, 1e3)
            return mu + torch.log(k)
        return mu                                           # gumbel / rev_gumbel / rev_exp / gaussian

    def _L(self):
        diag = F.softplus(self.gen_L_diag_raw).clamp(0.1, 1e3)
        L = torch.zeros(self.dim, self.dim, device=self.device, dtype=diag.dtype)
        L = L.index_put((self._tril_idx[0], self._tril_idx[1]), self.gen_L_offdiag)
        return L + torch.diag_embed(diag)

    # ------------------------------------------------------------------
    # Generator log-density  log f_T(t),  t: (..., dim) -> (...)
    # ------------------------------------------------------------------

    def log_f_T(self, t):
        g = self.generator
        mu = self._loc()
        if g in ('gumbel', 'rev_gumbel', 'rev_exp'):
            b = F.softplus(self.gen_log_scale).clamp(0.1, 1e3)
            w = (t - mu) / b
            if g == 'gumbel':                                   # max-Gumbel
                logf = -torch.log(b) - w - torch.exp(torch.clamp(-w, max=30.0))
            elif g == 'rev_gumbel':                             # min-Gumbel
                logf = -torch.log(b) + w - torch.exp(torch.clamp(w, max=30.0))
            else:                                               # reverse exp, support t <= mu
                logf = torch.where(t <= mu, -torch.log(b) + w,
                                   torch.full_like(w, -1e30))
            return logf.sum(dim=-1)
        if g == 'log_gamma':                                    # T = log Gamma(k, 1) + mu
            k = F.softplus(self.gen_log_shape).clamp(0.1, 1e3)
            u = t - mu
            logf = -torch.lgamma(k) + k * u - torch.exp(torch.clamp(u, max=30.0))
            return logf.sum(dim=-1)
        # gaussian
        L = self._L()
        r = (t - mu).reshape(-1, self.dim).transpose(0, 1)      # (dim, M)
        y = torch.linalg.solve_triangular(L, r, upper=False)    # (dim, M)
        quad = (y * y).sum(dim=0)
        logdet = torch.log(torch.diagonal(L)).sum()
        logf = -0.5 * self.dim * _LOG2PI - logdet - 0.5 * quad
        return logf.reshape(t.shape[:-1])

    # ------------------------------------------------------------------
    # log \int f_T(z + s 1) ds     data: (batch, dim) -> (batch,)
    # ------------------------------------------------------------------

    def log_integral_f_T(self, data):
        # The integrand  phi(s) = sum_j log f_{T_j}(z_j + s)  is log-concave (hence
        # unimodal) for every generator here.  Two-stage grid:
        #  (1) a wide generator-mode-anchored coarse scan to bracket the per-row
        #      peak (needed e.g. for the reverse-exponential generator, whose hard
        #      support t <= mu pins the mode at s* = min_j(mu_j - z_j));
        #  (2) a Newton step + Laplace curvature estimate at the peak, then a
        #      *per-row* fine trapezoid grid whose half-width tracks the local
        #      standard deviation sd = 1/sqrt(-phi'').  A fixed half-width /
        #      point count (the previous behaviour) silently under-resolves a
        #      narrow peak -- large log_gamma shape k, or an ill-conditioned
        #      gaussian Sigma, especially in high dimension -- which biases the
        #      trapezoid estimate differently for each family and corrupts both
        #      the MLE and the AIC log-likelihood.
        z = data
        zmin = z.min(dim=1).values
        zmax = z.max(dim=1).values
        with torch.no_grad():
            mode = self._gen_mode()
            anchor = (mode[None, :] - z).mean(dim=1)                    # (B,)
            span = (mode.max() - mode.min()).clamp(min=0.0)
            wide = 20.0 + (zmax - zmin) + span                          # (B,)
            u = torch.linspace(0.0, 1.0, 400, device=z.device)
            s_coarse = anchor[None, :] + (2.0 * u[:, None] - 1.0) * wide[None, :]  # (400, B)
            lg_c = self.log_f_T(z[None, :, :] + s_coarse[:, :, None])
            lg_c = torch.nan_to_num(lg_c, nan=-1e30, neginf=-1e30, posinf=-1e30)
            s_hat = torch.gather(s_coarse, 0, lg_c.argmax(0, keepdim=True)).squeeze(0)

            hfd = 1e-2

            def _phi(s):
                return torch.nan_to_num(self.log_f_T(z + s[:, None]),
                                        nan=-1e30, neginf=-1e30, posinf=-1e30)

            # Newton refinement of the per-row peak location (skipped on rows whose
            # finite-difference stencil straddles a hard-support wall, e.g. rev_exp).
            for _ in range(2):
                p0, pp, pm = _phi(s_hat), _phi(s_hat + hfd), _phi(s_hat - hfd)
                d1 = (pp - pm) / (2.0 * hfd)
                d2 = (pp - 2.0 * p0 + pm) / (hfd * hfd)
                wall = (p0 < -1e20) | (pp < -1e20) | (pm < -1e20)
                step = torch.where((d2 < -1e-8) & ~wall, d1 / d2, torch.zeros_like(d1))
                step = torch.maximum(torch.minimum(step, wide), -wide)
                s_hat = s_hat - step
            p0, pp, pm = _phi(s_hat), _phi(s_hat + hfd), _phi(s_hat - hfd)
            d2 = (pp - 2.0 * p0 + pm) / (hfd * hfd)
            wall = (p0 < -1e20) | (pp < -1e20) | (pm < -1e20)
            sd = (-d2).clamp_min(1e-8).rsqrt()                          # (B,)
            half_b = (10.0 * sd).clamp(1.5, self.s_half_width)          # (B,)
            half_b = torch.where(wall, torch.full_like(half_b, self.s_half_width), half_b)
        s_hat = s_hat.detach()
        half_b = half_b.detach()

        n_grid = self.num_integration_points
        unit = torch.linspace(-1.0, 1.0, n_grid, device=z.device)       # (G,)
        s_row = s_hat[None, :] + unit[:, None] * half_b[None, :]         # (G, B)
        t = z[None, :, :] + s_row[:, :, None]                           # (G, B, dim)
        logint = self.log_f_T(t)                                        # (G, B)
        logint = torch.nan_to_num(logint, nan=-1e30, neginf=-1e30, posinf=-1e30)
        m = logint.max(dim=0, keepdim=True).values.clamp(min=-1e29)
        stable = torch.exp(logint - m)                                  # (G, B)
        step_b = 2.0 * half_b / (n_grid - 1)                            # (B,)
        integ = step_b * (stable.sum(dim=0) - 0.5 * (stable[0] + stable[-1]))
        out = m.squeeze(0) + torch.log(integ + 1e-30)

        if getattr(self, 'debug_integral', False):
            with torch.no_grad():
                lap = 0.5 * _LOG2PI + torch.log(sd.clamp_min(1e-8)) + _phi(s_hat)
                ok = torch.isfinite(out) & torch.isfinite(lap) & ~wall
                d = float((out[ok] - lap[ok]).abs().max()) if ok.any() else float('nan')
                print(f'[log_integral_f_T] gen={self.generator}  max|trapz-Laplace|={d:.3g}  '
                      f'sd[min,med,max]=({float(sd.min()):.3g},'
                      f'{float(sd.median()):.3g},{float(sd.max()):.3g})')
        return torch.nan_to_num(out, nan=-1e30, neginf=-1e30, posinf=-1e30)

    def log_prob_T_mGPD_std(self, data):
        return self.log_integral_f_T(data) - data.max(dim=1).values

    # ------------------------------------------------------------------
    # Observational-scale log density
    # ------------------------------------------------------------------

    def _log_prob_chunk(self, x_data, very_negative):
        sigma = self.get_sigma().unsqueeze(0)
        gamma = self.get_gamma().unsqueeze(0)
        inside_raw = sigma + gamma * x_data
        valid = (inside_raw > 0).all(dim=1)
        out = torch.full((x_data.shape[0],), float(very_negative),
                         device=x_data.device, dtype=x_data.dtype)
        if valid.any():
            z = self.data_transform.inverse_transform(x_data[valid])
            lp = self.log_prob_T_mGPD_std(z) - torch.log(inside_raw[valid]).sum(dim=1)
            out = out.masked_scatter(valid, lp)
        return out

    def log_prob(self, x_data, very_negative=-1e30, chunk=256):
        return torch.cat([self._log_prob_chunk(x_data[i:i + chunk], very_negative)
                          for i in range(0, x_data.shape[0], chunk)])

    # ------------------------------------------------------------------
    # Training objective
    # ------------------------------------------------------------------

    def loss_components(self, x_data, very_negative=-1e30):
        sigma = self.get_sigma().unsqueeze(0)
        gamma = self.get_gamma().unsqueeze(0)
        inside_raw = sigma + gamma * x_data
        valid = (inside_raw > 0).all(dim=1)

        support_loss = torch.relu(-inside_raw).pow(2).sum(dim=1).mean()
        support_loss_weighted = self.penalty_lambda * support_loss

        # weak L2 anchor on the (only weakly identified) generator locations,
        # keeps them from drifting far enough to de-stabilise the integration grid
        loc_reg = 1e-3 * self.gen_loc_free.pow(2).sum()

        if valid.any():
            z = self.data_transform.inverse_transform(x_data[valid])
            lp = self.log_prob_T_mGPD_std(z) - torch.log(inside_raw[valid]).sum(dim=1)
            lp = lp[torch.isfinite(lp)]
            # leaky floor: rows the generator can barely represent (lp << FLOOR)
            # are capped so they don't dominate the batch mean, but keep a small
            # gradient so the optimiser can still widen the generator to cover
            # them -- avoids the "excluded row gets no gradient" dead end.
            FLOOR = -100.0
            lp = torch.where(lp > FLOOR, lp, FLOOR + 0.05 * (lp - FLOOR))
            nll = -lp.mean() if lp.numel() > 0 else inside_raw.sum() * 0.0
        else:
            nll = inside_raw.sum() * 0.0

        total_loss = nll + support_loss_weighted + loc_reg
        return {
            'nll': nll,
            'support_loss': support_loss,
            'support_loss_weighted': support_loss_weighted,
            'loc_reg': loc_reg,
            'total_loss': total_loss,
            'num_valid_rows': valid.sum(),
        }

    def forward(self, x_data, very_negative=-1e30, return_components=False):
        comps = self.loss_components(x_data, very_negative=very_negative)
        if return_components:
            return comps
        return comps['total_loss']

    def n_params(self):
        d = self.dim
        # sigma, gamma -- free only when fix_margin=False (DataTransform then
        # stores them as nn.Parameter rather than plain constant tensors)
        margins = 2 * d if isinstance(self.data_transform.log_sigma, nn.Parameter) else 0
        if self.generator == 'gaussian':
            gen = d + d * (d - 1) // 2 + (d - 1)           # L diag + L offdiag + mean
        else:
            gen = d + (d - 1)                              # scale/shape + location
        return margins + gen

    # ------------------------------------------------------------------
    # Sampling
    # ------------------------------------------------------------------

    @torch.no_grad()
    def sample(self, n_samples=1):
        d, dev = self.dim, self.device
        mu = self._loc()
        g = self.generator
        if g in ('gumbel', 'rev_gumbel', 'rev_exp'):
            b = F.softplus(self.gen_log_scale).clamp(0.1, 1e3)
            U = torch.rand(n_samples, d, device=dev).clamp(1e-12, 1 - 1e-12)
            if g == 'gumbel':
                samples_T = mu + b * (-torch.log(-torch.log(U)))
            elif g == 'rev_gumbel':
                samples_T = mu + b * torch.log(-torch.log(U))
            else:
                samples_T = mu - b * (-torch.log(U))
        elif g == 'log_gamma':
            k = F.softplus(self.gen_log_shape).clamp(0.1, 1e3).cpu().numpy()
            G = np.random.gamma(shape=k[None, :], scale=1.0, size=(n_samples, d))
            samples_T = mu + torch.log(
                torch.tensor(G, dtype=mu.dtype, device=dev).clamp_min(1e-30))
        else:
            L = self._L()
            Z = torch.randn(n_samples, d, device=dev)
            samples_T = mu + Z @ L.transpose(0, 1)
        samples_T = samples_T - samples_T.max(dim=1, keepdim=True).values
        samples_E = torch.empty(n_samples, 1, device=dev).exponential_(1.0)
        samples_y = samples_E + samples_T
        samples_x = self.data_transform.forward_transform(samples_y)
        return samples_x, samples_y, samples_T


# ---------------------------------------------------------------------------
# SGD maximum-likelihood fit of one parametric-generator mGPD
# ---------------------------------------------------------------------------

def fit_parametric_mgpd(generator, train_arr, val_arr, device,
                        epochs=300, batch_size=256, seed=1234, patience=40,
                        verbose=True, fix_margin=False, stop_at_epoch=None):
    """Fit one parametric-generator mGPD by stochastic-gradient MLE.

    Args:
        generator (str): one of GENERATORS.
        train_arr, val_arr (ndarray): n x d training / validation exceedances
            (threshold u = 0).
        device: torch device.
        epochs, batch_size, seed, patience: optimisation controls.
        verbose (bool): print progress every 50 epochs.
        fix_margin (bool): if True, sigma/gamma are held at (1, 0) and only the
            generator is fitted -- for data already on the standardized mGPD
            scale (the parametric analogue of a two-stage GPDFlow fit).
        stop_at_epoch (int | None): if set, hard-stop after this many epochs
            regardless of the patience rule (``epochs`` -- hence the cosine
            ``T_max`` -- is left unchanged).  Used for the final full-data refit
            at the median of the per-fold CV early-stop epochs.

    Returns:
        (model, best_val, hist): the fitted ParametricMGPD (loaded with the
        best-validation-NLL state), that best validation NLL, and a history dict
        (keys ``train_batch``, ``val``, ``best_epoch``).
    """
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    dim = train_arr.shape[1]
    model = ParametricMGPD(dim, generator, device, fix_margin=fix_margin).to(device)

    # marginal init from a per-margin GPD fit on the positive excesses (cf. cell 21);
    # skipped when fix_margin=True (data already standardized, sigma/gamma pinned).
    if not fix_margin:
        init_sigma, init_gamma = np.ones(dim), np.zeros(dim)
        for j in range(dim):
            col = train_arr[train_arr[:, j] > 0, j]
            if len(col) > 10:
                c, _, sc = stats.genpareto.fit(col, floc=0)
                init_gamma[j] = float(np.clip(c, -0.7, 0.7))
                init_sigma[j] = float(max(sc, 1e-3))
        with torch.no_grad():
            model.data_transform.log_sigma.copy_(
                torch.log(torch.tensor(init_sigma, dtype=torch.float32, device=device)))
            model.data_transform.theta.copy_(
                torch.atanh(torch.tensor(init_gamma, dtype=torch.float32, device=device)))

    train_data = torch.tensor(np.asarray(train_arr), dtype=torch.float32, device=device)
    val_data = torch.tensor(np.asarray(val_arr), dtype=torch.float32, device=device)
    loader = DataLoader(train_data, batch_size=batch_size, shuffle=True, drop_last=False)

    gen_groups = [{'params': [p], 'lr': 5e-3 if 'offdiag' in n else 1e-2}
                  for n, p in model.named_parameters() if n.startswith('gen_')]
    marg_groups = ([] if fix_margin else
                   [{'params': [model.data_transform.log_sigma], 'lr': 1e-2},
                    {'params': [model.data_transform.theta], 'lr': 1e-2}])
    optimizer = torch.optim.Adam(marg_groups + gen_groups)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, epochs)

    best_val, best_state, bad, best_epoch = float('inf'), None, 0, 0
    hist = {'train_batch': [], 'val': []}
    for epoch in range(epochs):
        model.train()
        last = float('nan')
        for xb in loader:
            loss = model.loss_components(xb)['total_loss']
            if not torch.isfinite(loss):
                bad = patience
                break
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            last = float(loss.detach())
        scheduler.step()

        model.eval()
        with torch.no_grad():
            vlp = model.log_prob(val_data)
            # Score every finite validation row, flooring (not dropping) the rows the
            # generator represents poorly.  Dropping them lets a family early-stop on a
            # state that fits fewer rows better, making `best_val` non-comparable across
            # generators (and vs the AIC cell, which floors rather than drops).
            VAL_FLOOR = -50.0
            vf = vlp[torch.isfinite(vlp)]
            vf = torch.where(vf > VAL_FLOOR, vf, torch.full_like(vf, VAL_FLOOR))
            vnll = float(-vf.mean()) if vf.numel() else float('inf')
        hist['val'].append(vnll); hist['train_batch'].append(last)

        if vnll < best_val - 1e-4:
            best_val, bad = vnll, 0
            best_epoch = epoch + 1
            best_state = {k: v.detach().cpu().clone()
                          for k, v in model.state_dict().items()}
        else:
            bad += 1
        if verbose and (epoch % 50 == 0 or epoch == epochs - 1 or bad >= patience):
            print(f'  [{generator:10s}] ep {epoch:3d}  batch loss {last:10.3f}'
                  f'   val NLL {vnll:9.3f}   (best {best_val:9.3f})')
        if bad >= patience:
            break
        if stop_at_epoch is not None and (epoch + 1) >= stop_at_epoch:
            if verbose:
                print(f'  [{generator:10s}] stopping at fixed epoch {epoch + 1} '
                      f'(stop_at_epoch).')
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    hist['best_epoch'] = best_epoch
    return model, best_val, hist


# ---------------------------------------------------------------------------
# Model-selection helper
# ---------------------------------------------------------------------------

def support_mask(model, x_data):
    """Boolean (n,) mask of rows inside the fitted GPD marginal support of model,
    i.e. sigma_j + gamma_j x_j > 0 for every margin j."""
    with torch.no_grad():
        s = model.get_sigma().unsqueeze(0)
        gm = model.get_gamma().unsqueeze(0)
    return ((s + gm * x_data) > 0).all(dim=1)
