"""5-fold cross-validation utilities for the GPDFlow FX application.

This module factors the three inline training loops of
``Exchange_rate_application.ipynb`` -- the likelihood-based GPDFlow (cell 19), the
two-stage MMD-based GPDFlow (cell 32) and the probit-copula RealNVP baseline
(cell 42) -- into reusable ``fit_*`` functions (mirroring
``parametric_mGPD.fit_parametric_mgpd``: reset RNGs, train with the *unchanged*
early-stopping rule, restore the best-checkpoint state, return
``(model, best_val, hist)``), plus ``run_cv_*`` drivers that

  * partition the exceedance slice ``samples`` into 5 folds (``make_folds``),
  * for fold k, train on the other 4 folds and early-stop on fold k,
  * draw ``samples_multiplier * |fold k|`` simulated observations from the fold-k
    model, and
  * return the 5 out-of-fold simulated blocks **separately** (``*_folds`` lists,
    raw scale).  The drivers compute no tail-dependence summary -- form ``chi`` /
    ``omega`` in the notebook (fold-mean of the per-fold estimates) so the
    evaluation threshold stays a free knob.

The parametric-mGPD driver (``run_cv_parametric_mgpd``) follows the same recipe
for a fixed AIC-selected generator family.

The likelihood driver also returns the fold-averaged marginal parameter
estimates for the marginal-parameter scatter plot.
"""

import math
import random

import numpy as np
import pandas as pd
import torch
import normflows as nf
import scipy.stats as stats
from scipy.stats import norm
from torch.utils.data import DataLoader

from GPDFlow_S import (
    S_mGPD_NF,
    ConditionalMaskedAffineFlow,
    ActNormConditional,
    ConditionalNormalizingFlow,
)
from parametric_mGPD import fit_parametric_mgpd
import Common_Functions as cf


# ---------------------------------------------------------------------------
# Fold construction
# ---------------------------------------------------------------------------

def make_folds(n, n_splits=5, seed=1234):
    """Partition ``range(n)`` into ``n_splits`` cross-validation folds.

    Uses the notebook's canonical ``np.random.default_rng(seed).permutation``
    idiom, then cuts the permutation into ``n_splits`` near-equal contiguous
    chunks.  Fold k holds out chunk k (``val_idx``) and trains on the union of
    the others (``train_idx``).

    Returns:
        list of ``(train_idx, val_idx)`` int64 arrays.  The ``val_idx`` blocks
        are disjoint and their union is exactly ``range(n)``, so the 5 simulated
        blocks concatenate back to ``n`` rows.
    """
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    chunks = np.array_split(perm, n_splits)
    folds = []
    for k in range(n_splits):
        val_idx = np.sort(chunks[k]).astype(np.int64)
        train_idx = np.sort(
            np.concatenate([chunks[j] for j in range(n_splits) if j != k])
        ).astype(np.int64)
        folds.append((train_idx, val_idx))
    return folds


# ---------------------------------------------------------------------------
# Model / transform factories (fresh module per fold)
# ---------------------------------------------------------------------------

def build_conditional_flow(dim, device, *, num_layers=4, latent_size_factor=2,
                           embed_dim_mult=3, seed=0):
    """Rebuild the exact S-representation flow of notebook cells 19 / 32.

    ``torch.manual_seed(seed)`` is called immediately before construction (as in
    both training cells) so the ``nn.Embedding`` / MLP init is identical every
    fold.
    """
    latent_size = dim - 1
    embed_dim = embed_dim_mult * latent_size
    base_s = nf.distributions.DiagGaussian(latent_size)

    torch.manual_seed(seed)
    b = torch.Tensor([1 if i % 2 == 0 else 0 for i in range(latent_size)])
    flows = []
    for i in range(num_layers):
        mask = b if i % 2 == 0 else (1 - b)
        s_net = nf.nets.MLP(
            [latent_size + embed_dim, latent_size_factor * latent_size, latent_size],
            init_zeros=True, output_fn='tanh',
        )
        t_net = nf.nets.MLP(
            [latent_size + embed_dim, latent_size_factor * latent_size, latent_size],
            init_zeros=True, output_fn='tanh',
        )
        flows.append(ConditionalMaskedAffineFlow(mask, s_net, t_net))
        flows.append(ActNormConditional(latent_size))

    return ConditionalNormalizingFlow(
        d=dim, embed_dim=embed_dim, flows=flows, base=base_s
    ).to(device)


def build_baseline_flow(dim, device, *, num_layers=4, hidden_mult=2):
    """Rebuild the unconditional probit-copula RealNVP of notebook cell 42.

    The caller is responsible for seeding torch/numpy/random beforehand (cell 42
    seeds once and does not reseed per layer).
    """
    nf_base = nf.distributions.DiagGaussian(dim)
    hidden = hidden_mult * dim
    nf_b = torch.tensor([1 if i % 2 == 0 else 0 for i in range(dim)], dtype=torch.float32)
    nf_flows = []
    for i in range(num_layers):
        s = nf.nets.MLP([dim, hidden, dim], init_zeros=True, output_fn='tanh')
        t = nf.nets.MLP([dim, hidden, dim], init_zeros=False)
        mask = nf_b if i % 2 == 0 else (1 - nf_b)
        nf_flows.append(nf.flows.MaskedAffineFlow(mask, t, s))
        nf_flows.append(nf.flows.ActNorm(dim))
    return nf.NormalizingFlow(nf_base, nf_flows).to(device)


def build_probit_copula(train_samples, dim, n_ref=None):
    """Reproduce the hybrid marginal model of notebook cell 41, fit on
    ``train_samples`` (shifted-slice coordinates).

    Per margin j: empirical body CDF for ``x <= 0``, univariate-GPD tail
    (``scipy.stats.genpareto.fit`` on the positive excesses) for ``x > 0``.

    Returns:
        ``(to_gaussian, from_gaussian, marginal_cdf, marginal_ppf)`` -- closures
        over the fitted per-margin parameters.
    """
    train_samples = np.asarray(train_samples, dtype=float)
    if n_ref is None:
        n_ref = len(train_samples)

    body_sorted = [np.sort(train_samples[:, j]) for j in range(dim)]
    p0 = np.array([(train_samples[:, j] > 0).mean() for j in range(dim)])
    gpd_sigma = np.empty(dim)
    gpd_gamma = np.empty(dim)
    for j in range(dim):
        col = train_samples[train_samples[:, j] > 0, j]
        c, _loc, scale = stats.genpareto.fit(col, floc=0)
        gpd_gamma[j] = c
        gpd_sigma[j] = scale

    def marginal_cdf(x, j):
        x = np.asarray(x, dtype=float)
        n = len(body_sorted[j])
        body = np.searchsorted(body_sorted[j], x, side='right') / (n + 1)
        tail = (1.0 - p0[j]) + p0[j] * stats.genpareto.cdf(
            np.clip(x, 0.0, None), c=gpd_gamma[j], scale=gpd_sigma[j]
        )
        return np.where(x > 0.0, tail, body)

    def marginal_ppf(u, j):
        u = np.asarray(u, dtype=float)
        cut = 1.0 - p0[j]
        body = np.quantile(body_sorted[j], np.clip(u, 0.0, 1.0))
        tail_u = np.clip((u - cut) / p0[j], 0.0, 1.0 - 1e-9)
        tail = stats.genpareto.ppf(tail_u, c=gpd_gamma[j], scale=gpd_sigma[j])
        return np.where(u > cut, tail, body)

    def to_gaussian(X):
        X = np.asarray(X, dtype=float)
        U = np.column_stack([marginal_cdf(X[:, j], j) for j in range(dim)])
        U = np.clip(U, 1.0 / (2 * n_ref), 1.0 - 1.0 / (2 * n_ref))
        return norm.ppf(U)

    def from_gaussian(W):
        W = np.asarray(W, dtype=float)
        U = norm.cdf(W)
        return np.column_stack([marginal_ppf(U[:, j], j) for j in range(dim)])

    return to_gaussian, from_gaussian, marginal_cdf, marginal_ppf


# ---------------------------------------------------------------------------
# Epoch-statistics helpers (shared by the two S_mGPD_NF fitters)
# ---------------------------------------------------------------------------

_LIK_KEYS = ('nll', 'support_loss', 'support_loss_weighted',
             'marg_loss', 'marg_loss_weighted', 'total_loss')
_MMD_KEYS = ('nll', 'support_loss', 'support_loss_weighted',
             'marg_loss', 'marg_loss_weighted',
             'mmd_loss', 'mmd_loss_weighted', 'total_loss')


def _init_stats(keys):
    d = {k: 0.0 for k in keys}
    d['count'] = 0
    return d


def _update_stats(stats_dict, comps, bsz, keys):
    for k in keys:
        stats_dict[k] += comps[k].item() * bsz
    stats_dict['count'] += bsz


def _finalize_stats(stats_dict):
    if stats_dict['count'] == 0:
        return stats_dict
    return {k: (v if k == 'count' else v / stats_dict['count'])
            for k, v in stats_dict.items()}


def _reset_seeds(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)


# ---------------------------------------------------------------------------
# 1. Likelihood-based GPDFlow  (from notebook cell 19)
# ---------------------------------------------------------------------------

def fit_gpdflow_likelihood(train_arr, val_arr, device, *, epochs=500,
                           batch_size=256, patience=100, min_delta=0.0,
                           penalty_lambda=1e6, marg_lambda=500.0,
                           num_layers=4, latent_size_factor=2,
                           seed=1234, flow_seed=0, grad_clip=1.0,
                           init_from_gpd=False, stop_at_epoch=None,
                           verbose=True, detect_anomaly=False):
    """Fit the likelihood-based GPDFlow on ``train_arr``, early-stopping on the
    validation NLL over ``val_arr`` (the exact rule of cell 19: patience 100 on
    ``val_stats['nll']``, *not* the total weighted loss).

    ``stop_at_epoch`` (optional): hard-stop after this many epochs regardless of
    the patience rule, while leaving ``epochs`` (hence the ``CosineAnnealingLR``
    ``T_max``) unchanged.  Used for the final full-data refit, which stops at the
    median of the per-fold CV early-stop epochs but keeps the CV LR schedule.

    Returns:
        ``(model, best_val, hist)`` where ``best_val`` is the best epoch-mean
        validation NLL and ``hist`` holds the per-epoch trajectories.
    """
    torch.autograd.set_detect_anomaly(bool(detect_anomaly))
    dim = train_arr.shape[1]

    flow_model = build_conditional_flow(
        dim, device, num_layers=num_layers,
        latent_size_factor=latent_size_factor, seed=flow_seed,
    )
    model = S_mGPD_NF(
        dim=dim, flow=flow_model, device=device,
        penalty_lambda=penalty_lambda, fix_margin=False, marg_lambda=marg_lambda,
    )
    optimizer = torch.optim.Adam([
        {'params': model.flow_model.parameters(), 'lr': 1e-3},
        {'params': model.data_transform.log_sigma, 'lr': 1e-2},
        {'params': model.data_transform.theta, 'lr': 1e-2},
        {'params': [model.log_pi_raw], 'lr': 1e-3},
    ])
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, epochs)

    _reset_seeds(seed)

    # Per-margin scipy GPD fit. Cell 19 computes and prints these but does NOT
    # copy them into the model (the DataTransform starts at sigma=1, gamma=0);
    # that stays the default (init_from_gpd=False) so the pre-CV numbers still
    # hold.  With init_from_gpd=True the fitted (sigma, gamma) seed the margins.
    init_sigmas, init_gammas = np.ones(dim), np.zeros(dim)
    for j in range(dim):
        col = train_arr[train_arr[:, j] > 0, j]
        if len(col) > 10:
            c, _loc, sc = stats.genpareto.fit(col, floc=0)
            init_gammas[j] = float(c)
            init_sigmas[j] = float(sc)
    if verbose:
        print(f'  scipy margin init  gamma~{np.round(init_gammas, 3)}')
        print(f'  scipy margin init  sigma~{np.round(init_sigmas, 3)}')
    if init_from_gpd:
        with torch.no_grad():
            model.data_transform.log_sigma.copy_(torch.log(
                torch.tensor(np.clip(init_sigmas, 1e-3, None),
                             dtype=torch.float32, device=device)))
            model.data_transform.theta.copy_(torch.atanh(
                torch.tensor(np.clip(init_gammas, -0.7, 0.7),
                             dtype=torch.float32, device=device)))

    train_data = torch.tensor(np.asarray(train_arr), dtype=torch.float32, device=device)
    val_data = torch.tensor(np.asarray(val_arr), dtype=torch.float32, device=device)
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True, drop_last=False)
    val_loader = DataLoader(val_data, batch_size=batch_size, shuffle=False, drop_last=False)

    best_val = float('inf')
    best_state = None
    best_epoch = 0
    bad_epochs = 0
    hist = {'train_total': [], 'val_total': [], 'val_nll': [], 'loss_history': []}

    for epoch in range(epochs):
        model.train()
        tr = _init_stats(_LIK_KEYS)
        for x_data in train_loader:
            comps = model(x_data, return_components=True)
            batch_loss = comps['total_loss']
            if torch.isnan(batch_loss) or torch.isinf(batch_loss):
                if verbose:
                    print(f'  [epoch {epoch}] NaN/Inf in TRAIN loss -- stopping.')
                bad_epochs = patience
                break
            optimizer.zero_grad(set_to_none=True)
            batch_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip)
            optimizer.step()
            _update_stats(tr, comps, x_data.shape[0], _LIK_KEYS)
            hist['loss_history'].append(batch_loss.item())
        tr = _finalize_stats(tr)
        hist['train_total'].append(tr['total_loss'] if tr['count'] else float('inf'))

        model.eval()
        va = _init_stats(_LIK_KEYS)
        with torch.no_grad():
            for x_val in val_loader:
                comps = model(x_val, return_components=True)
                if torch.isnan(comps['total_loss']) or torch.isinf(comps['total_loss']):
                    if verbose:
                        print(f'  [epoch {epoch}] NaN/Inf in VAL loss -- stopping.')
                    bad_epochs = patience
                    break
                _update_stats(va, comps, x_val.shape[0], _LIK_KEYS)
        va = _finalize_stats(va)
        hist['val_total'].append(va['total_loss'] if va['count'] else float('inf'))
        hist['val_nll'].append(va['nll'] if va['count'] else float('inf'))

        scheduler.step()

        if verbose and ((epoch + 1) % 25 == 0 or epoch == 0):
            print(f"  epoch {epoch + 1:4d}/{epochs} | train={hist['train_total'][-1]:.4f} "
                  f"| val={hist['val_total'][-1]:.4f} | val_nll={hist['val_nll'][-1]:.4f}")

        val_nll = va['nll']
        if (best_val - val_nll) > min_delta:
            best_val = val_nll
            best_epoch = epoch + 1
            bad_epochs = 0
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        else:
            bad_epochs += 1
        if bad_epochs >= patience:
            if verbose:
                print(f'  early stopping at epoch {epoch + 1} (best val NLL @ {best_epoch}).')
            break
        if stop_at_epoch is not None and (epoch + 1) >= stop_at_epoch:
            if verbose:
                print(f'  stopping at fixed epoch {epoch + 1} (stop_at_epoch).')
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    hist['best_epoch'] = best_epoch
    return model, best_val, hist


# ---------------------------------------------------------------------------
# 2. MMD-based (two-stage) GPDFlow  (from notebook cell 32)
# ---------------------------------------------------------------------------

def fit_gpdflow_mmd(train_arr, val_arr, device, *, epochs=500, batch_size=400,
                    patience=100, warmup_epochs=50, min_delta=1e-3,
                    penalty_lambda=1e4, mmd_lambda=25.0, mmd_samples_per_face=128,
                    mmd_censor_c=0.0, mmd_censor_tau=0.15,
                    mmd_indicator_alpha=10.0, mmd_indicator_rho=5.0,
                    num_layers=4, latent_size_factor=2,
                    seed=1234, flow_seed=0, grad_clip=1.0, stop_at_epoch=None,
                    verbose=True, detect_anomaly=False):
    """Fit the MMD-based GPDFlow (``fix_margin=True``, ``use_likelihood=False``)
    on standardized-scale ``train_arr``, early-stopping on the validation total
    loss over ``val_arr`` with a warm-up (the exact rule of cell 32: the NLL term
    is always 0 here, so checkpoint selection uses ``val_stats['total_loss']``).

    Every training step passes the full ``train_arr`` tensor as ``x_data_aux``
    (the marginal term is evaluated on it), matching cell 32.

    Returns:
        ``(model, best_val, hist)``.
    """
    torch.autograd.set_detect_anomaly(bool(detect_anomaly))
    dim = train_arr.shape[1]

    flow_model = build_conditional_flow(
        dim, device, num_layers=num_layers,
        latent_size_factor=latent_size_factor, seed=flow_seed,
    )
    model = S_mGPD_NF(
        dim=dim, flow=flow_model, device=device,
        penalty_lambda=penalty_lambda, fix_margin=True, marg_lambda=0.0,
        use_likelihood=False, mmd_lambda=mmd_lambda,
        mmd_samples_per_face=mmd_samples_per_face,
        mmd_censor_c=mmd_censor_c, mmd_censor_tau=mmd_censor_tau,
        mmd_indicator_alpha=mmd_indicator_alpha, mmd_indicator_rho=mmd_indicator_rho,
    )
    optimizer = torch.optim.Adam([
        {'params': model.flow_model.parameters(), 'lr': 1e-3},
        {'params': [model.log_pi_raw], 'lr': 1e-3},
    ])
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, epochs)

    _reset_seeds(seed)

    train_data = torch.tensor(np.asarray(train_arr), dtype=torch.float32, device=device)
    val_data = torch.tensor(np.asarray(val_arr), dtype=torch.float32, device=device)
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True, drop_last=False)
    val_loader = DataLoader(val_data, batch_size=batch_size, shuffle=False, drop_last=False)

    best_val = float('inf')
    best_state = None
    best_epoch = 0
    bad_epochs = 0
    hist = {'train_total': [], 'val_total': []}

    for epoch in range(epochs):
        model.train()
        tr = _init_stats(_MMD_KEYS)
        for x_data in train_loader:
            comps = model(x_data, return_components=True, x_data_aux=train_data)
            batch_loss = comps['total_loss']
            if torch.isnan(batch_loss) or torch.isinf(batch_loss):
                if verbose:
                    print(f'  [epoch {epoch}] NaN/Inf in TRAIN loss -- stopping.')
                bad_epochs = patience
                break
            optimizer.zero_grad(set_to_none=True)
            batch_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip)
            optimizer.step()
            _update_stats(tr, comps, x_data.shape[0], _MMD_KEYS)
        tr = _finalize_stats(tr)
        hist['train_total'].append(tr['total_loss'] if tr['count'] else float('inf'))

        model.eval()
        va = _init_stats(_MMD_KEYS)
        with torch.no_grad():
            for x_val in val_loader:
                comps = model(x_val, return_components=True)
                if torch.isnan(comps['total_loss']) or torch.isinf(comps['total_loss']):
                    if verbose:
                        print(f'  [epoch {epoch}] NaN/Inf in VAL loss -- stopping.')
                    bad_epochs = patience
                    break
                _update_stats(va, comps, x_val.shape[0], _MMD_KEYS)
        va = _finalize_stats(va)
        val_epoch_loss = va['total_loss'] if va['count'] else float('inf')
        hist['val_total'].append(val_epoch_loss)

        scheduler.step()

        if verbose and ((epoch + 1) % 25 == 0 or epoch == 0):
            print(f"  epoch {epoch + 1:4d}/{epochs} | train={hist['train_total'][-1]:.4f} "
                  f"| val={val_epoch_loss:.4f}")

        past_warmup = (epoch + 1) > warmup_epochs
        if past_warmup and best_val == float('inf'):
            improved = True
        else:
            improved = past_warmup and ((best_val - val_epoch_loss) > min_delta)

        if improved:
            best_val = val_epoch_loss
            best_epoch = epoch + 1
            bad_epochs = 0
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        elif past_warmup:
            bad_epochs += 1
        if bad_epochs >= patience:
            if verbose:
                print(f'  early stopping at epoch {epoch + 1} (best val loss @ {best_epoch}).')
            break
        if stop_at_epoch is not None and (epoch + 1) >= stop_at_epoch:
            if verbose:
                print(f'  stopping at fixed epoch {epoch + 1} (stop_at_epoch).')
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    hist['best_epoch'] = best_epoch
    return model, best_val, hist


# ---------------------------------------------------------------------------
# 3. Probit-copula RealNVP baseline  (from notebook cell 42)
# ---------------------------------------------------------------------------

def fit_flow_baseline(W_train, W_val, device, *, epochs=400, batch_size=256,
                      lr=1e-3, weight_decay=1e-4, sigma_jitter=0.10,
                      num_layers=4, hidden_mult=2, warmup=20, patience=40,
                      min_delta=0.0, seed=0, grad_clip=1.0, stop_at_epoch=None,
                      verbose=True):
    """Fit the unconditional RealNVP on pseudo-Gaussian data by maximum
    likelihood, selecting / early-stopping on the held-out **NLL** after a
    warm-up.

    ``stop_at_epoch`` (if set) hard-stops after that many epochs regardless of
    the patience rule, while leaving ``epochs`` (hence the ``CosineAnnealingLR``
    ``T_max``) unchanged -- used for the final full-data refit at the median of
    the per-fold CV early-stop epochs.

    Returns:
        ``(flow, best_val, hist)`` where ``best_val`` is the best held-out NLL.
    """
    _reset_seeds(seed)
    dim = W_train.shape[1]
    flow = build_baseline_flow(dim, device, num_layers=num_layers, hidden_mult=hidden_mult)

    Wt = torch.tensor(np.asarray(W_train), dtype=torch.float32, device=device)
    Wv = torch.tensor(np.asarray(W_val), dtype=torch.float32, device=device)
    loader = DataLoader(Wt, batch_size=batch_size, shuffle=True, drop_last=False)

    opt = torch.optim.Adam(flow.parameters(), lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)

    best_val = float('inf')
    best_state = None
    best_epoch = 0
    bad = 0
    hist = {'train_nll': [], 'val_nll': []}

    for epoch in range(epochs):
        flow.train()
        running, count = 0.0, 0
        for wb in loader:
            opt.zero_grad(set_to_none=True)
            loss = flow.forward_kld(wb + sigma_jitter * torch.randn_like(wb))
            if torch.isnan(loss) or torch.isinf(loss):
                if verbose:
                    print(f'  [epoch {epoch}] NaN/Inf in train loss -- stopping.')
                bad = patience
                break
            loss.backward()
            torch.nn.utils.clip_grad_norm_(flow.parameters(), max_norm=grad_clip)
            opt.step()
            running += loss.item() * wb.shape[0]
            count += wb.shape[0]
        sched.step()
        train_nll = running / max(count, 1)

        flow.eval()
        with torch.no_grad():
            val_nll = flow.forward_kld(Wv).item()
        hist['train_nll'].append(train_nll)
        hist['val_nll'].append(val_nll)

        past_warmup = (epoch + 1) > warmup
        if past_warmup and (best_val - val_nll) > min_delta:
            best_val = val_nll
            best_epoch = epoch + 1
            bad = 0
            best_state = {k: v.detach().cpu().clone()
                          for k, v in flow.state_dict().items()}
        elif past_warmup:
            bad += 1

        if verbose and ((epoch + 1) % 25 == 0 or epoch == 0):
            print(f'  epoch {epoch + 1:4d}/{epochs} | train_nll={train_nll:.3f} '
                  f'| val_nll={val_nll:.3f}')

        if past_warmup and bad >= patience:
            if verbose:
                print(f'  early stopping at epoch {epoch + 1} (best val NLL @ {best_epoch}).')
            break
        if stop_at_epoch is not None and (epoch + 1) >= stop_at_epoch:
            if verbose:
                print(f'  stopping at fixed epoch {epoch + 1} (stop_at_epoch).')
            break

    if best_state is not None:
        flow.load_state_dict(best_state)
    flow.eval()
    hist['best_epoch'] = best_epoch
    return flow, best_val, hist


# ---------------------------------------------------------------------------
# g_std standardization helpers (two-stage MMD path, from notebook cells 28/36)
# ---------------------------------------------------------------------------

_Z_FLOOR = -30.0


def _gpd_standardize(X, sig, gam):
    """Forward g_std (CLAUDE.md Sec 3.1), applied to threshold-shifted coords.

    ``X`` is ``(n, d)`` in shifted-slice coordinates (threshold already at 0).
    Margins that fall below the fitted GPD support map to ``-inf``.
    """
    gam0 = np.isclose(gam, 0.0)
    gam_safe = np.where(gam0, 1.0, gam)
    inside = 1.0 + gam_safe * X / sig
    inside_ok = inside > 0
    Z = np.where(gam0, X / sig,
                 np.log(np.where(inside_ok, inside, 1.0)) / gam_safe)
    Z = np.where((~gam0) & (~inside_ok), -np.inf, Z)
    return Z


def _gpd_unstandardize(Z, sig, gam, thres_row):
    """Inverse g_std back to the raw observational scale (cell 36).

    ``gam * Z`` is clamped at 15 before ``exp`` (mirrors
    ``DataTransform.forward_transform``) to guard gamma<0 samples near the
    upper endpoint.
    """
    gam0 = np.isclose(gam, 0.0)
    gam_safe = np.where(gam0, 1.0, gam)
    t = np.clip(gam_safe * Z, None, 15.0)
    X = np.where(gam0, sig * Z, sig * (np.exp(t) - 1.0) / gam_safe)
    return thres_row.reshape(1, -1) + X


# ---------------------------------------------------------------------------
# CV drivers
# ---------------------------------------------------------------------------

def run_cv_gpdflow_likelihood(samples, folds, device, thres, *, epochs=500,
                              samples_multiplier=1, init_from_gpd=False,
                              verbose=True, **fit_kw):
    """5-fold CV for the likelihood-based GPDFlow.

    For each fold: train on the other 4 folds, early-stop on the held-out fold,
    then draw ``samples_multiplier * |fold|`` observations from the fold model.
    The out-of-fold simulated blocks are kept **separate** (one per fold) -- they
    are NOT concatenated, because a pool of draws from 5 separately fitted flows
    is a 5-component mixture that is not itself an mGPD.

    ``init_from_gpd`` (default ``False``): if ``True``, warm-start each fold's
    GPDFlow marginal parameters from a per-margin univariate GPD fit of that
    fold's training exceedances (forwarded to ``fit_gpdflow_likelihood``);
    ``False`` keeps the cold start at ``sigma = 1``, ``gamma = 0``.

    No tail-dependence summary is computed here: keep the per-fold draws and
    compute ``chi`` / ``omega`` in the notebook, so the threshold can be changed
    without refitting.

    Returns a dict with keys: ``models``, ``best_vals``, ``hist``,
    ``sigma_hat_folds`` / ``gamma_hat_folds`` (n_folds x dim),
    ``sigma_hat`` / ``gamma_hat`` (fold means),
    ``samples_simu_folds`` (list of n_folds raw-scale ``(n_k, dim)`` arrays).
    """
    samples = np.asarray(samples, dtype=np.float32)
    thres = np.asarray(thres, dtype=np.float64).reshape(1, -1)
    dim = samples.shape[1]

    models, best_vals, hists = [], [], []
    sigma_folds, gamma_folds, sim_blocks = [], [], []

    for k, (tr_idx, va_idx) in enumerate(folds):
        if verbose:
            print(f'[likelihood] fold {k + 1}/{len(folds)} '
                  f'(train {len(tr_idx)}, val {len(va_idx)})')
        model, bv, h = fit_gpdflow_likelihood(
            samples[tr_idx], samples[va_idx], device,
            epochs=epochs, init_from_gpd=init_from_gpd, verbose=verbose, **fit_kw,
        )
        models.append(model)
        best_vals.append(bv)
        hists.append(h)
        sigma_folds.append(model.get_sigma().detach().cpu().numpy())
        gamma_folds.append(model.get_gamma().detach().cpu().numpy())

        n_k = int(samples_multiplier * len(va_idx))
        with torch.no_grad():
            xs, _, _ = model.sample(n_k)
        sim_blocks.append(xs.detach().cpu().numpy() + thres)   # raw scale, per fold

    sigma_hat_folds = np.vstack(sigma_folds)
    gamma_hat_folds = np.vstack(gamma_folds)

    return {
        'models': models, 'best_vals': best_vals, 'hist': hists,
        'sigma_hat_folds': sigma_hat_folds, 'gamma_hat_folds': gamma_hat_folds,
        'sigma_hat': sigma_hat_folds.mean(axis=0),
        'gamma_hat': gamma_hat_folds.mean(axis=0),
        'samples_simu_folds': sim_blocks,
    }


def run_cv_gpdflow_mmd(samples, folds, device, thres, *, epochs=500,
                       samples_multiplier=1,
                       gamma_clip=0.5, min_tail_count=10, verbose=True, **mmd_kw):
    """5-fold CV for the two-stage MMD-based GPDFlow.

    Per fold the Stage-1 univariate GPD is **refit on the fold's training rows
    only** (no leakage), then used to standardize both the train and val rows,
    re-threshold at 0, and finally map the fold model's standardized draws back
    to the raw scale.  The raw out-of-fold blocks are kept **separate** (one per
    fold), not concatenated.  No tail-dependence summary is computed here --
    compute ``chi`` / ``omega`` in the notebook from ``samples_simu_2s_folds``.

    Returns a dict with keys: ``models``, ``best_vals``, ``hist``,
    ``dropped_rows_train`` / ``dropped_rows_val`` (per fold),
    ``uni_sigma_2s_folds`` / ``uni_gamma_2s_folds`` (n_folds x dim),
    ``samples_simu_2s_folds`` (list of n_folds raw-scale arrays), ``dim_2s``.
    """
    samples = np.asarray(samples, dtype=np.float64)
    thres = np.asarray(thres, dtype=np.float64).reshape(1, -1)
    dim = samples.shape[1]

    models, best_vals, hists = [], [], []
    sig_folds, gam_folds = [], []
    dropped_tr, dropped_va, raw_blocks = [], [], []

    for k, (tr_idx, va_idx) in enumerate(folds):
        tr, va = samples[tr_idx], samples[va_idx]

        # Stage 1: per-margin GPD refit on fold-train positive excesses only.
        sig_k = np.ones(dim)
        gam_k = np.zeros(dim)
        for j in range(dim):
            col = tr[tr[:, j] > 0, j]
            if len(col) > min_tail_count:
                c, _loc, sc = stats.genpareto.fit(col, floc=0)
                gam_k[j] = float(np.clip(c, -gamma_clip, gamma_clip))
                sig_k[j] = float(max(sc, 1e-3))
        sig_folds.append(sig_k.copy())
        gam_folds.append(gam_k.copy())

        std_tr = _gpd_standardize(tr, sig_k, gam_k)
        std_va = _gpd_standardize(va, sig_k, gam_k)

        keep_tr = np.any(std_tr > 0, axis=1)
        keep_va = np.any(std_va > 0, axis=1)
        dropped_tr.append(int((~keep_tr).sum()))
        dropped_va.append(int((~keep_va).sum()))

        train_2s = np.clip(std_tr[keep_tr], _Z_FLOOR, None).astype(np.float32)
        val_2s = np.clip(std_va[keep_va], _Z_FLOOR, None).astype(np.float32)

        assert train_2s.shape[1] == dim, (train_2s.shape, dim)

        if verbose:
            print(f'[mmd] fold {k + 1}/{len(folds)} '
                  f'(train {train_2s.shape[0]}/{len(tr_idx)}, '
                  f'val {val_2s.shape[0]}/{len(va_idx)})')

        model, bv, h = fit_gpdflow_mmd(
            train_2s, val_2s, device, epochs=epochs, verbose=verbose, **mmd_kw,
        )
        models.append(model)
        best_vals.append(bv)
        hists.append(h)

        n_k = int(samples_multiplier * len(va_idx))
        with torch.no_grad():
            _, ys, _ = model.sample(n_k)          # standardized-scale draws (fix_margin)
        z_k = ys.detach().cpu().numpy()
        raw_blocks.append(_gpd_unstandardize(z_k, sig_k, gam_k, thres[0]))

    return {
        'models': models, 'best_vals': best_vals, 'hist': hists,
        'dropped_rows_train': dropped_tr, 'dropped_rows_val': dropped_va,
        'uni_sigma_2s_folds': np.vstack(sig_folds),
        'uni_gamma_2s_folds': np.vstack(gam_folds),
        'samples_simu_2s_folds': raw_blocks,
        'dim_2s': dim,
    }


def run_cv_flow_baseline(samples, folds, device, thres, *, epochs=400,
                         samples_multiplier=1, max_draw_rounds=20,
                         verbose=True, **bl_kw):
    """5-fold CV for the probit-copula RealNVP baseline.

    Per fold the hybrid marginal model is refit on the fold's training rows, the
    flow is trained on the pseudo-Gaussian train data and early-stopped on the
    held-out NLL, then ``samples_multiplier * |fold|`` finite draws are mapped
    back to the raw scale.  The 5 out-of-fold blocks are kept **separate** (one
    per fold), NOT concatenated -- as with the GPDFlow drivers, a pool of draws
    from 5 separately fitted flows is a mixture, not the target distribution.
    No tail-dependence summary is computed here -- form ``chi`` / ``omega`` in
    the notebook as the fold-mean of the per-fold estimates.

    Returns a dict with keys: ``models``, ``best_vals``, ``hist``,
    ``samples_simu_bl_folds`` (list of n_folds raw-scale ``(n_k, dim)`` arrays),
    ``last_fold`` (artifacts for the notebook's diagnostic cells).
    """
    samples = np.asarray(samples, dtype=np.float64)
    thres = np.asarray(thres, dtype=np.float64).reshape(1, -1)
    dim = samples.shape[1]

    models, best_vals, hists, sim_blocks = [], [], [], []
    last_fold = None

    for k, (tr_idx, va_idx) in enumerate(folds):
        to_g, from_g, _cdf, _ppf = build_probit_copula(samples[tr_idx], dim)
        W_train_k = to_g(samples[tr_idx])
        W_val_k = to_g(samples[va_idx])

        if verbose:
            print(f'[baseline] fold {k + 1}/{len(folds)} '
                  f'(train {len(tr_idx)}, val {len(va_idx)})')

        flow, bv, h = fit_flow_baseline(
            W_train_k, W_val_k, device, epochs=epochs, verbose=verbose, **bl_kw,
        )
        models.append(flow)
        best_vals.append(bv)
        hists.append(h)

        n_k = int(samples_multiplier * len(va_idx))
        collected = []
        have = 0
        flow.eval()
        for _ in range(max_draw_rounds):
            with torch.no_grad():
                Ws, _ = flow.sample(max(2 * (n_k - have), 256))
            Ws = Ws.detach().cpu().numpy()
            Ws = Ws[np.isfinite(Ws).all(axis=1)]
            collected.append(Ws)
            have += Ws.shape[0]
            if have >= n_k:
                break
        W_star_k = np.vstack(collected)[:n_k]
        sim_blocks.append(from_g(W_star_k) + thres)   # raw scale, per fold
        last_fold = {'W_train': W_train_k, 'W_val': W_val_k,
                     'W_star': W_star_k, 'from_gaussian': from_g}

    return {
        'models': models, 'best_vals': best_vals, 'hist': hists,
        'samples_simu_bl_folds': sim_blocks,
        'last_fold': last_fold,
    }


def run_cv_parametric_mgpd(samples, folds, device, thres, generator, *,
                           epochs=300, samples_multiplier=1, seed=1234,
                           patience=40, verbose=True):
    """5-fold CV for the traditional parametric mGPD, for one fixed generator.

    ``generator`` is the AIC-selected family (notebook cell 50).  Per fold: train
    on the other 4 folds with the margins learned jointly (``fix_margin=False``,
    like the likelihood driver -- no per-fold Stage-1 refit), early-stop on the
    held-out fold, then draw ``samples_multiplier * |fold|`` observations.  The 5
    out-of-fold blocks are kept **separate**; form ``chi`` / ``omega`` in the
    notebook as the fold-mean of the per-fold estimates.

    Returns a dict with keys: ``models``, ``best_vals``, ``hist``,
    ``samples_simu_par_folds`` (list of n_folds raw-scale ``(n_k, dim)`` arrays).
    """
    samples = np.asarray(samples, dtype=np.float64)
    thres = np.asarray(thres, dtype=np.float64).reshape(1, -1)

    models, best_vals, hists, sim_blocks = [], [], [], []

    for k, (tr_idx, va_idx) in enumerate(folds):
        if verbose:
            print(f'[parametric:{generator}] fold {k + 1}/{len(folds)} '
                  f'(train {len(tr_idx)}, val {len(va_idx)})')
        model, bv, h = fit_parametric_mgpd(
            generator, samples[tr_idx], samples[va_idx], device,
            epochs=epochs, seed=seed, patience=patience, verbose=verbose,
        )
        models.append(model)
        best_vals.append(bv)
        hists.append(h)

        n_k = int(samples_multiplier * len(va_idx))
        model.eval()
        with torch.no_grad():
            xs, _, _ = model.sample(n_k)
        xs = xs.detach().cpu().numpy()
        xs = xs[np.isfinite(xs).all(axis=1)]
        sim_blocks.append(xs + thres)   # raw scale, per fold

    return {
        'models': models, 'best_vals': best_vals, 'hist': hists,
        'samples_simu_par_folds': sim_blocks,
    }


# ---------------------------------------------------------------------------
# Monte Carlo simulation study (Simulation_1S.ipynb)
# ---------------------------------------------------------------------------

# Default 25-dim marginal-parameter vectors of the Simulation 1 DGP
# (Simulation_1S.ipynb cell 5).  Sliced to ``:dim`` inside ``simulate_1_mgpd``.
_SIM1_SIGMA_FULL = np.array([0.5, 1.2, 1.0, 1.5, 0.8] * 10)   # shape (50,)
_SIM1_GAMMA_FULL = np.array([-0.1, 0.2, 0.0, 0.1, 0.3] * 10)  # shape (50,)

# Default Burr XII marginal parameters of the Simulation 2 DGP
# (Simulation_2S.ipynb cell "In[21]").  Tiled to ``:dim`` inside
# ``simulate_2_tcopula_burr``.  ``k`` is scipy's ``burr12`` shape ``d``.
_SIM2_C_BASE = np.array([2.0, 3.0, 2.5, 1.5, 3.0])
_SIM2_K_BASE = np.array([1.5, 1.0, 2.0, 3.0, 3.5])


def simulate_1_mgpd(*, n_samples=1000, dim=25, alpha=1.0, beta=0.0,
                    sigma_full=None, gamma_full=None, seed=42):
    """One S-representation mGPD sample of the Simulation 1 DGP.

    Independent reverse-exponential generator ``T_j = log(U_j) / alpha - beta``
    with ``U_j ~ Unif(0, 1)`` i.i.d., standardised mGPD ``Z = E + T - max(T)``
    (``E ~ Exp(1)``), then the marginal forward transform

        X_ij = sigma_j * Z_ij                               if gamma_j == 0
        X_ij = (sigma_j / gamma_j) * (exp(gamma_j Z_ij) - 1) otherwise.

    This is the same construction as ``Simulation_1S.ipynb`` cell 5 and as
    ``cf.sim_revexp_T_mgpd`` (whose scale is ``1 / alpha``).  It draws every
    random number from a single ``np.random.default_rng(seed)`` stream instead of
    the cell's ``np.random.seed`` + ``torch.manual_seed`` +
    ``torch.Tensor.exponential_`` mix, so it does **not** bit-reproduce the
    notebook's ``samples_origin`` -- the DGP and its theoretical
    ``chi = 1 - 1/(1 + 2 alpha)`` are identical.  Cell 5 is left unchanged.

    Args:
        n_samples: number of observations to draw.
        dim: dimension (``sigma_full`` / ``gamma_full`` are sliced to ``:dim``).
        alpha, beta: reverse-exponential generator parameters (equal across
            margins).
        sigma_full, gamma_full: marginal scale / shape vectors, length >= dim.
            Default to the Simulation 1 25-dim vectors.
        seed: RNG seed for this replicate.

    Returns:
        ``X`` -- raw-scale array, shape ``(n_samples, dim)``.
    """
    if sigma_full is None:
        sigma_full = _SIM1_SIGMA_FULL
    if gamma_full is None:
        gamma_full = _SIM1_GAMMA_FULL
    sigma = np.asarray(sigma_full, dtype=np.float64)[:dim]
    gamma = np.asarray(gamma_full, dtype=np.float64)[:dim]

    rng = np.random.default_rng(seed)
    U = rng.uniform(0.0, 1.0, size=(n_samples, dim))
    T = np.log(U) / alpha - beta
    S = T - T.max(axis=1, keepdims=True)          # S = T - max(T)
    E = rng.exponential(1.0, size=n_samples)[:, None]
    Z = E + S                                     # standardised mGPD

    gam0 = (gamma == 0.0)
    gam_safe = np.where(gam0, 1.0, gamma)
    X = np.where(gam0, sigma * Z,
                 (sigma / gam_safe) * (np.exp(gam_safe * Z) - 1.0))
    return X


def cv_pairwise_chi(fold_blocks, p):
    """Fold-mean rank-based pairwise-chi matrix from a list of draw blocks.

    ``fold_blocks`` is a list of raw-scale ``(n_k, dim)`` arrays (one per CV
    fold, e.g. ``run_cv_gpdflow_likelihood(...)['samples_simu_folds']``).  Returns
    the average over folds of ``cf.pairwise_chi_from_data(block, p=p)`` -- the
    folds are kept separate because a pool of draws from 5 separately fitted
    flows is a mixture, not an mGPD.
    """
    return np.mean([cf.pairwise_chi_from_data(b, p=p) for b in fold_blocks],
                   axis=0)


_SIM1_MODEL_DEFAULTS = dict(
    epochs=500, samples_multiplier=30, patience=100,
    num_layers=4, latent_size_factor=3,
    penalty_lambda=1e4, marg_lambda=200, seed=1234,
)


def run_simulation_study_1(n_sim=20, device=None, *, base_seed=42,
                           sim_kwargs=None, model_kwargs=None,
                           init_from_gpd=False,
                           n_splits=5, fold_seed=1234, chi_q=0.95, verbose=True):
    """Repeat the Simulation 1 DGP + 5-fold-CV GPDFlow fit ``n_sim`` times.

    Each replicate: draw a fresh sample from ``simulate_1_mgpd`` (per-replicate
    seed ``base_seed + r``), partition it into ``n_splits`` folds, run
    ``run_cv_gpdflow_likelihood`` (the likelihood / "S representation" path), and
    keep the **aggregated CV estimators**:

      * ``sigma_hat`` / ``gamma_hat`` -- the fold-mean marginal parameters
        (``cv['sigma_hat']`` / ``cv['gamma_hat']``);
      * ``chi_hat`` -- the fold-mean rank-based pairwise-chi matrix at level
        ``chi_q`` (``cv_pairwise_chi(cv['samples_simu_folds'], chi_q)``),
        stored as its upper-triangle vector.

    Defaults reproduce the current ``Simulation_1S.ipynb`` setting.  Pass a
    partial ``sim_kwargs`` / ``model_kwargs`` dict to override individual knobs
    (missing keys fall back to the defaults).

    Args:
        n_sim: number of simulation replicates.
        device: torch device (default: cuda if available else cpu).
        base_seed: replicate ``r`` uses DGP seed ``base_seed + r``.
        sim_kwargs: overrides for ``simulate_1_mgpd`` (``n_samples``, ``dim``,
            ``alpha``, ``beta``, ``sigma_full``, ``gamma_full``).  ``seed`` is
            ignored -- it is set per replicate.
        model_kwargs: overrides for the CV fit -- ``epochs``,
            ``samples_multiplier``, ``patience`` go to
            ``run_cv_gpdflow_likelihood``; the rest (``num_layers``,
            ``latent_size_factor``, ``penalty_lambda``, ``marg_lambda``,
            ``seed``, ...) are forwarded to ``fit_gpdflow_likelihood``.
        init_from_gpd: if True, warm-start each fold's GPDFlow marginal
            parameters from a per-margin univariate GPD fit
            (``scipy.stats.genpareto.fit(col, floc=0)`` on that fold's training
            exceedances), instead of the default cold start at
            ``sigma = 1``, ``gamma = 0``.  Forwarded via
            ``run_cv_gpdflow_likelihood(init_from_gpd=...)``.  May also be passed
            inside ``model_kwargs``; the explicit argument is the fallback.
        n_splits, fold_seed: passed to ``make_folds``.
        chi_q: exceedance level for the pairwise-chi estimate.
        verbose: forwarded to the CV driver.

    Returns:
        dict with keys ``sigma_hat`` / ``gamma_hat`` ``(n_sim, dim)``,
        ``chi_hat_pairs`` ``(n_sim, dim*(dim-1)/2)``,
        ``chi_hat_offdiag_mean`` ``(n_sim,)``,
        ``sigma_true`` / ``gamma_true`` ``(dim,)``, ``chi_true`` (float scalar),
        ``pair_idx`` (``np.triu_indices(dim, 1)``),
        ``best_epochs`` / ``best_vals`` ``(n_sim, n_splits)``,
        ``rep_seeds`` (list), ``config`` (dict).
    """
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    sim_cfg = dict(n_samples=1000, dim=25, alpha=1.0, beta=0.0)
    if sim_kwargs:
        sim_cfg.update(sim_kwargs)
    sim_cfg.pop('seed', None)   # set per replicate

    model_cfg = dict(_SIM1_MODEL_DEFAULTS)
    if model_kwargs:
        model_cfg.update(model_kwargs)
    epochs = model_cfg.pop('epochs')
    samples_multiplier = model_cfg.pop('samples_multiplier')
    # Marginal warm-start toggle: honour it if passed inside model_kwargs too,
    # explicit arg as fallback; pop it so it is not forwarded twice via **fit_kw.
    init_from_gpd = bool(model_cfg.pop('init_from_gpd', init_from_gpd))
    fit_kw = model_cfg   # patience + architecture/penalty knobs -> fit_gpdflow_likelihood

    dim = sim_cfg['dim']
    alpha = sim_cfg['alpha']
    sigma_full = sim_cfg.get('sigma_full')
    gamma_full = sim_cfg.get('gamma_full')
    sigma_true = (np.asarray(sigma_full, dtype=np.float64)[:dim]
                  if sigma_full is not None else _SIM1_SIGMA_FULL[:dim])
    gamma_true = (np.asarray(gamma_full, dtype=np.float64)[:dim]
                  if gamma_full is not None else _SIM1_GAMMA_FULL[:dim])
    chi_true = float(cf.chi_theorical(alpha, alpha))
    iu = np.triu_indices(dim, k=1)
    n_pairs = iu[0].size

    sigma_hat = np.empty((n_sim, dim))
    gamma_hat = np.empty((n_sim, dim))
    chi_hat_pairs = np.empty((n_sim, n_pairs))
    best_epochs = np.empty((n_sim, n_splits), dtype=int)
    best_vals = np.empty((n_sim, n_splits))
    rep_seeds = []
    thres = np.zeros((1, dim))

    for r in range(n_sim):
        rep_seed = base_seed + r
        rep_seeds.append(rep_seed)
        if verbose:
            print(f'\n===== simulation replicate {r + 1}/{n_sim} '
                  f'(DGP seed {rep_seed}) =====')

        samples = simulate_1_mgpd(seed=rep_seed, **sim_cfg)
        folds = make_folds(len(samples), n_splits=n_splits, seed=fold_seed)

        cv = run_cv_gpdflow_likelihood(
            samples, folds, device, thres,
            epochs=epochs, samples_multiplier=samples_multiplier,
            init_from_gpd=init_from_gpd,
            verbose=verbose, **fit_kw,
        )

        sigma_hat[r] = cv['sigma_hat']
        gamma_hat[r] = cv['gamma_hat']
        chi_mat = cv_pairwise_chi(cv['samples_simu_folds'], chi_q)
        chi_hat_pairs[r] = chi_mat[iu]
        best_epochs[r] = [h['best_epoch'] for h in cv['hist']]
        best_vals[r] = cv['best_vals']

    return {
        'sigma_hat': sigma_hat,
        'gamma_hat': gamma_hat,
        'chi_hat_pairs': chi_hat_pairs,
        'chi_hat_offdiag_mean': chi_hat_pairs.mean(axis=1),
        'sigma_true': sigma_true,
        'gamma_true': gamma_true,
        'chi_true': chi_true,
        'pair_idx': iu,
        'best_epochs': best_epochs,
        'best_vals': best_vals,
        'rep_seeds': rep_seeds,
        'config': {
            'n_sim': n_sim, 'base_seed': base_seed,
            'sim_kwargs': sim_cfg, 'model_kwargs': dict(_SIM1_MODEL_DEFAULTS,
                                                        **(model_kwargs or {})),
            'init_from_gpd': init_from_gpd,
            'n_splits': n_splits, 'fold_seed': fold_seed, 'chi_q': chi_q,
        },
    }


def summarize_simulation_study_1(results):
    """Monte Carlo bias and variance of the ``run_simulation_study_1`` estimators.

    For each estimator (``sigma`` / ``gamma`` per margin, ``chi`` per pair) with
    ``R`` replicate estimates ``{theta_hat_r}`` and truth ``theta``:

        bias         = mean_r(theta_hat_r) - theta
        variance     = var_r(theta_hat_r, ddof=1)          # unbiased sample var
        variance_pop  = var_r(theta_hat_r, ddof=0)          # = (R-1)/R * variance
        mse          = mean_r((theta_hat_r - theta) ** 2)   # empirical MSE
        rmse         = sqrt(mse)

    The bias-variance decomposition ``mse == bias**2 + variance_pop`` is an exact
    algebraic identity (expand ``(theta_hat_r - theta_bar + theta_bar - theta)**2``
    and sum).  With the *unbiased* ``variance`` it holds only up to the finite-R
    term ``bias**2 + variance - mse == variance / R``.  ``summary`` carries
    ``*_decomp_err`` = ``bias**2 + variance_pop - mse`` (numerically ~0) so the
    identity can be checked directly.

    Args:
        results: the dict returned by ``run_simulation_study_1``.

    Returns:
        ``(summary_df, summary)``.  ``summary`` holds the per-margin / per-pair
        arrays ``{name}_mean_est``, ``{name}_bias``, ``{name}_variance``,
        ``{name}_variance_pop``, ``{name}_mse``, ``{name}_rmse``,
        ``{name}_decomp_err`` for ``name`` in ``sigma`` / ``gamma`` / ``chi``.
        ``summary_df`` is a 3-row pandas DataFrame with columns ``mean_est``,
        ``mean_bias``, ``mean_abs_bias``, ``mean_variance``, ``mean_mse``,
        ``bias2_plus_var``, ``rmse`` (each averaged over margins / pairs;
        ``bias2_plus_var`` uses ``variance_pop`` so it equals ``mean_mse``).
    """
    R = np.asarray(results['sigma_hat']).shape[0]
    summary = {}
    rows = {}
    specs = [
        ('sigma', results['sigma_hat'], results['sigma_true']),
        ('gamma', results['gamma_hat'], results['gamma_true']),
        ('chi', results['chi_hat_pairs'], results['chi_true']),
    ]
    for name, hat, true in specs:
        hat = np.asarray(hat, dtype=np.float64)
        true = np.asarray(true, dtype=np.float64)
        mean_est = hat.mean(axis=0)
        bias = mean_est - true
        variance = hat.var(axis=0, ddof=1)
        variance_pop = hat.var(axis=0, ddof=0)
        mse = np.mean((hat - true) ** 2, axis=0)
        decomp_err = bias ** 2 + variance_pop - mse          # ~0 exactly
        rmse = np.sqrt(mse)

        summary[f'{name}_mean_est'] = mean_est
        summary[f'{name}_bias'] = bias
        summary[f'{name}_variance'] = variance
        summary[f'{name}_variance_pop'] = variance_pop
        summary[f'{name}_mse'] = mse
        summary[f'{name}_rmse'] = rmse
        summary[f'{name}_decomp_err'] = decomp_err

        rows[name] = {
            'mean_est': float(np.mean(mean_est)),
            'mean_bias': float(np.mean(bias)),
            'mean_abs_bias': float(np.mean(np.abs(bias))),
            'mean_variance': float(np.mean(variance)),
            'mean_mse': float(np.mean(mse)),
            'bias2_plus_var': float(np.mean(bias ** 2 + variance_pop)),
            'rmse': float(np.sqrt(np.mean(mse))),
        }

    summary['n_rep'] = R
    summary_df = pd.DataFrame.from_dict(rows, orient='index')[
        ['mean_est', 'mean_bias', 'mean_abs_bias', 'mean_variance',
         'mean_mse', 'bias2_plus_var', 'rmse']
    ]
    return summary_df, summary


# ---------------------------------------------------------------------------
# Monte Carlo simulation study 2 (Simulation_2S.ipynb):
# AR(1) Student-t copula + Burr XII margins, POT threshold exceedances.
# ---------------------------------------------------------------------------

_SIM2_MODEL_DEFAULTS = dict(_SIM1_MODEL_DEFAULTS)   # same flow / penalty defaults


def simulate_2_tcopula_burr(*, n_samples=30000, dim=20, rho=0.7, nu=5.0,
                            corr=None, c_full=None, k_full=None, seed=42):
    """One realisation of the Simulation 2 DGP: t-copula + Burr XII margins.

    Same construction as ``Simulation_2S.ipynb`` cell "In[21]", factored out.
    Draw ``Y ~ multivariate_t(0, R, nu)`` with ``R`` an AR(1) correlation matrix
    (``R[j, k] = rho ** |j - k|``, or the explicit ``corr``), map to uniforms with
    the t CDF, then apply Burr XII inverse CDFs
    ``F_j^{-1}(u) = burr12.ppf(u; c_j, d=k_j)``.

    Args:
        n_samples: number of observations to draw.
        dim: dimension.  ``c_full`` / ``k_full`` are tiled / sliced to ``:dim``.
        rho: AR(1) base correlation (ignored if ``corr`` is given).
        nu: degrees of freedom of the t-copula.
        corr: optional explicit ``(dim, dim)`` correlation matrix.
        c_full, k_full: Burr XII ``c`` / ``k`` base vectors (any length >= 1;
            tiled to ``dim``).  Default to the Simulation 2 5-vectors.
        seed: RNG seed for this replicate.

    Returns:
        ``X`` -- raw-scale array, shape ``(n_samples, dim)``.
    """
    c_base = _SIM2_C_BASE if c_full is None else np.asarray(c_full, dtype=np.float64)
    k_base = _SIM2_K_BASE if k_full is None else np.asarray(k_full, dtype=np.float64)
    c = np.tile(c_base, dim // len(c_base) + 1)[:dim]
    k = np.tile(k_base, dim // len(k_base) + 1)[:dim]

    R = cf.ar1_correlation(dim, rho) if corr is None else np.asarray(corr, dtype=np.float64)

    Y = stats.multivariate_t(loc=np.zeros(dim), shape=R, df=nu).rvs(
        size=n_samples, random_state=seed)
    Y = np.asarray(Y, dtype=np.float64).reshape(n_samples, dim)
    U = stats.t.cdf(Y, df=nu)
    X = np.column_stack([stats.burr12.ppf(U[:, j], c=c[j], d=k[j])
                         for j in range(dim)])
    return X


def burr12_excess_gpd_params(c, k, q):
    """Exact conditional-excess GPD parameters of a Burr XII margin at level ``q``.

    Reciprocal-hazard / von Mises parameterisation: for ``a(x) = S(x) / f(x)``
    (the reciprocal hazard function), the excess over a threshold ``u`` is
    generalized Pareto with scale ``sigma(u) = a(u)`` and shape
    ``gamma(u) = a'(u)`` -- no GPD fit / estimation involved.  For Burr XII
    (``S(x) = (1 + x^c)^{-k}``) this has the closed form

        u      = burr12.ppf(q; c, d=k)
        sigma  = (u ** (1 - c) + u) / (k * c)
        gamma  = ((1 - c) * u ** (-c) + 1) / (k * c)

    with ``gamma(u) -> 1 / (c * k)`` (the Burr XII tail index reciprocal) as
    ``q -> 1``.

    Args:
        c, k: Burr XII shape parameters (scalars or equal-length arrays; ``k`` is
            scipy's ``burr12`` shape ``d``).
        q: threshold probability level in ``(0, 1)``.

    Returns:
        ``(u, sigma, gamma)`` -- arrays broadcast to the shape of ``c`` / ``k``.
    """
    c = np.asarray(c, dtype=np.float64)
    k = np.asarray(k, dtype=np.float64)
    u = stats.burr12.ppf(q, c=c, d=k)
    sigma = (u ** (1.0 - c) + u) / (k * c)
    gamma = ((1.0 - c) * u ** (-c) + 1.0) / (k * c)
    return u, sigma, gamma


def run_simulation_study_2(n_sim=20, device=None, *, base_seed=42,
                           threshold=0.97, sim_kwargs=None, model_kwargs=None,
                           init_from_gpd=False,
                           n_splits=5, fold_seed=1234, chi_q=0.95, verbose=True):
    """Repeat the Simulation 2 DGP + threshold-exceedance 5-fold-CV fit ``n_sim`` times.

    Each replicate: draw a fresh sample from ``simulate_2_tcopula_burr``
    (per-replicate seed ``base_seed + r``); form the threshold-exceedance set --
    every row with at least one component above its per-margin empirical
    ``threshold``-quantile -- and subtract the threshold; partition the excess
    data into ``n_splits`` folds; run ``run_cv_gpdflow_likelihood``; and keep the
    aggregated CV estimators (fold-mean ``sigma_hat`` / ``gamma_hat`` and the
    fold-mean rank-based pairwise chi).

    The ``threshold`` argument is a scalar **quantile level** ``q``: per replicate
    the exceedance cutoff is ``np.quantile(samples_origin, q, axis=0)`` (matches
    ``Simulation_2S.ipynb`` cell "In[24]").

    RMSE targets (the ``*_true`` entries) are analytic, computed once:

      * ``chi_true`` -- per pair, the t-copula closed form
        ``cf.t_copula_pairwise_chi(R, nu)`` (the ``q -> 1`` limit);
      * ``sigma_true`` / ``gamma_true`` -- the exact conditional-excess GPD
        parameters of each Burr XII margin at level ``q`` via
        ``burr12_excess_gpd_params`` (reciprocal hazard, no fitting).

    Because the CV draws are exceedance-conditional, the pairwise chi is
    evaluated at the cluster-conditional level
    ``p_sim = 1 - (1 - chi_q) / p_cluster`` with
    ``p_cluster = n_exceed / n_origin`` (matches ``Simulation_2S.ipynb`` "In[35]").

    Args:
        n_sim: number of simulation replicates.
        device: torch device (default: cuda if available else cpu).
        base_seed: replicate ``r`` uses DGP seed ``base_seed + r``.
        threshold: scalar quantile level ``q`` for the POT exceedance cutoff and
            for the analytic marginal truth.
        sim_kwargs: overrides for ``simulate_2_tcopula_burr`` (``n_samples``,
            ``dim``, ``rho``, ``nu``, ``corr``, ``c_full``, ``k_full``).  ``seed``
            is ignored -- it is set per replicate.
        model_kwargs: overrides for the CV fit -- ``epochs``,
            ``samples_multiplier``, ``patience`` go to
            ``run_cv_gpdflow_likelihood``; the rest (``num_layers``,
            ``latent_size_factor``, ``penalty_lambda``, ``marg_lambda``,
            ``seed``, ...) are forwarded to ``fit_gpdflow_likelihood``.
        init_from_gpd: warm-start each fold's GPDFlow margins from a per-fold
            univariate GPD fit of the exceedance data (see
            ``run_simulation_study_1``).  May also be passed inside
            ``model_kwargs``.
        n_splits, fold_seed: passed to ``make_folds``.
        chi_q: target (unconditional) exceedance level for the pairwise-chi
            estimate; re-levelled per replicate via ``p_cluster``.
        verbose: forwarded to the CV driver.

    Returns:
        dict with the same keys as ``run_simulation_study_1`` --
        ``sigma_hat`` / ``gamma_hat`` ``(n_sim, dim)``,
        ``chi_hat_pairs`` ``(n_sim, n_pairs)``, ``chi_hat_offdiag_mean``
        ``(n_sim,)``, ``sigma_true`` / ``gamma_true`` ``(dim,)``,
        ``chi_true`` ``(n_pairs,)`` **vector** (per pair, not a scalar),
        ``pair_idx``, ``best_epochs`` / ``best_vals`` ``(n_sim, n_splits)``,
        ``rep_seeds``, ``config`` -- plus Simulation-2 extras ``u_true``
        ``(dim,)``, ``chi_true_mat`` ``(dim, dim)``, ``p_cluster`` / ``p_sim``
        ``(n_sim,)``, ``threshold``.
    """
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    sim_cfg = dict(n_samples=30000, dim=20, rho=0.7, nu=5.0)
    if sim_kwargs:
        sim_cfg.update(sim_kwargs)
    sim_cfg.pop('seed', None)   # set per replicate

    model_cfg = dict(_SIM2_MODEL_DEFAULTS)
    if model_kwargs:
        model_cfg.update(model_kwargs)
    epochs = model_cfg.pop('epochs')
    samples_multiplier = model_cfg.pop('samples_multiplier')
    init_from_gpd = bool(model_cfg.pop('init_from_gpd', init_from_gpd))
    fit_kw = model_cfg

    dim = sim_cfg['dim']
    rho = sim_cfg['rho']
    nu = sim_cfg['nu']
    corr = sim_cfg.get('corr')
    R = cf.ar1_correlation(dim, rho) if corr is None else np.asarray(corr, dtype=np.float64)

    c_base = np.asarray(sim_cfg.get('c_full', _SIM2_C_BASE), dtype=np.float64)
    k_base = np.asarray(sim_cfg.get('k_full', _SIM2_K_BASE), dtype=np.float64)
    c_vec = np.tile(c_base, dim // len(c_base) + 1)[:dim]
    k_vec = np.tile(k_base, dim // len(k_base) + 1)[:dim]

    iu = np.triu_indices(dim, k=1)
    n_pairs = iu[0].size

    chi_true_mat = np.asarray(cf.t_copula_pairwise_chi(R, nu), dtype=np.float64)
    chi_true = chi_true_mat[iu]
    u_true, sigma_true, gamma_true = burr12_excess_gpd_params(c_vec, k_vec, threshold)

    sigma_hat = np.empty((n_sim, dim))
    gamma_hat = np.empty((n_sim, dim))
    chi_hat_pairs = np.empty((n_sim, n_pairs))
    best_epochs = np.empty((n_sim, n_splits), dtype=int)
    best_vals = np.empty((n_sim, n_splits))
    p_cluster_arr = np.empty(n_sim)
    p_sim_arr = np.empty(n_sim)
    rep_seeds = []

    for r in range(n_sim):
        rep_seed = base_seed + r
        rep_seeds.append(rep_seed)
        if verbose:
            print(f'\n===== simulation replicate {r + 1}/{n_sim} '
                  f'(DGP seed {rep_seed}) =====')

        samples_origin = simulate_2_tcopula_burr(seed=rep_seed, **sim_cfg)
        thres_vec = np.quantile(samples_origin, threshold, axis=0)
        cond = (samples_origin > thres_vec).any(axis=1)
        samples = (samples_origin[cond] - thres_vec).astype(np.float32)

        p_cluster = samples.shape[0] / samples_origin.shape[0]
        p_sim_raw = 1.0 - (1.0 - chi_q) / p_cluster
        p_sim = float(np.clip(p_sim_raw, 0.5, 0.999))
        if verbose or p_sim != p_sim_raw:
            msg = (f'[replicate {r + 1}] p_cluster={p_cluster:.4f}  '
                   f'p_sim={p_sim:.4f}')
            if p_sim != p_sim_raw:
                msg += (f'  (clipped from {p_sim_raw:.4f}; p_cluster too small '
                        f'for chi_q={chi_q})')
            print(msg)

        folds = make_folds(len(samples), n_splits=n_splits, seed=fold_seed)

        cv = run_cv_gpdflow_likelihood(
            samples, folds, device, thres_vec.reshape(1, -1),
            epochs=epochs, samples_multiplier=samples_multiplier,
            init_from_gpd=init_from_gpd,
            verbose=verbose, **fit_kw,
        )

        sigma_hat[r] = cv['sigma_hat']
        gamma_hat[r] = cv['gamma_hat']
        chi_mat = cv_pairwise_chi(cv['samples_simu_folds'], p_sim)
        chi_hat_pairs[r] = chi_mat[iu]
        best_epochs[r] = [h['best_epoch'] for h in cv['hist']]
        best_vals[r] = cv['best_vals']
        p_cluster_arr[r] = p_cluster
        p_sim_arr[r] = p_sim

    return {
        'sigma_hat': sigma_hat,
        'gamma_hat': gamma_hat,
        'chi_hat_pairs': chi_hat_pairs,
        'chi_hat_offdiag_mean': chi_hat_pairs.mean(axis=1),
        'sigma_true': sigma_true,
        'gamma_true': gamma_true,
        'chi_true': chi_true,
        'chi_true_mat': chi_true_mat,
        'u_true': u_true,
        'pair_idx': iu,
        'best_epochs': best_epochs,
        'best_vals': best_vals,
        'p_cluster': p_cluster_arr,
        'p_sim': p_sim_arr,
        'threshold': threshold,
        'rep_seeds': rep_seeds,
        'config': {
            'n_sim': n_sim, 'base_seed': base_seed, 'threshold': threshold,
            'sim_kwargs': sim_cfg,
            'model_kwargs': dict(_SIM2_MODEL_DEFAULTS, **(model_kwargs or {})),
            'init_from_gpd': init_from_gpd,
            'n_splits': n_splits, 'fold_seed': fold_seed, 'chi_q': chi_q,
        },
    }


def summarize_simulation_study_2(results):
    """Monte Carlo bias / variance / RMSE of the ``run_simulation_study_2`` estimators.

    Identical computation to ``summarize_simulation_study_1`` -- it broadcasts a
    per-pair ``chi_true`` vector as-is -- exposed under a Scenario-2 name for
    symmetry.  See ``summarize_simulation_study_1`` for the returned structure.
    """
    return summarize_simulation_study_1(results)
