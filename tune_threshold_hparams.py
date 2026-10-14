"""Per-threshold hyperparameter tuning for Exchange_rate_threshold_sensitivity.ipynb.

For every threshold level q and model (likelihood GPDFlow, two-stage MMD GPDFlow,
probit-copula RealNVP, parametric mGPD) a small coordinate-wise grid search over
the Table 2 search spaces is run with 5-fold CV.  A config is scored by

    RMSE( fold-mean chi_hat(0.98),  empirical chi(0.98) )   over all pairs,

where chi_hat is computed from each fold model's draws at the re-levelled
cluster-conditional level q_clus = 1 - (1 - 0.98) / p_cluster (as in the notebook).

Hidden widths: Table 2 reports H for d = 20, i.e. H = factor * 19 for GPDFlow
(latent_size_factor) and H = factor * 20 for the baseline (hidden_mult), so the
search space is coded directly in factors {2, 3, 4}.

Every scored config is cached in RES_DIR/tuning/q{q}_{model}.json (resumable);
the winner of each (q, model) is merged into RES_DIR/best_configs.json, which
the notebook reads.  The median of the winner's per-fold early-stop epochs is
stored with it and used as ``stop_at_epoch`` of the notebook's full-data fit.

Usage:
    python tune_threshold_hparams.py --q 0.95 0.96 --models lik mmd bl par [--smoke]
"""

import argparse
import fcntl
import json
import os
import time

import numpy as np
import torch

import Common_Functions as cf
from cv_utils import (
    make_folds,
    cv_pairwise_chi,
    run_cv_gpdflow_likelihood,
    run_cv_gpdflow_mmd,
    run_cv_flow_baseline,
    run_cv_parametric_mgpd,
)

Q_EVAL = 0.98
SIM_MULT = 200
DIR_OUT = "/home/pgrad2/2448355h/My_PhD_Project/GPDFlow/Data_and_Model/"
ORIGIN_CACHE = DIR_OUT + "threshold_sensitivity/samples_origin.npy"

# --- Non-tuned settings: as in Exchange_rate_threshold_sensitivity.ipynb ------
BASE = {
    'lik': dict(penalty_lambda=1e4, seed=1234),
    'mmd': dict(penalty_lambda=1e4, mmd_samples_per_face=128,
                mmd_censor_c=0, mmd_censor_tau=0.2, seed=1234),
    'bl': dict(warmup=20, seed=0),
    'par': dict(seed=1234),
}
PATIENCE = dict(lik=100, mmd=100, bl=40)       # par: driver default (40)

# --- Search spaces (Table 2; widths as factors) -------------------------------
# Coordinate-wise sweep: start at START, sweep each key of SWEEP in order,
# keep the best value before moving on.
START = {
    'lik': dict(num_layers=3, latent_size_factor=3, marg_lambda=30),
    'mmd': dict(num_layers=3, latent_size_factor=2,
                mmd_indicator_alpha=5, mmd_indicator_rho=3),
    'bl': dict(num_layers=6, hidden_mult=4),
    'par': dict(generator='gaussian'),
}
SWEEP = {
    'lik': [('num_layers', [3, 4, 5, 6]),
            ('latent_size_factor', [2, 3, 4]),
            ('marg_lambda', [50, 100, 200, 300])],
    'mmd': [('num_layers', [3, 4, 5, 6]),
            ('latent_size_factor', [2, 3, 4]),
            ('mmd_indicator_alpha', [1, 3, 5, 10]),
            ('mmd_indicator_rho', [1, 3, 5])],
    'bl': [('num_layers', [3, 4, 5, 6]),
           ('hidden_mult', [2, 3, 4])],
    'par': [('generator', ['gaussian', 'gumbel', 'rev_gumbel', 'rev_exp', 'log_gamma'])],
}


def _cfg_key(cfg):
    return json.dumps(cfg, sort_keys=True)


def _load_json(path, default):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return default


def _merge_best(path, q, model, entry):
    """Merge one winner into best_configs.json (locked: several q may run at once)."""
    with open(path + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        best = _load_json(path, {})
        best.setdefault(str(q), {})[model] = entry
        tmp = path + '.tmp'
        with open(tmp, 'w') as f:
            json.dump(best, f, indent=1)
        os.replace(tmp, path)
        fcntl.flock(lk, fcntl.LOCK_UN)


def score_config(model, cfg, samples, folds, device, thres, q_clus, chi_emp, iu,
                 epochs):
    """5-fold CV of one config -> (rmse_chi of the fold-mean chi, fold epochs)."""
    if model == 'lik':
        cv = run_cv_gpdflow_likelihood(
            samples, folds, device, thres, epochs=epochs, samples_multiplier=SIM_MULT,
            patience=PATIENCE['lik'], **BASE['lik'], **cfg)
        blocks = cv['samples_simu_folds']
    elif model == 'mmd':
        cv = run_cv_gpdflow_mmd(
            samples, folds, device, thres, epochs=epochs, samples_multiplier=SIM_MULT,
            patience=PATIENCE['mmd'], **BASE['mmd'], **cfg)
        blocks = cv['samples_simu_2s_folds']
    elif model == 'bl':
        cv = run_cv_flow_baseline(
            samples, folds, device, thres, epochs=epochs, samples_multiplier=SIM_MULT,
            patience=PATIENCE['bl'], **BASE['bl'], **cfg)
        blocks = cv['samples_simu_bl_folds']
    else:
        cv = run_cv_parametric_mgpd(
            samples, folds, device, thres, cfg['generator'], epochs=epochs,
            samples_multiplier=SIM_MULT, **BASE['par'])
        blocks = cv['samples_simu_par_folds']

    fold_epochs = [int(h['best_epoch']) for h in cv['hist']]
    blocks = [b[np.isfinite(b).all(axis=1)] for b in blocks]
    del cv
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    if min(len(b) for b in blocks) < 100:
        return float('inf'), fold_epochs
    chi_cv = cv_pairwise_chi(blocks, q_clus)
    rmse = float(cf.rmse_chi(chi_cv[iu], chi_emp[iu]))
    return (rmse if np.isfinite(rmse) else float('inf')), fold_epochs


def tune(q, model, samples_origin, chi_emp, iu, device, res_dir, epochs):
    n_obs = samples_origin.shape[0]
    thres = np.quantile(samples_origin, q, axis=0)
    cond = np.any(samples_origin > thres, axis=1)
    samples = samples_origin[cond, :] - thres.reshape(1, -1)
    p_cluster = samples.shape[0] / n_obs
    q_clus = 1.0 - (1.0 - Q_EVAL) / p_cluster
    folds = make_folds(len(samples), n_splits=5, seed=1234)
    print(f'\n##### q = {q} | {model}: slice {samples.shape}, '
          f'p_cluster = {p_cluster:.4f}, q_clus = {q_clus:.4f}', flush=True)

    log_path = res_dir + f'tuning/q{q}_{model}.json'
    log = _load_json(log_path, [])
    scored = {_cfg_key(r['cfg']): r for r in log}

    def evaluate(cfg):
        key = _cfg_key(cfg)
        if key in scored:
            return scored[key]['rmse_chi']
        print(f'--- q = {q} | {model} | {cfg}', flush=True)
        t0 = time.time()
        try:
            rmse, fold_epochs = score_config(model, cfg, samples, folds, device,
                                             thres, q_clus, chi_emp, iu, epochs)
        except Exception as e:                       # e.g. degenerate generator
            print(f'    FAILED: {e!r}', flush=True)
            rmse, fold_epochs = float('inf'), []
        rec = dict(cfg=dict(cfg), rmse_chi=rmse, fold_epochs=fold_epochs,
                   median_epoch=int(np.median(fold_epochs)) if fold_epochs else None,
                   elapsed_sec=round(time.time() - t0, 1))
        print(f'    RMSE_chi(CV) = {rmse:.5f}  fold epochs = {fold_epochs}  '
              f'({rec["elapsed_sec"]} s)', flush=True)
        log.append(rec)
        scored[key] = rec
        with open(log_path, 'w') as f:
            json.dump(log, f, indent=1)
        return rmse

    cur = dict(START[model])
    for name, grid in SWEEP[model]:
        scores = {v: evaluate({**cur, name: v}) for v in grid}
        cur[name] = min(grid, key=lambda v: scores[v])   # ties -> earlier value

    win = scored[_cfg_key(cur)]
    print(f'===== q = {q} | {model}: best {cur}  RMSE_chi(CV) = {win["rmse_chi"]:.5f}',
          flush=True)
    _merge_best(res_dir + 'best_configs.json', q, model,
                dict(cfg=cur, median_epoch=win['median_epoch'],
                     cv_rmse_chi=win['rmse_chi'], fold_epochs=win['fold_epochs']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=float, nargs='+', default=[0.95, 0.96, 0.97, 0.975])
    ap.add_argument('--models', nargs='+', default=['lik', 'mmd', 'bl', 'par'],
                    choices=['lik', 'mmd', 'bl', 'par'])
    ap.add_argument('--smoke', action='store_true', help='5 epochs, separate output dir')
    args = ap.parse_args()

    epochs = 5 if args.smoke else 500
    res_dir = DIR_OUT + ('threshold_sensitivity_smoke/' if args.smoke
                         else 'threshold_sensitivity/')
    os.makedirs(res_dir + 'tuning/', exist_ok=True)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print('device:', device)

    samples_origin = np.load(ORIGIN_CACHE)
    iu = np.triu_indices(samples_origin.shape[1], k=1)
    chi_emp = cf.pairwise_chi_from_data(samples_origin, p=Q_EVAL)

    for q in args.q:
        for m in args.models:
            tune(q, m, samples_origin, chi_emp, iu, device, res_dir, epochs)


if __name__ == '__main__':
    main()
