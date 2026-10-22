"""Per-threshold hyperparameter tuning for Simulation_2S_best_config.ipynb.

Same NSF search spaces, non-tuned settings and coordinate-wise 5-fold-CV sweep
as tune_threshold_hparams.py, applied to the simulated data of the notebook
(one-factor Gumbel copula + Burr XII margins, n = 10000, d = 20).  A config is
scored by

    RMSE( fold-mean chi_hat(0.999),  reference chi(0.999) )   over all pairs,

where the reference chi_theo is the rank-based pairwise chi at p = 0.999 of a
2e6-row Monte-Carlo draw of the copula (notebook cell 5), and chi_hat is
computed from each fold model's draws at q_clus = 1 - (1 - 0.999) / p_cluster.

Every scored config is cached in RES_DIR/tuning/q{q}_{model}.json (resumable);
the winner of each (q, model) and the median of its per-fold early-stop epochs
are merged into RES_DIR/best_configs.json.

Usage:
    python tune_sim2s_hparams.py --q 0.95 0.96 0.97 --models lik mmd bl par [--smoke]
"""

import argparse
import fcntl
import os

import numpy as np
import pyvinecopulib as pv
import torch
from scipy.stats import burr12

import Common_Functions as cf
import tune_threshold_hparams as tts

Q_EVAL = 0.999
N_SAMPLES = 10000
N_CHI_TRUTH = 2_000_000
DIM = 20
RES_NAME = "sim2s_nsf"

# --- DGP: as in Simulation_2S_best_config.ipynb (cell 4) ----------------------
C_BASE = np.array([2.0, 3.0, 2.5, 1.5, 3.0])
K_BASE = np.array([1.5, 1.0, 2.0, 3.0, 3.5])
C_VALS = np.tile(C_BASE, DIM // len(C_BASE) + 1)[:DIM]
K_VALS = np.tile(K_BASE, DIM // len(K_BASE) + 1)[:DIM]
THETA = 1.0 + np.linspace(0, 7.0, DIM + 2)[1:-1]        # (dim,) in (1, 8)


def simulate_factor_copula(n, thetas, seed=42):
    rng = np.random.default_rng(seed)
    V = rng.random(n)
    U = np.empty((n, len(thetas)))
    for j, th in enumerate(thetas):
        bicop = pv.Bicop(pv.BicopFamily.gumbel, parameters=np.array([[th]]))
        U[:, j] = bicop.hinv1(np.column_stack([V, rng.random(n)]))   # U_j | V
    return U


def load_data(cache_dir):
    """(samples_origin, chi_theo), generated once and cached (locked: several q may start at once)."""
    path_x, path_chi = cache_dir + 'samples_origin.npy', cache_dir + 'chi_theo.npy'
    with open(cache_dir + 'data.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if not (os.path.exists(path_x) and os.path.exists(path_chi)):
            U = simulate_factor_copula(N_SAMPLES, THETA, seed=42)
            samples_origin = np.column_stack([
                burr12.ppf(U[:, j], c=C_VALS[j], d=K_VALS[j]) for j in range(DIM)])
            U_chi = simulate_factor_copula(N_CHI_TRUTH, THETA, seed=42)
            chi_theo = cf.pairwise_chi_from_data(U_chi, p=Q_EVAL)
            del U_chi
            np.save(path_x, samples_origin)
            np.save(path_chi, chi_theo)
        fcntl.flock(lk, fcntl.LOCK_UN)
    return np.load(path_x), np.load(path_chi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=float, nargs='+', default=[0.95, 0.96, 0.97])
    ap.add_argument('--models', nargs='+', default=['lik', 'mmd', 'bl', 'par'],
                    choices=['lik', 'mmd', 'bl', 'par'])
    ap.add_argument('--smoke', action='store_true', help='5 epochs, separate output dir')
    args = ap.parse_args()

    epochs = 5 if args.smoke else 500
    res_dir = tts.DIR_OUT + RES_NAME + ('_smoke/' if args.smoke else '/')
    os.makedirs(res_dir + 'tuning/', exist_ok=True)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print('device:', device)

    samples_origin, chi_theo = load_data(res_dir)
    iu = np.triu_indices(samples_origin.shape[1], k=1)
    chi_off = chi_theo[iu]
    print(f'samples_origin {samples_origin.shape} | chi_ref off-diagonal: min {chi_off.min():.3f} '
          f'median {np.median(chi_off):.3f} max {chi_off.max():.3f}', flush=True)

    for q in args.q:
        for m in args.models:
            tts.tune(q, m, samples_origin, chi_theo, iu, device, res_dir, epochs,
                     q_eval=Q_EVAL)


if __name__ == '__main__':
    main()
