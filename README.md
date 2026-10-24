# GPDFlow: Generative Multivariate Threshold Exceedance Modeling via Normalizing Flows

This repository contains the code to reproduce the results in the paper

**"GPDFlow: Generative Multivariate Threshold Exceedance Modeling via Normalizing Flows"**
by Chenglei Hu and Daniela Castro-Camilo.

GPDFlow models multivariate threshold exceedances by embedding a normalizing flow in the
multivariate generalized Pareto distribution (mGPD). The marginal GPD scale and shape
parameters are estimated explicitly, and the flow models the extremal dependence structure.
The model has an exact likelihood and fast sampling. It is implemented in PyTorch with the
`normflows` package, using Real NVP or neural spline flow (NSF) coupling layers.

## Reproducing the paper

| Paper section | Notebook | Description |
|---|---|---|
| Simulation study 1 | `Simulation_1S.ipynb` | Monte Carlo study of the likelihood-based GPDFlow on data from a reverse-exponential mGPD. Reports bias and RMSE of the marginal parameters and of pairwise χ against the known truth. |
| Simulation study 2 | `Simulation_2S_cv_study.ipynb` | One-factor Gumbel copula with Burr XII margins (d = 20). Compares the likelihood-based GPDFlow, the MMD-based GPDFlow, a parametric mGPD and a probit-copula NSF baseline on pairwise χ, marginal quantiles, joint and exactly-k exceedance probabilities, and quantiles of the maximum. |
| Application | `Exchange_rate_application.ipynb` | Daily GBP exchange rates from the Bank of England (downloaded inside the notebook). The notebook fits ARMA–GARCH filters through `rpy2`/`fGarch`, fits the four models to the threshold exceedances of the residuals, and compares pairwise χ/ω, simultaneous-exceedance probabilities and exceedance-breadth distributions. |

The hyperparameters hard-coded in `Simulation_2S_cv_study.ipynb` and
`Exchange_rate_application.ipynb` come from the tuning scripts below. You don't need to rerun those scripts to reproduce the paper.

### Supplementary notebooks

- `Exchange_rate_threshold_sensitivity.ipynb` refits the four exchange-rate models at
  threshold levels q ∈ {0.95, 0.96, 0.97, 0.975}, using the configurations tuned by
  `tune_threshold_hparams.py`.
- `threshold_sensitivity_check.ipynb` plots the tuned models' cross-validated χ RMSE
  against the threshold level, from the outputs of both tuning scripts.
- `Simulation_2S_best_config.ipynb` fits each model once on the Simulation 2 data. It is the
  single-run counterpart of `Simulation_2S_cv_study.ipynb`.

## Python modules

| File | Contents |
|---|---|
| `GPDFlow_S.py` | The GPDFlow model `S_mGPD_NF`, based on the S-representation S = Z − max(Z) of the standardized mGPD with one conditional flow per face. Also contains the marginal GPD transform `DataTransform`, the conditional Real NVP (`ConditionalMaskedAffineFlow`) and spline (`ConditionalSplineCoupling`) coupling layers, and both training objectives: the likelihood and the likelihood-free MMD loss. |
| `cv_utils.py` | Builders for the flows (`build_conditional_flow`, `build_baseline_flow`, `flow_type='realnvp'` or `'nsf'`) and the probit-copula baseline. Training loops with early stopping (`fit_gpdflow_likelihood`, `fit_gpdflow_mmd`, `fit_flow_baseline`). 5-fold cross-validation and full-data drivers. The data-generating processes and Monte Carlo drivers for both simulation studies (`run_simulation_study_1`, `run_simulation_study_2`), and the true/estimated exceedance-probability summaries. |
| `parametric_mGPD.py` | Parametric mGPD baseline (`ParametricMGPD`, `fit_parametric_mgpd`) with five T-representation generators: `gumbel`, `rev_gumbel`, `rev_exp`, `log_gamma` and `gaussian`. |
| `Common_Functions.py` | Shared utilities: simulating the reverse-exponential T-mGPD, theoretical and empirical pairwise χ(p)/ω(p), simultaneous-exceedance probabilities, and the error metrics (RMSE of γ, log σ and χ, TDE, LJPE, RMSLE). |
| `tune_threshold_hparams.py` | Per-threshold hyperparameter search for the exchange-rate models. It runs a coordinate-wise grid search scored by the 5-fold CV RMSE between the fitted and empirical χ(0.98), and writes `Data_and_Model/threshold_sensitivity_nsf/best_configs.json`. |
| `tune_sim2s_hparams.py` | The same search for Simulation 2, scored against the true χ(0.999). It writes `Data_and_Model/sim2s_nsf/best_configs.json`. |

Example usage of the tuning scripts:

```bash
python tune_threshold_hparams.py --q 0.95 0.96 0.97 0.975 --models lik mmd bl par
python tune_sim2s_hparams.py     --q 0.95 0.96 0.97       --models lik mmd bl par
```

Add `--smoke` for a quick test run. Here `lik` is the likelihood-based GPDFlow, `mmd` the
MMD-based GPDFlow, `bl` the probit-copula flow baseline and `par` the parametric mGPD.

## Directories

- `Data_and_Model/` holds the downloaded exchange-rate data (`BoE_exchange_rates.csv`),
  the tuning results, and the figures written by the notebooks (`Data_and_Model/Plots/`).

## Requirements

- Python 3 with `torch`, `normflows`, `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`,
  `pyvinecopulib` (for the copula data-generating processes) and `requests`.
- R with the `fGarch` package, accessed from Python through `rpy2`. Only the exchange-rate
  notebooks need it.
- A CUDA GPU is used automatically when one is available. Otherwise the code runs on the CPU.

## Notes

- The notebooks set the output folder through a `dir_out` variable with an absolute path.
  Change it to a location on your machine before running them.
- The manuscript notebooks have a `SMOKE` flag. Set it to `True` for a fast top-to-bottom
  check before the full run.
