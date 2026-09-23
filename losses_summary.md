# GPDFlow_S.py: Marginal Loss Functions

## Setup

For margin $j$, the marginal transform parameters are

$$
\sigma_j = \exp(\text{log\_sigma}_j), \qquad
\gamma_j = 0.5\tanh(\theta_j).
$$

Exceedances are observations $x_{ij}$ satisfying $x_{ij} > 0$ and $\sigma_j + \gamma_j x_{ij} > 0$
(the `_exceedance_mask` condition). Margin $j$ is *active* if it has at least
`quantile_min_exceedances` such observations.

The standardizing transform (`DataTransform.inverse_transform`) is

$$
y_{ij} = g_{\mathrm{std}}(x_{ij}; \sigma_j, \gamma_j) =
\begin{cases}
\dfrac{x_{ij}}{\sigma_j}, & \gamma_j = 0, \\[8pt]
\dfrac{1}{\gamma_j}\log\!\left(1 + \dfrac{\gamma_j x_{ij}}{\sigma_j}\right), & \gamma_j \neq 0.
\end{cases}
$$

## `mse_quantile_loss`: a differentiable Q–Q goodness-of-fit discrepancy

**Formal statement.** For each active margin $j$, sort the exceedance-standardized values
$y_{ij}=g_{\mathrm{std}}(x_{ij};\sigma_j,\gamma_j)$ to obtain order statistics
$y_{(1)}\le\cdots\le y_{(n_j)}$, and compute the corresponding unit-exponential quantiles under the
Weibull plotting-position convention,

$$
p_k=\frac{k}{n_j+1}, \qquad t_k=-\log(1-p_k), \qquad k=1,\dots,n_j.
$$

The per-margin loss is the mean squared error between empirical order statistics and theoretical
quantiles,

$$
L_j = \frac{1}{n_j}\sum_{k=1}^{n_j}\big(y_{(k)}-t_k\big)^2,
$$

and the total loss averages over the set $\mathcal{A}$ of active margins:

$$
\mathcal{L}_{\mathrm{mse\text{-}q}} = \text{mse\_quantile\_loss} = \frac{1}{|\mathcal{A}|}\sum_{j\in\mathcal{A}} L_j.
$$

**Statistical justification.** Under the model's own null hypothesis — that
$Z_j\mid Z_j>0\sim\mathrm{Exp}(1)$ exactly, a property inherited from the T-representation of the
mGPD — the plotting positions $p_k=k/(n_j+1)$ are the classical unbiased estimators of the CDF
value at the $k$-th order statistic of an i.i.d. sample (since $\mathbb{E}[U_{(k)}]=k/(n+1)$ for
uniform order statistics), and $t_k=F^{-1}(p_k)$ are the corresponding reference quantiles. $L_j$
is thus a mean squared quantile–quantile (Q–Q) discrepancy: it measures the deviation of the
empirical order statistics from the 45° reference line of a standard exponential Q–Q plot, the
diagnostic tool classically used to assess univariate GPD goodness of fit. Framing this diagnostic
as a differentiable loss makes it directly optimizable, rather than only usable as a post hoc check.

Its role is complementary to `loss_marg` below: whereas `loss_marg` is a pointwise likelihood
criterion sensitive to the density at each observation, $\mathcal{L}_{\mathrm{mse\text{-}q}}$ is a
rank-based, distribution-shape criterion sensitive to systematic departures across the whole
empirical distribution (e.g., tail miscalibration that a pointwise likelihood may under-penalize).
Because it operates on ranks rather than raw likelihood values, it is comparatively robust to a
small number of extreme or poorly fit observations dominating the marginal parameter gradient.

**Suggested manuscript language.**

> To further calibrate the marginal parameters against the exponentiality implied by the model's
> own standardization, we introduce a differentiable relaxation of the classical quantile–quantile
> (Q–Q) diagnostic used to assess univariate generalized Pareto fit. For each margin, we compare
> the order statistics of the standardized exceedances against their expected values under the
> unit-exponential null, using the standard Weibull plotting-position estimator $p_k=k/(n+1)$, and
> penalize the mean squared deviation (Eq. Y). This term provides a rank-based, distribution-shape
> criterion that complements the pointwise marginal likelihood $\mathcal{L}_{\mathrm{marg}}$,
> directly discouraging systematic tail miscalibration that a pointwise likelihood objective may
> not sufficiently penalize.

## `loss_marg`: an exact auxiliary marginal likelihood

**Formal statement.** For $x_{ij}>0$ with $\sigma_j+\gamma_j x_{ij}>0$, define

$$
\ell_{ij}(\sigma_j,\gamma_j) = g_{\mathrm{std}}(x_{ij};\sigma_j,\gamma_j) + \log\big(\sigma_j+\gamma_j x_{ij}\big),
$$

where the first term is $y_{ij}$ (as above) and the second uses
`inside_raw` $= \sigma_j + \gamma_j x_{ij}$. Averaging over the $N_j$ exceedances within each
active margin, then over the $d$ active margins:

$$
\mathcal{L}_{\mathrm{marg}} = \text{loss\_marg} = \frac{1}{|\mathcal{A}|} \sum_{j \in \mathcal{A}} \frac{1}{N_j}
\sum_{i \,:\, x_{ij} > 0} \ell_{ij}(\sigma_j,\gamma_j).
$$

**Statistical justification.** The T-representation of the mGPD implies an exact distributional
identity: each standardized margin, conditional on exceedance, is unit-exponential,
$Z_j = g_{\mathrm{std}}(X_j;\sigma_j,\gamma_j) \mid Z_j>0 \sim \mathrm{Exp}(1)$ — the classical
closed-form univariate margin of a multivariate generalized Pareto distribution. Applying the
change-of-variables formula to this known marginal law gives the univariate conditional exceedance
density in the original scale,

$$
f_{X_j}(x;\sigma_j,\gamma_j) = \exp\{-g_{\mathrm{std}}(x;\sigma_j,\gamma_j)\}\cdot\big(\sigma_j+\gamma_j x\big)^{-1}, \qquad x>0,\ \sigma_j+\gamma_j x>0,
$$

which is precisely the classical generalized Pareto density with scale $\sigma_j$ and shape
$\gamma_j$. Consequently $\ell_{ij}=-\log f_{X_j}(x_{ij};\sigma_j,\gamma_j)$, and
$\mathcal{L}_{\mathrm{marg}}$ is the exact negative log-likelihood of the univariate GPD margins
implied by the model.

This term is not an ad hoc regularizer but a **model-consistent auxiliary likelihood**: it exploits
an analytic invariant of the mGPD class (the univariate margins are known in closed form,
independent of the generator's dependence structure) to supply a direct, low-variance gradient for
$(\sigma_j,\gamma_j)$. This complements the joint negative log-likelihood, whose dependence on
$(\sigma_j,\gamma_j)$ is mediated only implicitly through a numerically approximated
one-dimensional integral over the generator density, and can therefore carry higher variance or
slower-converging gradient information for the marginal parameters, particularly early in training
or under coarse numerical integration.

**Suggested manuscript language.**

> Because every margin of the mGPD is, conditional on exceedance, exactly unit-exponential after
> standardization, the marginal parameters $(\sigma_j,\gamma_j)$ admit a closed-form univariate
> likelihood independent of the generator's dependence structure. We include this exact marginal
> log-likelihood, $\mathcal{L}_{\mathrm{marg}}$ (Eq. X), as an auxiliary training objective. Unlike
> the joint likelihood, whose dependence on $(\sigma_j,\gamma_j)$ is mediated through a numerically
> approximated integral over the flow-induced generator density, $\mathcal{L}_{\mathrm{marg}}$
> provides an unbiased, closed-form gradient for the marginal parameters, stabilizing their
> estimation independently of the fidelity of the numerical integration used for the joint density.

## `support_loss`

This term (in `loss_components`, GPDFlow_S.py:549-550) enforces the support constraint of the
marginal GPD: the standardizing transform requires $\sigma_j + \gamma_j x_{ij} > 0$ for every
margin (this is the domain of the GPD — when $\gamma_j < 0$ it corresponds to $x_j$ being bounded
above by $-\sigma_j/\gamma_j$). It is a soft (squared-hinge) penalty rather than a hard constraint,
so it can be optimized with gradient descent.

Let

$$
\text{inside}_{ij} = \sigma_j + \gamma_j x_{ij}.
$$

The violation of the support condition for observation $i$, margin $j$, is

$$
r_{ij} = \operatorname{ReLU}(-\text{inside}_{ij}) = \max\big(0,\, -(\sigma_j + \gamma_j x_{ij})\big),
$$

which is zero when $x_{ij}$ lies in the valid GPD support and positive (equal to the magnitude of
the violation) otherwise. The support loss is the mean over observations of the summed squared
violation across margins:

$$
\text{support\_loss} = \frac{1}{n} \sum_{i=1}^{n} \sum_{j=1}^{d} r_{ij}^{\,2}
= \frac{1}{n} \sum_{i=1}^{n} \sum_{j=1}^{d} \Big[\operatorname{ReLU}\big(-(\sigma_j+\gamma_j x_{ij})\big)\Big]^{2}.
$$

It enters the total loss with weight `penalty_lambda`:

$$
\text{support\_loss\_weighted} = \lambda_{\text{support}} \cdot \text{support\_loss}.
$$
