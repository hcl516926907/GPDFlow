import numpy as np
import torch
import normflows as nf
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
from GPDFlow import T_mGPD_NF

dir_out = "/home/pgrad2/2448355h/My_PhD_Project/01_Output/GPDFlow/"

def sim_revexp_T_mgpd(n, d, a, beta, sig, gamma, MGPD=True, std=False):
    """
    Sample from the T representation of mGPD with the reverse exponential generator, i.e. 
    f_T(t_1,...,t_d) = \prod_{j=1}^d exp((t_j + beta_j)/a_j)

    Args: 
        n (integer): number of samples 
        d (integer): dimension of the samples
        a (array): scale parameters of the reverse exponential distribution
        beta (array): location parameters of the reverse exponential distribution
        sig (array): marginal scale parameter of the mGPD
        gamma (array): marginal shape parameter of the mGPD
        MGDP (boolean): generate samlples on mGPD scale (i.e. samples are dependent on sig and gamma) if true
        std (boolean): generate samlples on standardized scale (i.e. sig = 1 and gamma = 0) if true
    
    Returns:
        samples of n*d shape or a dictionary that constains two samples, each has n*d shape.


    """
    E = np.random.exponential(1, n )[:,None]
    T_total = []
    for j in range(d):
        U = np.random.uniform(0, 1, n)
        T = a[j]*np.log(U) - beta[j]
        T_total.append(T)
    T_total = np.column_stack(T_total)

    S = T_total - np.max(T_total, axis=1)[:,None]
    Z = E + S 

    if std and not MGPD:
        return Z

    X = []
    for j in range(d):
        if gamma[j] != 0:
            X.append(sig[j] * (np.exp(Z[:, j] * gamma[j]) - 1) / gamma[j])
        else:
            X.append(sig[j] * Z[:, j])
    X = np.column_stack(X)

    if MGPD and not std:
        return X
    if std and MGPD:
        return {'X': X, 'Z': Z}
    

def chi_theorical(alpha1, alpha2):
    """
    Calculate the theorcial chi of the T representation of mGPD with the reverse exponential generator in dimension 2.
    Formula can be found in the supporting material of "Peaks over thresholds modeling with multivariate generalized Pareto distributions". 


    Args: 
        alpha1 (float): inverse of the first scale parameter in the reverse exponential distribution 
        alpha2 (float): inverse of the second scale parameter in the reverse exponential distribution 
    
    Returns:
        chi (float): theorical chi value 

    """

    if alpha1 <= 0 or alpha2 <= 0:
        raise ValueError("alpha1 and alpha2 must be greater than 0")
    
    # Identify alpha_(1) = max(alpha1, alpha2) and alpha_(2) = min(alpha1, alpha2)
    alpha_max = max(alpha1, alpha2)
    alpha_min = min(alpha1, alpha2)
    
    # Calculate the components of the formula
    term1 = (1 + 1 / alpha_max) / (1 + 1 / alpha_min)
    term1_power = term1 ** (1 + alpha_min)
    term2 = (alpha_max / alpha_min) * (1 / (1 + alpha1 + alpha2))
    
    # Final chi value
    chi = 1 - term1_power * term2
    return chi

def pairwise_chi_theorical(alpha_vector):
    """
    Calculate the pairwise theorcial chi of the T representation of mGPD with the reverse exponential generator 
    for all pairs.

    Args: 
        alpha_vector (array): inverse of the scale parameter in the reverse exponential distribution 
    
    Returns:
        chi_dict (dict): dictionary with pair names as keys and corresponding chi values as values
    """

    n = len(alpha_vector)
    chi_dict = {}
    
    for i in range(n-1):
        for j in range(i+1,n):
            chi_dict[f"{i+1}-{j+1}"] =  chi_theorical(alpha_vector[i], alpha_vector[j])
    
    return chi_dict


def ar1_correlation(d, rho):
    """
    AR(1) correlation matrix R[i, j] = rho ** |i - j|  (d × d, unit diagonal).

    Args:
        d (integer): dimension
        rho (float): base correlation in (-1, 1)

    Returns:
        R (ndarray): d × d correlation matrix
    """
    idx = np.arange(d)
    return rho ** np.abs(idx[:, None] - idx[None, :])


def t_copula_pairwise_chi(rho, nu):
    """
    Closed-form upper-tail dependence coefficient chi for a Student-t copula:

        chi = 2 * T_{nu+1}( -sqrt( (nu + 1) (1 - rho) / (1 + rho) ) )

    Args:
        rho (float or ndarray): scalar correlation, or a d × d correlation matrix
                                (e.g. from ar1_correlation)
        nu (float): degrees of freedom of the t-copula

    Returns:
        chi (float or ndarray): same shape as rho; if rho is a d × d matrix the
                                diagonal is set to 1.0
    """
    from scipy.stats import t as _t
    rho = np.asarray(rho, dtype=float)
    arg = -np.sqrt(np.clip((nu + 1.0) * (1.0 - rho) / (1.0 + rho), 0.0, None))
    chi = 2.0 * _t.cdf(arg, df=nu + 1.0)
    if chi.ndim == 2 and chi.shape[0] == chi.shape[1]:
        np.fill_diagonal(chi, 1.0)
    return chi


def pairwise_chi_from_data(X, p=0.95):
    """
    Compute the d×d matrix of pairwise empirical upper tail-dependence coefficients
    from an n×d data matrix.

    Uses the standard extreme-value estimator

        chi_hat_{ij}(p) = P_hat(U_i > p, U_j > p) / (1 - p)

    where U_{t,i} = rank(X_{t,i}) / n are the probability-integral-transformed
    (pseudo-uniform) observations.  The estimator is symmetric in i and j.

    Args:
        X (ndarray): n × d data matrix (raw or standardised scale)
        p (float):   quantile level in (0, 1), e.g. 0.95 or 0.99

    Returns:
        chi_matrix (ndarray): d × d symmetric matrix; diagonal entries = 1.0
    """
    X = np.asarray(X, dtype=float)
    n, d = X.shape

    # Rank-based probability-integral transform → U in (0, 1)
    from scipy.stats import rankdata
    U = np.column_stack([rankdata(X[:, j]) / n for j in range(d)])  # (n, d)

    exceed = U > p          # (n, d) bool
    n_marg = exceed.sum(0)  # (d,): marginal exceedance counts

    # n_joint[i, j] = #{t : U_{t,i} > p AND U_{t,j} > p}  (symmetric)
    n_joint = exceed.T.astype(float) @ exceed.astype(float)  # (d, d)

    # Denominator: expected count under independence = n * (1-p)
    denom = n * (1.0 - p)

    chi_matrix = n_joint / denom
    np.fill_diagonal(chi_matrix, 1.0)
    return chi_matrix


def pairwise_omega_from_data(X, p=0.95):
    """
    Compute the d×d matrix of pairwise empirical coefficients of tail correlation
    from an n×d data matrix.

    Uses the estimator

        omega_hat_{ij}(p) = P_hat(U_i > p or U_j > p) / (1 - p)

    where U_{t,i} = rank(X_{t,i}) / n are the probability-integral-transformed
    (pseudo-uniform) observations.  The estimator is symmetric in i and j.

    Args:
        X (ndarray): n × d data matrix (raw or standardised scale)
        p (float):   quantile level in (0, 1), e.g. 0.95 or 0.99

    Returns:
        omega_matrix (ndarray): d × d symmetric matrix; diagonal entries = 1.0
    """
    X = np.asarray(X, dtype=float)
    n, d = X.shape

    # Rank-based probability-integral transform → U in (0, 1)
    from scipy.stats import rankdata
    U = np.column_stack([rankdata(X[:, j]) / n for j in range(d)])  # (n, d)

    exceed = U > p          # (n, d) bool
    n_marg = exceed.sum(0)  # (d,): marginal exceedance counts

    # n_and[i, j] = #{t : U_{t,i} > p AND U_{t,j} > p}  (symmetric)
    n_and = exceed.T.astype(float) @ exceed.astype(float)  # (d, d)

    # Inclusion-exclusion: #{t : U_{t,i} > p OR U_{t,j} > p}
    n_or = n_marg[:, None] + n_marg[None, :] - n_and

    # Denominator: expected count under independence = n * (1-p)
    denom = n * (1.0 - p)

    omega_matrix = n_or / denom
    np.fill_diagonal(omega_matrix, 1.0)
    return omega_matrix


def pairwise_chi_from_data_dict(X, p=0.95):
    """
    Same as pairwise_chi_from_data but returns a dict keyed by "i-j" (1-indexed),
    matching the format used by pairwise_chi_theorical.

    Args:
        X (ndarray): n × d data matrix
        p (float):   quantile level

    Returns:
        chi_dict (dict): keys "i-j" for i < j, values chi_hat_{ij}(p)
    """
    chi_matrix = pairwise_chi_from_data(X, p)
    d = chi_matrix.shape[0]
    chi_dict = {}
    for i in range(d - 1):
        for j in range(i + 1, d):
            chi_dict[f"{i+1}-{j+1}"] = float(chi_matrix[i, j])
    return chi_dict


def empirical_upper_tail_dependence(X, p):
    """
    Calculate the empirical pairwise chi(p) of the data

    Args: 
        X (array): two-dimensional data matrix
        p (array): quantile level of the q.  
    
    Returns:
        chi_u (array): one-dimensional array that contains the chi(p) at p-quantile
    """

    # Compute the quantiles
    quantile = np.quantile(X, p, axis = 0)
    
    # Count joint exceedances
    joint_exceedance = np.sum((X[:,0] > quantile[0]) & (X[:,1] > quantile[1]))
    exceedance_X1 = np.sum(X[:,0] > quantile[0])
    
    # Calculate lambda_u
    if exceedance_X1 == 0:  # Avoid division by zero
        return 0.0
    chi_u = joint_exceedance / exceedance_X1
    return chi_u




def GPDFlow(dim):
    """
    A quick set up of a GPDFlow with the architecture in the simulation section.

    Args: 
        dim (integer): dimension of the model
        model (NN.module): GPDFlow with initialization
        n_monte_carlo: number of samples from the model for monte carlo estimations
        n_experiments: number of repeat of the simulation scenario  
    
    Returns:
        model (NN.module): GPDFlow
    
    """
    base = nf.distributions.DiagGaussian(dim)

    num_layers = 16
    torch.manual_seed(0)

    latent_size = dim
    b = torch.Tensor([1 if i % 2 == 0 else 0 for i in range(latent_size)])
    flows = []
    for i in range(num_layers):
        s = nf.nets.MLP([latent_size, 4 * latent_size, latent_size], init_zeros=True, output_fn='tanh')
        t = nf.nets.MLP([latent_size, 4 * latent_size, latent_size], init_zeros=True, output_fn='tanh')
        if i % 2 == 0:
            flows += [nf.flows.MaskedAffineFlow(b, t, s)]
        else:
            flows += [nf.flows.MaskedAffineFlow(1 - b, t, s)]
        flows += [nf.flows.ActNorm(latent_size)]


    f_T_model = nf.NormalizingFlow(base, flows)
    f_T_model = f_T_model.to(device)



    model = T_mGPD_NF(dim=dim, flow =f_T_model, device=device, s_min=-10,
                    s_max = 10, num_integration_points=1000, penalty_lambda=10000, fix_margin=False)
    return model


def pairwise_chi_empirical(dim, n_monte_carlo, n_experiments=100):
    """
    Load the esimated Calculate the all empirical pairwise chi(p) of the data simulated by GPDFlow

    Args: 
        dim (integer): dimension of the model
        model (NN.module): GPDFlow with initialization
        n_monte_carlo: number of samples from the model for monte carlo estimations
        n_experiments: number of repeat of the simulation scenario  
    
    Returns:
        sample_chi (dict): dictionary with the name of dimension pairs as the key and a list of
        empirical chi in all experiments as the value
    """

    # Compute lambda_u for each quantile
    model = GPDFlow(dim)
    sample_chi = {}
    for i in range(dim-1):
        for j in range(i+1,dim):
            chi_values = []
            for _ in range(n_experiments):  # Step 3: Repeat 200 times
                # load the estimated weight
                model.load_state_dict(torch.load(dir_out + f'model_{dim}D_100_{_}.pt', weights_only=True))
                model.eval()
                samples_obs, samples_std, samples_T= model.sample(n_monte_carlo)
                sampled_data = samples_obs.cpu().data.numpy()
                chi_values.append(empirical_upper_tail_dependence(sampled_data[:,[i,j]], 0.99))
            sample_chi[f"{i+1}-{j+1}"] =  chi_values
    
    return sample_chi


def marginal_parameter_monte_carlo(dim, n_experiments=100):
    """
    Load the estimated sigma and gamma of GPDFlow models in a simulation scenario, 
    and save them in two dictionaries.
    
    Args: 
        dim (integer): dimension of the model.
        n_experiments (integer): Number of repeat of the simulation
    
    Returns:
        sigma_dict(dict)， gamma_dict(dict): two dictionaries contains the estimated sigma and gamma of 
        GPDFlow models in a simulation scenario. 
    
    """
    sigma_dict = {f'{i+1}': [] for i in range(dim)}
    gamma_dict = {f'{i+1}': [] for i in range(dim)}
    model = GPDFlow(dim)
    for _ in range(n_experiments):  # Step 3: Repeat 200 times
        model.load_state_dict(torch.load(dir_out + f'model_{dim}D_100_{_}.pt', weights_only=True))
        model.eval()
        sigma_hat = model.data_transform.get_sigma().cpu().data.numpy()
        gamma_hat = model.data_transform.get_gamma().cpu().data.numpy()
        for i in range(dim):
            sigma_dict[f'{i+1}'].append(sigma_hat[i])
            gamma_dict[f'{i+1}'].append(gamma_hat[i])
    return sigma_dict, gamma_dict


def empirical_tail_dependence_measure(X, p, cond = 'and'):
    """
    Calculate the empirical chi(p) and omega(p) of the data

    Args: 
        X (array): two-dimensional data matrix
        p (array): quantile level of the q.  
        cond (string): an indicator of whether calculating chi(p) or omega(p). If cond == 'and',
        then return chi(p), otherwise omega(p).
    
    Returns:
        measure (array): one-dimensional array that contains the chi(p)/omega(p) at p-quantile
    """

    # Compute the quantiles
    quantile = np.quantile(X, p, axis = 0)
    
    # Count joint exceedances
    joint_exceedance = X[:,0] > quantile[0]
    if cond == 'and':
        joint_exceedance = np.all(X > quantile, axis=1)
    else:
        joint_exceedance = np.any(X > quantile, axis=1)
        
    exceedance_X1 = np.sum(X[:,0] > quantile[0])
    joint_exceedance = np.sum(joint_exceedance)
    # Calculate lambda_u
    if exceedance_X1 == 0:  # Avoid division by zero
        return 0.0
    measure = joint_exceedance / exceedance_X1
    return measure


def simultaneous_exceedance_probability(X, u, m):
    """
    Empirical probability that more than m out of d margins simultaneously
    exceed a common absolute threshold u: P(sum_j 1{X_j > u} > m).

    Args:
        X (array): n x d data matrix whose rows are checked against u.
        u (float or array): absolute threshold applied directly to the columns
        of X. A scalar is compared against every margin; an array of length d
        gives a per-margin threshold. No quantile transform is applied, so X
        and u must be on the same scale (pass simulated samples rescaled back
        to the raw data scale).
        m (int): exceedance-count threshold.

    Returns:
        prob (float): empirical P(sum_j 1{X_j > u} > m).
    """
    X = np.asarray(X, dtype=float)
    n, d = X.shape
    u = np.broadcast_to(np.asarray(u, dtype=float), (d,))
    exceed_count = np.sum(X > u, axis=1)
    return np.mean(exceed_count > m)


def common_threshold_grid(X_ref, n=100, q_lo=0.97, q_hi=0.995):
    """
    Grid of common absolute thresholds over the range where a single threshold
    lies between every margin's q_lo and q_hi empirical quantile:
        lo = max_j quantile(X_ref[:, j], q_lo)
        hi = min_j quantile(X_ref[:, j], q_hi)
    Returns np.linspace(lo, hi, n). Raises if lo >= hi.
    """
    X_ref = np.asarray(X_ref, dtype=float)
    lo = np.quantile(X_ref, q_lo, axis=0).max()
    hi = np.quantile(X_ref, q_hi, axis=0).min()
    if not lo < hi:
        raise ValueError(f"empty threshold grid: lo={lo:.4g} >= hi={hi:.4g}")
    return np.linspace(lo, hi, n)


def rmse_gamma(gamma_hat, gamma_true):
    """
    Compute RMSE for the shape parameter gamma across d margins.

    Args:
        gamma_hat  (array): estimated shape parameters, length d
        gamma_true (array): true shape parameters, length d

    Returns:
        rmse (float)
    """
    gamma_hat = np.asarray(gamma_hat)
    gamma_true = np.asarray(gamma_true)
    return np.sqrt(np.mean((gamma_hat - gamma_true) ** 2))


def rmse_log_sigma(sigma_hat, sigma_true):
    """
    Compute RMSE for the scale parameter sigma on the log scale across d margins.

    Args:
        sigma_hat  (array): estimated scale parameters, length d
        sigma_true (array): true scale parameters, length d

    Returns:
        rmse (float)
    """
    sigma_hat = np.asarray(sigma_hat, dtype=float)
    sigma_true = np.asarray(sigma_true, dtype=float)
    return np.sqrt(np.mean((np.log(sigma_hat) - np.log(sigma_true)) ** 2))


def tde_chi(chi_hat, chi_true):
    """
    Average pairwise absolute error for the upper tail-dependence coefficient chi.

    TDE = (2 / (d*(d-1))) * sum_{j<k} |chi_hat_jk - chi_true_jk|

    Args:
        chi_hat  (ndarray): d × d symmetric matrix of estimated chi values
                            (as returned by pairwise_chi_from_data)
        chi_true (ndarray): d × d symmetric matrix of true chi values

    Returns:
        tde (float)
    """
    chi_hat = np.asarray(chi_hat, dtype=float)
    chi_true = np.asarray(chi_true, dtype=float)
    d = chi_hat.shape[0]
    idx = np.triu_indices(d, k=1)
    return np.mean(np.abs(chi_hat[idx] - chi_true[idx]))


def rmse_chi(chi_hat, chi_true):
    """
    RMSE for the upper tail-dependence coefficient chi over all pairs.

    RMSE_chi = sqrt((2 / (d*(d-1))) * sum_{j<k} (chi_hat_jk - chi_true_jk)^2)
             = sqrt(mean((chi_hat - chi_true)^2))  over all pairs

    Args:
        chi_hat  (array): estimated pairwise chi values, length d*(d-1)/2
        chi_true (array): true pairwise chi values, same length

    Returns:
        rmse (float)
    """
    chi_hat = np.asarray(chi_hat)
    chi_true = np.asarray(chi_true)
    return np.sqrt(np.mean((chi_hat - chi_true) ** 2))


def ljpe(p_hat, p_true, eps=1e-12):
    """
    Log joint probability error: average absolute log-ratio of estimated
    to true exceedance probabilities.

    LJPE = (1/L) * sum_l |log(max(p_hat_l, eps)) - log(p_true_l)|

    Args:
        p_hat  (array): estimated exceedance probabilities, length L
        p_true (array): true exceedance probabilities, length L
        eps    (float): numerical floor applied to p_hat before log (default 1e-12)

    Returns:
        ljpe_val (float)
    """
    p_hat = np.asarray(p_hat, dtype=float)
    p_true = np.asarray(p_true, dtype=float)
    return np.mean(np.abs(np.log(np.maximum(p_hat, eps)) - np.log(p_true)))


def rmsle_p(p_hat, p_true, eps=1e-12):
    """
    Root mean squared log error for exceedance probabilities.

    RMSLE_p = sqrt((1/L) * sum_l (log(max(p_hat_l, eps)) - log(p_true_l))^2)

    Args:
        p_hat  (array): estimated exceedance probabilities, length L
        p_true (array): true exceedance probabilities, length L
        eps    (float): numerical floor applied to p_hat before log (default 1e-12)

    Returns:
        rmsle (float)
    """
    p_hat = np.asarray(p_hat, dtype=float)
    p_true = np.asarray(p_true, dtype=float)
    return np.sqrt(np.mean((np.log(np.maximum(p_hat, eps)) - np.log(p_true)) ** 2))


def empirical_survival(X, u, sign):
    """
    Calculate the empirical threshold exceedance probability 

    Args: 
        X (array): n*d-dimensional data matrix
        u (array): d-dimensional vector of the threshold  
        sign (string): a d-length string combination of ">" and "<",  which tells the way of threshold exceedance. For example,
        "><" means P(X_1 > u_1, X_2 < u_2)

    Returns:
        prob (float): empirical threshold exceedance probability
    """

    n_obs = X.shape[0]
    n_dim = X.shape[1]
    cond = np.array([True]*n_obs)
    sign = list(sign)
    for j in range(n_dim):
        s = sign.pop(0)
        if s == '>':
            cond &= X[:,j] > u[j]
        else:
            cond &= X[:,j] < u[j]
    # Count joint exceedances
    joint_exceedance = np.sum(cond)

    return joint_exceedance/n_obs