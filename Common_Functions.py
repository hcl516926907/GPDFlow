import numpy as np
import torch

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


def pairwise_chi_empirical(dim, model, n_monte_carlo, n_experiments):
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

    probs = np.linspace(0.80, 0.99, 100)
    # Compute lambda_u for each quantile
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


def marginal_parameter_monte_carlo(dim):
    """
    Load the estimated sigma and gamma of GPDFlow models in a simulation scenario, 
    and save them in two dictionaries.
    
    Args: 
        dim (integer): dimension of the model.
    
    Returns:
        sigma_dict(dict)， gamma_dict(dict): two dictionaries contains the estimated sigma and gamma of 
        GPDFlow models in a simulation scenario. 
    
    """
    sigma_dict = {f'{i+1}': [] for i in range(dim)}
    gamma_dict = {f'{i+1}': [] for i in range(dim)}

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