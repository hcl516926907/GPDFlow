import torch
import torch.nn as nn

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

    def __init__(self, dim, device):
        super().__init__()
        self.dim = dim
        # Store log_sigma so sigma = exp(log_sigma) > 0
        self.log_sigma = nn.Parameter(torch.zeros(dim, device=device))
        
        self.theta = nn.Parameter(torch.zeros(dim, device=device))
        
    def get_sigma(self):
        sigma = torch.exp(self.log_sigma)
        sigma = torch.clamp(sigma, min = 1e-6, max = 1e6)
        return sigma
    
    def get_gamma(self):
        gamma = -0.5 + torch.sigmoid(self.theta)  + 1e-6 
        gamma = gamma.sign() * torch.max(gamma.abs(), torch.tensor(1e-6, device=gamma.device))
        return gamma
        
    def forward_transform(self, y):
        """
        y -> x
        """
#         y = torch.clamp(y, min= -1e6, max = 1e6)
        
        sigma = self.get_sigma()
        gamma = self.get_gamma()
        
        # We'll do it elementwise, but in a vectorized way.
        # y and x are shape (batch_size, dim).
    
        x_out = torch.zeros_like(y)

        # For gamma_j = 0 => x_j = sigma_j * y_j
        # For gamma_j != 0 => x_j = sigma_j*(exp(gamma_j*y_j) - 1)/gamma_j
        # We can use torch.where(...) to handle each dimension separately.

        # Expand so shapes match for broadcasting
        sigma = sigma.unsqueeze(0)  # (1, dim)
        gamma = gamma.unsqueeze(0)  # (1, dim)

        # Exponential case: x_j = sigma_j*(exp(gamma_j*y_j) - 1)/gamma_j
        x_out = sigma * (torch.exp(torch.clamp(gamma * y, max=10)) - 1.0) / gamma
        
#         x_out = torch.clamp(x_out, min= -1e10, max = 1e10)
        return x_out

    def inverse_transform(self, x):
        """
        x -> y
        """
#         x = torch.clamp(x, min= -1e10, max = 1e10)
        
        sigma = self.get_sigma()
        gamma = self.get_gamma()

        y_out = torch.zeros_like(x)

        # For gamma_j = 0 => y_j = x_j / sigma_j
        # For gamma_j != 0 => y_j = (1/gamma_j)* log(1 + gamma_j*x_j / sigma_j)
        sigma = sigma.unsqueeze(0)  # (1, dim)
        gamma = gamma.unsqueeze(0)  # (1, dim)

#         gamma_is_zero = (gamma.abs() < 1e-6)

#         # Linear case: y_j = x_j / sigma_j
#         y_linear = x / sigma

        # Exponential (inverse) case: y_j = (1/gamma_j)* log(1 + gamma_j*x_j / sigma_j)
        inside = 1.0 + (gamma * x) / sigma
        inside = torch.clamp(inside, min=1e-12) 
        y_out = (1.0 / gamma) * torch.log(inside)

#         y_out = torch.clamp(y_out, min= -1e6, max = 1e6)
        return y_out

    def forward(self, z, reverse=False):
        """
        For normflows compatibility:
          - If reverse=False: we do y->x
          - If reverse=True: we do x->y
        """
        if not reverse:
            return self.forward_transform(z)
        else:
            return self.inverse_transform(z)



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
    def __init__(self, dim, flow, device, s_min , s_max, num_integration_points, penalty_lambda):
        super().__init__()
        self.dim = dim
        # Learnable data transformation
        self.data_transform = DataTransform(dim, device)
        # RealNVP flow model
        self.flow_model = flow
        self.device = device
        self.s_values = torch.linspace(s_min, s_max, num_integration_points, device=device)
        self.s_values = self.s_values.reshape(-1, 1, 1)
        self.s_min = s_min
        self.s_max = s_max
        self.num_integration_points = num_integration_points
        self.penalty_lambda = penalty_lambda
        
    def log_integral_f_T(self, data):
            # Expand batch_x to match s_values
        batch_size = data.shape[0]
        dim = data.shape[1]
        
        x_expanded = data.unsqueeze(0)  # Shape (1, effective_batch_size, dim)
        x_expanded = x_expanded.expand(self.num_integration_points, -1, -1)  # Shape (num_points, batch_size, dim)

        # Expand s_values to match batch size and dimension
        s_expanded = self.s_values.expand(-1, batch_size, 1)  # Shape (num_points, batch_size, 1)

        # Compute x + s for all s_values
        x_plus_s = x_expanded + s_expanded  # Broadcasting over the last dimension
        x_plus_s = x_plus_s.reshape(-1, dim)  # Flatten to (num_points * batch_size, dim)

        log_f_T = self.flow_model.log_prob(x_plus_s)  # Shape (num_points * batch_size,)
        log_integrand = log_f_T.reshape(self.num_integration_points, batch_size)
        
        # 1. log-sum-exp over 'num_points' dimension
        max_vals, _ = torch.max(log_integrand, dim=0, keepdim=True)  # shape (1, batch_size)
        stable_exp = torch.exp(log_integrand - max_vals)             # shape (num_points, batch_size)
        
        delta_s = (self.s_max - self.s_min) / (self.num_integration_points - 1)
        sum_exp = torch.trapz(stable_exp, dx=delta_s, dim=0)         # shape (batch_size,)
        # 2. Now put the max_vals back in:
        log_integral = max_vals.squeeze(0) + torch.log(sum_exp + 1e-40)  # shape (batch_size,)
        
        return log_integral
    
    def log_prob_T_mGPD_std(self, data):
        log_integral = self.log_integral_f_T(data)
        max_T = torch.max(data, dim=1)[0]
        
        log_prob = log_integral - max_T
        return log_prob
    
    
        

    def forward(self, x_data):
        """
        Return minus log p_data(x_data).
        We do:
          y = T_inv(x_data)  [since T_inv is x->y]
          log p_data(x_data) = log p_flow(y) + log|det(d y / d x)|
        """
        # 1) x->y
        y = self.data_transform.inverse_transform(x_data)

        # 2) log prob under flow
        log_prob_y = self.log_prob_T_mGPD_std(y)

        # 3) log | det dy/dx
        sigma = self.data_transform.get_sigma()
        gamma = self.data_transform.get_gamma()

        sigma = sigma.unsqueeze(0)  # (1, dim)
        gamma = gamma.unsqueeze(0)  # (1, dim)

        # Build inside = sigma + gamma*x
        inside = sigma + gamma * x_data

        # log_abs_detJ per sample = - sum_j log(inside_j)
        log_abs_detJ = -torch.sum(torch.log(inside.abs() + 1e-12), dim=1)
        
        # => log p(x_data) = log p(y) + log|det dy/dx|
        log_prob_x = log_prob_y + log_abs_detJ
        
        ######## Soft penalty: penalize negative (sigma_j + gamma_j*x_j)
        negative_part = torch.relu(-inside)  # = max(0, -(inside)) => positive if inside<0
        # penalty = sum_j [ negative_part^2 ] over j, average over batch
        penalty_per_sample = (negative_part ** 2).sum(dim=1)
        penalty = self.penalty_lambda * penalty_per_sample.mean()
        
        return -log_prob_x.mean() +  penalty 

    def sample(self, n_samples=1):
        """
        Sample from the learned distribution in x-space by:
          1) Sample y ~ flow
          2) Convert y->x using forward_transform
        """
        self.flow_model.eval()
        samples_T,_ = self.flow_model.sample(n_samples)
        self.flow_model.train()
        
        samples_T_max = torch.max(samples_T,axis=1,keepdim=True)[0]
        samples_T_1 = samples_T - samples_T_max
        
        samples_E  = torch.empty(n_samples, device=self.device)
        samples_E = samples_E.exponential_(1.0).unsqueeze(1)
        
        samples_y = samples_E + samples_T_1

        samples_x = self.data_transform.forward_transform(samples_y)  # y->x
        return samples_x, samples_y, samples_T
