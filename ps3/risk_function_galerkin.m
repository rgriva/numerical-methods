function risk_vector = risk_function_galerkin(a_galerkin, n_galerkin, kgrid, zgrid, state, P, alpha, mu, beta, delta)
% This computes the risk function! The output is a column vector of
% R(gamma, K0), where K0 are collocation points.
% gamma: projection coefficients
% K0: capital points
% kgrid: grid points for capital
% zgrid: grid points for log(TFP)
% state: current log(TFP) state index
% P: transition matrix
% alpha, mu, beta, delta: parameters

intercept = zeros(n_galerkin, 1);
risk_vector = ones(n_galerkin,1);

element_grid = linspace(kgrid(1), kgrid(end), n_galerkin + 1);
for i=1:length(risk_vector)
    % Integral limits
    lower_bound = element_grid(i);
    upper_bound = element_grid(i+1);
    
    % Defining the function to be integrated
    inner_handle = @(k) (risk_function_finel(intercept, a_galerkin, k, kgrid, n_galerkin, zgrid, state, P, alpha, mu, beta, delta).*a_galerkin(i).*k(:))';
    
    % Integrating using 'integral' function
    risk_vector(i) = integral(inner_handle, lower_bound, upper_bound);
end

% Function end
end
    
    


