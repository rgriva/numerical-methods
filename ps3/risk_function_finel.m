function risk_vector = risk_function_finel(intercept, a_finel, K0, kgrid, n_elements, zgrid, state, P, alpha, mu, beta, delta)
% This computes the risk function! The output is a column vector of
% R(intercept, a_finel, K0), where K0 are collocation points.
% gamma: projection coefficients
% K0: capital points
% kgrid: grid points for capital
% zgrid: grid points for log(TFP)
% state: current log(TFP) state index
% P: transition matrix
% alpha, mu, beta, delta: parameters
% The risk-vector is just the Euler Equation.

C0 = C_proj_finel(intercept, a_finel, K0, kgrid, n_elements);
K1 = (exp(zgrid(state))* K0(:).^alpha + (1 - delta)*K0(:) - C0);
C1 = C_proj_finel(intercept, a_finel, K1, kgrid, n_elements);

% Conditional measure from the Markovian Process
measure = P(state,:);

% TMS is column vector
tms = (C1./C0).^(-mu);

% Marginal productivity is a nk x nz matrix
pmg = (1 - delta + alpha * exp(zgrid).* K1.^(alpha -1));

risk_vector = beta*(tms .* pmg)*measure' - 1;
end