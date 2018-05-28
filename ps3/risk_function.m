function risk_vector = risk(gamma, K0, kgrid, alpha, mu, beta, delta)
% This computes the risk function! The output is a column vector of
% R(gamma, K0), where K0 are collocation points.

C0 = C_proj(gamma, K0, kgrid);
K1 = K0.^alpha + (1 - delta)*K0 - C0;
C1 = C_proj(gamma, K1, kgrid);

tms = (C1./C0).^(-mu);

risk_vector = beta * tms .* (1 - delta + alpha * K1.^(alpha-1));
risk_vector = risk_vector';

% Function end
end