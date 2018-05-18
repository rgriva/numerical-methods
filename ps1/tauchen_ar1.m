function [grid, P] = tauchen_ar1(mu, rho, sigma2, N, m)
% This function implements Tauchen (1986) discretization method. It takes as
% arguments the average mu of an AR1 process, error variance sigma2,
% persistence coefficient rho, the number N of points in a grid and the scaling
% parameter m.
% Outputs are the grid (state values) and the corresponding transition matrix.

sigma = sqrt(sigma2);

theta_N = m * sigma/sqrt((1 - rho^2)) ;
theta_1 = - m * sigma/((1 - rho^2)^0.5);
delta = (theta_N - theta_1)/(N -1);

% Shifting the grid to account for mu
grid = linspace(theta_1, theta_N, N) + mu;

% I use a vector trick to compute P faster
[J, I] = meshgrid(grid, grid);

u_bound = (1/sigma) * (J + delta/2 - (1-rho)*mu - rho*I);
l_bound = (1/sigma) * (J - delta/2 - (1-rho)*mu - rho*I);

P = normcdf(u_bound) - normcdf(l_bound);

P(:,1) = normcdf((1/sigma) * (grid(1) - (1 - rho)*mu - rho*I(:,1) + delta/2));
P(:,N) = 1 - normcdf((1/sigma) * (grid(N) - (1 - rho)*mu - rho*I(:,1) - delta/2));

grid = grid';
end
