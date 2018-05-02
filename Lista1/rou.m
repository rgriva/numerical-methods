function [grid, P] = rou(mu, rho, sigma2, N)
% This function implements Rouwenhorst discretization method. It takes as
% arguments the average mu of an AR1 process, error variance sigma2,
% persistence coefficient rho, the number of points in a grid N.
% Outputs are the grid (state values) and the corresponding transition matrix.

theta_N = sqrt(sigma2/(1 - rho^2))* sqrt(N - 1);
theta_1 = - theta_N;

grid = linspace(theta_1, theta_N, N) + mu;

p = (1 + rho)/2;

% Initializing the transition matrix
P_init = [p, 1-p;
          1-p, p];
      
% Looping over N points (I couldn't figure out a why to vectorize this)
for n = 3:N
    P_new = p * [P_init, zeros(n - 1, 1); zeros(n - 1, 1)', 0] + ...
            (1 - p)*[zeros(n - 1, 1), P_init; 0, zeros(n - 1, 1)'] + ...
            (1 - p)*[zeros(n - 1, 1)', 0; P_init, zeros(n - 1, 1)] + ...
            p * [0, zeros(n - 1, 1)'; zeros(n - 1, 1), P_init];
        
     % Updating variables
     P_init = P_new;
end

% Normalizing conditional measures
for n = 1:N
    P_new(n, :) = P_new(n, :)/sum(P_new(n, :));
end

P = P_new;
grid = grid';
    
