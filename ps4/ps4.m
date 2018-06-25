% Metodos Numericos - EPGE/FGV 2018
% Instructor: Cezar Santos
% Problem Set 4 - Raul Guarini Riva

clc; close all; clear all;

%% Calibration
beta = 0.96;
q = 0.96;
%rho = input('Input the value for rho: ');
%gamma = input('Input the value for gamma: ');
%sigma = input('Input the value for sigma: ');

rho = 0.9;
sigma = 0.01;
sigma2 = sigma^2;
gamma = 1.0001;
m = 3;      % Tauchen 'bandwidth' parameter
mu = 0;

u = @(c) (c.^(1-gamma) - 1)/(1 - gamma);

nz = 9;
nk = 500;

%% Tauchen method
[zgrid, P] = tauchen_ar1(mu, rho, sigma2, nz, m);
disp('Tauchen discretization done');

%% Solving the individual household problem

% Grid for assets
upper_bound = 10;
a_grid = linspace(0, upper_bound, nk)';
% Grid for endowment
e_grid = exp(zgrid);

% Initial guess for value function and iteration parameters
V0 = repmat(a_grid, 1, nz);
max_it = 1000;
tol = 1e-3;
dist = tol +1;
it = 1;

% Creating the arrays for brute force optimization
[endownment, a, new_a] = meshgrid(e_grid, a_grid, a_grid);

% Computing consumption on each possible state and choice
C = max(endownment + a - q*new_a, 0);

% Computing utility on all possible states and choices
U = u(C);

% Value function iteration
disp('Starting value function iteration...')
tic
while it < max_it & dist > tol
    EV = V0 * P';
    EV = permute(repmat(EV, 1, 1, nk), [3, 2, 1]);
    H = U + beta*EV;
    [TV, index] = max(H, [], 3);
    it = it + 1;        % Updating iterations
    dist = max(max(abs(TV - V0)));       % Computing error
    V0 = TV;
    
    if rem(it, 100) == 0
        disp('Current iteration:')
        disp(it)
        disp('Current norm:')
        disp(dist)
    end
end

V = TV;
g = a_grid(index);
disp('Fixed point found!')
toc

%% Computing the invariant distribution

% Initial distribution (starting uniformly)
pi0 = ones(nk, nz);
pi0 = pi0/length(pi0(:));

% This indicator function shows if a given level of assets is the best
% choice under a given state. For example, if indicator(2,3,6) == 1, it
% means that under the second value for current assets and third value for
% endownment, the sixth value of assets is the best choice
indicator = zeros(nk, nz, nk);
for ik = 1:nk
    indicator(:,:, ik) = (g == a_grid(ik));
end

next_pi = zeros(size(pi0));     % Allocating memory

inv_it = 1;         % Invariant distribution iterations counter
max_inv_it = 1000;
inv_tol = 1e-4;     % Invariant distribution convergence parameter
inv_error = inv_tol + 1;

disp(' ')
disp('Looking for a stationary distribution over the state space...')
tic
while inv_error > inv_tol & inv_it < max_inv_it
    % Computing next period distribution
    for ik = 1:nk
        effective_measure = pi0 .* indicator(:,:,ik);
        next_pi(ik, :) = sum(effective_measure*P);
    end
    
    % Updating everything
    inv_error = norm(next_pi - pi0);        % Updating error
    inv_it = inv_it + 1;                    % Updating iterations
    pi0 = next_pi;                          % Updating the loop
end
disp('Invariant distribution found!')
toc
disp('Robustness check --> Displaying the sum of all elements (should be 1):')
sum_elements = sum(sum(next_pi));
sum_elements
