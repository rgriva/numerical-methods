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
upper_bound = 100;
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
disp('Fixed point found')
toc
