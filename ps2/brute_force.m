% Metodos Numericos - EPGE/FGV 2018
% Instructor: Cezar Santos
% [Vectorized] brute force implementation value function iteration

% The goal of this script is to provide a benchmark in terms of
% performance that will be useful in comparing other approaches against it.

clc; close all; clear all;
tic;            % Starting timer

%% Item 2 - Calibration and Steady State
beta = 0.987;
mu = 2;
alpha = 1/3;
delta = 0.012;
rho = 0.95;
sigma = 0.007;
m = 3;      % Tauchen's scale parameter

kss = ((alpha * beta)/(1 - beta*(1 - delta)))^(1/(1 - alpha));

% Defining functional forms for preferences and technology
u = @(c) (c.^(1-mu) - 1)/(1-mu);
f = @(z, k) exp(z).*k.^(alpha);

%% Item 3 - Value Function Iteration
nk = 500;
nz = 7;
kgrid = linspace(0.75*kss, 1.25*kss, nk)';      % Column vector!!

% I use Tauchen's method, already implemented in the previous problem set
% in order to discretize productivity shock (via tauchen_ar1 function!)
[zgrid, P] = tauchen_ar1(0, rho, sigma^2, nz, m);
zgrid = zgrid';         % Row vector!
disp('Tauchen discretization done.')

% Defining parameters for value function iteration:
% In the matrix form, rows are different values of capital and columns are
% different values for z

V0 = repmat(sqrt(kgrid), 1, nz);     % Concave and increasing guess
max_it = 1000;
tol = 1e-3;

%% Iteration
norm = tol + 1;
it = 1;
tic;
[K, Z, new_K] = meshgrid(kgrid, zgrid, kgrid);
K = permute(K, [2, 1, 3]);
Z = permute(Z, [2, 1, 3]);
new_K = permute(new_K, [2, 1, 3]);

% Computing consumption on each possible state and choice
C = max(f(Z,K) + (1-delta)*K - new_K,0);
% All possible utilities
U = u(C);

disp('Starting value function iteration through the good and old brute force...')
while it < max_it & norm > tol
    EV = V0 * P';
    EV = permute(repmat(EV, 1, 1, nk), [3, 2, 1]);
    H = U + beta*EV;
    [TV, index] = max(H, [], 3);
    it = it + 1;        % Updating iterations
    norm = max(max(abs(TV - V0)));       % Computing error
    V0 = TV;
    
    if rem(it, 100) == 0
        disp('Current iteration:')
        disp(it)
        disp('Current norm:')
        disp(norm)
    end
end

V = TV;
g = kgrid(index);
toc

%% Plotting
figure; hold on
for i = 1:nz
    plot(g(:,i), 'DisplayName', strcat('iz = ', num2str(i)))
end
title('Policy Function')
legend('show')
hold off
grid on

figure; hold on;
for i = 1:nz
    plot(V(:,i), 'DisplayName', strcat('iz = ', num2str(i)))
end
title('Value Function')
legend('show')
hold off
grid on

