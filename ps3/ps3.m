% Metodos Numericos - EPGE/FGV 2018
% Instructor: Cezar Santos
% Problem Set 3 - Raul Guarini Riva

clc; close all; clear all;

%% Calibration and Steady State

% Same calibration as the previous problem set
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

nk = 500;
nz = 7;
kgrid = linspace(0.75*kss, 1.25*kss, nk)';      % Column vector!!

% I use Tauchen's method, already implemented in the previous problem set
% in order to discretize productivity shock (via tauchen_ar1 function!)
[zgrid, P] = tauchen_ar1(0, rho, sigma^2, nz, m);
zgrid = zgrid';         % Row vector!
disp('Tauchen discretization done.')
disp(' ')   % Gimme some space!

%% Problem 1 - Chebyshev Polynomials

tic;
% I will follow the steps from the slides
d = 6;         % Order of target polynomial

a = 2/(max(kgrid) - min(kgrid));
b = -1*((max(kgrid) + min(kgrid)) / (max(kgrid) - min(kgrid)));

C = zeros(nk, nz);

% The implementation cannot use d=7 directly. On my tests, these was
% unreliable. It would generate strange beheavior (negative consumptio). So
% I procede in a "multigrid" for gamma.

disp('Starting the projection on Chebyshev Polynomials...')
for iz = 1:nz
    gamma_optimal = ones(2,1);      % Just to start the loop
    state = iz;
    
    % Compute gamma_optimal for a fixed value of current z
    for id = 2:d
        
        % Finding roots and reverting them back to original scale
        [~, roots] = chebyshev_poly(id+1, 0);
        K0 = ((roots - b)/a)';
        
        % Function handle
        R = @(gamma) risk_function(gamma, K0, kgrid, zgrid, state, P, alpha, mu, beta, delta);
        
        % Update initial condition
        gamma0 = zeros(id+1, 1);
        gamma0(1:id) = gamma_optimal;
        
        % Solving the system
        options = optimset('Display','off');     % Turning off dialogs
        gamma_optimal = fsolve(R, gamma0, options);
    end
    
    % Interpolated policy function for consumption
    C(:, iz) = C_proj(gamma_optimal, kgrid, kgrid);
end
disp('Spectral projection done.')
toc

disp(' ')

% Computing the policy function for capital
g = exp(zgrid).*kgrid.^alpha + (1-delta)*kgrid - C;

% Computing Euler Errors
u_marginal = @(c) c.^(-mu);
u_marginal_inverse = @(u) u.^(-1/mu);
pmg = @(K_new, Z) alpha*exp(Z).*K_new.^(alpha - 1) + 1 - delta;

% Computing consumption tomorrow interpolating the policy function for
% comsumption found previously

next_C = zeros(nk, nz);

for iz = 1:nz
    next_C(:, iz) = interp1(kgrid, C(:, iz), g(:, iz));
end

E = u_marginal(next_C).* pmg(g, zgrid) * P';
EE = log10(abs(1 - u_marginal_inverse(beta*E)./C));

%% Plotting Item 1
disp('Press any key to plot results from Item 1')
pause;

set(0,'defaultAxesFontSize',16);
figure('position', [100,10,900, 1400]); 
subplot(3,1,1)
hold on
for i = 1:nz
    plot(kgrid, g(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title(sprintf('Policy Function for Capital Stock (using %d Chebyshev Polynomials)', d+1))
xlabel('Capital Stock')
legend('show', 'Location', 'southeast')
hold off
grid on

subplot(3,1,2)
hold on
for i = 1:nz
    plot(kgrid, C(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title(sprintf('Policy Function for Consumption (using %d Chebyshev Polynomials)', d+1))
xlabel('Capital Stock')
legend('show', 'Location', 'southeast')
hold off
grid on

% Plotting Euler Errors

subplot(3,1,3)
hold on
for i = 1:nz
    plot(kgrid, EE(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title(sprintf('Euler Errors (using %d Chebyshev Polynomials)', d+1))
xlabel('Capital Stock')
legend('show', 'Location', 'southeast')
hold off
grid on

%% Problem 2 - Finite Elements + Collocation
tic
C_finel = zeros(nk,nz);
n_elements = 50;

disp('Starting the projection on Finite Elements with Collocation...')
for iz = 1:nz
    state = iz;
    K0 = linspace(kgrid(1), kgrid(end), 2*n_elements)';     % This should be a column vector!
    
    % Initial condition
    a_finel0 = ones(n_elements, 1)/(n_elements);
    intercept = ones(n_elements, 1);
    
    % Joining parameters
    params0 = [intercept; a_finel0];
    
    % Solving the system
    options = optimset('Display','off');     % Turning off dialogs
    R_finel = @(params) risk_function_finel(params(1:n_elements), params(n_elements+1 :end), K0, kgrid, n_elements, zgrid, state, P, alpha, mu, beta, delta);
    params_optimal = fsolve(R_finel, params0, options);
    
    % Separating again
    a_finel_optimal = params_optimal(n_elements+1:end);
    intercept_optimal = params_optimal(1:n_elements);
    
    % Computing projection
    C_finel(:,iz) = C_proj_finel(intercept_optimal, a_finel_optimal, kgrid, kgrid, n_elements);
end
disp('Finite Elements projection done.')
toc

% Computing the policy function for capital
g_finel = exp(zgrid).*kgrid.^alpha + (1-delta)*kgrid - C_finel;

% Computing Euler Errors
next_C_finel = zeros(nk, nz);

for iz = 1:nz
    next_C_finel(:, iz) = interp1(kgrid, C_finel(:, iz), g_finel(:, iz));
end

E = u_marginal(next_C_finel).* pmg(g_finel, zgrid) * P';
EE_finel = log10(abs(1 - u_marginal_inverse(beta*E)./C_finel));

%% Plotting Results
disp(' ')
disp('Press any key to plot results from Item 2 - Collocation Points')
pause;

set(0,'defaultAxesFontSize',16);
figure('position', [100,10,900, 1400]); 
subplot(3,1,1)
hold on
for i = 1:nz
    plot(kgrid, g_finel(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title(sprintf('Policy Function for Capital Stock (using Collocation and %d FE)', n_elements))
xlabel('Capital Stock')
legend('show', 'Location', 'Best')
hold off
grid on

subplot(3,1,2)
hold on
for i = 1:nz
    plot(kgrid, C_finel(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title(sprintf('Policy Function for Consumption (using Collocation and %d FE)', n_elements))
xlabel('Capital Stock')
legend('show', 'Location', 'Best')
hold off
grid on

% Plotting Euler Errors
subplot(3,1,3)
hold on
for i = 1:nz
    plot(kgrid, EE_finel(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title(sprintf('Euler Errors (using Collocation and %d FE)', n_elements))
xlabel('Capital Stock')
legend('show', 'Location', 'Best')
hold off
grid on

figure
hold on
for i = 1:nz
    plot(kgrid, C(:,i) - C_finel(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title(sprintf('Difference of Consumption Functions (%d FE vs %d Chebyshev Polynomials)', [n_elements, d+1]))
xlabel('Capital Stock')
legend('show', 'Location', 'Best')
hold off
grid on

%% Problem 2 - Finite Elements + Galerkin
n_galerkin = 15;
intercept = zeros(n_galerkin, 1);
C_galerkin = zeros(nk, nz);
%intercept = zeros(n_galerkin, 1);
tic
for iz = 1:nz
    state = iz;
    
    % Initial condition
    a_galerkin0 = ones(n_galerkin, 1)/(n_galerkin);
    
    % Computing the optimal parameters
    options = optimset('Display','off');     % Turning off dialogs
    R_galerkin = @(a_galerkin) risk_function_galerkin(a_galerkin, n_galerkin, kgrid, zgrid, state, P, alpha, mu, beta, delta);
    [a_galerkin_optimal, value] = fsolve(R_galerkin, a_galerkin0, options);
    
    C_galerkin(:, iz) = C_proj_finel(intercept, a_galerkin_optimal, kgrid, kgrid, n_galerkin);
end
toc
plot(kgrid, C_galerkin)
