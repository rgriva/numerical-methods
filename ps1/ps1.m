% Metodos Numericos - EPGE/FGV 2018
% Instructor: Cezar Santos
% Problem Set 1 - Raul Guarini Riva
clc; close all; clear all;
tic;            % Starting timer

%% Item 1
% tauchen_ar1.m is a companion function that discretizes an AR(1) process in
% the form 
% [theta_{t + 1} = (1 - rho)*mu + rho * theta_{t} + e_{t+1}]. 
% The %error term e is a Gaussian iid process with mean zero and variance
% sigma2.

% Calibration
mu = 0;
rho = 0.99;
sigma = 0.007;
sigma2 = sigma^2;
N = 9;
m = 3;      % Scaling parameter

[grid_tauchen, P_tauchen] = tauchen_ar1(mu, rho, sigma2, N, m);
disp('Tauchen discretization done.')

%% Item 2
% The same previous comment applies here for the Rouwenhorst method.
[grid_rou, P_rou] = rou(mu, rho, sigma2, N);
disp('Rouwenhorst discretization done.')

%% Item 3
burn = 100;
periods = 10000;

% Drawing the gaussian shock sequence
shocks = normrnd(0, sigma, [periods + burn, 1]);

% Simulating continuous process
% theta_0 = mu;
cont_path = zeros(periods + burn, 1);      % Just to allocate memory
cont_path(1) = mu;

for t=2:(periods+burn)
    cont_path(t) = (1-rho)*mu + rho*cont_path(t-1) + shocks(t);
end

% Deleting the burn periods, so initial condition is stochastic
cont_path = cont_path(burn+1:end);
disp('Continuous path simulated.')

% Mapping the normal shocks on uniform ones
uni_shocks = normcdf(shocks, 0, sigma);

% Simulating the Tauchen path
tauchen_path = zeros(periods+burn, 1);      % Allocating memory again
tauchen_path(1) = grid_tauchen(floor(N/2));     % Starting "in the middle" of the grid
state = floor(N/2);     % This is the index of the states, with N possible values

for t = 2:(periods + burn)
    cum_sum = cumsum(P_tauchen(state, :));      % Cumulative sum of each adequate line
    ind = sum(uni_shocks(t) >= cum_sum);
    state = ind + 1;        % Update the state
    tauchen_path(t) = grid_tauchen(state);
end

tauchen_path = tauchen_path(burn+1:end);
disp('Tauchen path simulated.')

% Simulating Rouwenhorst path
rou_path = zeros(periods+burn, 1);      % Allocating memory again
rou_path(1) = grid_rou(floor(N/2));     % Starting "in the middle" of the grid
state = floor(N/2);     % This is the index of the states, with N possible values

for t = 2:(periods + burn)
    cum_sum = cumsum(P_rou(state, :));      % Cumulative sum of each adequate line
    ind = sum(uni_shocks(t) >= cum_sum);
    state = ind + 1;        % Update the state
    rou_path(t) = grid_rou(state);
end

rou_path = rou_path(burn+1:end);
disp('Rouwenhorst path simulated.')



% Computing MSE
MSE_tauchen = sum((tauchen_path - cont_path).^2)/periods;
MSE_rou = sum((rou_path - cont_path).^2)/periods;
disp(' ')
disp('Mean Square Error of Tauchen Discretization:')
disp(MSE_tauchen)
disp('Mean Square Error of Rouwenhorst Discretization:')
disp(MSE_rou)

%% Item 4
% Estimating using Tauchen
model_tauchen = polyfit(tauchen_path(1:end -1), tauchen_path(2:end), 1);
rho_hat_tauchen = model_tauchen(1);
mu_hat_tauchen = model_tauchen(2)/(1 - rho_hat_tauchen);
disp(' ')
disp('Estimated rho (Tauchen):')
disp(rho_hat_tauchen)

% Estimating using Rouwenhorst
model_rou = polyfit(rou_path(1:end -1), rou_path(2:end), 1);
rho_hat_rou = model_rou(1);
mu_hat_rou = model_rou(2)/(1 - rho_hat_rou);
disp('Estimated rho (Rouwenhorst):')
disp(rho_hat_rou)

toc

%% Plotting stuff
%  Plotting stuff
figure;

subplot(2,1,1)
plot(cont_path, 'blue')
hold on
plot(tauchen_path, 'red')
legend('Continuous', 'Tauchen')
title(sprintf('Tauchen Process vs Continuous Process (N = %d, rho = %.2f)', N, rho), 'FontSize', 14)
xlabel('Time', 'FontSize', 14)
grid on
hold off

subplot(2,1,2)
plot(cont_path, 'blue')
hold on
plot(rou_path, 'green')
legend('Continuous', 'Rouwenhorst')
title(sprintf('Rouwenhorst Process vs Continuous Process (N = %d, rho = %.2f)', N, rho), 'FontSize', 14)
xlabel('Time', 'FontSize', 14)
grid on
hold off