% Metodos Numericos - EPGE/FGV 2018
% Instructor: Cezar Santos
% Problem Set 2 - Raul Guarini Riva

clc; close all; clear all;

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
disp(' ')   % Gimme some space!

% Defining parameters for value function iteration:
% In the matrix form, rows are different values of capital and columns are
% different values for z

V0 = repmat(sqrt(kgrid), 1, nz);     % Concave and increasing guess
max_it = 1000;
tol = 1e-3;

% Real iteration
disp('Solving the functional equation with VFinder_Iterated...')
tic     % Starting timer
[V_iterated, g_iterated] = VFinder_Iterated(u, f, delta, beta, V0, P, kgrid, zgrid, max_it, tol);
toc
time_iterated = toc;
%% Plotting Results
set(0,'defaultAxesFontSize',16);
figure('position', [100,10,1000, 600]); 
subplot(2,1,2)
hold on
for i = 1:nz
    plot(kgrid, g_iterated(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title('Policy Function')
xlabel('Capital Stock')
legend('show')
hold off
grid on

subplot(2,1,1)
hold on
for i = 1:nz
    plot(kgrid, V_iterated(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title('Value Function')
xlabel('Capital Stock')
legend('show')
hold off
grid on

%% Item 4
% Percentage of times full optimization is done:
pct = 1;      % 10% of iterations!


disp('Solving the functional equation with VFinder_Accelerated...')
tic
[V_accelerated, g_accelerated] = VFinder_Accelerated(u, f, delta, beta, V0, P, kgrid, zgrid, max_it, tol, pct);
toc
time_accelerated = toc;

%% Plot
% Comparing results
figure('position', [100,10,1100, 400]);
subplot(1,2,1)
hold on
for i = 1:nz
    plot(kgrid, V_iterated(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title('VFinder\_Iterated')
xlabel('Capital Stock')
legend('show')
grid on
hold off

subplot(1,2,2)
hold on
for i = 1:nz
    plot(kgrid, V_accelerated(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title('VFinder\_Accelerated')
xlabel('Capital Stock')
legend('show')
grid on
hold off