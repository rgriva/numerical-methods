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
title('Policy Function (via VFinder\_Iterated)')
xlabel('Capital Stock')
legend('show')
hold off
grid on

subplot(2,1,1)
hold on
for i = 1:nz
    plot(kgrid, V_iterated(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title('Value Function (via VFinder\_Iterated)')
xlabel('Capital Stock')
legend('show')
hold off
grid on

%% Item 4 - Accelerator
% Percentage of times full optimization is done:
pct = 0.1;      % 10% of iterations!

disp(' ')
disp('Solving the functional equation with VFinder_Accelerated...')
tic
[V_accelerated, g_accelerated] = VFinder_Accelerated(u, f, delta, beta, V0, P, kgrid, zgrid, max_it, tol, pct);
toc
time_accelerated = toc;

%% Plotting
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

%% Item 5 - Multigrid
% These are the grid sizes that will be used in the multigrid scheme

disp(' ')
disp('To start multigrid method, press any key.')
pause;
disp(' ')
disp('Starting multigrid method...')
nk_multi = [100, 250, 500];

tic
for i = 1:length(nk_multi)
    if i == 1
        kgrid_multi = linspace(0.75*kss, 1.25*kss, nk_multi(i))';
        V0 = repmat(sqrt(kgrid_multi), 1, nz);
    end
    
    [V_multi, g_multi] = VFinder_Accelerated(u, f, delta, beta, V0, P, kgrid_multi, zgrid, max_it, tol, 0.1);
    
    % Defining new grid and redefine V0
    if i < length(nk_multi)
        next_kgrid_multi = linspace(0.75*kss, 1.25*kss, nk_multi(i+1))';
        V0 = zeros(nk_multi(i+1), nz);
        for iz = 1:nz
            % Linear interpolation is used along each fixed iz. This is
            % fast and the value function doesn't have so much curvature
            % near the steady state, as the previous plots has shown
            V0(:,iz) = interp1(kgrid_multi, V_multi(:,iz), next_kgrid_multi);
        end
        
        kgrid_multi = next_kgrid_multi;
    end
end
toc

 %% Plotting Results
set(0,'defaultAxesFontSize',16);
figure('position', [100,10,1000, 600]); 
subplot(2,1,2)
hold on
for i = 1:nz
    plot(kgrid, g_multi(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title('Policy Function (via VFinder\_Accelerated + Multigrid)')
xlabel('Capital Stock')
legend('show')
hold off
grid on

subplot(2,1,1)
hold on
for i = 1:nz
    plot(kgrid, V_iterated(:,i), 'DisplayName', strcat('iz ={ }', num2str(i)))
end
title('Value Function (via VFinder\_Accelerated + Multigrid)')
xlabel('Capital Stock')
legend('show')
hold off
grid on

%% Item 6 - Endogenous Grid

% I start using polyfit to numerically invert the "Cash-at-Hand" LHS that
% will appear when inverting the Euler Equation.

[Z, K_extended] = meshgrid(zgrid, linspace(0.5*kss, 1.5*kss, 10000));
m = exp(Z).*K_extended.^alpha + (1-delta)*K_extended; 

[Z, K_new] = meshgrid(zgrid, kgrid);        % each row is a value for K and each column a value for log(z)
K = K_new;      % Just for keeping a good notation

% Creating a polynomial representation
 deg = 4;        % polynomial degree
 fitter = zeros(nz, deg+1);         % This holds the estimates coefficients
 for iz = 1:nz
     fitter(iz, :) = polyfit(m(:, iz), K_extended(:,1), deg);
 end

% The fitter object can be understood as a numerical inverse. For a fixed
% z, I informed the data generated in m(k,z) as the independent variable to
% the fitter and kgrid as the dependent variable.

% ATTENTION: Each ROW of the fitter represents a given log(z) shock value.

u_marginal = @(c) c.^(-mu);
u_marginal_inverse = @(u) u.^(-1/mu);
pmg = @(K_new, Z) alpha*exp(Z).*K_new.^(alpha -1) + 1 - delta;

exo_grid = linspace(0.75*kss, 1.25*kss, nk)';
endo_grid = zeros(nk, nz);      % We will have a different endogenous grid for each value of z;
g_endogenous = zeros(nk, nz);    % Allocating memory for the policy function
dist = tol + 1;
it = 0;

c0 = repmat(kgrid, 1, nz)/max(kgrid);      % Increasing guess for consumption function
c = zeros(nk, nz);       % Allocating memory for comsumption policy function

while it < max_it && dist > tol
    E = u_marginal(c0).*pmg(K_new, Z) * P';
    rhs = K_new + u_marginal_inverse(beta * E);
    
    % Computing the endogenous grid and recovering policy function for
    % capital
    
    for iz = 1:nz
        endo_grid(:, iz) = polyval(fitter(iz,:), rhs(:, iz));
        g_endogenous(:, iz) = interp1(endo_grid(:, iz), kgrid, kgrid, 'linear');
    end
    g_endogenous = fillmissing(g_endogenous, 'linear');
    
    c = exp(Z).*K.^(alpha) + (1-delta)*K - g_endogenous;
    
    it = it+1;
    dist = norm(c - c0);
    c0 = c;
end

approx = norm(g_endogenous - g_iterated)

% Computing Euler errors


        








