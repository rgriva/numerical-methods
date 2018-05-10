function [V, g] = VFinder_Accelerated(u, f, delta, beta, V0, P, kgrid, zgrid, max_it, tol, pct)
% [V, g] = VFINDER_ACCELERATED(u, f, delta, beta, V0, P, kgrid, zgrid, max_it, tol, pct) finds the Value Function 
% and Policy Function of the standard Social
% Planner's problem. It exploits concavity of V (due to the concavity of
% the utility function) and monotonicity. It uses the value function
% iteration method. It also implements the Accelerator method. For the
% Standard approach, see VFINDER_ITERATED.
% Parameters are:
% u: function handle for utility
% f: function handle for technology
% delta: depreciation
% beta: discounting factor
% V0: initial guess for iteration
% P: transition matrix for the Markovian Shocks
% kgrid: grid values for capital
% zgrid: grid values for z shock (log of TFP)
% max_it = maximum number of iterations
% tol = convergence criterium
% pct = percentage of iterations in which we perform some sort of
% optimization