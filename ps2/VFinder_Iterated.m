function [V, g] = VFinder_Iterated(u, f, delta, beta, V0, P, kgrid, zgrid, max_it, tol)
% [V, g] = VFINDER_ITERATED(u, f, delta, beta, V0, P, kgrid, zgrid, max_it, tol) finds the Value Function 
% and Policy Function of the standard Social
% Planner's problem. It exploits concavity of V (due to the concavity of
% the utility function) and monotonicity. It uses the value function
% iteration method.
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

nk = length(kgrid);
nz = length(zgrid);

g = repmat(kgrid, 1, nz);
norm = tol + 1;
it = 1;

% Trying to vectorize calculations a little bit
% Every row is for a given K, every column for a given Z and every
% 'diagonal' for a new_k, which is what we choose
[K, Z, new_K] = meshgrid(kgrid, zgrid, kgrid);
K = permute(K, [2, 1, 3]);
Z = permute(Z, [2, 1, 3]);
new_K = permute(new_K, [2, 1, 3]);

% Computing consumption on each possible state and choice
C = max(f(Z,K) + (1-delta)*K - new_K,0);
% All possible utilities
U = u(C);

% Starting iterations

while it < max_it & norm > tol
    % Procede with another iteration
    EV = V0 * P';
    EV = permute(repmat(EV, 1, 1, nk), [3, 2, 1]);
    if it == 1;
    % Brute force for the first time
        H = U + beta*EV;
        [TV, index] = max(H, [], 3);
        it = it +1;
    
    % For other iterations than the first, we exploit monotonicity and concavity
    else        
        for iz = 1:nz
            for ik = 1:nk
                pointer = index(max(ik-1, 1),iz);
                H = U(ik, iz, pointer) + beta*EV(ik, iz, pointer);
                H_next = U(ik, iz, pointer+1) + beta*EV(ik, iz, pointer+1);
                while H_next > H & pointer < (nk-1)
                    pointer = pointer + 1;
                    H = H_next;
                    H_next = U(ik, iz, pointer+1) + beta*EV(ik, iz, pointer+1);
                end
                index(ik, iz) = pointer;
                TV(ik, iz) = U(ik, iz, pointer) + beta*EV(ik, iz, pointer);
            end
        end
    end
    
    norm = max(max(abs(TV - V0)));       % Computing error
    V0 = TV;        % Updating the loop
    it = it + 1;
    
    if rem(it, 100) == 0
        disp('Current iteration:')
        disp(it)
        disp('Current norm:')
        disp(norm)
    end
end

g = kgrid(index);
V = TV;
disp('Fixed point found!!!')
            
            
        

