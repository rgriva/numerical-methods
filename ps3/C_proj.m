function proj = C_proj(gamma, K, kgrid)
% proj = C_PROJ(gamma, K) computes the projection of the comsumption
% function on a point K of kgrid using the vector gamma of projection coefficients.
% The order of projection is the size of gamma (for order d, choose a d+1
% dimensional gamma)

a = 2/(max(kgrid) - min(kgrid));
b = -1*((max(kgrid) + min(kgrid)) / (max(kgrid) - min(kgrid)));

% Remapping into -1 and 1
point = a*K + b;

T = zeros(length(gamma), length(K));
for i=1:length(gamma)
    T(i,:) = chebyshev_poly(i-1, point);
end

proj = gamma'*T;

%End function
end
