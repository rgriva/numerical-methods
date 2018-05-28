function proj = C_proj(gamma, K, kgrid)
% proj = C_PROJ(gamma, K) computes the projection of the comsumption
% function on a point K of kgrid using the vector gamma of projection coefficients.
% The order of projection is the size of gamma (for order d, choose a d+1
% dimensional gamma)

point = (K-min(kgrid))/(max(kgrid) - min(kgrid));
T = zeros(length(gamma),1);
T(1) = 1;
T(2) = point;

for i = 3:length(gamma)
    T(i) = 2*point*T(i-1) - T(i-2);
end

proj = gamma'*T;

