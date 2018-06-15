function proj = C_proj_finel(intercept, a_finel, K, kgrid, n_elements)
% proj = C_PROJ_FINEL(intercept, a_finel, K, kgrid, n_elements) 
% computes the projection of the comsumption 
% function on a point K using the vector a of projection coefficients and
% intercepts, each with the size as the number os elements.
% It uses the companion function FIND_ELEMENT

% if length(a_finel) ~= n_elements
%     error('Parameter lenght is different from the number of finite elements.')
% end

element_index = find_element(K, kgrid, n_elements);
proj = intercept(element_index) + a_finel(element_index).*K(:);
end
            



