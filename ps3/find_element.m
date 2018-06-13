function element_index = find_element(K, kgrid, n_elements)
% FIND_ELEMENT(K, kgrid, n_elements) divides the space of kgrid in
% n_elements and throws back the correspoding index of the element that K
% belongs to. This also works with vectors.
% It uses the 'Divide and Conquer' approach due to performance concerns.

endpoint = kgrid(end);
startpoint = kgrid(1);
element_grid = linspace(startpoint, endpoint, n_elements + 1);
element_index = zeros(size(K));

for i = 1:length(K)
    % Divide and Conquer
    lower_bound_index = 1;
    upper_bound_index = n_elements + 1;
    current_K = K(i);
    guess_index = round((n_elements + 1)/2);
    guess = element_grid(guess_index);
    while (upper_bound_index - lower_bound_index) > 1
        if current_K >= guess
            % update lower bound
            lower_bound_index = guess_index;
            % update guess
            guess_index = round((lower_bound_index + upper_bound_index)/2);
            guess = element_grid(guess_index);
        else
            % update upper bound
            upper_bound_index = guess_index;
            % update guess
            guess_index = round((lower_bound_index + upper_bound_index)/2);
            guess = element_grid(guess_index);
        end
    end
    
    element_index(i) = guess_index-1;
end
