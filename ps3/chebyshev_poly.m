function [image, roots] = chebyshev_poly(order, x)
% [image, roots] = CHEBYSHEV_POLY(order, x) computes Chebyshev's polynomial 
% of given order at point x. 
% They are defined only on [-1,1]. It returns the image as the first output
% and the roots as a second output.

if length(x) == 1 && (x > 1 || x < -1)
    error('Point in the domain (second argument) is out of bound')
end

if mod(order, 1) ~= 0
    error('Order is not an integer')
end


image = cos(order * acos(x));


i = 1:order;
roots = -1*cos((2*i - 1)/(2*order) * pi);
end