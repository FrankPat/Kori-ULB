function Y = heaviside(X)

%HEAVISIDE    Step function.
%    HEAVISIDE(X) is 0 for X < 0 and 1 for X > 0.
%    The value HEAVISIDE(0) is 0.5 by default. function Y = heaviside(X)
% HEAVISIDE Step function

    Y = zeros(size(X))+0.5;
    Y(X < 0) = 0;
    Y(X > 0) = 1;
    
end


