function [m]=MscaledPoly(x,p)

% Kori-ULB
% Polynomial function for Plume melt parameterization underneath ice
% shelves

    m=zeros(size(x));
    for i=1:12
        m=m+p(i)*x.^(i-1);
    end
    
end


