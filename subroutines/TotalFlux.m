function [flux]=TotalFlux(H,ux,uy,delta)

% Kori-ULB
% Calculation of ice sheet flux using staggered grids

    ux1=circshift(ux,[0 1]); % ux(i,j-1)
    uy1=circshift(uy,[1 0]); % uy(i-1,j)
    flux=H.*sqrt((0.5*(ux+ux1)).^2+(0.5*(uy+uy1)).^2)*delta;
    
end


