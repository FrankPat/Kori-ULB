function [alpha_local] = ComputeAlphaLocal(HB, shMASK, ctr)

% Kori-ULB
% Compute local basal slope following Favier et al., 2019 (Appendix B)
% 
% Inputs:
%   HB       - Ice draft depth (negative downwards) 
%   shMASK   - Ice shelf mask indicating shelf regions 
%
% Output:
%   alpha_local - Local slope angle (radians)

    HB(shMASK==0)=NaN;

    % Compute slope components using a central difference approach
    slope_x = CheckSlope1D(HB, ctr, 1); % Compute along x-direction
    slope_y = CheckSlope1D(HB, ctr, 2); % Compute along y-direction

    % Compute local slope magnitude
    alpha_local = atan(sqrt(slope_x.^2 + slope_y.^2));

end


