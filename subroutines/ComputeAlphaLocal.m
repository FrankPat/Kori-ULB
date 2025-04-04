function [alpha_local] = ComputeAlphaLocal(HB, shMASK, ctr)
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
slope_x = CheckSlopeOneDimension(HB, ctr, 1); % Compute along x-direction
slope_y = CheckSlopeOneDimension(HB, ctr, 2); % Compute along y-direction

% Compute local slope magnitude
alpha_local = atan(sqrt(slope_x.^2 + slope_y.^2));

end

function slope = CheckSlopeOneDimension(HB, ctr, dim)
% Compute the basal slope in one dimension
%
% Inputs:
%   HB  - Ice draft depth (negative downwards) [matrix]
%   ctr.delta - Grid spacing in the given direction [scalar]
%   dim - Dimension to compute slope (1 for x, 2 for y) [scalar]
%
% Output:
%   slope - Slope along the specified dimension [matrix]

% Shift the matrix in both positive and negative directions
HB_plus = circshift(HB, -1, dim);
HB_minus = circshift(HB, 1, dim);

% Compute slopes using a centered difference method
slope_both = (HB_minus - HB_plus) ./ sqrt((2 * ctr.delta)^2);
slope_right = (HB - HB_plus) ./ sqrt(ctr.delta^2);
slope_left = (HB_minus - HB) ./ sqrt(ctr.delta^2);

% Combine slope calculations, prioritizing valid values
slope = slope_both;
slope(isnan(slope)) = slope_right(isnan(slope));
slope(isnan(slope)) = slope_left(isnan(slope));

% Set remaining NaNs to zero
slope(isnan(slope)) = 0;
end
