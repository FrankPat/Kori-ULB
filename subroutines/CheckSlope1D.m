function slope = CheckSlope1D(HB, ctr, dim)

% Kori-ULB
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


