function [M_hat]=ComputeMhat(X_hat)
 
% Kori-ULB
% Compute M_hat for plume parameterisation.
% This function computes the M_hat (or M0) for the plume parameterisation.
% This is Equation 26 in Lazeroms et al. 2019.
%
% Parameters
%     X_hat : scalar or array
%         Coordinate describing distance from plume origin.
% Returns
%     M_hat : scalar or array
%     Dimensionless melt rate, to be multiplied with the Mterm in Eq 28a.
% here M_hat = M0(X_tilde), Eq 26 in Lazeroms et al. 2019

    M_hat=1./(2*sqrt(2)).*(3*(1-X_hat).^(4/3)-1).*sqrt(1-(1-X_hat).^(4/3));
    
end


