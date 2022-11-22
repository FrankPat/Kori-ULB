function [Mterm] = ComputeMterm(T_in,S_in,Tf,c_rho_1,c_tau,gamma,E0, ...
    thermal_forcing_factor,sina,par)

% Kori-ULB
% Compute M-term for plume parameterisation.
% This function computes the M-term for the plume parameterisation.
% This is the beginning of Equation 28(a) in Lazeroms et al. 2019.
%
% Parameters
%     T_in : scalar (or array?)
%         Ambient temperature in degrees C.
%     S_in : scalar (or array?)
%         Ambient salinity in psu.
%     Tf : scalar (or array?)
%         Freezing temperature
%     c_rho_1 : scalar
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     c_tau : scalar (or array?)
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     gamma: scalar
%         Effective thermal Stanton number. Can be modulated for tuning.
%     E0: scalar
%         Entrainment coefficient. Can be modulated for tuning.
%     alpha: scalar or array
%         Slope angle in rad (must be positive).
%     thermal_forcing_factor : scalar (or array?)
%         Factor to be multiplied to T0-Tf in the end of Eq. 28a - 
%         either thermal forcing or thermal forcing average.
% Returns
%     Mterm : scalar or array
%         Term to be multiplied with M_hat in Eq 28a.

    Mterm=sqrt((par.beta_coeff_lazero.*S_in.*par.g)./(par.lambda3.* ...
        (par.Latent./par.cp0).^3)).*sqrt(max(0,(1-c_rho_1.*gamma)./ ...
        (par.Cd + E0.*sina))).*((gamma.*E0.*sina)./(gamma+c_tau+ ...
        E0.*sina)).^(3/2).*(T_in-Tf).*thermal_forcing_factor;

end


