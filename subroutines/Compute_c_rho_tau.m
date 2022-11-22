function [c_rho_1,c_rho_2,c_tau]=Compute_c_rho_tau(gamma,S_in,par)

% Kori-ULB
% Compute the c-constants for plume parameterisation.
% This function computes c_rho_1, c_rho_2 and c_tau for the plume 
% parameterisation. They are constants given in Table 1 of 
% Lazeroms et al. 2019.
% 
% Parameters
%     gamma: scalar
%         Effective thermal Stanton number. Can be modulated for tuning.
%     S_in : scalar (or array?)
%         Ambient salinity in psu.
% Returns
%     c_rho_1 : scalar (or array?)
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     c_rho_2 : scalar
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     c_tau : scalar (or array?)
%         Constant given in Table 1 in Lazeroms et al. 2019.
    
    c_rho_1=(par.Latent*par.alpha_coeff_lazero)./(par.cp0*gamma* ...
        par.beta_coeff_lazero.* S_in);
    c_rho_2=-par.lambda1*par.alpha_coeff_lazero./par.beta_coeff_lazero;
    c_tau=c_rho_2./c_rho_1;

end


