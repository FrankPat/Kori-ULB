function [Melt]=PlumeMelt2019(T_in, S_in, ice_draft_depth, zGL, sina, par, ctr)

% Kori-ULB
% Apply the plume parametrization.
% This function computes the basal melt based on a plume parametrization
% (see Lazeroms et al. 2018 and Lazeroms et al. 2019).
%
% Parameters
%   T_in : scalar (or array?)
%       Ambient temperature in degrees C.
%   S_in : scalar (or array?)
%       Ambient salinity in psu.
%   ice_draft_depth : scalar or array
%       Depth of the ice draft in m (depth is negative!).
%   zGL: scalar or array
%       Depth of the grounding line where the source of the plume is in m (depth is negative!).
%   alpha: scalar or array
%       Slope angle in rad (must be positive).
%   gamma: scalar
%       Effective thermal Stanton number. Can be modulated for tuning.
%   E0: scalar
%       Entrainment coefficient. Can be modulated for tuning.
%   picop: Boolean
%       Option defining which Mterm function to use.
% Returns
%   melt_rate : scalar or array
%       Melt rate in m ice per second.

    E0=par.Eo;
    gamma=ctr.gammaTplume;

    [c_rho_1,c_rho_2,c_tau]=Compute_c_rho_tau(gamma,S_in,par);

    % freezing temperature at the grounding line
    Tf=par.lambda1*S_in+par.lambda2+par.lambda3*zGL; % Ocean freezing point
    thermal_forcing=T_in-Tf;

    x_hat=ComputeXhat(ice_draft_depth,zGL,T_in,Tf,E0,c_tau,gamma,sina,par);
    M_hat=ComputeMhat(x_hat);
    Mterm=ComputeMterm(T_in,S_in,Tf,c_rho_1,c_tau,gamma,E0, ...
        thermal_forcing,sina,par);
    Melt=Mterm.*M_hat.*par.secperyear; % m ice per year

end


