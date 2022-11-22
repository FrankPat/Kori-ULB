function [x_hat]=ComputeXhat(ice_draft_depth,zGL,T_in,Tf,E0,c_tau, ...
    gamma,sina,par)

% Kori-ULB
% Compute x_hat (or x_tilda) for plume parameterisation.
% This function computes x_hat (or x_tilda) for the plume parameterisation.
% It is a dimensionless coordinate describing distance from plume origin.
% This is Equation 28(b) in Lazeroms et al. 2019.
%
%  Parameters
%     ice_draft_depth : scalar (or array?)
%         Depth of the ice draft in m (depth is negative!).
%     zGL: scalar (or array?)
%         Depth of the grounding line where the source of the plume
%         is in m (depth is negative!).
%     T_in : scalar (or array?)
%         Ambient temperature in degrees C.
%     Tf : scalar (or array?)
%         Freezing temperature
%     E0: scalar
%         Entrainment coefficient. Can be modulated for tuning.
%     c_tau : scalar (or array?)
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     alpha: scalar or array
%         Slope angle in rad (must be positive).
%     gamma: scalar
%         Effective thermal Stanton number. Can be modulated for tuning.
% Returns
%     x_hat : scalar or array
%         Dimensionless coordinate describing distance from plume origin.
%         Has to be between 0 and 1
%     x_tilda in Eq 28b in Lazeroms et al. 2019

    x_hat=par.lambda3*(ice_draft_depth-zGL)./((T_in-Tf).*(1+ ...
        par.C_eps_lazero.*((E0*sina)./(gamma+c_tau+E0*sina)).^(3/4)));
    % all of this derivation is only valid for x in [0,1]
    x_hat=max(min(x_hat,1),0);
    x_hat(zGL>=0)=0;
    
end


