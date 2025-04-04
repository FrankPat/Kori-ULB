function [damage]=NyeDamage(par,ctr,dudx,dvdy,dudy,dvdx,eta,H,HAF,MASK)

% Kori-ULB
% compute total crevasse depth d, either 0, the surface crevasse depth or
% the total, whichever is the larger, given the vertically integrated
% 'damaged' stretching stress t
%
% The crevasse depths are given by the Nye zero-stress rule, used recently
% in e.g Nick et al., 2010 JoG (eq.4 and 6)-our t is their Rxx
%
% After Sun et al. (2017)

    eps=1e-8;
    dw=zeros(ctr.imax,ctr.jmax);
    [lambda0,~]=PrincipalStrain(dudx,dvdy,dudy,dvdx);
    lambda0=2*lambda0.*eta; % convert strain to stress
    % basal crevasses
    db=(lambda0./(par.rho*par.g*(H+eps))-max(HAF,0))*par.rho/(par.rhow-par.rho);
    db(MASK==1)=0; % Javi - no  basal crevasses for grounded ice
    % surface crevasses
    ds=lambda0./(par.rho*par.g*(H+eps))+par.rhow*dw/par.rho;
    damage=max(0,min(max(ds,ds+db),H*par.dlim));
    
end


