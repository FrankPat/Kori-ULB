function ds=SurfaceDamage(ctr,par,dudx,dvdy,dudy,dvdx,eta,H)

% Kori-ULB
% Surface Damage functions

% ctr.SFdam=0: No damage
% ctr.SFdam=1: Nye damage function      (following Sun et al., 2017, Nick et al., 2010)
% ctr.SFdam=2: Weertman damage function (following Lai et al., 2020)
% ctr.SFdam=3: Kachuck based damage (removing floating contr.)
% ctr.SFdam=4: Lai damage function      (following Lai et al., 2020)

    eps=1e-8; % avoid zero values
    [lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); % 1st/2nd principal strain
    % convert strain to stress -- note that eta is the vertically integrated viscosity
    % hence, the ice thickness is considered in there
    tau1=2*lambda1.*eta./(H+eps);

    dw=zeros(ctr.imax,ctr.jmax); % Water depth in the surface crevasse (Sun2017, Nick2010) -- TO DO!
    if ctr.SFdam==1
        ds=tau1./(par.rho*par.g)+par.rhow.*dw/par.rho;
    elseif ctr.SFdam==2
        ds=pi*0.5*tau1./(par.rho*par.g)+par.rhow*dw/par.rho;
    elseif ctr.SFdam==3
        alpha=lambda2./lambda1;
        ds=ds.*(2+alpha);
        %    ds=tau1.*(2+alpha)./(par.rho*par.g*(H+eps))+par.rhow*dw/par.rho;
    elseif ctr.SFdam==4
        F=1.122;
        f=1.068;
        ds=(tau1.*pi*F)./(6*par.rho*par.g*f);
    else
        ds=zeros(ctr.imax,ctr.jmax);
    end
    ds=max(0,min(ds,H));

end


