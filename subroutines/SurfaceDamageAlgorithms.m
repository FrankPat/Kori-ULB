function [ds]=SurfaceDamageAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,eta,H)

    % Kori-ULB
    % Damage functions

    % ctr.srfdamage=0: No damage
    % ctr.srfdamage=1: Nye damage function      (following Sun et al., 2017)
    % ctr.srfdamage=2: Weertman damage function (following Lai et al., 2020)
    % ctr.srfdamage=3: Kachuck based damage (removing floating contr.)   
    % ctr.srfdamage=4: Lai damage function      (following Lai et al., 2020)

    % Initialize to zeros
    ds=zeros(ctr.imax,ctr.jmax);

    if ctr.srfdamage~=0
        eps=1e-8; % avoid zero values
        [lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); % 1st principal strain
        % convert strain to stress
	% note that eta is the vertically integrated viscosity
	% hence, the ice thickness is considered in there
        tau1=2*lambda1.*eta;
	
	if ctr.srfdamage==1 || ctr.srfdamage==2
	    % No water thickness (TO DO!)
	    dw=zeros(ctr.imax,ctr.jmax);
            ds=tau1./(par.rho*par.g*(H+eps))+par.rhow*dw/par.rho;
            if ctr.srfdamage==2
                ds=ds.*pi*0.5;
	    end
	end

        if ctr.srfdamage==3
            alpha=lambda2./lambda1;
	    ds=tau1.*(2+alpha)./(par.rho*par.g*(H+eps))+par.rhow*dw/par.rho;
        end

	if ctr.srfdamage==4
            F=1.122;
            f=1.068;
            ds=(tau1.*pi*F)./(6*par.rho*par.g*f*(H+eps));
	end
    end
    % Limit to damage limit
    ds=max(0,min(ds,H*par.dlim));   
end
