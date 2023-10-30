function [db]=BasalDamageAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,eta,H,HAF)
    % Damage functions
    % ctr.bsldamage=0: No damage
    % ctr.bsldamage=1: Nye damage function      (following Sun et al., 2017: 10.5194/tc-11-2543-2017)
    % ctr.bsldamage=2: Weertman damage function (following Lai et al., 2020: 10.1038/s41586-020-2627-8)
    %                  Is the factor pi/2 valid for basal crevasses?
    % ctr.bsldamage=3: Kachuck damage function  (following Kachuck et al., 2022: 10.1017/jog.2022.12 )

    % Initialize to zeros
    db=zeros(ctr.imax,ctr.jmax);

    if ctr.bsldamage~=0
        eps=1e-8; % avoid zero values
        [lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); % 1st/2nd principal strain
        % convert strain to stress
	% note that eta is the vertically integrated viscosity
	% hence, the ice thickness is considered in there
        tau1=2*lambda1.*eta;
	
	if ctr.bsldamage==1 || ctr.bsldamage==2
            db=(par.rho/(par.rhow-par.rho))*((tau1./(par.rho*par.g*(H+eps)))-max(HAF,0));
            if ctr.srfdamage==2
                db=db.*pi*0.5;
	    end
	end

	if ctr.bsldamage==3
            alpha=lambda2./lambda1;
            db=(par.rho/(par.rhow-par.rho))*((2+alpha).*tau1./(par.rho*par.g*(H+eps));
	end
    end
    % Limit to damage limit
    db=max(0,min(db,H*par.dlim));   
end
