function [db]=BasalDamageAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,eta,H,HAF)
    % Damage functions
    % ctr.bsldamage=0: No damage
    % ctr.bsldamage=1: Nye damage function      (following Sun et al., 2017: 10.5194/tc-11-2543-2017)
    % ctr.bsldamage=2: Weertman damage function (following Lai et al., 2020: 10.1038/s41586-020-2627-8)
    %                  Is the factor pi/2 valid for basal crevasses?
    % ctr.bsldamage=3: Kachuck damage function  (following Kachuck et al., 2022: 10.1017/jog.2022.12 )

    % Initialize to zeros
    db=zeros(ctr.imax,ctr.jmax);

    eps=1e-8; % avoid zero values
    [lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); % 1st/2nd principal strain
    % convert strain to stress
    % note that eta is the vertically integrated viscosity
    % hence, the ice thickness is considered in there
    tau1=2*lambda1.*eta./(H+eps);

    % CWR: crevasse width ratio
    % CW:  crevasse width [m]
    CR  = 20.0;
    CWR = CR./ctr.delta;

    if ctr.bsldamage~=0
        %eps=1e-8; % avoid zero values
        %[lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); % 1st/2nd principal strain
        % convert strain to stress
	% note that eta is the vertically integrated viscosity
	% hence, the ice thickness is considered in there
        %tau1=2*lambda1.*eta./(H+eps);

	% in order for crevasses to open
        % EffStr=dudx.^2+dvdy.^2+dudx.*dvdy+0.25*(dudy+dvdx).^2;
        % yield srength from strain
        % tau_c=2.*eta.*EffStr./(H+eps);
        %if tau1 < par.tauice
        %    tau1=0.0;
        %end

	if ctr.bsldamage==1 || ctr.bsldamage==2 || ctr.bsldamage==3
            db=(par.rho/(par.rhow-par.rho))*((tau1./(par.rho*par.g))-max(HAF,0));
            if ctr.bsldamage==2
                db=(par.rho/(par.rhow-par.rho))*((pi*0.5*tau1./(par.rho*par.g))-max(HAF,0));
	    end
	    if ctr.bsldamage==3
		alpha=lambda2./lambda1;
                db=db.*(2+alpha);
            end
	end

	%if ctr.bsldamage==3
        %    alpha=lambda2./lambda1;
        %    db=(par.rho/(par.rhow-par.rho))*(tau1.*(2+alpha)./(par.rho*par.g*(H+eps)));
	%end
    
	%db=max(0,min(db,H.*par.dlim));

    end
    % Limit to damage limit
    %db=max(0,min(db,H*par.dlim));   

    db=max(0.0,db.*CWR);

    for i=1:ctr.imax
            for j=1:ctr.jmax
                if tau1(i,j)<par.tauice
                    db(i,j)=0.0;
        end
    end

end
