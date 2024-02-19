function [ds]=SurfaceDamageAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,eta,H,MASK)

    % Kori-ULB
    % Damage functions

    % ctr.srfdamage=0: No damage
    % ctr.srfdamage=1: Nye damage function      (following Sun et al., 2017)
    % ctr.srfdamage=2: Weertman damage function (following Lai et al., 2020)
    % ctr.srfdamage=3: Kachuck based damage (removing floating contr.)   
    % ctr.srfdamage=4: Lai damage function      (following Lai et al., 2020)

    % Initialize to zeros
    ds=zeros(ctr.imax,ctr.jmax);

    eps=1e-8; % avoid zero values
    [lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); % 1st principal strain
    % convert strain to stress
    % note that eta is the vertically integrated viscosity
    % hence, the ice thickness is considered in there
    tau1=2*lambda1.*eta./(H+eps);
  
    % CWR: crevasse width ratio
    % CW:  crevasse width [m]
    CR  = 20.0;
    CWR = CR./ctr.delta;

    if ctr.srfdamage~=0
        %eps=1e-8; % avoid zero values
        %[lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); % 1st principal strain
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
        % first principal stress 150 kPa
        %if tau1./(H+eps)<150000
        %        tau1=0.0
        %end

        dw=zeros(ctr.imax,ctr.jmax);

	if ctr.srfdamage==1 || ctr.srfdamage==2 || ctr.srfdamage==3
	    % No water thickness (TO DO!)
	    %dw=zeros(ctr.imax,ctr.jmax);
            ds=tau1./(par.rho*par.g)+par.rhow*dw/par.rho;
            if ctr.srfdamage==2
                ds=pi*0.5*tau1./(par.rho*par.g)+par.rhow*dw/par.rho;
	    end
	    if ctr.srfdamage==3
		alpha=lambda2./lambda1;
                ds=ds.*(2+alpha);
            end
	end

        %if ctr.srfdamage==3
        %    alpha=lambda2./lambda1;
	%    ds=tau1.*(2+alpha)./(par.rho*par.g*(H+eps))+par.rhow*dw/par.rho;
        %end

	if ctr.srfdamage==4
            F=1.122;
            f=1.068;
            ds=(tau1.*pi*F)./(6*par.rho*par.g*f);
	end
    
        %if tau1 < par.tauice
        %    ds=0.0;
        %end

        % ds cannot be bigger than ice thickness
        %ds=max(0,min(ds,H.*par.dlim));

    end
    
    % Limit to damage limit
    %ds=max(0,min(ds,H*par.dlim));
    % jablasco: damage can only occur if yield strength of ice is reached

    ds=max(0.0,ds.*CWR);

    % Damage can only form if Effecive Stress is > yield-strength ice
    for i=1:ctr.imax
    	    for j=1:ctr.jmax
                if tau1(i,j)<par.tauice
                    ds(i,j)=0.0;
    	end
    end

end
