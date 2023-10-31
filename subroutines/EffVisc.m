function [eta,dudx,dvdy,dudy,dvdx]=EffVisc(A,uxssa,uyssa,H,par,MASK, ...
    glMASK,shelftune,ctr)

% Kori-ULB
% Effective viscosity of the SSA solution. On the borders of the domain, a
% fixed value is determined based on a fixed ice thickness (Hshelf) that
% determines eta1

    dudx=(uxssa-circshift(uxssa,[0 1]))/ctr.delta;
    dudx(:,1)=dudx(:,2);
    dudx(:,ctr.jmax)=dudx(:,ctr.jmax-1);
    dvdy=(uyssa-circshift(uyssa,[1 0]))/ctr.delta;
    dvdy(1,:)=dvdy(2,:);
    dvdy(ctr.imax,:)=dvdy(ctr.imax-1,:);
    dudy=(circshift(uxssa,[-1 1])+circshift(uxssa,[-1 0])- ...
        circshift(uxssa,[1 1])-circshift(uxssa,[1 0]))/(4*ctr.delta);
    dudy(1,:)=dudy(2,:);
    dudy(ctr.imax,:)=dudy(ctr.imax-1,:);
    dudy(:,1)=dudy(:,2);
    dudy(:,ctr.jmax)=dudy(:,ctr.jmax-1);
    dvdx=(circshift(uyssa,[0 -1])+circshift(uyssa,[1 -1])- ...
        circshift(uyssa,[0 1])-circshift(uyssa,[1 1]))/(4*ctr.delta);
    dvdx(1,:)=dvdx(2,:);
    dvdx(ctr.imax,:)=dvdx(ctr.imax-1,:);
    dvdx(:,1)=dvdx(:,2);
    dvdx(:,ctr.jmax)=dvdx(:,ctr.jmax-1);
    
    if ctr.mismip>=1
        dvdx(:,1)=-dvdx(:,2);
        dvdx(:,ctr.jmax)=dvdx(:,ctr.jmax-1);
        dvdy(1,:)=dvdy(3,:);
        dudy(1,:)=dudy(3,:);
        if ctr.mismip==1
            dvdy(ctr.imax,:)=dvdy(ctr.imax-3,:);
            dudy(ctr.imax,:)=dudy(ctr.imax-3,:);
        else
            dvdy(ctr.imax,:)=dvdy(ctr.imax-1,:);
            dudy(ctr.imax,:)=dudy(ctr.imax-1,:);
        end
    end
    EffStr=dudx.^2+dvdy.^2+dudx.*dvdy+0.25*(dudy+dvdx).^2;
    EffStr=max(EffStr,1e-12);
    eta=0.5*H.*A.^(-1./par.n).*EffStr.^((1-par.n)/(2*par.n));
    eta(MASK==0)=eta(MASK==0)./shelftune(MASK==0);  %VL: 2D shelftune

    % Bassis et al., (2021) regularization
    % DOI: 0.1126/science.abf6271
    if crt.bassis==1
        d_grain=5e-3; % Diffussion creep. Grain size [m]
        eta_diff=H.*A.^(-1./par.n)/(2.*d_grain.^2);
        tau_c=2.*eta.*EffStr./H; % yield strength from strain
        if tau_c > par.tauice
            % Once reached the failure, the effective stress decreases with increasing strain rate.
            tau_y=max(tau_c-(tau_c-par.taulim).*EffStr./(H.*par.strcrit),par.taulim);
            eta_plas=H.*tau_y./(2.*EffStr);
        else
            % If failure not reached, the ice yield strength value doesn't change
            eta_plas = H.*par.tauice./(2.*EffStr);
        end
        % numerical convergence value
        eta_min=1e6;
	% regularized viscosity
        eta=eta_min+1./(1./eta_diff+1./eta_disl+1./eta_plas);
    end

    if ctr.shelf==1 || ctr.schoof>0
        MASKb=ones(ctr.imax,ctr.jmax); % use constant eta on edges of ice shelf
        if ctr.mismip==0
            MASKb(3:ctr.imax-2,3:ctr.jmax-2)=0;
        elseif ctr.mismip==1
            MASKb(:,1:ctr.jmax-2)=0;
        elseif ctr.mismip==2
            MASKb(1:ctr.imax-2,1:ctr.jmax-2)=0;
        end
        MASKb(MASK==1)=0;
        % remove outliers (10%) - especially important for basins where
        % grounded parts may exist on the domain boundary
        eta(MASKb==1)=trimmean(eta(MASKb==1),10);
        
        % Instead of calculating effective viscosity on the sea ice,
        % keep constant viscosity. Need to further check how to deal
        % with this. May be quoted
        eta(glMASK==6)=1e7;
    end
    
end


