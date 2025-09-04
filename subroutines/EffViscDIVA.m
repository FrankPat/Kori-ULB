function [eta,etaD,dudx,dvdy,dudy,dvdx]=EffViscDIVA(A3d,betax,betay,ubx,uby, ...
    etaD,H,damage,uxssa,uyssa,zeta,MASK,glMASK,shelftune,ctr,par)

% Kori-ULB
% Effective viscosity of the DIVA solution. On the borders of the domain, a
% fixed value is determined either based on the mean viscosity of the ice
% shelf or a specific value for the ocean.

    eps=1e-8;
    [dudx,dudy,dvdx,dvdy]=VelocityGradients(uxssa,uyssa,ctr);
    
    % Basal stress to compute full 3D viscosity etaD(x,y,z).
    % betax and betax are defined on the velocity grid.
    taubx=betax.*ubx;
    tauby=betay.*uby;
    % Stagger tau values to the ice thickness grid (eta) to divide by eta_diva.
    taubx=0.5*(taubx+circshift(taubx,[0 1]));
    tauby=0.5*(tauby+circshift(tauby,[1 0]));
    
    taubxD=repmat(taubx,[1,1,ctr.kmax]);
    taubyD=repmat(tauby,[1,1,ctr.kmax]);
    zetaD=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);

    % Eq. 36 Lipscomb et al. (2019). Add a small number to prevent divide
    % by zero
    dudzD=taubxD.*zetaD./max(etaD,eps);
    dvdzD=taubyD.*zetaD./max(etaD,eps);
    
    % Strain rates including vertical derivatives.
    dudxD=repmat(dudx,[1,1,ctr.kmax]);
    dvdyD=repmat(dvdy,[1,1,ctr.kmax]);
    dudyD=repmat(dudy,[1,1,ctr.kmax]);
    dvdxD=repmat(dvdx,[1,1,ctr.kmax]);
    EffStrD=dudxD.^2+dvdyD.^2+dudxD.*dvdyD+0.25*(dudyD+dvdxD).^2+ ...
      0.25*(dudzD.^2+dvdzD.^2);
    % Obtain 2D strain rates by vertically averaging the 3D DIVA strain rates.
    EffStrD=max(EffStrD,1e-12);
    
    % Include damage in the calculation of ice viscosity
    scale_eta=((H-min(min(H.*par.damlim,damage),H-eps))./(H+eps)).^(-par.n);
    Astar=A3d.*repmat(scale_eta,[1,1,ctr.kmax]);
    etaD=0.5*Astar.^(-1./par.n).*EffStrD.^((1-par.n)/(2*par.n));
    
    % Vertical integration of effective viscosity
    eta=zeros(ctr.imax,ctr.jmax);
    for k=2:ctr.kmax
        eta=eta+0.5*(etaD(:,:,k)+etaD(:,:,k-1))*(zeta(k)-zeta(k-1));
    end
    eta=eta.*H;
    eta(MASK==0)=eta(MASK==0)./shelftune(MASK==0);  %VL: 2D shelftune
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
        eta(MASKb==1)=trimmean(eta(MASKb==1),10);
        
        % Instead of calculating effective viscosity on the sea ice,
        % keep constant viscosity. Need to further check how to deal
        % with this. May be quoted
        eta(glMASK==6)=ctr.OceanVisc; % Default: 1.0e7. Daniel: 8.0e9
        % 1e8 is approx for ice shelf of 200m thick
    end

end


