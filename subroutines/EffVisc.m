function [eta,dudx,dvdy,dudy,dvdx]=EffVisc(A,uxssa,uyssa,H,par,eta1,MASK,ctr)

% Kori-ULB
% Effective viscosity of the SSA solution. On the borders of the domain, a
% fixed value is determined based on a fixed ice thickness (Hshelf) that
% determines eta1

    dudx=(uxssa-circshift(uxssa,[0 1]))/ctr.delta;
    dvdy=(uyssa-circshift(uyssa,[1 0]))/ctr.delta;
    dudy=(circshift(uxssa,[-1 1])+circshift(uxssa,[-1 0])- ...
        circshift(uxssa,[1 1])-circshift(uxssa,[1 0]))/(4*ctr.delta);
    dvdx=(circshift(uyssa,[0 -1])+circshift(uyssa,[1 -1])- ...
        circshift(uyssa,[0 1])-circshift(uyssa,[1 1]))/(4*ctr.delta);
    EffStr=dudx.^2+dvdy.^2+dudx.*dvdy+0.25*(dudy+dvdx).^2;
    eta=0.5*H.*A.^(-1./par.n).*EffStr.^((1-par.n)/(2*par.n));
    eta(eta<=0)=NaN;
    eta(MASK==0)=eta(MASK==0)./ctr.shelftune(MASK==0);  %VL: 2D shelftune
    eta(isnan(eta))=1e7;
    
    % NEED TO CHECK THE VISCOSITY ON EDGES OF ICE SHELF/SEA ICE!!!
    % NOW TAKEN AS A CONSTANT VALUE FOR STABILITY AS MOSTLY SEA ICE
    
    if ctr.mismip>=1
        MASKb=ones(ctr.imax,ctr.jmax); % use constant eta on edges of ice shelf
        eta(:,1)=eta(:,3); % ice divide - symmetric
        eta(1,:)=eta(3,:); % periodic BC
        if ctr.mismip==1
            eta(ctr.imax,:)=eta(ctr.imax-2,:); % periodic BC
            MASKb(1:ctr.imax,1:ctr.jmax-2)=0;
        else
            MASKb(1:ctr.imax-1,1:ctr.jmax-2)=0;
        end
        eta(MASKb==1 & MASK==0)=eta1(MASKb==1 & MASK==0);
        % boundary conditions on strain components
        dudx(:,1)=-dudx(:,2);
        dudx(:,ctr.jmax)=dudx(:,ctr.jmax-1);
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
    else
        MASKb=ones(ctr.imax,ctr.jmax); % use constant eta on edges of ice shelf
        MASKb(3:ctr.imax-2,3:ctr.jmax-2)=0;
        eta(MASKb==1 & MASK==0)=eta1(MASKb==1 & MASK==0);
    end
    eta=min(max(eta,1e5),1e15);
    
end


