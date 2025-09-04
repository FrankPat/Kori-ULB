function [dudx,dudy,dvdx,dvdy]=VelocityGradients(uxssa,uyssa,ctr)

% Kori-ULB
% Horizontal velocity gradients for calculation of effective strain.

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

end


