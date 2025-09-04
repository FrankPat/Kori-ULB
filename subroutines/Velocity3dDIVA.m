function [ut,vt,wt]=Velocity3dDIVA(ubx,uby,etaD,H,gradsx,gradsy, ...
    gradHx,gradHy,Mb,Bmelt,beta2,zeta,ctr)

% Kori-ULB
% 3D velocity field from 3D viscosity field based on DIVA velocities.

    u=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    v=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    F1=zeros(ctr.imax,ctr.jmax,ctr.kmax);

    % horizontal velocities in 3D
    for k=ctr.kmax-1:-1:1
        etalayer=max(1e-8,(etaD(:,:,k)+etaD(:,:,k+1)));
        F1(:,:,k)=F1(:,:,k+1)+H.*(zeta(k)+zeta(k+1))*(zeta(k+1)-zeta(k))./ ...
            etalayer;
    end
    F1x=0.5*(F1+circshift(F1,[0 -1 0]));    % F1-velocity on u grid
    F1y=0.5*(F1+circshift(F1,[-1 0 0]));    % F1-velocity on v grid
    betax=0.5*(beta2+circshift(beta2,[0 -1]));
    betay=0.5*(beta2+circshift(beta2,[-1 0]));
    for k=1:ctr.kmax
        u(:,:,k)=ubx.*(1+betax.*F1x(:,:,k)); % u on u-grid
        v(:,:,k)=uby.*(1+betay.*F1y(:,:,k)); % v on v-grid
    end
    
    % Horizontal velocities on H-grid
    ut=0.5*(u+circshift(u,[0 1]));
    vt=0.5*(v+circshift(v,[1 0]));
    
    % Vertical velocity from mass conservation (on H-grid) based on
    % horizontal velocities on u and v grids
    wt=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    utjm=circshift(u,[0 1 0]);
    vtim=circshift(v,[1 0 0]);
    
    wt(:,:,ctr.kmax)=ut(:,:,ctr.kmax).*(gradsx-gradHx)+vt(:,:,ctr.kmax).* ...
        (gradsy-gradHy);
    for k=ctr.kmax-1:-1:1
        horu=0.25*(circshift(H,[0 -1])+2*H+circshift(H,[0 1])).* ...
            (u(:,:,k)+u(:,:,k+1)-utjm(:,:,k)-utjm(:,:,k+1))/(2*ctr.delta);
        horv=0.25*(circshift(H,[-1 0])+2*H+circshift(H,[1 0])).* ...
            (v(:,:,k)+v(:,:,k+1)-vtim(:,:,k)-vtim(:,:,k+1))/(2*ctr.delta);
        veru=(u(:,:,k)-u(:,:,k+1))/(zeta(k)-zeta(k+1)).* ...
            (gradsx-zeta(k)*gradHx);
        verv=(v(:,:,k)-v(:,:,k+1))/(zeta(k)-zeta(k+1)).* ...
            (gradsy-zeta(k)*gradHy);
        wt(:,:,k)=(horu+horv+veru+verv)*(zeta(k)-zeta(k+1))+wt(:,:,k+1);
    end
    
%     wt(:,:,1)=ut(:,:,1).*gradsx+vt(:,:,1).*gradsy-Mb+Bmelt;
%     for k=2:ctr.kmax
%         horu=0.25*(circshift(H,[0 -1])+2*H+circshift(H,[0 1])).* ...
%             (u(:,:,k)+u(:,:,k-1)-utjm(:,:,k)-utjm(:,:,k-1))/(2*ctr.delta);
%         horv=0.25*(circshift(H,[-1 0])+2*H+circshift(H,[1 0])).* ...
%             (v(:,:,k)+v(:,:,k-1)-vtim(:,:,k)-vtim(:,:,k-1))/(2*ctr.delta);
%         veru=(u(:,:,k)-u(:,:,k-1))/(zeta(k)-zeta(k-1)).* ...
%             (gradsx-zeta(k)*gradHx);
%         verv=(v(:,:,k)-v(:,:,k-1))/(zeta(k)-zeta(k-1)).* ...
%             (gradsy-zeta(k)*gradHy);
%         wt(:,:,k)=(horu+horv+veru+verv)*(zeta(k)-zeta(k-1))+wt(:,:,k-1);
%     end
    
end


