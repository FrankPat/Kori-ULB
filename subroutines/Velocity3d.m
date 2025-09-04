function [ut,vt,wt,pl]=Velocity3d(udx,udy,ubx,uby,pxy,zeta,Mb,Bmelt,gradsx, ...
    gradsy,gradHx,gradHy,ctr)

% Kori-ULB
% 3D velocity field from 2D horizontal velocities derived with SIA or
% hybrid model. Uses shape functions.

    % horizontal velocities on H grid
    uxdt=0.5*(udx+circshift(udx,[0 1]));
    uydt=0.5*(udy+circshift(udy,[1 0]));
    uxbt=0.5*(ubx+circshift(ubx,[0 1]));
    uybt=0.5*(uby+circshift(uby,[1 0]));
    pl=repmat(pxy,[1,1,ctr.kmax]);
    zl=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    ushllib=(pl+2)./(pl+1).*(1-zl.^(pl+1));
    ut=repmat(uxdt,[1,1,ctr.kmax]).*ushllib+repmat(uxbt,[1,1,ctr.kmax]);
    vt=repmat(uydt,[1,1,ctr.kmax]).*ushllib+repmat(uybt,[1,1,ctr.kmax]);

    % vertical velocity according to Pattyn (2010) with Lliboutry shape
    % function
    wshllib=1-(pl+2).*zl./(pl+1)+1./(pl+1).*zl.^(pl+2);
    wt=repmat(-Mb,[1,1,ctr.kmax]).*wshllib-Bmelt+ut.* ...
        (repmat(gradsx,[1,1,ctr.kmax])-zl.*repmat(gradHx,[1,1,ctr.kmax])) ...
        +vt.*(repmat(gradsy,[1,1,ctr.kmax])-zl.*repmat(gradHy,[1,1,ctr.kmax]));
    
end


