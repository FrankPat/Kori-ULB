function [tmp,ctr]=Temperature3d(tmp,Mb,Ts,pxy,par, ...
    ctr,dt,gradsx,gradsy,gradHx,gradHy,udx,udy,ub,ubx,uby,zeta,gradxy,H, ...
    dzc,dzm,dzp,G,taudxy,A,DeltaT,MASK,Bmelt,cnt)

% Kori-ULB
% 3d englacial temperature calculation in ice sheet and ice shelves

    tmp(:,:,1)=Ts+par.T0;
    
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
    
    % horizontal advection
    dTdxm=(circshift(tmp,[0 -1])-tmp)/ctr.delta;
    dTdxp=(tmp-circshift(tmp,[0 1]))/ctr.delta;
    advecx=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    advecx(ut>0)=ut(ut>0).*dTdxp(ut>0)*dt;
    advecx(ut<0)=ut(ut<0).*dTdxm(ut<0)*dt;
    dTdym=(circshift(tmp,[-1 0])-tmp)/ctr.delta;
    dTdyp=(tmp-circshift(tmp,[1 0]))/ctr.delta;
    advecy=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    advecy(vt>0)=vt(vt>0).*dTdyp(vt>0)*dt;
    advecy(vt<0)=vt(vt<0).*dTdym(vt<0)*dt;
    
    % vertical velocity according to Pattyn (2010) with Lliboutry shape
    % function

    wshllib=1-(pl+2).*zl./(pl+1)+1./(pl+1).*zl.^(pl+2);
    w=repmat(-max(Mb,1e-5),[1,1,ctr.kmax]).*wshllib-Bmelt+ut.* ...
        (repmat(gradsx,[1,1,ctr.kmax])-zl.*repmat(gradHx,[1,1,ctr.kmax])) ...
        +vt.*(repmat(gradsy,[1,1,ctr.kmax])-zl.*repmat(gradHy,[1,1,ctr.kmax]));
        
    % Internal heating
    repz=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    dudz=repmat(2*A.*taudxy.^par.n.*H.*(pxy+2)/(par.n+2),[1,1,ctr.kmax]).* ...
        repz.^pl;
    fric=par.rho*par.g*par.kdif*dt*repz.*dudz.*repmat(sqrt(gradxy), ...
        [1,1,ctr.kmax])/par.K;
    repmask=repmat(MASK,[1,1,ctr.kmax]);
    fric(repmask==0)=0; % no frictional heat in ice shelves
    extraterm=max(min(fric-advecx-advecy,5),-5); %10

    % Temperature solution
    repH2=repmat(H,[1,1,ctr.kmax]);
    repH=max(repH2,1e-8);
    Tp=par.pmp*repH2.*repz;
    atp=(2*par.kdif*par.secperyear./(repH.*dzm)-w)*dt./(repH.*dzc);
    btp=1+2*par.kdif*par.secperyear*dt./((repH.^2).*dzp.*dzm);
    ctp=(2*par.kdif*par.secperyear./(repH.*dzp)+w)*dt./(repH.*dzc);
    ftp=ones(ctr.imax,ctr.jmax,ctr.kmax);
    gtp=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    % basal BC with strain heating (correction with ud)
    gtp(:,:,ctr.kmax)=(G+taudxy.*ub/par.secperyear).*H.*dzm(:,:,ctr.kmax)/par.K;
    
    % ice shelves basal boundary condition
    mask=zeros(ctr.imax,ctr.jmax);
    mask(MASK==0)=1;
    repmask=repmat(mask,[1,1,ctr.kmax]);
    repmask(:,:,1:ctr.kmax-1)=0;
    TBshelf=par.T0+repmat(min(par.Toi+ctr.meltfactor*DeltaT(cnt)- ...
        0.12e-3*par.rho*H/par.rhow,0),[1,1,ctr.kmax]);
    ftp(repmask==1)=0;
    gtp(repmask==1)=TBshelf(repmask==1);
    
    % Temperature solution
    for k=ctr.kmax-1:-1:2
        ftp(:,:,k)=atp(:,:,k)./(btp(:,:,k)-ctp(:,:,k).*ftp(:,:,k+1));
        gtp(:,:,k)=(tmp(:,:,k)+extraterm(:,:,k)+ ...
            ctp(:,:,k).*gtp(:,:,k+1))./(btp(:,:,k)-ctp(:,:,k).*ftp(:,:,k+1));
    end
    for k=2:ctr.kmax
        tmp(:,:,k)=tmp(:,:,k-1).*ftp(:,:,k)+gtp(:,:,k);
    end
    c3=ceil(par.rhom/par.rho)*10+ceil(par.rhow/par.rho);
    ctr.runmode(DeltaT(cnt)==c3)=5;
    
    % Correction for unstable temperature profiles
    % find an anomaly where temperature decreases with depth
    [ipos,jpos]=find((tmp(:,:,ctr.kmax-1)-tmp(:,:,ctr.kmax))>3 | ...
        (tmp(:,:,ctr.kmax-2)-tmp(:,:,ctr.kmax-1))>3 | ...
        (tmp(:,:,ctr.kmax-3)-tmp(:,:,ctr.kmax-2))>3);
    nMASK=zeros(size(MASK));
    nMASK(sub2ind(size(nMASK),ipos,jpos))=1;
    % use analytical solution at the domain boundary when using basins
    if ctr.basin==1
        nMASK(1,:)=1;
        nMASK(ctr.imax,:)=1;
        nMASK(:,1)=1;
        nMASK(:,ctr.jmax)=1;
    end
    % apply linear temperature profile when MASK=0
    nMASK(nMASK==1 & MASK==0)=2;
    Tgrad=-G/par.K;
    l=sqrt(2*par.kdif*(H+1e-8)./max(Mb,1e-8)*par.secperyear);
    nMASKz=repmat(nMASK,[1,1,ctr.kmax]);
    repl=repmat(l,[1,1,ctr.kmax]);
    repTs=repmat(Ts,[1,1,ctr.kmax]);
    repTgrad=repmat(Tgrad,[1,1,ctr.kmax]);
    tmpb=repTs+sqrt(pi)*0.5*repl.*repTgrad.* ...
        (erf((1-repz).*repH./repl)-erf(repH./repl))+par.T0;
    tmp(nMASKz==1)=tmpb(nMASKz==1);
    Tshelf=repTs+par.T0-(repTs+par.T0-TBshelf).*repz;
    tmp(nMASKz==2)=Tshelf(nMASKz==2); 

    % Correction for pmp
    tmp(tmp>par.T0-Tp)=par.T0-Tp(tmp>par.T0-Tp);

    mintemp=min(min(Ts))+min(DeltaT)+par.T0-5;
    tmp(tmp<mintemp)=mintemp;
    
    if ctr.mismip>0
        tmp(:,1,:)=tmp(:,3,:);
        tmp(:,end,:)=tmp(:,end-1,:);
        tmp(1,:,:)=tmp(3,:,:);
        tmp(end,:,:)=tmp(end-2,:,:);
    end
    
end


