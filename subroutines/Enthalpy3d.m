function [E,Epmp,wat,CTSm,CTSp,Dbw,Dfw,Ht,tmp]= ...
    Enthalpy3d(par,ctr,E,Mb,Ts,G,A,H, ...
    pxy,dt,gradsx,gradsy,gradHx,gradHy,gradxy,taudxy, ...
    udx,udy,ub,ubx,uby,zeta,dzc,dzm,dzp,DeltaT,DeltaTo,MASK, ...
    Bmelt,CTSm,CTSp,Hw,Ht,Dbw,Dfw,etaD,beta2,cnt)

% Kori-ULB
% 3d englacial enthalpy calculation in ice sheet and ice shelves
 
    % surface BC
    E(:,:,1)=par.cp*(Ts+par.T0-par.Tref);

    % Enthalpy at pressure melting point
    repz=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]); 
    repH=max(repmat(H,[1,1,ctr.kmax]),1e-8);
    repH2=repmat(H,[1,1,ctr.kmax]);
    Tp=par.pmp*repH2.*repz;
    Tpmp=par.T0-Tp;
    Epmp=par.cp*(Tpmp-par.Tref);

    % ice diffusivity
    kdif=(par.Kc./par.rho)+(((par.K0-par.Kc)/par.rho)*heaviside(E-Epmp));
        
    % Drain excess englacial water computed at previous time step
    if ctr.drain>0
        DFlux=(par.rhof./par.rho).*Dfw.*par.Latent; 
        Bmelt=Bmelt+Dbw;
    else
        DFlux=zeros(size(E));
    end
    
    % horizontal and vertical velocities on H grid
    if ctr.SSA==3
        [ut,vt,wt]=Velocity3dDIVA(ubx,uby,etaD,H,gradsx,gradsy, ...
            gradHx,gradHy,Mb,Bmelt,beta2,zeta,ctr);
    else
        [ut,vt,wt,pl]=Velocity3d(udx,udy,ubx,uby,pxy,zeta,Mb,Bmelt,gradsx, ...
            gradsy,gradHx,gradHy,ctr);
    end
    
    % horizontal advection
    dEdxm=(circshift(E,[0 -1])-E)/ctr.delta;
    dEdxp=(E-circshift(E,[0 1]))/ctr.delta;
    advecx=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    advecx(ut>0)=ut(ut>0).*dEdxp(ut>0)*dt;
    advecx(ut<0)=ut(ut<0).*dEdxm(ut<0)*dt;
    dEdym=(circshift(E,[-1 0])-E)/ctr.delta;
    dEdyp=(E-circshift(E,[1 0]))/ctr.delta;
    advecy=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    advecy(vt>0)=vt(vt>0).*dEdyp(vt>0)*dt;
    advecy(vt<0)=vt(vt<0).*dEdym(vt<0)*dt;
        
    % Internal/strain heating [CHECK THIS BECAUSE CHANGED]
    repz=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    if ctr.SSA==3
        dudz=zeros(ctr.imax,ctr.jmax,ctr.kmax);
        for k=2:ctr.kmax
            dudz(:,:,k)=(ut(:,:,k)-ut(:,:,k-1))/(zeta(k)-zeta(k-1)).* ...
                gradsx+(vt(:,:,k)-vt(:,:,k-1))/(zeta(k)- ...
                zeta(k-1)).*gradsy;
        end
        fric=par.g*dt*repz.*dudz;
    else
        dudz=repmat(2*A.*taudxy.^par.n.*H.*(pxy+2)/ ...
            (par.n+2),[1,1,ctr.kmax]).*repz.^pl;
        fric=par.g*dt*repz.*dudz.*repmat(sqrt(gradxy),[1,1,ctr.kmax]);
    end
    repmask=repmat(MASK,[1,1,ctr.kmax]);
    fric(repmask==0)=0; % no frictional heat in ice shelves
    extraterm=max(min((fric-advecx-advecy-DFlux),5*par.cp),-5*par.cp); 
   
    % central difference for diffusion & upstream difference for advection
    atp=dt*((2*kdif*par.secperyear./(repH.^2.*dzm.*dzc))-wt./(dzm.*repH)); % E(k-1)
    btp=1+dt*((2*kdif*par.secperyear./(repH.^2.*dzp.*dzm))-wt./(dzm.*repH)); % E(k)
    ctp=dt*((2*kdif*par.secperyear./(repH.^2.*dzp.*dzc))); % E(k+1)

    % correction for discontinuity at the CTS
    Kb=kdif;
    Ka=kdif;
    Kb(:,:,2:end-1)=hamean(kdif(:,:,1:end-2),kdif(:,:,2:end-1));
    Ka(:,:,2:end-1)=hamean(kdif(:,:,3:end),kdif(:,:,2:end-1));
    
    % Correction CTS: upstream differences for advection
    atp(CTSp==1)=dt*((2*(Kb(CTSp==1))*par.secperyear./(repH(CTSp==1).^2.* ...
        dzm(CTSp==1).*dzc(CTSp==1)))-wt(CTSp==1)./(dzm(CTSp==1).*repH(CTSp==1))); % E(k-1)
    btp(CTSp==1)=1+dt*(((Ka(CTSp==1)+Kb(CTSp==1))*par.secperyear./ ...
        (repH(CTSp==1).^2.*dzp(CTSp==1).*dzm(CTSp==1)))-wt(CTSp==1)./ ...
        (dzm(CTSp==1).*repH(CTSp==1))); % E(k)
    ctp(CTSp==1)=dt*((2*(Ka(CTSp==1))*par.secperyear./(repH(CTSp==1).^2.* ...
        dzp(CTSp==1).*dzc(CTSp==1)))); % E(k+1)
    atp(CTSm==1)=dt*((2*(Kb(CTSm==1))*par.secperyear./(repH(CTSm==1).^2.* ...
        dzm(CTSm==1).*dzc(CTSm==1)))-wt(CTSm==1)./(dzm(CTSm==1).*repH(CTSm==1))); % E(k-1)
    btp(CTSm==1)=1+dt*(((Ka(CTSm==1)+Kb(CTSm==1))*par.secperyear./ ...
        (repH(CTSm==1).^2.*dzp(CTSm==1).*dzm(CTSm==1)))-wt(CTSm==1)./ ...
        (dzm(CTSm==1).*repH(CTSm==1))); % E(k)
    ctp(CTSm==1)=dt*((2*(Ka(CTSm==1))*par.secperyear./(repH(CTSm==1).^2.* ...
        dzp(CTSm==1).*dzc(CTSm==1)))); % E(k+1)

    ftp=ones(ctr.imax,ctr.jmax,ctr.kmax);
    gtp=zeros(ctr.imax,ctr.jmax,ctr.kmax);

    % BASAL CONDITIONS
    % Cold base and dry or about to get frozen
    Cld=(E(:,:,ctr.kmax)<Epmp(:,:,ctr.kmax) & Hw==0);
    % Cold base and wet
    Blw=(E(:,:,ctr.kmax)<Epmp(:,:,ctr.kmax) & Hw>0);
    % Temperate base
    Abv=((E(:,:,ctr.kmax)>=Epmp(:,:,ctr.kmax)) & Ht==0 & Hw>0);
    % Temperate layer
    Ptv=((E(:,:,ctr.kmax)>=Epmp(:,:,ctr.kmax)) & Ht>0 & Hw>0);

    % Otherwise
    Cld(Cld==0 & Blw==0 & Abv==0 & Ptv==0 & Hw==0)=1;
    Abv(Cld==0 & Blw==0 & Abv==0 & Ptv==0 & Hw>0)=1;

    % BASAL BOUNDARY CONDITIONS
    % basal BC for cold ice + dry bed or bed about to get frozen
    gtp1=((G+taudxy.*ub/par.secperyear).*H.*dzm(:,:,ctr.kmax))*par.cp/par.K; 
    ftp1=ftp(:,:,ctr.kmax);
    % basal BC for temperate base with no positive layer of temperate ice
    gtp3=Epmp(:,:,ctr.kmax);
    ftp3=zeros(ctr.imax,ctr.jmax); % dirichlet bc
    % basal BC for temperate base with positive thickness of temperate ice  
    gtp4=zeros(ctr.imax,ctr.jmax);
    ftp4=ones(ctr.imax,ctr.jmax); % insulated bc
    % mixed bc
    gtp(:,:,ctr.kmax)=Cld.*gtp1 + Blw.*gtp3 + Abv.*gtp3 + Ptv.*gtp4;
    ftp(:,:,ctr.kmax)=Cld.*ftp1 + Blw.*ftp3 + Abv.*ftp3 + Ptv.*ftp4;

    % ice shelves basal boundary condition
    mask=zeros(ctr.imax,ctr.jmax);
    mask(MASK==0)=1;
    repmask=repmat(mask,[1,1,ctr.kmax]);
    repmask(:,:,1:ctr.kmax-1)=0;
    EBshelf=par.cp*((par.T0+repmat(min(par.Toi+DeltaTo(cnt)- ...
              0.12e-3*par.rho*H/par.rhow,0),[1,1,ctr.kmax]))-par.Tref);
    ftp(repmask==1)=0;
    gtp(repmask==1)=EBshelf(repmask==1);

    % Enthalpy solution
    for k=ctr.kmax-1:-1:2
        ftp(:,:,k)=atp(:,:,k)./(btp(:,:,k)-ctp(:,:,k).*ftp(:,:,k+1));
        gtp(:,:,k)=(E(:,:,k)+extraterm(:,:,k)+ ...
            ctp(:,:,k).*gtp(:,:,k+1))./(btp(:,:,k)-ctp(:,:,k).*ftp(:,:,k+1));
    end
    for k=2:ctr.kmax
        E(:,:,k)=E(:,:,k-1).*ftp(:,:,k)+gtp(:,:,k); 
    end
    c3=ceil(par.rhom/par.rho)*10+ceil(par.rhow/par.rho);
    ctr.runmode(DeltaT(cnt)==c3)=5;
        
    % OR: function that calculate CTS position + temperate layer thickness
    [CTSm,CTSp,Ht,E,idx]=CalculateCTS(ctr,E,Epmp,MASK,H,zeta);
    
    % Temperature field  
    Tp=par.pmp*repH2.*repz;
    Tpmp=par.T0-Tp;
    tmp=min((E/par.cp)+par.Tref,Tpmp);
    tmp(E>=Epmp)=Tpmp(E>=Epmp); % temperate ice

    % Water Content 
    wat=max((E-Epmp)./par.Latent,0); % Aschwanden et al, 2012
    wat(E<Epmp)=0; % no water content in cold ice
    wat(tmp<Tpmp)=0; % avoid truncature errors with Epmp
    
    % ---------------------------------------------------------------------
    % Recalculate enthalpy profiles in cold ice
    % Keep prior enthalpy guess for temperate ice
    % ENTM scheme in Greve & Blatter, 2015, 2016
    % ---------------------------------------------------------------------

    E0=E;     % prior guess enthalpy
    tmp0=tmp; % prior guess temperature
    wat0=wat; % prior guess water content

    % surface BC
    E(:,:,1)=par.cp*(Ts+par.T0-par.Tref);
    
    % ice diffusivity
    kdif=(par.Kc./par.rho)+(((par.K0-par.Kc)/par.rho)*heaviside(E-Epmp));
    Kc=par.K./par.cp; kdif(E<Epmp)=Kc./par.rho;
    kdif(E>=Epmp)=par.K0./par.rho;
    
    % horizontal advection
    dEdxm=(circshift(E,[0 -1])-E)/ctr.delta;
    dEdxp=(E-circshift(E,[0 1]))/ctr.delta;
    advecx=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    advecx(ut>0)=ut(ut>0).*dEdxp(ut>0)*dt;
    advecx(ut<0)=ut(ut<0).*dEdxm(ut<0)*dt;
    dEdym=(circshift(E,[-1 0])-E)/ctr.delta;
    dEdyp=(E-circshift(E,[1 0]))/ctr.delta;
    advecy=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    advecy(vt>0)=vt(vt>0).*dEdyp(vt>0)*dt;
    advecy(vt<0)=vt(vt<0).*dEdym(vt<0)*dt;
    extraterm=max(min((fric-advecx-advecy-DFlux),5*par.cp),-5*par.cp);

    % Temperature/Enthalpy solution
    repH2=repmat(H,[1,1,ctr.kmax]); repH=max(repH2,1e-8);
    % centred difference for diffusion & upstream difference for advection
    atp=dt*((2*kdif*par.secperyear./(repH.^2.*dzm.*dzc))-wt./(dzm.*repH)); % E(k-1)
    btp=1+dt*((2*kdif*par.secperyear./(repH.^2.*dzp.*dzm))-wt./(dzm.*repH)); % E(k)
    ctp=dt*((2*kdif*par.secperyear./(repH.^2.*dzp.*dzc))); % E(k+1)

    % basal BC at the CTS
    ftp=ones(ctr.imax,ctr.jmax,ctr.kmax);
    gtp=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    gtp(CTSm==1)=Epmp(CTSm==1);
    ftp(CTSm==1)=0;

    % Recalculate enthalpy profiles for ice columns possessing a CTS
    maskE=Ht>0|CTSm(:,:,ctr.kmax)==1; % base is temperate
    for k=ctr.kmax-1:-1:2
        mask_calc=(k<idx);
        ftp(:,:,k)=mask_calc.*(atp(:,:,k)./(btp(:,:,k)-ctp(:,:,k).*ftp(:,:,k+1)))+...
           ~mask_calc.*0;
        gtp(:,:,k)=mask_calc.*((E(:,:,k)+extraterm(:,:,k)+ ...
          ctp(:,:,k).*gtp(:,:,k+1))./(btp(:,:,k)-ctp(:,:,k).*ftp(:,:,k+1)))+...
           ~mask_calc.*Epmp(:,:,k);
    end
    for k=2:ctr.kmax
        % mask
        mask_cold=(k<idx)&maskE==1;
        mask_temp=(k>idx)&maskE==1;
        mask_CTS=(k==idx);

        % update enthalpy profiles in cold ice above the CTS (mask_cold)
        % keep prior guess enthalpy profiles in temperate ice (mask_temp)
        E(:,:,k)=mask_cold.*(E(:,:,k-1).*ftp(:,:,k)+gtp(:,:,k))+mask_CTS.*Epmp(:,:,k)+mask_temp.*E0(:,:,k);
        tmp(:,:,k)=mask_cold.*((E(:,:,k)./par.cp)+par.Tref)+mask_CTS.*Tpmp(:,:,k)+mask_temp.*tmp0(:,:,k);
        wat(:,:,k)=(mask_cold|mask_CTS).*0+mask_temp.*wat0(:,:,k);
    end

    % Keep prior guess enthalpy profiles for cold ice columns
    E(repmat(maskE,[1 1 ctr.kmax])==0)=E0(repmat(maskE,[1 1 ctr.kmax])==0);
    tmp(repmat(maskE,[1 1 ctr.kmax])==0)=tmp0(repmat(maskE,[1 1 ctr.kmax])==0);
    wat(repmat(maskE,[1 1 ctr.kmax])==0)=0;
    
    % use analytical solution at the domain boundary when using basins
    nMASK=zeros(size(MASK));
    if ctr.basin==1
        nMASK(1,:)=1;
        nMASK(ctr.imax,:)=1;
        nMASK(:,1)=1;
        nMASK(:,ctr.jmax)=1;
    end
    if ctr.mismip>0
        E(:,1,:)=E(:,3,:);
        E(:,end,:)=E(:,end-1,:);
        E(1,:,:)=E(3,:,:);
        E(end,:,:)=E(end-2,:,:);
    end

    % Correction for unstable temperature profiles
    % find an anomaly where temperature decreases with depth
    [ipos,jpos]=find((tmp(:,:,ctr.kmax-1)-tmp(:,:,ctr.kmax))>3 | ...
        (tmp(:,:,ctr.kmax-2)-tmp(:,:,ctr.kmax-1))>3 | ...
        (tmp(:,:,ctr.kmax-3)-tmp(:,:,ctr.kmax-2))>3);
    nMASK(sub2ind(size(nMASK),ipos,jpos))=1;

    % find an anomaly at the CTS or in the temperate layer
    [ipos,jpos]=find(Ptv==1 & E(:,:,ctr.kmax)<Epmp(:,:,ctr.kmax));
    nMASK(sub2ind(size(nMASK),ipos,jpos))=1;
    % apply linear temperature profile when MASK=0 & unstable
    TBshelf=par.T0+repmat(min(par.Toi+DeltaTo(cnt)- ...
        0.12e-3*par.rho*H/par.rhow,0),[1,1,ctr.kmax]);
    nMASK(nMASK==1 & MASK==0)=2;
    nMASKz=repmat(nMASK,[1,1,ctr.kmax]);
    repTs=repmat(Ts,[1,1,ctr.kmax]);
    Tgrad=-G/par.K; 
    repTgrad=repmat(Tgrad,[1,1,ctr.kmax]);
    l=sqrt(2*par.kdif*(H+1e-8)./max(Mb,1e-8)*par.secperyear); 
    repl=repmat(l,[1,1,ctr.kmax]);        
    tmpb=repTs+sqrt(pi)*0.5*repl.*repTgrad.*(erf((1-repz).*repH./repl)-erf(repH./repl))+par.T0;
    tmp(nMASKz==1)=tmpb(nMASKz==1);
    Tshelf=repTs+par.T0-(repTs+par.T0-TBshelf).*repz;
    tmp(nMASKz==2)=Tshelf(nMASKz==2);  

    % Correction for pmp
    tmp(tmp>par.T0-Tp)=par.T0-Tp(tmp>par.T0-Tp);
    mintemp=min(min(Ts))+min(DeltaT)+par.T0-5;
    tmp(tmp<mintemp)=mintemp;
    % Convert corrected temperature to enthalpy
    E=(nMASKz>=1).*((tmp-par.Tref)*par.cp) + (nMASKz<1).*E;
    % Update water content
    wat=(nMASKz>=1).*max(0,(E-Epmp)./par.Latent) + (nMASKz<1).*wat;
    % Enthalpy cannot exceed enthalpy of liquid water
    Eliq=Epmp+par.Latent;
    E(E>=Eliq)=Eliq(E>=Eliq);
    E(E<((mintemp-par.Tref)*par.cp))=(mintemp-par.Tref)*par.cp;
    
    if ctr.drain>0    
        % any water exceeding the prescribed threshold of 001 (= 1%) is 
        % drained instantaneously to the bed (Greve, 1997; Greve & Blatter, 2016)
        wat_max=0.01;
        Dw = max(wat-wat_max,0.0);

        % ensure fraction drained does not exceed the difference (PISM)
        dz=repH.*dzm; 
        Dfw=zeros(size(wat));
        Dfw(wat > wat_max) = Dw(wat > wat_max).*dt;
        Dfw(wat > wat_max) = min(Dfw(wat > wat_max), wat(wat > wat_max) - wat_max);
        Dw(wat > wat_max) = Dfw(wat > wat_max) .* dz(wat > wat_max);

        % Drain to bed (Wang et al, 2020)
        Dbw=sum(Dw,3);
        Dbw=Dbw/dt; % m/an

        % Correction on water content
        wat(wat>wat_max)=wat_max;
    end
    Dbw(Dbw>1)=1; % correction for unrealistic high values
    Dbw(Dbw<0)=0; 

    % get enthalpy again (to be consistent with new water content)
    E=(tmp<Tpmp).*((tmp-par.Tref)*par.cp) + (tmp>=Tpmp).*((wat.*par.Latent)+Epmp);

end


