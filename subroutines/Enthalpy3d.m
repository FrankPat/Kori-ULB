function [E,Epmp,wat,CTSm,CTSp,Bmelt,Dbw,Dfw,Hw,Ht,tmp]= ...
    Enthalpy3d(par,ctr,E,Mb,Ts,G,A,H, ...
    pxy,dt,gradsx,gradsy,gradHx,gradHy,gradxy,taudxy, ...
    udx,udy,ub,ubx,uby,zeta,dzc,dzm,dzp,DeltaT,MASK, ...
    Bmelt,Epmp,CTSm,CTSp,Hw,Ht,Dbw,Dfw,cnt)

% Kori-ULB
% 3d englacial enthalpy calculation in ice sheet and ice shelves
 
    % ice diffusivity
    kdif=par.kdif*(1.+(1e-5-1.)*heaviside(E-Epmp));
    repz=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]); 
    repH=max(repmat(H,[1,1,ctr.kmax]),1e-8);
    repH2=repmat(H,[1,1,ctr.kmax]);

    % surface BC
    E(:,:,1)=par.cp*(Ts+par.T0-par.Tref);
    
    % Drain excess englacial water computed at previous time step
    if ctr.drain>0
        DFlux=(par.rhof./par.rho).*Dfw.*par.Latent; 
        Bmelt=Bmelt+Dbw;
    else
        DFlux=zeros(size(E));
    end
    
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
    
    % vertical velocity according to Pattyn (2010) with Lliboutry shape
    % function
    wshllib=1-(pl+2).*zl./(pl+1)+1./(pl+1).*zl.^(pl+2);
    w=repmat(-Mb,[1,1,ctr.kmax]).*wshllib-max(Bmelt,-5e-2)+ut.* ...
        (repmat(gradsx,[1,1,ctr.kmax])-zl.*repmat(gradHx,[1,1,ctr.kmax])) ...
        +vt.*(repmat(gradsy,[1,1,ctr.kmax])-zl.*repmat(gradHy,[1,1,ctr.kmax]));
    
    % Internal/strain heating [CHECK THIS BECAUSE CHANGED]
    dudz=repmat(2*A.*taudxy.^par.n.*H.*(pxy+2)/(par.n+2),[1,1,ctr.kmax]).*repz.^pl; 
    fric=par.rho.*par.g.*repmat(sqrt(gradxy),[1,1,ctr.kmax]).*repz.*dudz.*dt./par.rho;
    repmask=repmat(MASK,[1,1,ctr.kmax]);
    fric(repmask==0)=0; % no frictional heat in ice shelves
    extraterm=max(min((fric-advecx-advecy-DFlux),5*par.cp),-5*par.cp); 
   
    % central difference for diffusion & upstream difference for advection
    atp=dt*((2*kdif*par.secperyear./(repH.^2.*dzm.*dzc))-w./(dzm.*repH)); % E(k-1)
    btp=1+dt*((2*kdif*par.secperyear./(repH.^2.*dzp.*dzm))-w./(dzm.*repH)); % E(k)
    ctp=dt*((2*kdif*par.secperyear./(repH.^2.*dzp.*dzc))); % E(k+1)

    % correction for discontinuity at the CTS
    Kb=kdif;
    Ka=kdif;
    Kb(:,:,2:end-1)=hamean(kdif(:,:,1:end-2),kdif(:,:,2:end-1));
    Ka(:,:,2:end-1)=hamean(kdif(:,:,3:end),kdif(:,:,2:end-1));
    
    % Correction CTS: upstream differences for advection
    atp(CTSp==1)=dt*((2*(Kb(CTSp==1))*par.secperyear./(repH(CTSp==1).^2.* ...
        dzm(CTSp==1).*dzc(CTSp==1)))-w(CTSp==1)./(dzm(CTSp==1).*repH(CTSp==1))); % E(k-1)
    btp(CTSp==1)=1+dt*(((Ka(CTSp==1)+Kb(CTSp==1))*par.secperyear./ ...
        (repH(CTSp==1).^2.*dzp(CTSp==1).*dzm(CTSp==1)))-w(CTSp==1)./ ...
        (dzm(CTSp==1).*repH(CTSp==1))); % E(k)
    ctp(CTSp==1)=dt*((2*(Ka(CTSp==1))*par.secperyear./(repH(CTSp==1).^2.* ...
        dzp(CTSp==1).*dzc(CTSp==1)))); % E(k+1)
    atp(CTSm==1)=dt*((2*(Kb(CTSm==1))*par.secperyear./(repH(CTSm==1).^2.* ...
        dzm(CTSm==1).*dzc(CTSm==1)))-w(CTSm==1)./(dzm(CTSm==1).*repH(CTSm==1))); % E(k-1)
    btp(CTSm==1)=1+dt*(((Ka(CTSm==1)+Kb(CTSm==1))*par.secperyear./ ...
        (repH(CTSm==1).^2.*dzp(CTSm==1).*dzm(CTSm==1)))-w(CTSm==1)./ ...
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
    EBshelf=par.cp*((par.T0+repmat(min(par.Toi+ctr.meltfactor*DeltaT(cnt)- ...
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
    
    % CTS position
    CTS=zeros(ctr.imax,ctr.jmax,ctr.kmax); CTSm=CTS; CTSp=CTS;
    % Temperate layer thickness
    Ht=zeros(ctr.imax,ctr.jmax);
    Hx=zeros(ctr.imax,ctr.jmax);

    % only keep the first value for the CTS where the condition apply
    % --> assumption that only one CTS exists for a given ice column
    first_one = false(size(CTS, 1), size(CTS, 2));
    second_one = false(size(CTS, 1), size(CTS, 2));

    for k=ctr.kmax-1 :-1:2
        CTSm(:,:,k)=E(:,:,k) >= Epmp(:,:,k) & E(:,:,k-1) < Epmp(:,:,k-1) & ~first_one; % first grid point where E>=Epmp
        first_one = first_one | CTSm(:,:,k);
        CTSp(:,:,k)=E(:,:,k) < Epmp(:,:,k) & E(:,:,k+1) >= Epmp(:,:,k+1) & ~second_one; % last grid point where E<Epmp
        second_one = second_one | CTSp(:,:,k);
    end
    for i=1:ctr.imax
        for j=1:ctr.jmax
            if MASK(i,j)>0 && H(i,j)>0
                % Find index of the CTS
                idx = find(CTSm(i,j,:)==1);
                if E(i,j,ctr.kmax)>=Epmp(i,j,ctr.kmax) && E(i,j,ctr.kmax-1)>=Epmp(i,j,ctr.kmax-1)
                    % temperate layer thickness & CTS
                    if idx~=0 
                        Ht(i,j) = H(i,j).*(1-zeta(idx)); 
                    end
                elseif E(i,j,ctr.kmax)>=Epmp(i,j,ctr.kmax) && E(i,j,ctr.kmax-1)<Epmp(i,j,ctr.kmax-1)
                    CTSp(i,j,:)=0;
                    CTSm(i,j,:)=0;
                    Ht(i,j)=0;
                    % Correction when CTS is at the bed-ice interface
                    CTSm(i,j,ctr.kmax)=1;
                    E(i,j,ctr.kmax)=Epmp(i,j,ctr.kmax);
                else % no CTS if cold layer of ice above the bed
                     CTSp(i,j,:)=0;
                     CTSm(i,j,:)=0;
                     Ht(i,j)=0;
                end
                if Ht(i,j)>0  % positive thickness of temperate ice 
                    for k=2:ctr.kmax 
                         % Find if ice is cold below the CTS while it should be temperate
                         if k>idx && E(i,j,k)<Epmp(i,j,k)
                         E(i,j,k)=Epmp(i,j,k);
                         end
                         % Find if temperate conditions are present above the CTS while it should be cold
                         if k<idx && E(i,j,k)>Epmp(i,j,k)
                         E(i,j,k)=Epmp(i,j,k);
                         end
                         % Find anomaly in the enthalpy field in cold ice around the CTS
                         if (k<(idx+5) && k>(idx-5)) && abs(E(i,j,k-1)-E(i,j,k))>=par.cp
                            Hx(i,j)=1;
                         end
                    end
                end
            end
        end
    end
    % Temperature field  
    Tp=par.pmp*repH2.*repz;
    Tpmp=par.T0-Tp;
    tmp=min((E/par.cp)+par.Tref,Tpmp);
    tmp(E>=Epmp)=Tpmp(E>=Epmp); % temperate ice
    
    % Enthalpy correction
    Epmp=par.cp*(Tpmp-par.Tref);

    % Water Content 
    wat=max((E-Epmp)./par.Latent,0); % Aschwanden et al, 2012
    wat(E<Epmp)=0; % no water content in cold ice
    wat(tmp<Tpmp)=0; % avoid truncature errors with Epmp
    
    % Correction for unstable temperature profiles
    % find an anomaly where temperature decreases with depth
    [ipos,jpos]=find((tmp(:,:,ctr.kmax-1)-tmp(:,:,ctr.kmax))>3 | ...
        (tmp(:,:,ctr.kmax-2)-tmp(:,:,ctr.kmax-1))>3 | ...
        (tmp(:,:,ctr.kmax-3)-tmp(:,:,ctr.kmax-2))>3);
    nMASK(sub2ind(size(nMASK),ipos,jpos))=1;

    %OR: find an anomaly at the CTS or in the temperate layer
    [ipos,jpos]=find(Ptv==1 & E(:,:,ctr.kmax)<Epmp(:,:,ctr.kmax));
    nMASK(sub2ind(size(nMASK),ipos,jpos))=1;
    % apply linear temperature profile when MASK=0 & unstable
    TBshelf=par.T0+repmat(min(par.Toi+ctr.meltfactor*DeltaT(cnt)- ...
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
        Dfw(wat > 0.01) = Dw(wat > 0.01).*dt;
        Dfw(wat > 0.01) = min(Dfw(wat > 0.01), wat(wat > 0.01) - 0.01);
        Dw(wat > 0.01) = Dfw(wat > 0.01) .* dz(wat > 0.01);

        % Drain to bed (Wang et al, 2020)
        Dbw=sum(Dw,3);
        Dbw=Dbw/dt; % m/an

        % Correction on water content
        wat(wat>wat_max)=wat_max;
    end
    Dbw(Dbw>1)=1; % correction for unrealistic high values

    % YELMO code: get enthalpy again (to be consistent with new water content)
    E=(tmp<Tpmp).*((tmp-par.Tref)*par.cp) + (tmp>=Tpmp).*((wat.*par.Latent)+Epmp);

    % Basal melting 
    % UPDATE BASAL CONDITIONS
    % Cold base and dry
    Cld=(E(:,:,ctr.kmax)<Epmp(:,:,ctr.kmax) & Hw==0);
    % Cold base and wet
    Blw=(E(:,:,ctr.kmax)<Epmp(:,:,ctr.kmax) & Hw>0);
    % Temperate base
    Abv=((E(:,:,ctr.kmax)>=Epmp(:,:,ctr.kmax)) & Ht==0 & Hw>0);
    % Temperate layer
    Ptv=((E(:,:,ctr.kmax)>=Epmp(:,:,ctr.kmax)) & Ht>0 & Hw>0);
    repE=E(:,:,ctr.kmax);
    repEpmp=Epmp(:,:,ctr.kmax); 
    repE(Blw==1)=repEpmp(Blw==1);
    E(:,:,ctr.kmax)=repE;
    repE(Abv==1)=repEpmp(Abv==1);
    E(:,:,ctr.kmax)=repE;
    % Basal Heat sources
    Fb=ub.*taudxy/par.secperyear; % frictional heating
    BasalHeat=G+Fb; % [W/m2]=[J s-1 m-2]
    % Enthalpy gradient at the base
    dEdz=(E(:,:,ctr.kmax)-E(:,:,ctr.kmax-1))./(max(H,1e-8).*dzm(:,:,ctr.kmax)); % [J kg-1 m-1]
    % Basal (non advective) Heat Flux
    Qq = (par.K*dEdz/par.cp).*(Ptv==0) + (par.K*1e-5*dEdz/par.cp).*(Ptv==1);
    % Basal melting (m/an)
    % Aschwanden et al. (2012)
    Bmelt=min(1,max(-1,((BasalHeat-Qq).*par.secperyear./((1-wat(:,:,ctr.kmax))*par.Latent*par.rho)))); 
    Bmelt(MASK~=1)=0;
    Bmelt(H==0)=0; 
    Bmelt(Cld==1)=0; 
    Bmelt(Bmelt<0 & Hw<=0)=0; 

    % Water layer thickness 
    % with a constant drainage rate of 1 mm/a
    % more stable if only applied when considerable Hw
    Hw=max(0,min(Hw+(Bmelt+Dbw-(par.Cdr.*(Hw>=par.Wmax)))*dt,par.Wmax));
    Hw(Hw<0)=0;    
    Hw(MASK~=1)=0; % only for the grounded ice sheet
end


