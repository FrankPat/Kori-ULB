function [tmp,E,Epmp,wat]=InitTemp3d(G,taudxy,ub,ud,par,H,Mb,zeta,ctr,Ts,MASK,DeltaT)

% Kori-ULB
% Temperature field initialization with analytical solution for temperature
% profiles; different analytical solution for ice shelves

    uvel=ub+par.udfrac*ud;
    Tgrad=-(G+taudxy.*uvel/par.secperyear)/par.K;
    l=sqrt(2*par.kdif*(H+1e-8)./max(Mb,1e-8)*par.secperyear);
    repl=repmat(l,[1,1,ctr.kmax]);
    repH2=repmat(H,[1,1,ctr.kmax]);
    repH=max(repH2,1e-8);
    repTs=repmat(Ts,[1,1,ctr.kmax]);
    repTgrad=repmat(Tgrad,[1,1,ctr.kmax]);
    repz=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    tmp=repTs+sqrt(pi)*0.5*repl.*repTgrad.* ...
        (erf((1-repz).*repH./repl)-erf(repH./repl))+par.T0;
    Tp=par.pmp*repH2.*repz;

    % Ice shelf and open ocean
    repmask=repmat(MASK,[1,1,ctr.kmax]);
    repM=repmat(max(Mb,1e-5),[1,1,ctr.kmax]);
    TBshelf=repmat(min(par.Toi+ctr.meltfactor*DeltaT- ...
        0.12e-3*par.rho*H/par.rhow,0),[1,1,ctr.kmax]);
    c2=exp(repM.*repH/(par.kdif*par.secperyear));
    c1=exp(repM.*repz.*repH/(par.kdif*par.secperyear));
    shelftmp=((repTs-TBshelf).*c1+TBshelf-repTs.*c2)./(1-c2)+par.T0;
    shelftmp(repH<=1)=TBshelf(repH<=1)+par.T0;
    tmp(repmask==0)=shelftmp(repmask==0);

    % Correction for pmp
    tmp(tmp>par.T0-Tp)=par.T0-Tp(tmp>par.T0-Tp);
    
    % Enthalpy
    if ctr.Enthalpy==1
        E=par.cp.*(tmp-par.Tref);
        Epmp=par.cp.*(par.T0-Tp-par.Tref);
        % Water content
        wat=max(0,(E-Epmp)./par.Latent);
    else
        [E,Epmp,wat]=deal(false);
    end
end


