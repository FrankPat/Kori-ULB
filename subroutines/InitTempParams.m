function [tmp,Tb,zeta,dzc,dzp,dzm,E,Epmp,wat]=InitTempParams(ctr,par,tmp,Ts,H,E,wat)

% Kori-ULB
% Initialization of vertical derivatives for temperature profiles in the 3d
% temperature calculation

    zeta=CalculateZeta(ctr.kmax,0.015);
    dzc=circshift(zeta,[0 -1])-circshift(zeta,[0 1]);
    dzp=circshift(zeta,[0 -1])-zeta;
    dzm=zeta-circshift(zeta,[0 1]);
    dzc=repmat(reshape(dzc,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    dzp=repmat(reshape(dzp,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    dzm=repmat(reshape(dzm,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);

    if islogical(tmp)==1
        tmp=repmat(Ts+par.T0,[1,1,ctr.kmax]);
    end
    repH=repmat(H,[1,1,ctr.kmax]);
    repz=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    Tp=par.pmp*repH.*repz;
    % Correction for pmp
    tmp(tmp>par.T0-Tp)=par.T0-Tp(tmp>par.T0-Tp);
    Tb=tmp(:,:,ctr.kmax)-par.T0;
    if ctr.Enthalpy==1
        wat=zeros(size(tmp));
        E=par.cp*(tmp-par.Tref);
        Epmp=par.cp*(Tp+par.T0-par.Tref);
    else
        Epmp=false;
    end
end


