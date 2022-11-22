function [tmp,Tb,zeta,dzc,dzp,dzm]=InitTempParams(ctr,par,tmp,Ts)

% Kori-ULB
% Initialization of vertical derivatives for temperature profiles in the 3d
% temperature calculation

    if islogical(tmp)==1
        tmp=repmat(Ts+par.T0,[1,1,ctr.kmax]);
        Tb=Ts;
    else
        Tb=tmp(:,:,ctr.kmax)-par.T0;
    end
    zeta=CalculateZeta(ctr.kmax,0.015);
    dzc=circshift(zeta,[0 -1])-circshift(zeta,[0 1]);
    dzp=circshift(zeta,[0 -1])-zeta;
    dzm=zeta-circshift(zeta,[0 1]);
    dzc=repmat(reshape(dzc,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    dzp=repmat(reshape(dzp,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    dzm=repmat(reshape(dzm,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);

end


