function [A,Ax,Ay,Ad,Tbc]=ThermoCoupling(ctr,par,Tb,H,bMASK,bMASKm, ...
    bMASKx,bMASKy,wat)

% Kori-ULB
% Thermomechanical coupling using Arrhenius relationship

    if ctr.Tcalc>0
        Tbc=Tb+par.pmp*H;
    else
        Tbc=false;
    end
    A=zeros(ctr.imax,ctr.jmax)+ctr.Ao;
    if ctr.Tcalc==2
        A=0.5*par.atune*((Tbc<-6.5)*par.a1+(Tbc>=-6.5)*par.a2).* ...
            exp(((Tbc<-6.5)*par.Q1+(Tbc>=-6.5)*par.Q2)./par.R.* ...
            (1./(par.T0-par.pmp*H)-1./(Tb+par.T0)));
        if wat
            A=A.*(1+181.25*min(wat(:,:,ctr.kmax),0.03)); % Adding effect of water content with enthalpy
        end
        [Ax,Ay,Ad]=StaggeredA(A);
    else
        Ax=A;
        Ay=A;
        Ad=A;
    end
    if ctr.basin==1
        A(bMASK==1)=par.A0; % A on h-grid
        Ax(bMASKx==1)=par.A0;
        Ay(bMASKy==1)=par.A0;
        Ad(bMASKm==1)=par.A0; % A on d-grid
    end

end


