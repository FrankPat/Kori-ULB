function [A,Ax,Ay,Ad]=ThermoCoupling(ctr,par,Tb,Tbc,H,bMASK,bMASKm,bMASKx,bMASKy)

% Kori-ULB
% Thermomechanical coupling using Arrhenius relationship

    A=zeros(ctr.imax,ctr.jmax)+ctr.Ao;
    if ctr.Tcalc==2
        A=0.5*par.atune*par.a1*exp(par.Q1./par.R*(1./(par.T0+par.pmp*H)-1./ ...
            (Tbc+par.T0)));
        A(Tb>-6.5)=0.5*par.atune*par.a2*exp(par.Q2./par.R*(1./(par.T0+ ...
            par.pmp*H(Tb>-6.5))-1./(Tbc(Tb>-6.5)+par.T0)));
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


