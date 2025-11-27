function [A,Ax,Ay,Ad,A3d,Tbc]=ThermoCoupling(ctr,par,Tb,H,bMASK,bMASKm, ...
    bMASKx,bMASKy,tmp,zeta,wat)

% Kori-ULB
% Thermomechanical coupling using Arrhenius relationship

    if ctr.Tcalc>0
        Tbc=Tb+par.pmp*H;
    else
        Tbc=false;
    end
    A=zeros(ctr.imax,ctr.jmax)+ctr.Ao;
    if ctr.Tcalc==2
        if ctr.SSA==3
            A3d=zeros(ctr.imax,ctr.jmax,ctr.kmax);
            repz=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
            repH=repmat(H+1e-8,[1,1,ctr.kmax]);
            Tstar=min(tmp+par.pmp*repH.*repz,par.T0);
            A3d(Tstar<=par.T0-10)=par.atune*par.a1D*exp(-par.Q1D./ ...
                (par.R*(Tstar(Tstar<=par.T0-10))));
            A3d(Tstar>par.T0-10)=par.atune*par.a2D*exp(-par.Q2D./ ...
                (par.R*(Tstar(Tstar>par.T0-10))));
            if ctr.Enthalpy>0
                A3d=A3d.*(1+181.25*min(wat,0.03)); % Adding effect of water content with enthalpy
            end
        else
            A3d=false;
        end
        A=0.5*par.atune*((Tbc<-6.5)*par.a1+(Tbc>=-6.5)*par.a2).* ...
            exp(((Tbc<-6.5)*par.Q1+(Tbc>=-6.5)*par.Q2)./par.R.* ...
            (1./(par.T0-par.pmp*H)-1./(Tb+par.T0)));
        if ctr.Enthalpy>0
            A=A.*(1+181.25*min(wat(:,:,ctr.kmax),0.03)); % Adding effect of water content with enthalpy
        end
        [Ax,Ay,Ad]=StaggeredA(A);
    else
        Ax=A;
        Ay=A;
        Ad=A;
        if ctr.SSA==3
            A3d=zeros(ctr.imax,ctr.jmax,ctr.kmax)+ctr.Ao;
        else
            A3d=false;
        end
    end
    if ctr.basin==1
        A(bMASK==1)=par.A0; % A on h-grid
        Ax(bMASKx==1)=par.A0;
        Ay(bMASKy==1)=par.A0;
        Ad(bMASKm==1)=par.A0; % A on d-grid
    end

end


