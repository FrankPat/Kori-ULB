function [FMB,FMR]=VerticalFaceMeltVio(ctr,par,SLR,B,Melt,MASK,glMASK,he)

% Kori-ULB
% Calculation of ice shelf frontal melt

if ctr.FrontalMelt==1

    Dz=SLR-B;
    Dz(MASK==0)=par.rho*he(MASK==0)/par.rhow;
    FMB=max(Dz.*Melt/ctr.delta,0);
    FMB(glMASK~=5)=0;
    FMR=FMB.*ctr.delta./he; % retreat rate due to vertical face melt

else

    FMR=zeros(ctr.imax,ctr.jmax);
    FMB=FMR;

end

end
        
