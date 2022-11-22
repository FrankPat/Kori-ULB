function [FMB]=VerticalFaceMelt(ctr,par,SLR,B,Melt,MASK,glMASK,he)

% Kori-ULB
% Calculation of ice shelf frontal melt

    Dz=SLR-B;
    Dz(MASK==0)=par.rho*he(MASK==0)/par.rhow;
    FMB=max(Dz.*Melt/ctr.delta,0);
    FMB(glMASK~=5)=0;
    
end

        
