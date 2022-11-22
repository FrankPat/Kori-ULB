function [HAF,MASK,HB,sn]=Floatation(par,B,SLR,H,MASK)

% Kori-ULB
% Determination of MASK and height of the bottom of ice shelves based on
% floatation. Height above buoyancy is also returned.

    H(MASK==0 & H<=par.SeaIceThickness)=0;
    HAF=B-SLR+H*par.rho/par.rhow;
    MASK(HAF<0)=0;
    MASK(HAF>=0)=1;
    HB=max(SLR-par.rho/par.rhow*H,B);
    sn=HB+H;
    
end


