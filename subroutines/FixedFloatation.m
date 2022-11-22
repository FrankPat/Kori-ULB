function [sn,HB]=FixedFloatation(par,B,SLR,H,MASK)

% Kori-ULB
% height of the bottom of ice shelves based on floatation for a fixed MASK

    HB=B;
    HB(MASK==0)=SLR(MASK==0)-par.rho/par.rhow*H(MASK==0);
    HB(HB<B)=B(HB<B);
    sn=H+HB;
    
end


