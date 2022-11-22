function zeta=CalculateZeta(kmax,min_layer)

% Kori-ULB
% Calculation of the vertical scaled coordinate system for temperature
% calculation. The vertical layers are unevenly distributed with smaller
% layers near the bottom of the ice sheet.

    if (kmax-1)*min_layer > 1. % if number of layers is too high for min_layer
        zeta=0:(1/(kmax-1)):1;
    else
        xz=[1 kmax-1 kmax];
        yz=[0 1.-min_layer 1];
        p=polyfit(xz,yz,2);
        nxz=1:1:kmax;
        zeta=p(1)*nxz.^2+p(2)*nxz+p(3);
    end
    zeta(1)=0;
    zeta(kmax)=1;
    
end


