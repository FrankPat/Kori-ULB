function [T0o,S0o,zb]=OceanVarBox(numsh,To,So,ShelfN,HB)

% Kori-ULB
% Temperature and salinity in boxes of the PICO model

    T0o=zeros(numsh,1); % mean box T0
    S0o=zeros(numsh,1); % mean box S0
    zb=zeros(numsh,1); % mean box HB

    for k=1:numsh
        T0o(k)=mean(To(ShelfN==k));
        S0o(k)=mean(So(ShelfN==k));
        zb(k)=mean(HB(ShelfN==k));
    end
end


