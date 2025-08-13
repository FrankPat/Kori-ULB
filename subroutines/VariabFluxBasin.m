function [flux]=VariabFluxBasin(variab,H,MASK,seaice,ZB)

% Kori-ULB
% Flux parameter per basin

    for i=1:max(ZB(:))
        flux(i)=sum(variab(ZB==i & (MASK==1 | (MASK==0 & H>seaice))));
    end
    
end


