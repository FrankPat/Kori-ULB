function flux=VariabFlux(variab,H,MASK,seaice)

% Kori-ULB
% Calculate mass balance component fluxes

    flux=sum(variab(MASK==1 | (MASK==0 & H>seaice)));
    
end


