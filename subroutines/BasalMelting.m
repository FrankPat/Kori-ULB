function Bmelt=BasalMelting(ctr,par,G,taudxy,ub,H,tmp,dzm,MASK)

% Kori-ULB
% Basal melting underneath the grounded ice sheet based on geothermal heat
% flux and basal sliding
    
    Bgrad=(G+taudxy.*ub/par.secperyear).*H/par.K;
    Tgrad=(tmp(:,:,ctr.kmax)-tmp(:,:,ctr.kmax-1))./dzm(:,:,ctr.kmax);
    Bmelt=min(1,max(1e-8,(Bgrad-Tgrad)*par.K*par.secperyear./ ...
        (par.rho*par.Latent*H)));
    Bmelt(MASK==0)=0;
    
end


