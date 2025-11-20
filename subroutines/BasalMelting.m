function [Bmelt,Hw]=BasalMelting(ctr,par,G,taudxy,ub,H,tmp,dzm,MASK, ...
    Hw,Ht,E,Epmp,wat,Dbw)

% Kori-ULB
% Basal melting underneath the grounded ice sheet based on geothermal heat
% flux and basal sliding
    
    if ctr.Enthalpy>0
        % UPDATE BASAL CONDITIONS
        % Cold base and dry
        Cld=(E(:,:,ctr.kmax)<Epmp(:,:,ctr.kmax) & Hw==0);
        % Cold base and wet
        Blw=(E(:,:,ctr.kmax)<Epmp(:,:,ctr.kmax) & Hw>0);
        % Temperate base
        Abv=((E(:,:,ctr.kmax)>=Epmp(:,:,ctr.kmax)) & Ht==0 & Hw>0);
        % Temperate layer
        Ptv=((E(:,:,ctr.kmax)>=Epmp(:,:,ctr.kmax)) & Ht>0 & Hw>0);
        repE=E(:,:,ctr.kmax);
        repEpmp=Epmp(:,:,ctr.kmax); 
        repE(Blw==1)=repEpmp(Blw==1);
        E(:,:,ctr.kmax)=repE;
        repE(Abv==1)=repEpmp(Abv==1);
        E(:,:,ctr.kmax)=repE;
        % Basal Heat sources
        Fb=ub.*taudxy/par.secperyear; % frictional heating
        BasalHeat=G+Fb; % [W/m2]=[J s-1 m-2]
        % Enthalpy gradient at the base
        dEdz=(E(:,:,ctr.kmax)-E(:,:,ctr.kmax-1))./(max(H,1e-8).* ...
            dzm(:,:,ctr.kmax)); % [J kg-1 m-1]
        % Basal (non advective) Heat Flux
        Qq = (par.K*dEdz/par.cp).*(Ptv==0) + (par.K0*dEdz).*(Ptv==1);
        % Basal melting (m/a) based on Aschwanden et al. (2012)
        Bmelt=min(1,max(-1,((BasalHeat-Qq).*par.secperyear./ ...
            ((1-wat(:,:,ctr.kmax))*par.Latent*par.rho))));
        Bmelt(MASK~=1|H==0)=0; % only for frounded ice
        Bmelt(Cld==1)=0; % no basal melting if cold basal ice
        Bmelt(Bmelt<0 & Hw<=0)=0; % no refreezing in the absence of water

        % Water layer thickness 
        % with a constant drainage rate of 1 mm/a
        % more stable if only applied when considerable Hw
        dt=ctr.dt.*ctr.intT;
        Hw=max(0,min(Hw+(Bmelt+Dbw-(par.Cdr.*(Hw>=par.Wmax)))*dt,par.Wmax));
        Hw(MASK~=1)=0; % only for the grounded ice sheet
    else
        Bgrad=(G+taudxy.*ub/par.secperyear).*H/par.K;
        Tgrad=(tmp(:,:,ctr.kmax)-tmp(:,:,ctr.kmax-1))./dzm(:,:,ctr.kmax);
        Bmelt=min(1,max(1e-8,(Bgrad-Tgrad)*par.K*par.secperyear./ ...
            (par.rho*par.Latent*H)));
        Bmelt(MASK==0)=0;
        Hw=false;
    end
    
end


