function mbcomp=BasinFlux(ctr,par,acc,Smelt,runoff,rain,Mb,Pr, ...
    H,Hn,Bmelt,Melt,CMB,FMB,MASK,bMASK,mbcomp,B,Bn,SLR,ZB)

    % (1) surface mass balance
    % (2) dhdt
    % (3) accumulation
    % (4) surface melt
    % (5) runoff
    % (6) rain
    % (7) basal melt grounded ice
    % (8) sub-shelf melt
    % (9) calving flux
    % (10) ice shelf frontal melt
    % (11) dynamical component
    
    if ctr.basin==1
        MASK=MASK+bMASK; % MASK=2 for outside basin domain and not accounted for
    end
    mbcomp(1,:)=VariabFluxBasin(max(-H,Mb),H,MASK,par.SeaIceThickness,ZB); % smb
    for i=1:max(ZB(:))
        mbcomp(2,i)=(sum(Hn(Hn>par.SeaIceThickness & ZB==i))-sum(H(H>par.SeaIceThickness & ZB==i)))/ctr.dt; % dhdt
    end
    if ctr.PDDcalc==1
        mbcomp(3,:)=VariabFluxBasin(acc,H,MASK,par.SeaIceThickness,ZB); % acc
        mbcomp(4,:)=VariabFluxBasin(min(H+acc,Smelt),H,MASK,par.SeaIceThickness,ZB); % Smelt
        mbcomp(5,:)=VariabFluxBasin(min(H+acc,runoff),H,MASK,par.SeaIceThickness,ZB); % runoff
        mbcomp(6,:)=VariabFluxBasin(rain,H,MASK,par.SeaIceThickness,ZB); % rain
    else
        mbcomp(3,:)=VariabFluxBasin(Pr,H,MASK,par.SeaIceThickness,ZB); % acc
        mbcomp(5,:)=VariabFluxBasin(min(H+Pr,runoff),H,MASK,par.SeaIceThickness,ZB); % runoff
    end
    if ctr.Tcalc>=1
        mbcomp(7,:)=VariabFluxBasin(min(H+max(-H,Mb),Bmelt),H,MASK,par.SeaIceThickness,ZB); % Bmelt
    end
    if ctr.meltfunc>=1 && ctr.glMASKexist==1
        mbcomp(8,:)=VariabFluxBasin(min(H+max(-H,Mb),Melt),H,MASK,par.SeaIceThickness,ZB);
    end
    if ctr.calving>=1 && ctr.shelf==1
        mbcomp(9,:)=VariabFluxBasin(min(H+max(-H,Mb),CMB),H,MASK,par.SeaIceThickness,ZB); % CMB
        mbcomp(10,:)=VariabFluxBasin(min(H+max(-H,Mb),FMB),H,MASK,par.SeaIceThickness,ZB); % FMB
    end
    mbcomp(11,:)=mbcomp(1,:)-mbcomp(2,:)-mbcomp(7,:)-mbcomp(8,:);

    % Grounded ice sheet components
    % (12) surface mass balance
    % (13) dhdt
    % (14) accumulation
    % (15) surface melt
    % (16) runoff
    % (17) rain
    % (18) basal melt grounded ice
    % (19) grounding line flux
    % (20) Rate of VAF
    % (21) dynamical VAF component
    
    MASK(MASK==0)=NaN;
    mbcomp(12,:)=VariabFluxBasin(max(-H,Mb),H,MASK,par.SeaIceThickness,ZB); % smb
    for i=1:max(ZB(:))
        mbcomp(13,i)=(sum(Hn(MASK==1 & ZB==i))-sum(H(MASK==1 & ZB==i)))/ctr.dt; % dhdt
    end
    if ctr.PDDcalc==1
        mbcomp(14,:)=VariabFluxBasin(acc,H,MASK,par.SeaIceThickness,ZB); % acc
        mbcomp(15,:)=VariabFluxBasin(min(H+acc,Smelt),H,MASK,par.SeaIceThickness,ZB); % Smelt
        mbcomp(16,:)=VariabFluxBasin(min(H+acc,runoff),H,MASK,par.SeaIceThickness,ZB); % runoff
        mbcomp(17,:)=VariabFluxBasin(rain,H,MASK,par.SeaIceThickness,ZB); % rain
    else
        mbcomp(14,:)=VariabFluxBasin(Pr,H,MASK,par.SeaIceThickness,ZB); % acc
        mbcomp(16,:)=VariabFluxBasin(min(H+Pr,runoff),H,MASK,par.SeaIceThickness,ZB); % runoff
    end
    if ctr.Tcalc>=1
        mbcomp(18,:)=VariabFluxBasin(min(H+max(-H,Mb),Bmelt),H,MASK,par.SeaIceThickness,ZB); % Bmelt
    end
    mbcomp(19,:)=mbcomp(12,:)-mbcomp(13,:)-mbcomp(18,:);
    
    VAF=max(0,H+min(B-SLR,0)*(par.rhow/par.rho));
    VAFn=max(0,Hn+min(Bn-SLR,0)*(par.rhow/par.rho));
    for i=1:max(ZB(:))
        mbcomp(20,i)=(sum(sum(VAFn(ZB==i)-VAF(ZB==i))))/ctr.dt;
    end
    mbcomp(21,:)=mbcomp(12,:)-mbcomp(20,:)-mbcomp(18,:);
    
    % convert all components to Gt/a (water equivalent)
    mbcomp=mbcomp*ctr.delta^2*par.rho/1e12;
end

