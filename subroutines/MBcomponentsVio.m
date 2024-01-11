function mbcomp=MBcomponents(ctr,par,acc,Smelt,runoff,rain,Mb,Pr, ...
    H,Hn,Bmelt,Melt,CMB,FMB,MASK,bMASK,mbcomp,B,Bn,SLR)

% Kori-ULB
% Mass balance components of the ice sheet and ice shelf

% Global ice sheet components (ice sheet + ice shelf)
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
    mbcomp(1)=VariabFlux(max(-H,Mb),H,MASK,par.SeaIceThickness); % smb
    mbcomp(2)=(sum(Hn(Hn>par.SeaIceThickness))-sum(H(H>par.SeaIceThickness)))/ctr.dt;
    if ctr.PDDcalc==1
        mbcomp(3)=VariabFlux(acc,H,MASK,par.SeaIceThickness); % acc
        mbcomp(4)=VariabFlux(min(H+acc,Smelt),H,MASK,par.SeaIceThickness); % Smelt
        mbcomp(5)=VariabFlux(min(H+acc,runoff),H,MASK,par.SeaIceThickness); % runoff
        mbcomp(6)=VariabFlux(rain,H,MASK,par.SeaIceThickness); % rain
    else
        mbcomp(3)=VariabFlux(Pr,H,MASK,par.SeaIceThickness); % acc
        mbcomp(5)=VariabFlux(min(H+Pr,runoff),H,MASK,par.SeaIceThickness); % runoff
    end
    if ctr.Tcalc>=1
        mbcomp(7)=VariabFlux(min(H+max(-H,Mb),Bmelt),H,MASK,par.SeaIceThickness); % Bmelt
    end
    if ctr.meltfunc>=1 && ctr.glMASKexist==1
        mbcomp(8)=VariabFlux(min(H+max(-H,Mb),Melt),H,MASK,par.SeaIceThickness);
    end
    if ctr.calving>=1 && ctr.shelf==1
        mbcomp(9)=VariabFlux(min(H+max(-H,Mb),CMB),H,MASK,par.SeaIceThickness); % CMB
        mbcomp(10)=VariabFlux(min(H+max(-H,Mb),FMB),H,MASK,par.SeaIceThickness); % FMB
    end
    mbcomp(11)=mbcomp(1)-mbcomp(2)-mbcomp(7)-mbcomp(8);
    
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
    mbcomp(12)=VariabFlux(max(-H,Mb),H,MASK,par.SeaIceThickness); % smb
    mbcomp(13)=(sum(Hn(MASK==1))-sum(H(MASK==1)))/ctr.dt;
    if ctr.PDDcalc==1
        mbcomp(14)=VariabFlux(acc,H,MASK,par.SeaIceThickness); % acc
        mbcomp(15)=VariabFlux(min(H+acc,Smelt),H,MASK,par.SeaIceThickness); % Smelt
        mbcomp(16)=VariabFlux(min(H+acc,runoff),H,MASK,par.SeaIceThickness); % runoff
        mbcomp(17)=VariabFlux(rain,H,MASK,par.SeaIceThickness); % rain
    else
        mbcomp(14)=VariabFlux(Pr,H,MASK,par.SeaIceThickness); % acc
        mbcomp(16)=VariabFlux(min(H+Pr,runoff),H,MASK,par.SeaIceThickness); % runoff
    end
    if ctr.Tcalc>=1
        mbcomp(18)=VariabFlux(min(H+max(-H,Mb),Bmelt),H,MASK,par.SeaIceThickness); % Bmelt
    end
    mbcomp(19)=mbcomp(12)-mbcomp(13)-mbcomp(18);
    VAF=max(0,H+min(B-SLR,0)*par.rhow/par.rho);
    VAFn=max(0,Hn+min(Bn-SLR,0)*par.rhow/par.rho);
    mbcomp(20)=(sum(sum(VAFn-VAF)))/ctr.dt;
    mbcomp(21)=mbcomp(12)-mbcomp(20)-mbcomp(18);
    % convert all components to Gt/a (water equivalent)
    mbcomp=mbcomp*ctr.delta^2*par.rho/1e12;
    
end


