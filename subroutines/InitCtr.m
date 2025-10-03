function [ctr,fc]=InitCtr(ctr,fc,default)

% Kori-ULB
% Initialization of ctr structure
% When not defined a priori, default values are assumed

    ctr.plotH(any(ismember(fields(ctr),'plotH'))==0)=0;
    ctr.runmode(any(ismember(fields(ctr),'runmode'))==0)=0;
    ctr.restart(any(ismember(fields(ctr),'restart'))==0)=0;
    ctr.diagnostic(any(ismember(fields(ctr),'diagnostic'))==0)=0;
    ctr.inverse(any(ismember(fields(ctr),'inverse'))==0)=0;
    ctr.timeslice(any(ismember(fields(ctr),'timeslice'))==0)=0;
    ctr.meltfac(any(ismember(fields(ctr),'meltfac'))==0)=0;
    ctr.meltfunc(any(ismember(fields(ctr),'meltfunc'))==0)=0;
    ctr.slidfac(any(ismember(fields(ctr),'slidfac'))==0)=1; % default no sliding perturbation
    ctr.schoof(any(ismember(fields(ctr),'schoof'))==0)=0;
    ctr.shelf(any(ismember(fields(ctr),'shelf'))==0)=0;
    ctr.SSA(any(ismember(fields(ctr),'SSA'))==0)=0;
    ctr.MbType(any(ismember(fields(ctr),'MbType'))==0)=0;
    ctr.MbConst(any(ismember(fields(ctr),'MbConst'))==0)=0;
    ctr.TsConst(any(ismember(fields(ctr),'TsConst'))==0)=0;
    ctr.TsType(any(ismember(fields(ctr),'TsType'))==0)=0;
    ctr.Tcalc(any(ismember(fields(ctr),'Tcalc'))==0)=0;
    ctr.Tinit(any(ismember(fields(ctr),'Tinit'))==0)=0;
    ctr.Enthalpy(any(ismember(fields(ctr),'Enthalpy'))==0)=0;
    ctr.drain(any(ismember(fields(ctr),'drain'))==0)=1;
    ctr.BedAdj(any(ismember(fields(ctr),'BedAdj'))==0)=0;
    ctr.m(any(ismember(fields(ctr),'m'))==0)=default.m; % default linear sliding
    ctr.p(any(ismember(fields(ctr),'p'))==0)=0;
    ctr.kmax(any(ismember(fields(ctr),'kmax'))==0)=default.kmax; % default number of z-levels    
    ctr.subwaterflow(any(ismember(fields(ctr),'subwaterflow'))==0)=0;
    ctr.SlidAdjust(any(ismember(fields(ctr),'SlidAdjust'))==0)=0;
    ctr.calving(any(ismember(fields(ctr),'calving'))==0)=0;
    ctr.LimitFront(any(ismember(fields(ctr),'LimitFront'))==0)=0;
    ctr.FrontalMelt(any(ismember(fields(ctr),'FrontalMelt'))==0)=0;
    ctr.CR(any(ismember(fields(ctr),'CR'))==0)=0;
    ctr.WV(any(ismember(fields(ctr),'WV'))==0)=0;
    ctr.HydroFrac(any(ismember(fields(ctr),'HydroFrac'))==0)=0;
    ctr.GeoidCalc(any(ismember(fields(ctr),'GeoidCalc'))==0)=0;
    ctr.starttime(any(ismember(fields(ctr),'starttime'))==0)=0;
    ctr.NumCheck(any(ismember(fields(ctr),'NumCheck'))==0)=0;
    ctr.shelfBC(any(ismember(fields(ctr),'shelfBC'))==0)=0;
    ctr.mismip(any(ismember(fields(ctr),'mismip'))==0)=0;
    ctr.basin(any(ismember(fields(ctr),'basin'))==0)=0;
    ctr.damage(any(ismember(fields(ctr),'damage'))==0)=0;
    ctr.PDDcalc(any(ismember(fields(ctr),'PDDcalc'))==0)=0;
    ctr.monthly(any(ismember(fields(ctr),'monthly'))==0)=0;
    ctr.runoffcorr(any(ismember(fields(ctr),'runoffcorr'))==0)=0;
    ctr.intT(any(ismember(fields(ctr),'intT'))==0)=default.intT;
    ctr.Hinv(any(ismember(fields(ctr),'Hinv'))==0)=default.Hinv;
    ctr.Tinv(any(ismember(fields(ctr),'Tinv'))==0)=default.Tinv;
    ctr.stopoptim(any(ismember(fields(ctr),'stopoptim'))==0)=default.stopoptim;
    ctr.HinvMelt(any(ismember(fields(ctr),'HinvMelt'))==0)=default.HinvMelt;
    ctr.TinvMelt(any(ismember(fields(ctr),'TinvMelt'))==0)=default.TinvMelt;
    ctr.radnorm(any(ismember(fields(ctr),'radnorm'))==0)=default.radnorm;
    ctr.snapshot(any(ismember(fields(ctr),'snapshot'))==0)=default.snapshot;
    ctr.BetaIter(any(ismember(fields(ctr),'BetaIter'))==0)=default.BetaIter;
    ctr.shelftune(any(ismember(fields(ctr),'shelftune'))==0)=default.shelftune;
    ctr.meltfactor(any(ismember(fields(ctr),'meltfactor'))==0)=default.meltfactor;
    ctr.Ao(any(ismember(fields(ctr),'Ao'))==0)=default.Ao;
    ctr.u0(any(ismember(fields(ctr),'u0'))==0)=default.u0;
    ctr.plotGL(any(ismember(fields(ctr),'plotGL'))==0)=default.plotGL;
    ctr.upstream(any(ismember(fields(ctr),'upstream'))==0)=default.upstream;
    ctr.ItSolv(any(ismember(fields(ctr),'ItSolv'))==0)=default.ItSolv;
    ctr.Asin(any(ismember(fields(ctr),'Asin'))==0)=default.Asin;
    ctr.taulim(any(ismember(fields(ctr),'taulim'))==0)=default.taulim;
    ctr.TRdam(any(ismember(fields(ctr),'TRdam'))==0)=1;
    ctr.THdam(any(ismember(fields(ctr),'THdam'))==0)=0;
    ctr.HLdam(any(ismember(fields(ctr),'HLdam'))==0)=1;
    if ctr.damage==1
        ctr.SFdam(any(ismember(fields(ctr),'SFdam'))==0)=1;
        ctr.BSdam(any(ismember(fields(ctr),'BSdam'))==0)=1;
    else
        ctr.SFdam(any(ismember(fields(ctr),'SFdam'))==0)=0;
        ctr.BSdam(any(ismember(fields(ctr),'BSdam'))==0)=0;
    end
    ctr.damlim(any(ismember(fields(ctr),'damlim'))==0)=default.damlim;
    ctr.OceanVisc(any(ismember(fields(ctr),'OceanVisc'))==0)=default.OceanVisc;
    if any(ismember(fields(ctr),'gammaT'))==0
        if ctr.meltfunc==1
            ctr.gammaT=default.gammaTlin;
        elseif ctr.meltfunc==2
            ctr.gammaT=default.gammaTquad;
        elseif ctr.meltfunc==9
            ctr.gammaT=default.gamma0_quad;
        elseif ctr.meltfunc==91
            ctr.gammaT=default.gamma0_quad_slope;
        else
            ctr.gammaT=default.gammaTpico;
        end
    end
    ctr.C(any(ismember(fields(ctr),'C'))==0)=default.Cpico;
    ctr.gammaTplume(any(ismember(fields(ctr),'gammaTplume'))==0)=default.gammaTplume;
    ctr.M0(any(ismember(fields(ctr),'M0'))==0)=default.M0picop;
    ctr.Hcrit(any(ismember(fields(ctr),'Hcrit'))==0)=default.Hcrit;
    ctr.YearlyMeans(any(ismember(fields(ctr),'YearlyMeans'))==0)=0;
    ctr.SnapList(any(ismember(fields(ctr),'SnapList'))==0)=0;
    
    if any(ismember(fields(fc),'DeltaT'))==0
        fc.DeltaT=zeros(ctr.nsteps,1);
    end
    if any(ismember(fields(fc),'DeltaSL'))==0
        fc.DeltaSL=zeros(ctr.nsteps,1);
    end
    if any(ismember(fields(fc),'butfac'))==0
        fc.butfac=ones(ctr.nsteps,1);
    end
    ctr.FreqHydro(any(ismember(fields(ctr),'FreqHydro'))==0)=1;
end


