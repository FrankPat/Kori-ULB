function [Melt_mean,Bmelt_mean,Ts_mean,Mb_mean,To_mean,So_mean,TF_mean, ...
    CR_mean,FMR_mean,fluxmx_mean,fluxmy_mean,Smelt_mean,runoff_mean, ...
    rain_mean,acc_mean]=YearlyMeans(Melt,Bmelt,Ts,Mb,To,So,TF,CR,FMR, ...
    fluxmx,fluxmy,Smelt,runoff,rain,acc,cnt,ctr,Melt_mean,Bmelt_mean, ...
    Ts_mean,Mb_mean,To_mean,So_mean,TF_mean,CR_mean,FMR_mean,fluxmx_mean, ...
            fluxmy_mean,Smelt_mean,runoff_mean,rain_mean,acc_mean)

% Kori-ULB
% Calculation of yearly mean values        

    if mod(cnt,1/ctr.dt)==1 % Initialise yearly-mean fields
        Melt_mean=Melt.*ctr.dt;
        Bmelt_mean=Bmelt.*ctr.dt;
        Ts_mean=Ts.*ctr.dt;
        Mb_mean=Mb.*ctr.dt;
        To_mean=To.*ctr.dt;
        So_mean=So.*ctr.dt;
        TF_mean=TF.*ctr.dt;
        CR_mean=CR.*ctr.dt;
        FMR_mean=FMR.*ctr.dt;
        fluxmx_mean=fluxmx.*ctr.dt;
        fluxmy_mean=fluxmy.*ctr.dt;
        Smelt_mean=Smelt.*ctr.dt;
        runoff_mean=runoff.*ctr.dt;
        rain_mean=rain.*ctr.dt;
        acc_mean=acc.*ctr.dt;
    else
        Melt_mean=Melt_mean+Melt.*ctr.dt;
        Bmelt_mean=Bmelt_mean+Bmelt.*ctr.dt;
        Ts_mean=Ts_mean+Ts.*ctr.dt;
        Mb_mean=Mb_mean+Mb.*ctr.dt;
        To_mean=To_mean+To.*ctr.dt;
        So_mean=So_mean+So.*ctr.dt;
        TF_mean=TF_mean+TF.*ctr.dt;
        CR_mean=CR_mean+CR.*ctr.dt;
        FMR_mean=FMR_mean+FMR.*ctr.dt;
        fluxmx_mean=fluxmx_mean+fluxmx.*ctr.dt;
        fluxmy_mean=fluxmy_mean+fluxmy.*ctr.dt;
        Smelt_mean=Smelt_mean+Smelt.*ctr.dt;
        runoff_mean=runoff_mean+runoff.*ctr.dt;
        rain_mean=rain_mean+rain.*ctr.dt;
        acc_mean=acc_mean+acc.*ctr.dt;
    end

end


