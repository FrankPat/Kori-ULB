function [Melt_mean,Bmelt_mean,Ts_mean,Mb_mean,To_mean,So_mean,TF_mean,CMB_mean,FMB_mean,fluxmx_mean,fluxmy_mean]=YearlyMeans(Melt,Bmelt,Ts,Mb,To,So,TF,CMB,FMB,fluxmx,fluxmy,cnt,ctr,Melt_mean,Bmelt_mean,Ts_mean,Mb_mean,To_mean,So_mean,TF_mean,CMB_mean,FMB_mean,fluxmx_mean,fluxmy_mean)

if cnt==1 % for first snapshot
    Melt_mean=Melt;
    Bmelt_mean=Bmelt;
    Ts_mean=Ts;
    Mb_mean=Mb;
    To_mean=To;
    So_mean=So;
    TF_mean=TF;
    CMB_mean=CMB;
    FMB_mean=FMB;
    fluxmx_mean=fluxmx;
    fluxmy_mean=fluxmy;
elseif mod(cnt,1/ctr.dt)==2 % Initialise yearly-mean fields
    Melt_mean=Melt.*ctr.dt;
    Bmelt_mean=Bmelt.*ctr.dt;
    Ts_mean=Ts.*ctr.dt;
    Mb_mean=Mb.*ctr.dt;
    To_mean=To.*ctr.dt;
    So_mean=So.*ctr.dt;
    TF_mean=TF.*ctr.dt;
    CMB_mean=CMB.*ctr.dt;
    FMB_mean=FMB.*ctr.dt;
    fluxmx_mean=fluxmx.*ctr.dt;
    fluxmy_mean=fluxmy.*ctr.dt;
else
    Melt_mean=Melt_mean+Melt.*ctr.dt;
    Bmelt_mean=Bmelt_mean+Bmelt.*ctr.dt;
    Ts_mean=Ts_mean+Ts.*ctr.dt;
    Mb_mean=Mb_mean+Mb.*ctr.dt;
    To_mean=To_mean+To.*ctr.dt;
    So_mean=So_mean+So.*ctr.dt;
    TF_mean=TF_mean+TF.*ctr.dt;
    CMB_mean=CMB_mean+CMB.*ctr.dt;
    FMB_mean=FMB_mean+FMB.*ctr.dt;
    fluxmx_mean=fluxmx_mean+fluxmx.*ctr.dt;
    fluxmy_mean=fluxmy_mean+fluxmy.*ctr.dt;
end

end