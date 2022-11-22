function [Tsf,Mbf,Prf,Evpf,runofff,cnt_atm,snp_atm,Mb_update]= ...
    ATMupdate(fc,ctr,time,cnt,Ts0,Mb0,Pr0,Evp0,runoff0,cnt_atm, ...
    snp_atm,Mb_update,Tsf,Mbf,Prf,Evpf,runofff)

% Kori-ULB
% Atmospheric forcing (precipitation, runoff, evaporation) based on
% sequential input files from (regional) climate models

    if fc.forcingATM==1
        if time(cnt)==fc.atm_Tinit % Initialise at first snapshot year
            cnt_atm=1;
            snp_atm=1;
        end
        if cnt_atm==1 % snapshot
            if time(cnt)<=fc.atm_Tend
                if any(ismember(fields(fc),'atm_Ts_fname'))
                    load([fc.atm_Ts_fname,num2str(snp_atm,'%03i')]);
                    Tsf=double(Ts); % make sure matrix is double
                else
                    Tsf=Ts0+fc.DeltaT(cnt); % simplified forcing otherwise
                end
                if any(ismember(fields(fc),'atm_Mb_fname'))
                    load([fc.atm_Mb_fname,num2str(snp_atm,'%03i')]);
                    Mbf=double(Mb); % make sure matrix is double
                    Mb_update=1;
                else % If Mb forcing does not exist, check for Mb component forcing
                    Mb_update=0;
                    Mbf=Mb0;
                    if any(ismember(fields(fc),'atm_Pr_fname'))
                        load([fc.atm_Pr_fname,num2str(snp_atm,'%03i')]);
                        Prf=double(Pr); % make sure matrix is double
                    else
                        Prf=Pr0;
                    end
                    if any(ismember(fields(fc),'atm_Evp_fname'))
                        load([fc.atm_Evp_fname,num2str(snp_atm,'%03i')]);
                        Evpf=double(Evp); % make sure matrix is double
                    else
                        Evpf=Evp0;
                    end
                    if any(ismember(fields(fc),'atm_runoff_fname')) && ctr.PDDcalc==0
                        load([fc.atm_runoff_fname,num2str(snp_atm,'%03i')]);
                        runofff=double(runoff); % make sure matrix is double
                    else
                        runofff=runoff0;
                    end
                end
                snp_atm=snp_atm+1;
            else % If forcing shorter than simulation time
                if snp_atm==fc.atm_snapshots+1
                    snp_atm=fc.atm_snapshots+1-fc.atm_nrep; % Re-initialise snapshot
                end
                if any(ismember(fields(fc),'atm_Ts_fname'))
                    load([fc.atm_Ts_fname,num2str(snp_atm,'%03i')]);
                    Tsf=double(Ts); % make sure matrix is double
                else
                    Tsf=Ts0+fc.DeltaT(cnt); % simplified forcing otherwise
                end
                if any(ismember(fields(fc),'atm_Mb_fname'))
                    load([fc.atm_Mb_fname,num2str(snp_atm,'%03i')]);
                    Mbf=double(Mb); % make sure matrix is double
                    Mb_update=1;
                else % If direct Mb forcing does not exist, check for Mb components forcing
                    Mb_update=0;
                    Mbf=Mb0;
                    if any(ismember(fields(fc),'atm_Pr_fname'))
                        load([fc.atm_Pr_fname,num2str(snp_atm,'%03i')]);
                        Prf=double(Pr); % make sure matrix is double
                    else
                        Prf=Pr0;
                    end
                    if any(ismember(fields(fc),'atm_Evp_fname'))
                        load([fc.atm_Evp_fname,num2str(snp_atm,'%03i')]);
                        Evpf=double(Evp); % make sure matrix is double
                    else
                        Evpf=Evp0;
                    end
                    if any(ismember(fields(fc),'atm_runoff_fname')) && ctr.PDDcalc==0
                        load([fc.atm_runoff_fname,num2str(snp_atm,'%03i')]);
                        runofff=double(runoff); % make sure matrix is double
                    else
                        runofff=runoff0;
                    end
                end
                snp_atm=snp_atm+1;
            end
        end
        cnt_atm=cnt_atm+1;
        cnt_atm(cnt_atm>fc.atm_cnt)=1;
    else
        Tsf=Ts0+fc.DeltaT(cnt); % simplified forcing otherwise
        Mbf=Mb0;
        Prf=Pr0;
        Evpf=Evp0;
        runofff=runoff0;
    end

end


