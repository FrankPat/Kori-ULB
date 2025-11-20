function [Tof,Sof,TFf,cnt_ocn,snp_ocn]=OCEANupdate(fc,time,cnt,So0, ...
    To0,Tof,Sof,TFf,cnt_ocn,snp_ocn)

% Kori-ULB
% Update ocean forcing based on external forcing data

    if fc.forcingOCEAN==1
        if time(cnt)==fc.ocn_Tinit % Initialise at first snapshot year
            cnt_ocn=1;
            snp_ocn=1;
        end
        if cnt_ocn==1 % snapshot
            if time(cnt)<=fc.ocn_Tend
                if any(ismember(fields(fc),'ocn_TF_fname'))
                    load([fc.ocn_TF_fname,num2str(snp_ocn,'%03i')]);
                    TFf=double(TF); % make sure matrix is double
                    Tof=To0+fc.DeltaTo(cnt); 
                    % simplified forcing otherwise
                    Sof=So0;
                else % If TF forcing does not exist, check for To or So forcing
                    TFf=false;
                    if any(ismember(fields(fc),'ocn_To_fname'))
                        load([fc.ocn_To_fname,num2str(snp_ocn,'%03i')]);
                        Tof=double(To); % make sure matrix is double
                    else
                        Tof=To0+fc.DeltaTo(cnt); 
                        % simplified forcing otherwise
                    end
                    if any(ismember(fields(fc),'ocn_So_fname'))
                        load([fc.ocn_So_fname,num2str(snp_ocn,'%03i')]);
                        Sof=double(So); % make sure matrix is double
                    else
                        Sof=So0;
                    end
                end
                snp_ocn=snp_ocn+1;
            else % If forcing shorter than simulation time
                if snp_ocn==fc.ocn_snapshots+1
                    snp_ocn=fc.ocn_snapshots+1-fc.ocn_nrep; 
                    % Re-initialise snapshot
                end
                if any(ismember(fields(fc),'ocn_TF_fname'))
                    load([fc.ocn_TF_fname,num2str(snp_ocn,'%03i')]);
                    TFf=double(TF); % make sure matrix is double
                    Tof=To0+fc.DeltaTo(cnt); 
                    % simplified forcing otherwise
                    Sof=So0;
                else % If TF forcing does not exist, check for To or So forcing
                    TFf=false;
                    if any(ismember(fields(fc),'ocn_To_fname'))
                        load([fc.ocn_To_fname,num2str(snp_ocn,'%03i')]);
                        Tof=double(To); % make sure matrix is double
                    else
                        Tof=To0+fc.DeltaTo(cnt); 
                        % simplified forcing otherwise
                    end
                    if any(ismember(fields(fc),'ocn_So_fname'))
                        load([fc.ocn_So_fname,num2str(snp_ocn,'%03i')]);
                        Sof=double(So); % make sure matrix is double
                    else
                        Sof=So0;
                    end
                end
                snp_ocn=snp_ocn+1;
            end
        end
        cnt_ocn=cnt_ocn+1;
        cnt_ocn(cnt_ocn>fc.ocn_cnt)=1;
    else
        TFf=false;
        Tof=To0+fc.DeltaTo(cnt); % simplified forcing otherwise
        Sof=So0;
    end

end


