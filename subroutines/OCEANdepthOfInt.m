function [To,So,TF]=OCEANdepthOfInt(fc,ctr,par,Tof,Sof,TFf,H,HB,B,MASK,glMASK,ShelfN,numsh)

% Kori-ULB
% Interpolate ocean variables to depth of interest (bottom of ice shelf)

    if any(ismember(fields(fc),'z')) % Depth-profile of To and So data or TF
        if ctr.meltfunc==3 || ctr.meltfunc==4 
            % compute average depth of continental shelf for each ice shelf
            % front (Burgard et al. 2021 -- PROTECT)
            front_bot_depth_avg=zeros(ctr.imax,ctr.jmax)-500;
            for b=1:numsh
                avg_B=mean(B(ShelfN==b & glMASK==5));
                if isnan(avg_B)==1
                    avg_B=mean(B(ShelfN==b)); % some small shelves have front 
                                              % not well define, then average 
                                              % over shelf area
                end
                front_bot_depth_avg(ShelfN==b)=avg_B;
            end
            depth_of_int=front_bot_depth_avg;
            To=InterpToDepthOfInt(Tof,fc.z,depth_of_int,ctr);
            So=InterpToDepthOfInt(Sof,fc.z,depth_of_int,ctr);
            if islogical(TFf)==0
                TF=InterpToDepthOfInt(TFf,fc.z,depth_of_int,ctr);
            else
                TF=TFf;
            end
        else
            % compute deepest entrance depth of continental shelf for each 
            % ice shelf (Burgard et al. 2021 -- PROTECT)
            front_bot_depth_max=zeros(ctr.imax,ctr.jmax)-2500;
            for b=1:numsh
                max_B=min(B(ShelfN==b & glMASK==5));
                if isempty(max_B)==1
                    max_B=min(B(ShelfN==b)); % some small shelves have front 
                                             % not well define, then average 
                                             % over shelf area
                end
                front_bot_depth_max(ShelfN==b)=max_B;
            end
            depth_of_int=HB; % interpolate at ice draft
            depth_of_int(HB<front_bot_depth_max)= ...
                front_bot_depth_max(HB<front_bot_depth_max);
                % Limit depth of int to the deepest entrance depth 
                % of continental shelf
            if islogical(TFf)==0
                TF=InterpToDepthOfInt(TFf,fc.z,depth_of_int,ctr);
                To=Tof;
                So=Sof;
            else
                TF=TFf;
                To=InterpToDepthOfInt(Tof,fc.z,depth_of_int,ctr);
                So=InterpToDepthOfInt(Sof,fc.z,depth_of_int,ctr);
            end
        end
    else
        if (ctr.meltfunc==3 || ctr.meltfunc==4) && ...
                any(ismember(fields(fc),'PICO_basins')) 
            % PICO method of updating values by basin (Kreutzer et al. 2021)
            [contshelfMASK]=ContinentalShelfMASK(MASK,H,B,par);
                % update continental shelf mask
            SO=MASK*0;
            TO=MASK*0;
            for b=1:19
                % compute average of edge grid cells within a specific
                % PICO basin and assign that value to So or To field
                avg_SO=mean(mean(Sof(fc.PICO_basins==b & contshelfMASK==1)));
                avg_TO=mean(mean(Tof(fc.PICO_basins==b & contshelfMASK==1)));
                SO(fc.PICO_basins==b)=avg_SO;
                TO(fc.PICO_basins==b)=avg_TO;
            end
            To=TO;
            So=SO;
            TF=TFf;
        else
            TF=TFf;
            To=Tof;
            So=Sof;
        end
    end

end


