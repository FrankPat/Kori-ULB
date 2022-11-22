function [angnorm]=NormCalc(ig,jg,ih,jh,glMASK,ctr)

% Kori-ULB
% Function normcalc (to calculate the direction of the grounding line)

nwin = round(ctr.radnorm*1e3/ctr.delta);  %radnorm in [km]
nwin = nwin + 2; %make sure to have all grid cells within radius radnorm

zavx = 0;
zavy = 0;

for j=jh-nwin:jh+nwin
    jj=max(1,min(j,ctr.jmax));
    for i=ih-nwin:ih+nwin
        ii=max(1,min(i,ctr.imax));
        if glMASK(ii,jj)==0 || glMASK(ii,jj)>2
            zdx = (j - 0.5*(jg+jh));
            zdy = (i - 0.5*(ig+ih));
            zdist = sqrt(zdx^2 + zdy^2)*ctr.delta/1e3; 
            if zdist>0 && zdist<=ctr.radnorm
                zavx = zavx + zdx; 
                zavy = zavy + zdy;
            end
        end
    end
end
angnorm = atan2(zavy,zavx);

end


