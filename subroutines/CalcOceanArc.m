function [arcocn,distocn,distgl]=CalcOceanArc(io,jo,MASK,imax,jmax, ...
    delta,narc,tanarc,ifi,idir)

% Kori-ULB
% Calculation of the angle of ice shelves to open ocean (Pollard and
% DeConto)
    
    arcocn=0;
    distocn=1e6;
    distgl=1e7;
    for mm=1:narc
        idirmm=idir(mm);
        tanarcmm=tanarc(mm);
        calcdist=0;
        if ifi(mm)==1
            ii=io;
            if idirmm==-1
                ttmax=io-1;
            else
                ttmax=imax-io-1;
            end
            for tt=1:ttmax
                ii=ii+idirmm;
                jj=round((ii-io)*tanarcmm)+jo;
                if jj<1 || jj>jmax
                    arcocn=arcocn+360./narc;
                    break
                end
                if MASK(ii,jj)==1
                    zdist=delta*sqrt((ii-io)^2+(jj-jo)^2);
                    distgl=min(distgl,zdist);
                    break
                elseif MASK(ii,jj)==0 && calcdist==0
                    zdist=delta*sqrt((ii-io)^2+(jj-jo)^2);
                    distocn=min(distocn,zdist);
                    calcdist=1;
                end
            end
        else
            jj=jo;
            if idirmm==-1
                ttmax=jo-1;
            else
                ttmax=jmax-jo-1;
            end
            for tt=1:ttmax
                jj=jj+idirmm;
                ii=round((jj-jo)/tanarcmm)+io;
                if ii<1 || ii>imax
                    arcocn=arcocn+360./narc;
                    break
                end
                if MASK(ii,jj)==1
                    zdist=delta*sqrt((ii-io)^2+(jj-jo)^2);
                    distgl=min(distgl,zdist);
                    break
                elseif MASK(ii,jj)==0 && calcdist==0
                    zdist=delta*sqrt((ii-io)^2+(jj-jo)^2);
                    distocn=min(distocn,zdist);
                    calcdist=1;
                end
            end
        end
        if tt==ttmax
            arcocn=arcocn+360./narc;
        end
    end
end


