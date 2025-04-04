function [H,Hn]=IceShelfContinuity(ctr,row,col,H,Hn,glMASK)

% Kori-ULB
% This algorithm systematically looks for adjacent points and ensures
% ice thickness continuity along the growth direction of the ice shelf.
% The idea is consistent with the concept that the difference between
% calving rate and front velocity gives the ice shelf terminus migration.
    
    for r=1:length(row)
        i = row(r);
        j = col(r);
        im = i;
        jm = j;
        ip = i;
        jp = j;
        for k=1:ctr.imax
            im = im - k;
            ip = ip + k;
            jp = jp + k;
            jm = jm - k;
            % Ensure indices within limits.
            if (im < 2) || (jm < 2) || (ip > ctr.imax-1) || (jp > ctr.jmax-1)
                break;
            end
            if glMASK(im,j)==4 || glMASK(im,j)==3 %&& glMASK_old(im,j)==4 instead??
                H(i,j)=H(im,j); % Necessary to add this?
                Hn(i,j)=Hn(im,j);
                break; 
            end
            if glMASK(ip,j)==4 || glMASK(ip,j)==3 %&& glMASK(ip,j)==4
                H(i,j)=H(ip,j);
                Hn(i,j)=Hn(ip,j);
                break;
            end
            if glMASK(i,jm)==4 || glMASK(i,jm)==3 %&& glMASK(i,jm)==4
                H(i,j)=H(i,jm);
                Hn(i,j)=Hn(i,jm);
                break;
            end
            if glMASK(i,jp)==4 || glMASK(i,jp)==3 %&& glMASK(i,jp)==4
                H(i,j)=H(i,jp);
                Hn(i,j)=Hn(i,jp);
                break;
            end
            if glMASK(im,jm)==4 || glMASK(im,jm)==3 %&& glMASK(im,jm)==4
                H(i,j)=H(im,jm);
                Hn(i,j)=Hn(im,jm);
                break;
            end
            if glMASK(im,jp)==4 || glMASK(im,jp)==3 %&& glMASK(im,jp)==4
                H(i,j)=H(im,jp);
                Hn(i,j)=Hn(im,jp);
                break;
            end
            if glMASK(ip,jm)==4 || glMASK(ip,jm)==3 %&& glMASK(ip,jm)==4
                H(i,j)=H(ip,jm);
                Hn(i,j)=Hn(ip,jm);
                break;
            end
            if glMASK(ip,jp)==4 || glMASK(ip,jp)==3 %&& glMASK(ip,jp)==4
                H(i,j)=H(ip,jp);
                Hn(i,j)=Hn(ip,jp);
                break;
            end
        end
    end
end


