function [bMASKm,bMASKx,bMASKy]=StaggeredBMASK(ctr,bMASK)

% Kori-ULB
% MASK for basin calculations on staggered grid

    if ctr.basin==1
        bMASKm=round((bMASK+circshift(bMASK,[-1 0])+circshift(bMASK,[0 -1])+ ...
            circshift(bMASK,[-1 -1]))/4); 
        bMASKx=round((bMASK+circshift(bMASK,[0 -1]))/2);
        bMASKy=round((bMASK+circshift(bMASK,[-1 0]))/2);
    else
        bMASKm=false;
        bMASKx=false;
        bMASKy=false;
    end
    
end


