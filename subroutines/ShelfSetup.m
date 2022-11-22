function [ShelfN,numsh,shMASK,MASKlk]=ShelfSetup(MASK,glMASK,H,ctr,par)

% Kori-ULB
% Delimit ice shelves and number them separately (ShelfN). Also returns a
% shelf MASK and identifies lakes (locally floating conditions within
% grounded ice sheet

    %LZ: find all 'lakes'
    [~,MASKlk]=bwboundaries(MASK,4);
    MASKlk(MASK==1)=0;
    MASKlk(MASKlk~=0)=1;

    shMASK=zeros(ctr.imax,ctr.jmax); % Mask of ice shelf
    shMASK(MASKlk==0 & (glMASK==3 | glMASK==4 | glMASK==5))=1;
    shMASK(H<5*par.SeaIceThickness)=0;
    
    %LZ: new function for calculation of ShelfN
    ShelfN=bwlabel(shMASK);
    numsh=max(ShelfN(:));

end


