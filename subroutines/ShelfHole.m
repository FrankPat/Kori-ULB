function [MASKHole,HoleCT]=ShelfHole(MASK,H,par)

% Kori-ULB
% Checking for eventual holes in ice shelves and returning a MASK that
% identifies these

    MASK(MASK==0 & H>par.SeaIceThickness)=3;
    [~,MASKHole]=bwboundaries(MASK>0,4);
    MASKHole(MASK==1 | H>par.SeaIceThickness)=0;
    MASKHole(MASKHole~=0)=1;
    H1=circshift(MASKHole,[-1 0]); % i+1,j
    H2=circshift(MASKHole,[1 0]); % i-1,j
    H3=circshift(MASKHole,[0 -1]); % i,j+1
    H4=circshift(MASKHole,[0 1]); % i,j-1
    HoleCT=H1==1 | H2==1 | H3==1 | H4==1;
end


