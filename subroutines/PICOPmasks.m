function [MASKpicop,MASKb]=PICOPmasks(MASK,MASKlk,H,ctr,par)

% Kori-ULB
% Defines the MASK for PICOP domain (ice shelf + grounding line + first
% ocean grid point) and the boundary (grounding line and first ocean
% grid point)

    MASKpicop=zeros(ctr.imax,ctr.jmax); % initialize
    MASKb=zeros(ctr.imax,ctr.jmax); % ice shelf boundary points
    MASK=MASK+MASKlk; % FP: make sure that lakes are removed
    % grounding line (=2)
    MASK1=circshift(MASK,[0 -1]);
    MASK2=circshift(MASK,[0 1]);
    MASK3=circshift(MASK,[-1 0]);
    MASK4=circshift(MASK,[1 0]);
    MASKpicop(MASK==1 & (MASK1==0 | MASK2==0 | MASK3==0 | MASK4==0))=1;
    MASKb(MASK==1 & (MASK1==0 | MASK2==0 | MASK3==0 | MASK4==0))=1;
    % ice shelves (=1)
    MASKpicop(MASK==0 & H>par.SeaIceThickness)=1;
    % adjacent floating grid cell
    MASK1=circshift(MASKpicop,[0 -1]);
    MASK2=circshift(MASKpicop,[0 1]);
    MASK3=circshift(MASKpicop,[-1 0]);
    MASK4=circshift(MASKpicop,[1 0]);
    MASKb(MASKpicop==0 & MASK==0 & (MASK1==1 | MASK2==1 | MASK3==1 | MASK4==1))=1;
    MASKpicop(MASKpicop==0 & MASK==0 & (MASK1==1 | MASK2==1 | MASK3==1 | MASK4==1))=1;
    % remove ice shelves near border of domain
%     MASKpicop(:,1:3)=0;
%     MASKpicop(:,ctr.jmax-2:ctr.jmax)=0;
%     MASKpicop(1:3,:)=0;
%     MASKpicop(ctr.imax-2:ctr.imax,:)=0;

end


