function [arcocn,distocn,distgl]=OceanArc(MASK,H,MASKlk,ctr,par)

% Kori-ULB
% Calculation of the angle of ice shelves to open ocean (Pollard and
% DeConto)

    arcocn=zeros(ctr.imax,ctr.jmax);
    distocn=zeros(ctr.imax,ctr.jmax)+1000e3;
    distgl=zeros(ctr.imax,ctr.jmax)+10000e3;
    narc=72; % number of directions (every 5 degree)
    angarc=zeros(1,narc);
    tanarc=zeros(1,narc);
    ifi=zeros(1,narc);
    idir=zeros(1,narc);
    MASK(MASK==0&H>par.SeaIceThickness)=3;  %mark shelves in MASK
    MASK1=circshift(MASK,[0 -1]); % MASK(i,j+1)
    MASK2=circshift(MASK,[0 1]); % MASK(i,j-1)
    MASK3=circshift(MASK,[-1 0]); % MASK(i+1,j)
    MASK4=circshift(MASK,[1 0]); % MASK(i-1,j)
    MASK5=circshift(MASK,[-1 -1]); % MASK(i+1,j+1)
    MASK6=circshift(MASK,[-1 1]); % MASK(i+1,j-1)
    MASK7=circshift(MASK,[1 -1]); % MASK(i-1,j+1)
    MASK8=circshift(MASK,[1 1]); % MASK(i-1,j-1)
    ifdo=zeros(ctr.imax,ctr.jmax);
    ifdo(MASK~=1&MASKlk==0&(MASK1~=0|MASK2~=0|MASK3~=0|MASK4~=0|MASK5~=0|MASK6~=0|MASK7~=0|MASK8~=0))=1;
    [jo,io] = meshgrid(1:ctr.imax,1:ctr.jmax);
    distocn(MASK==0)=0;

    for mm=1:narc
        angarc(mm)=-pi+(mm-0.5)*(2*pi)/narc;
        tanarc(mm)=tan(angarc(mm));
        if abs(angarc(mm))>=0.75*pi
            ifi(mm)=1;
            idir(mm)=-1;
        elseif angarc(mm)>=-0.75*pi && angarc(mm)<=-0.25*pi
            ifi(mm)=0;
            idir(mm)=-1;
        elseif abs(angarc(mm))<=0.25*pi
            ifi(mm)=1;
            idir(mm)=1;
        else
            ifi(mm)=0;
            idir(mm)=1;
        end
    end

    [arcocn(ifdo==1),distocn(ifdo==1),distgl(ifdo==1)]=arrayfun(@(io,jo) CalcOceanArc(io,jo,MASK,ctr.imax,ctr.jmax,ctr.delta,narc,tanarc,ifi,idir),io(ifdo==1),jo(ifdo==1));

    arcocn(MASK~=1&ifdo==0&MASKlk==0)=360;
    distocn(MASK==0)=0;
    distgl(MASK~=1&ifdo==0)=0;
end


