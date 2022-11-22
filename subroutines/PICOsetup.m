function [Bk,Ak,Bmax]=PICOsetup(glMASK,H,ctr,par,ShelfN,numsh,shMASK)

% Kori-ULB
% Setup of the PICO model to determine ice shelf area
% (Ak) and the ocean boxes underneath each shelf (Bk)

    % Add first ocean point in glMASK
    MASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
    MASK2=circshift(glMASK,[0 1]); % glMASK(i,j-1)
    MASK3=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
    MASK4=circshift(glMASK,[1 0]); % glMASK(i-1,j)
    glMASK(glMASK==6 & (MASK1==5 | MASK2==5 | ...
        MASK3==5 | MASK4==5))=7; % First floating point (=7)
        
    % relative distance to GL: gMASK
    gMASK=ones(ctr.imax,ctr.jmax)-2;
    gMASK(glMASK==2)=0; % start at GL position
    gMASK(glMASK==3 | glMASK==4 | glMASK==5 | glMASK==7)=-2;
    for i=0:1000
        sg=sum(sum(gMASK==-2));
        if i>0 && (sg-sg0)==0
            break;
        end
        MASK1=circshift(gMASK,[0 -1]); % glMASK(i,j+1)
        MASK2=circshift(gMASK,[0 1]); % glMASK(i,j-1)
        MASK3=circshift(gMASK,[-1 0]); % glMASK(i+1,j)
        MASK4=circshift(gMASK,[1 0]); % glMASK(i-1,j)
        gMASK(gMASK==-2 & (MASK1==i | MASK2==i | MASK3==i | MASK4==i))=i+1;
        sg0=sg;
    end

    % relative distance from front: fMASK
    fMASK=ones(ctr.imax,ctr.jmax)-2;
    fMASK(glMASK==5)=0; % start at front position
    fMASK(glMASK==2 | glMASK==3 | glMASK==4)=-2;
    for i=0:1000
        sg=sum(sum(fMASK==-2));
        if i>0 && (sg-sg0)==0
            break;
        end
        MASK1=circshift(fMASK,[0 -1]); % glMASK(i,j+1)
        MASK2=circshift(fMASK,[0 1]); % glMASK(i,j-1)
        MASK3=circshift(fMASK,[-1 0]); % glMASK(i+1,j)
        MASK4=circshift(fMASK,[1 0]); % glMASK(i-1,j)
        fMASK(fMASK==-2 & (MASK1==i | MASK2==i | MASK3==i | MASK4==i))=i+1;
        sg0=sg;
    end

    gMASK(gMASK<0 | H<5*par.SeaIceThickness)=0;
    fMASK(fMASK<0 | H<5*par.SeaIceThickness)=0;
    gMASK=gMASK*ctr.delta;
    fMASK=fMASK*ctr.delta;

    % define boxes and relative distance r
    dmax=max(max(gMASK));
    rd=gMASK./(gMASK+fMASK); % relative distance between GL and front
    
    % define number of boxes for different ice shelves
    nB=1+round(sqrt(gMASK/dmax)*(par.nbox-1));
    nD=zeros(ctr.imax,ctr.jmax);
    for i=1:numsh
        nD(ShelfN==i)=max(nB(ShelfN==i));
    end

    Bk=zeros(ctr.imax,ctr.jmax);
    for k=1:par.nbox
        LL=1-sqrt(abs((nD-k+1)./nD));
        UL=1-sqrt(abs((nD-k)./nD));
        Bk(rd>=LL & rd<=UL)=k;
        Bk(Bk>nD)=nD(Bk>nD);
    end
    Bk(glMASK==3)=1; % make sure that boxes near GL are always Box 1
    Bmax=zeros(numsh,1);
    for i=1:numsh
        Bmax(i)=max(Bk(ShelfN==i));
    end
    % Calculate surface of each box (Ak)
    Ak=shMASK; % size of each box within particular ice shelf
    for i=1:numsh
        Ak(ShelfN==i)=Ak(ShelfN==i)*sum(sum(ShelfN==i));
    end
    Bk(Ak==1)=1;
    Ak(nD>0)=Ak(nD>0)*ctr.delta^2./nD(nD>0);
end


