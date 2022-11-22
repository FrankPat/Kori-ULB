function [glMASK,MASK]=GroundingLineMask(MASK,bMASK,H,ctr)

% Kori-ULB
% Definition of glMASK when grounding lines appear in the domain
%   glMASK=1: grounded
%   glMASK=2: grounding line
%   glMASK=3: first floating grid point
%   glMASK=4: ice shelf
%   glMASK=5: calving front
%   glMASK=§: open ocean

    glMASK=MASK;
    MASK1=circshift(MASK,[0 -1]); % MASK(i,j+1)
    MASK2=circshift(MASK,[0 1]); % MASK(i,j-1)
    MASK3=circshift(MASK,[-1 0]); % MASK(i+1,j)
    MASK4=circshift(MASK,[1 0]); % MASK(i-1,j)

    glMASK(MASK==1 & (MASK1==0 | MASK2==0 | MASK3==0 | MASK4==0))=2;

    MASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
    MASK2=circshift(glMASK,[0 1]); % glMASK(i,j-1)
    MASK3=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
    MASK4=circshift(glMASK,[1 0]); % glMASK(i-1,j)

    % adjacent floating grid cell (=3)
    glMASK(MASK==0 & (MASK1==2 | MASK2==2 | MASK3==2 | MASK4==2))=3;

    if ctr.shelf==1
        glMASK(glMASK==0)=4; % ice shelf (=4)
        if ctr.calving>=1
            glMASK(glMASK==4 & H<5)=6;  % sea ice/ocean (=6)
            MASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
            MASK2=circshift(glMASK,[0 1]); % glMASK(i,j-1)
            MASK3=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
            MASK4=circshift(glMASK,[1 0]); % glMASK(i-1,j)
            glMASK((glMASK==3 | glMASK==4) & (MASK1==6 | MASK2==6 | ...
                MASK3==6 | MASK4==6))=5; % Calving front (=5)
        end
    end
    if ctr.mismip>=1
        glMASK(:,1)=glMASK(:,3);
        glMASK(1,:)=glMASK(3,:);
        if ctr.mismip==1
            glMASK(ctr.imax,:)=glMASK(ctr.imax-2,:);
        else
            glMASK(:,ctr.jmax)=glMASK(:,ctr.jmax-1);
            glMASK(ctr.imax,:)=glMASK(ctr.imax-1,:);
        end
    end
    if ctr.basin==1
        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(1,:)=1;
        MASKb(ctr.imax,:)=1;
        MASKb(:,1)=1;
        MASKb(:,ctr.jmax)=1;
        glMASK(MASKb==1 & MASK==1)=1;
        glMASK(MASKb==1 & MASK==0)=6;
        glMASK(bMASK==1)=1;
    end

end


