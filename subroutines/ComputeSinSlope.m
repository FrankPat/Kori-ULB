function [sina]=ComputeSinSlope(HB,vx,vy,ctr)

% Kori-ULB
% Calculate local basal slop of ice shelves in Plume model

    s1=circshift(HB,[-1 0]); % i+1,j
    s2=circshift(HB,[1 0]); % i-1,j
    s3=circshift(HB,[0 -1]); % i,j+1
    s4=circshift(HB,[0 1]); % i,j-1
    MV=ones(ctr.imax,ctr.jmax); % u>=1, v>=1
    MV(vx>=0 & vy<0)=2;
    MV(vx<0 & vy>=0)=3;
    MV(vx<0 & vy<0)=4;
    sina=sin(sqrt(((s4-HB)/ctr.delta).^2+((s2-HB)/ctr.delta).^2));
    sina(MV==2)=sin(sqrt(((s4(MV==2)-HB(MV==2))/ctr.delta).^2+ ...
        ((s1(MV==2)-HB(MV==2))/ctr.delta).^2));
    sina(MV==3)=sin(sqrt(((s3(MV==3)-HB(MV==3))/ctr.delta).^2+ ...
        ((s2(MV==3)-HB(MV==3))/ctr.delta).^2));
    sina(MV==4)=sin(sqrt(((s3(MV==4)-HB(MV==4))/ctr.delta).^2+ ...
        ((s1(MV==4)-HB(MV==4))/ctr.delta).^2));
    sina=min(sina,0.05); % limit on slope of ice shelf
    
end


