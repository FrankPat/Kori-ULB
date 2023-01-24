function [gradm,gradmx,gradmy,gradxy,gradsx,gradsy,gradHx,gradHy, ...
    Hm,Hmx,Hmy,Bmx,Bmy,signx,signy]=StaggeredGrid(sn,H,B,ctr)

% Kori-ULB
% Variables on staggered grids

    sn1=circshift(sn,[-1 0]); % sn(i+1,j)
    sn2=circshift(sn,[0 -1]); % sn(i,j+1)
    sn3=circshift(sn,[-1 -1]); % sn(i+1,j+1)
    sn4=circshift(sn,[0 1]); % sn(i,j-1)
    sn5=circshift(sn,[1 0]); % sn(i-1,j)

    Hm=(H+circshift(H,[-1 0])+circshift(H,[0 -1])+circshift(H,[-1 -1]))/4.; 
    gradm=((sn2+sn3-sn-sn1)/(2*ctr.delta)).^2+((sn1+sn3-sn-sn2)/(2*ctr.delta)).^2;
    gradxy=((sn2-sn4)/(2*ctr.delta)).^2+((sn1-sn5)/(2*ctr.delta)).^2;
    
    H1=circshift(H,[-1 0]); % sn(i+1,j)
    H2=circshift(H,[0 -1]); % sn(i,j+1)
    H4=circshift(H,[0 1]); % sn(i,j-1)
    H5=circshift(H,[1 0]); % sn(i-1,j)
    gradsx=(sn2-sn4)/(2*ctr.delta);
    gradsy=(sn1-sn5)/(2*ctr.delta);
    gradHx=(H2-H4)/(2*ctr.delta);
    gradHy=(H1-H5)/(2*ctr.delta);
    
    gradmx=(sn2-sn)/ctr.delta;
    Hmx=(H+circshift(H,[0 -1]))/2.;
    Bmx=(B+circshift(B,[0 -1]))/2.;
    gradmy=(sn1-sn)/ctr.delta;
    Hmy=(H+circshift(H,[-1 0]))/2.;
    Bmy=(B+circshift(B,[-1 0]))/2.;
    
    signx=sign(-gradmx);
    signy=sign(-gradmy);
    gradmx=gradmx.^2;
    gradmy=gradmy.^2;
    
    if ctr.mismip>0
        gradxy(:,1)=gradxy(:,3);
        gradxy(:,end)=gradxy(:,end-1);
        gradxy(1,:)=gradxy(3,:);
        gradxy(end,:)=gradxy(end-3,:);
    end

end


