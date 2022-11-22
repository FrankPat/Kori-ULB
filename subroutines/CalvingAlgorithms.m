function [CMB,LSF,he]=CalvingAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,glMASK,H,A, ...
    uxssa,uyssa,arcocn,B,runoff,MASK,MASKo,Ho,bMASK,LSF,node,nodes,VM)

% Kori-ULB
% Calving functions (UPDATE NEEDED)

    div=dudx+dvdy; % flow divergence
    [he,fi]=DefineEdgeThickness(ctr,par,glMASK,H); % Pollard 2015   %VL: add par
    Hshelf=H;
    Hshelf(fi<1)=he(fi<1);
    
    if ctr.calving==1 || ctr.calving==2
        dive=A.*(0.25*par.rho*par.g*Hshelf*(1.-par.rho/par.rhow)).^par.n;
        div(fi<1)=dive(fi<1); % divergence for shelf edge grid points
        div(div>dive)=dive(div>dive); % divergence not larger than maximum value
        d_s=2./(par.rho*par.g)*(max(0,div)./A).^(1/par.n); % dry surface crevasse depth
        d_b=2.*par.rho/((par.rhow-par.rho)*(par.rho*par.g))* ...
            (max(0,div)./A).^(1/par.n); % basal crevasse depth
        uxs1=circshift(uxssa,[0 1]); % ux(i+1,j)
        uys1=circshift(uyssa,[1 0]); % uy(i,j+1)
        usH=sqrt((0.5*(uxssa+uxs1)).^2+(0.5*(uyssa+uys1)).^2);
        d_a=Hshelf.*min(1,max(0,log(usH/1600))/log(1.2)); % additional crevasse deepening
        hc=150*max(0,min(1,(arcocn-70)./20)); % PD16: hc=200*max(0,min(1,(arcocn-40)./20));
        d_t=Hshelf.*max(0,min(1,(hc-Hshelf)/50)); % thin floating ice
        if ctr.HydroFrac==1
            RPDD=runoff; % annual surface melt + rainfall available after refreezing
            RPDD(runoff<1.5)=0; % R calculated as in DP16
            RPDD(runoff>3)=runoff(runoff>3).^2;
            RPDD(runoff<3 & runoff>1.5)=4*1.5*(runoff(runoff<3 & runoff>1.5)-1.5);
            d_w=100*RPDD;
        else
            d_w=zeros(ctr.imax,ctr.jmax);
        end
        % Calving
        ratio=(d_s+d_b+d_a+d_t+d_w)./he; 
        ratio_c=0.75; % Pollard(2015)
        CMB=3000*max(0,min(1,(ratio-ratio_c)./(1-ratio_c))).*Hshelf/ctr.delta;
        CMB(glMASK<4)=0;
    end
    if ctr.calving==2 % combination of 1 and 4 (making sure ice shelves don't get bigger)
        CMB(MASK==0 & MASKo==0 & Ho<par.SeaIceThickness)=50; % higher value than initial (LZ)
    end
    if ctr.calving==3
        wc=min(1,he/200.); % Pollard (2012)
        CMB=(1.-wc)*30+wc*3e5.*max(div,0).*he/ctr.delta; % Pollard (2012)
        if par.ArcOcean==1
            CMB=CMB.*max(0,min(1,(arcocn-70)./20)); % Vio: arcocean applied as in PD12
        end
        CMB(glMASK<5)=0;
    end
    if ctr.calving==4 % calving front kept at observed position (if ice is floating)
        CMB=zeros(ctr.imax,ctr.jmax);
        CMB(MASK==0 & MASKo==0 & Ho<par.SeaIceThickness)=50; % higher value than initial (LZ)
    end
    if ctr.calving>2 % Find holes in shelves + contour (where CMB applies)
        [MASKHole,HoleCT]=ShelfHole(MASK,H,par);
        CMB(MASKHole==1 | HoleCT==1)=0; % No calving in enclosed shelves holes
    end
    CMB(glMASK>=5 & B<-2700)=50; % Calving rule to prevent unrealistic areas of
                                 % thin ice extending seaward beyond continental
                                 % shelf (Vio)
    %VL: ice thinner than 5m can't persist --> equivalent to glMASK
    CMB(glMASK>=5 & H<5)=50;

    if ctr.calving==5 % LSF function calving (not operational)
        CMB=0;
        % determine eigenvalues of strain tensor
        StrTen=zeros(2,2);
        StrEig1=zeros(ctr.imax,ctr.jmax);
        StrEig2=zeros(ctr.imax,ctr.jmax);
        for i=1:ctr.imax
            for j=1:ctr.jmax
                StrTen(1,1)=dudx(i,j);
                StrTen(2,2)=dvdy(i,j);
                StrTen(1,2)=0.5*(dudy(i,j)+dvdx(i,j));
                StrTen(2,1)=StrTen(1,2);
                StrVal=eig(StrTen);
                StrEig1(i,j)=StrVal(1);
                StrEig2(i,j)=StrVal(2);
            end
        end
        EffStr=0.5*(max(0,StrEig1).^2+max(0,StrEig2).^2);
        sigmalim=0.8e5;
        q=sqrt(3)*A.^(-1/par.n).*EffStr.^(0.5/par.n)/sigmalim;
%             EigCalv=1e5*StrEig1.*StrEig2;
        FrontMelt=1000;
        wx=uxssa.*(1-q)-sign(uxssa).*FrontMelt; % w = u - c
        wy=uyssa.*(1-q)-sign(uyssa).*FrontMelt;
        if ctr.basin==1
            wx(bMASK==1)=0;
            wy(bMASK==1)=0;
        end
        wx=zeros(ctr.imax,ctr.jmax);
        wy=zeros(ctr.imax,ctr.jmax);
        LSF=LSFfunction(LSF,ctr,wx,wy,node,nodes,VM,MASK);
    end
end


