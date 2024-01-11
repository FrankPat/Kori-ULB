function [CMB,LSF,CR]=CalvingAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,glMASK,H,A, ...
    uxssa,uyssa,arcocn,B,runoff,MASK,MASKo,Ho,bMASK,LSF,node,nodes,VM,cnt,ux,uy,Melt,he,fi,FMR)

ux1=circshift(ux,[0 1]); % ux(i,j-1)
uy1=circshift(uy,[1 0]); % uy(i-1,j)

MAGV=max(0.0000000001,sqrt((0.5*(ux+ux1)).^2+(0.5*(uy+uy1)).^2));

XUV=-(0.5*(ux+ux1))./MAGV;
YUV=-(0.5*(uy+uy1))./MAGV;

uxh=(0.5*(ux+ux1));
uyh=(0.5*(uy+uy1));

if ctr.calving>=1 % LSF function calving. Generate a calving rate, CR

    div=dudx+dvdy; % flow divergence
    Hshelf=H;
    Hshelf(fi<1)=he(fi<1);

    if ctr.calving==1 % Direct, constant imposition of calving rate

        CR=zeros(size(LSF));
        CR(:,:)=ctr.CR;

    end

    if ctr.calving==2 % Direct, constant imposition of change in front positon. **ctr.WV=0 will fix calving front position to be unmoving**

        MAGV=sqrt((0.5*(ux+ux1)).^2+(0.5*(uy+uy1)).^2);
        CR=MAGV-ctr.WV;

    end

    if ctr.calving==3  % Pollard (2012) 

        wc=min(1,he/200.);
        CR=(1.-wc)*par.Calve_PD_MinThick+wc*par.Calve_PD_MaxCalveRate.*max(div,0);

    end

    if ctr.calving==4 % Pollard 2015

        dive=A.*(0.25*par.rho*par.g*Hshelf*(1.-par.rho/par.rhow)).^par.n;
        div(fi<1)=dive(fi<1); % divergence for shelf edge grid points
        div(div>dive)=dive(div>dive); % divergence not larger than maximum value
        d_s=2./(par.rho*par.g)*(max(0,div)./A).^(1/par.n); % dry surface crevasse depth
        d_b=2.*par.rho/((par.rhow-par.rho)*(par.rho*par.g))*(max(0,div)./A).^(1/par.n); % basal crevasse depth
        uxs1=circshift(uxssa,[0 1]); % ux(i+1,j)
        uys1=circshift(uyssa,[1 0]); % uy(i,j+1)
        usH=sqrt((0.5*(uxssa+uxs1)).^2+(0.5*(uyssa+uys1)).^2);
        d_a=Hshelf.*min(1,max(0,log(usH/par.Ucrit1))/log(par.Ucrit2/par.Ucrit1)); % additional crevasse deepening
        if par.ArcOcean==1
            hc=ctr.Hcrit*max(0,min(1,(arcocn-40)./20)); % PD16: hc=200*max(0,min(1,(arcocn-70)./20));
        else
            hc=zeros(ctr.imax,ctr.jmax)+ctr.Hcrit;
        end
        d_t=Hshelf.*max(0,min(1,(hc-Hshelf)./50)); % thin floating ice
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
        % ratio_c=0.75; % Pollard(2015)
        CR=par.Calve_P_MaxCalveRate*max(0,min(1,(ratio-par.Calve_P_CritCrevasse)./(1-par.Calve_P_CritCrevasse)));

    end

    if ctr.calving==5

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
        EffStr=max(0.0000000001,sqrt(0.5*(max(0,StrEig1).^2+max(0,StrEig2).^2)));
        tauVM=sqrt(3)*(EffStr./A).^(1./par.n);
        q=tauVM./ctr.taulim;
        CR=MAGV.*q;
    end

    if ctr.calving==6 % CalveMip thickness dependent calve rate

        CR= max(0,1+(ctr.Hcrit-H)/ctr.Hcrit).*MAGV;

    end

    if ctr.calving==7 % CalveMip Periodic forcing, ctr.CR_AMP is max rate of front position change

        Wv=-ctr.CR_AMP*sind(cnt*360/ctr.nsteps);
        MAGV=sqrt(ux.^2+uy.^2);
        CR=MAGV-Wv;

    end

    if ctr.calving==8 %CalveMip Periodic forcing, ctr.CR_AMP is max rate of front position change
        if cnt <10000
            Wv=-ctr.CR_AMP*sind(cnt*360/ctr.nsteps);
            MAGV=sqrt(ux.^2+uy.^2);
            CR=MAGV-Wv;
        else

            CR=zeros(size(LSF));
        end
    end

    if ctr.LimitFront==1
        CR(CR<MAGV & (MASKo==3 & circshift(MASKo,[-1 0])==0 | MASKo==3 & circshift(MASKo,[1 0])==0 ...
            | MASKo==3 & circshift(MASKo,[0 -1])==0 | MASKo==3 & circshift(MASKo,[0 1])==0)) ...
            =MAGV(CR<MAGV & (MASKo==3 & circshift(MASKo,[-1 0])==0 | MASKo==3 & circshift(MASKo,[1 0])==0 ...
            | MASKo==3 & circshift(MASKo,[0 -1])==0 | MASKo==3 & circshift(MASKo,[0 1])==0));
    end

    if ctr.calving==4 || ctr.calving==6 || ctr.LimitFront==1
        % extrapolate calving front CR value in the open ocean
        CR(glMASK==6)=NaN;
        [rows, cols] = find(~isnan(CR)); % Find the indices of non-NaN elements in CR
        values = CR(~isnan(CR)); % Find the corresponding values for non-NaN elements
        [nanRows, nanCols] = find(isnan(CR)); % Find the indices of NaN elements in CR
        % Interpolate the NaN elements in CR based on the non-NaN elements using 'nearest' method
        interpValues = griddata(cols, rows, values, nanCols, nanRows, 'nearest');
        % Replace the NaN elements in CR with the interpolated values
        CR(isnan(CR)) = interpValues;
    end

    % Diagnostic calculation of calving rate as a surface balance term over
    % the entire grid cell
    CMB=CR.*Hshelf/ctr.delta;
    CMB(glMASK~=5)=0;

    if  ctr.shelf==1 && ctr.FrontalMelt==1 %Add front melt to calve rate. Just uses vertical melt in calving front cell.....may need to improve this in future
        CR(glMASK==5)=CR(glMASK==5)+FMR(glMASK==5);
    end

    CRx=CR.*XUV;
    CRy=CR.*YUV;

    wx=uxh+CRx;
    wy=uyh+CRy;

    LSF=LSFfunction(LSF,ctr,wx,wy,node,nodes,VM,MASK); %Advect calving front position

    %     if ctr.CalveGround==1 %Allow calving of grounded ice
    %         LSF(glMASK<3)=1;
    %     end

    if ctr.CalveCirc==1 %Impose maximum calving front extent from a field, ctr.CF_Boundary
        load (ctr.CF_Boundary,'CIRC')
        LSF(CIRC<1)=-1;
    end

    if floor(cnt/par.LSFReset)==ceil(cnt/par.LSFReset)  % Reset LSF field for stability
        LSF(LSF<-1)=-1;
        LSF(LSF>1)=1;
    end

end
end