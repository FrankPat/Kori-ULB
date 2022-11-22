function [flw,Wd]=SubWaterFlux(ctr,par,H,HB,MASK,Bmelt)

% Kori-ULB
% Steady-state subglacial water flow model
% Method due to LeBrocq (GDS Warner model)
% Simplified implementation for ice thickness classes as a function of
% the distance over hich hydraulic gradient coupling takes place
% Water flux and depth calculated on h-grid

    flux=zeros(ctr.imax,ctr.jmax)-1;
    pot=par.rho*par.g*H+par.rhof*par.g*HB;
    pot(MASK==0)=0;
    [pot]=PotentialFilling(pot,ctr);
    % Correction of ice thickness due to filling of 'pot'
    H=(pot-par.rhof*par.g*HB)/(par.rho*par.g);
    % Hydraulic potential gradients on staggered d-grid
    gdsx=-par.rho*par.g*(circshift(H,[0 -1])-circshift(H,[0 1]))/ ...
        (2*ctr.delta)-par.rhof*par.g*(circshift(HB,[0 -1])- ...
        circshift(HB,[0 1]))/(2*ctr.delta);
    gdsy=-par.rho*par.g*(circshift(H,[-1 0])-circshift(H,[1 0]))/ ...
        (2*ctr.delta)-par.rhof*par.g*(circshift(HB,[-1 0])- ...
        circshift(HB,[1 0]))/(2*ctr.delta);

    % define multiplier for GDS
    Hlev=max(mean(H(MASK==1)),10);
    scale=Hlev*par.longcoupwater*2.; % scale of GDS window
    width=2.*scale;
    scale(width<=ctr.delta)=ctr.delta/2+1;
    maxlevel=2*round(width/ctr.delta-0.5)+1;
    frb=(maxlevel-1)/2;
    [ew,ns]=meshgrid(1:maxlevel,1:maxlevel); % size of GDS window

    kegdsx0=zeros(ctr.imax+2*frb,ctr.jmax+2*frb);
    kegdsy0=zeros(ctr.imax+2*frb,ctr.jmax+2*frb);
    kegdsx0(frb+1:frb+ctr.imax,frb+1:frb+ctr.jmax)=gdsx;
    kegdsy0(frb+1:frb+ctr.imax,frb+1:frb+ctr.jmax)=gdsy;

    % Define mulipliers as a function of ice thickness (scale) and convolute
    dist=sqrt((ctr.delta*(ew-frb-1)).^2+(ctr.delta*(ns-frb-1)).^2)/scale;
    multiplier=max(0.,1.-dist/2.);
    multiplier=multiplier/sum(sum((multiplier)));
    kegx=xcorr2(kegdsx0,multiplier);
    kegy=xcorr2(kegdsy0,multiplier);
    kegdsx=kegx(2*frb+1:ctr.imax+2*frb,2*frb+1:ctr.jmax+2*frb);
    kegdsy=kegy(2*frb+1:ctr.imax+2*frb,2*frb+1:ctr.jmax+2*frb);

    % Recursive algorithm
    gdsmag=abs(kegdsx)+abs(kegdsy);
    for i=1:ctr.imax-1
        for j=1:ctr.jmax-1
            if MASK(i,j)==1
                flux=DpareaWarGds(i,j,MASK,Bmelt, ...
                    kegdsx,kegdsy,par,ctr,flux,gdsmag,1);
            end
        end
    end

    corfac=zeros(ctr.imax,ctr.jmax)+1;
    denom=sqrt(kegdsx.^2+kegdsy.^2);
    corfac(gdsmag>0)=gdsmag(gdsmag>0)./denom(gdsmag>0);
    flw=min(max(flux./corfac/ctr.delta,0),1e5);
    flw(MASK==0)=0;
    Wd=min(par.Wdmax,max(par.Wdmin,(12*par.waterviscosity*flw/ ...
        mean(gdsmag(MASK==1))).^(1/3)));
end


