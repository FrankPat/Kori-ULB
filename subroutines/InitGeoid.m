function [Pg0,geoide,frg]=InitGeoid(par,ctr,MASK,H0,B0,B,SLR)

% Kori-ULB
% Initialization of geoid fields

    Pg0=zeros(ctr.imax,ctr.jmax);
    Pg0(MASK==1)=par.rho*H0(MASK==1)*ctr.delta^2.;
    Pg0(MASK==0)=(-par.rhow)*(B(MASK==0)-SLR(MASK==0))*ctr.delta^2.;
    Pg0=Pg0+par.rhom*B0*ctr.delta^2.; % addition of bed change (if applicable)
    [gkx,gky]=meshgrid(-par.geoidist:ctr.delta:par.geoidist,-par.geoidist:ctr.delta:par.geoidist);
    frg=round(par.geoidist/ctr.delta-0.5);
    dist=sqrt(gkx.^2+gky.^2);
    geoide=par.Re./(2*par.Me*sin(dist/(2*par.Re)));
    geoidmax=par.Re./(2*par.Me*sin(10/(2*par.Re)));
    geoide(geoide>geoidmax)=0;
    geoide(geoide<1.5e-18)=0;
    
end


