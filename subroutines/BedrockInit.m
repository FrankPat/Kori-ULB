function [B0,load0,frb,kei,Ll]=BedrockInit(ctr,par,H,B,SLR,MASK)

% Kori-ULB
% Initialization of bedrock loads (considering initial bedrock in isostatic
% equilibrium) for constant and variable flexural rigidity

    if ctr.Dbexist==0
        % Constant flexural rigidity
        Ll=(par.FlexRigid/par.rhom/par.g)^0.25;
        frb=round(6*Ll/ctr.delta);
        P=zeros(ctr.imax,ctr.jmax);
        P(MASK==1)=par.rho*par.g*H(MASK==1)*ctr.delta^2.;
        P(MASK==0)=(-par.rhow)*par.g*(B(MASK==0)-SLR(MASK==0))*ctr.delta^2.;
        P0=NaN(ctr.imax+2*frb,ctr.jmax+2*frb);
        P0(frb+1:frb+ctr.imax,frb+1:frb+ctr.jmax)=P;
        P0(isnan(P0))=P(ctr.imax,ctr.jmax);
        [kerx,kery]=meshgrid(-frb*ctr.delta/Ll:ctr.delta/Ll:frb*ctr.delta/Ll, ...
            -frb*ctr.delta/Ll:ctr.delta/Ll:frb*ctr.delta/Ll);
        ker=sqrt((kerx.^2+kery.^2));
        ker(ker<1e-8)=1e-8;
        kei=imag(besselk(0,ker*(1+1i)/sqrt(2)));
        load0=-xcorr2(P0,kei)*Ll^2./(2*pi*par.FlexRigid);
        load0=load0(2*frb+1:ctr.imax+2*frb,2*frb+1:ctr.jmax+2*frb);
        B0=B;
    else
        frb=false;
        kei=false;
        Ll=false;
        % Varying flexural rigidity (Kevin)
        load0=zeros(ctr.imax,ctr.jmax);
        load0(MASK==1)=par.rho*par.g*H(MASK==1);
        load0(MASK==0)=(-par.rhow*par.g)*(B(MASK==0)-SLR(MASK==0));
        B0=B;
    end
end


