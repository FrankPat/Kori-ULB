function [SLR,dSLR]=CalculateGeoid(par,ctr,DeltaSL,SLC,SLR,SLR0,Bn,Hn,Pg0,frg,geoide)

% Kori-ULB
% Calculation of geoid change due to mass changes of the ice sheet. The
% effect is stored in SLR (spatial variability of global sea level)

    [~,MASKn,~,~]=Floatation(par,Bn,SLR,Hn,zeros(ctr.imax,ctr.jmax));
    Pg1=zeros(ctr.imax,ctr.jmax);
    Pg1(MASKn==1)=par.rho*Hn(MASKn==1)*ctr.delta^2.;
    Pg1(MASKn==0)=(-par.rhow)*(Bn(MASKn==0)-SLR(MASKn==0)+DeltaSL)* ...
        ctr.delta^2.; % Takes into account previous geoid change but not
                      % DeltaSL (external SL forcing) as ocean mass
    Pg1=Pg1+par.rhom*Bn*ctr.delta^2.; % Takes into account land deformation
    DPg=Pg1-Pg0;
    % Fingerprint
    MaP=zeros(ctr.imax+2*frg,ctr.jmax+2*frg);
    MaP(frg+1:frg+ctr.imax,frg+1:frg+ctr.jmax)=MASKn;
    loadH=conv2fft(DPg,geoide); % FFT convolution
    % Determine local sea-level change
    Vd=sum(loadH(MaP==0))*ctr.delta^2; % Volume change in domain without mass addition
    loadG=loadH(frg+1:ctr.imax+frg,frg+1:ctr.jmax+frg);
    dSLR=loadG+SLC-Vd/par.Aoc; % previously VAF instead of SLC
    SLR=SLR0+dSLR+DeltaSL; % SLR = Initial SLR + Geoid perturbation + External SL forcing
        
end


