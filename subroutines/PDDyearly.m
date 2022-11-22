function [runoff,acc,rain,Smelt]=PDDyearly(ctr,par,Ts,Pr,lat,MASK,sn)

% Kori-ULB
% Scheme for yearly PDD model caculation

    % Yearly amplitude of the signal (applies to the Antarctic ice sheet)
    if islogical(lat)==1
        Ta=min(20+max(sn,0)/300.,30); % seasonal peak-to-peak amplitude in Ts
    else
        Ta=zeros(ctr.imax,ctr.jmax);
        Ta(MASK==1)=-11.06+0.003204*sn(MASK==1)-0.3584*lat(MASK==1);
        Ta(MASK==0)=-45.06-0.04598*sn(MASK==0)-0.969*lat(MASK==0);
        Ta(Ta<15)=15;
    end
    
    % Define idealized yearly cycle Tm:
    t=repmat(reshape(linspace(1,365.25,par.PDDsteps),1,1,par.PDDsteps), ...
        [ctr.imax,ctr.jmax,1]);
    Tm=repmat(Ts,[1,1,par.PDDsteps])-0.5*repmat(Ta,[1,1,par.PDDsteps]).* ...
        sin(2*pi*t/365.25);
    % PDD calculation from Calov and Greve (2005)
    nT=(Tm-par.PDDth)/(sqrt(2)*par.Tsigma);
    T=(par.Tsigma/sqrt(2*pi)*exp(-(nT.^2))+0.5*(Tm-par.PDDth).* ...
        erfc(-nT))+par.PDDth;

    % Rain fraction calculation
    nP=(Tm-par.Tsnow)/(sqrt(2)*par.Psigma);
    P=(par.Psigma/sqrt(2*pi)*exp(-(nP.^2))+0.5*(Tm-par.Tsnow).* ...
        erfc(-nP))+par.Tsnow;

    PDD=sum(T,3)*365.25/par.PDDsteps;
    R=sum(1-max(0,min((par.Train-P)/(par.Train-par.Tsnow),1)),3)/par.PDDsteps;
    
    % Melt model
    acc=max(0,Pr.*(1-R)); % Snow Accumulation Rate. R is the fraction of 
                          % precipitation that falls as rain
    rain=max(0,Pr.*R);
    ESM=min(max(acc,0),par.snowfac.*PDD); 
    % Effective snow melt (cannot exceed the amount of snow)
    PDDr=max(PDD-ESM./par.snowfac,0); % remaining PDDs after snow melt
    EIM=par.icefac.*PDDr; % Effective ice melt
    Smelt=ESM+EIM; % total surface melt

    % Refreezing model
    W_r=ESM+rain; % Available water mass - maximum amount of water that 
                  % can be retained provided enough energy is available
    P_r=max(-(min(Ts,0))*par.d_ice*par.K/(par.kdif*par.rho*par.Latent),0); 
    % Potential retention mass - maximum amount of liquid water that can 
    % be retained provided enough mass is available
    E_r=min(Pr,min(P_r,W_r)); 
    % Effective Refreezing - amount of water that is effectively 
    % retained at the end of the melting season
    % Er is limited by the annual precipitation
    runoff=EIM+W_r-E_r; % ice melt + rain and snowmelt liquid water beyond 
                        % the snowpack saturation

end


