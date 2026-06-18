function [runoff,acc,rain,Smelt]=PDDmonthly(Ts_yc,Pr_yc,par,ctr,sn)

% Kori-ULB
% Scheme for monthly PDD calculation

    Tm=Ts_yc; % Ts yearly cycle
    Ts=mean(Ts_yc,3); % Yearly mean
    Pr=Pr_yc; % Pr yearly cycle
    par.PDDsteps=size(Tm,3);

    % PDD calculation from Calov and Greve (2005)
    nT=(Tm-par.PDDth)/(sqrt(2)*par.Tsigma);
    T=(par.Tsigma/sqrt(2*pi)*exp(-(nT.^2))+0.5*(Tm-par.PDDth).*erfc(-nT)) ...
        +par.PDDth;
    PDD=T*365.25/par.PDDsteps;

    % % Rain fraction calculation
    % nP=(Tm-par.Tsnow)/(sqrt(2)*par.Psigma);
    % P=(par.Psigma/sqrt(2*pi)*exp(-(nP.^2))+0.5*(Tm-par.Tsnow).*erfc(-nP)) ...
    %     +par.Tsnow;
    % R=(1-max(0,min((par.Train-P)/(par.Train-par.Tsnow),1)));

    %% Rain fraction calculation following Hantel et al. (2000) %%%
    Trs_loc = par.Trs - par.gamma_rs .* (sn ./ 1000);       % If par.gamma_rs>0, transition temperature to depend on elevation
    Sfrac = 0.5 .* (1 - tanh(par.Sp .* (Tm - Trs_loc)));   % snow fraction
    R     = 1 - Sfrac;
    
    % Melt model
    snow_depth=zeros(ctr.imax,ctr.jmax);
    nyear=2; % Sufficient to initialize snow depth

    % Loop over timesteps through one year
    for y=1:nyear % run through several years to equilibrium (Tsai 2020)
        Smelt=zeros(ctr.imax,ctr.jmax); % Annual ice  and snow melt
        acc=zeros(ctr.imax,ctr.jmax);
        rain=zeros(ctr.imax,ctr.jmax);
        ESM=zeros(ctr.imax,ctr.jmax);
        EIM=zeros(ctr.imax,ctr.jmax);
        for nstep=1:par.PDDsteps
            PR=Pr(:,:,nstep)/par.PDDsteps; % Total Precipitation during nstep
            acc_nstep=max(0,PR.*(1-R(:,:,nstep))); 
            % Total Snow Accumulation during nstep. R is the fraction 
            % of precipitation that falls as rain
            rain_nstep=max(0,PR.*R(:,:,nstep)); % Total Rain during nstep
            snow_depth=max(0,snow_depth+acc_nstep);
            ESM_nstep=min(max(snow_depth,0),par.snowfac.*PDD(:,:,nstep)); 
            % Effective snow melt (cannot exceed the amount of snow) 
            % in m ice equivalent
            PDDr=max(PDD(:,:,nstep)-ESM_nstep./par.snowfac,0); 
            % remaining PDDs after snow melt
            snow_depth=max(0,snow_depth-ESM_nstep);
            EIM_nstep=par.icefac.*PDDr; % Effective ice melt
            Smelt_nstep=ESM_nstep+EIM_nstep;

            % Increment annual quantities
            ESM=ESM+ESM_nstep;
            EIM=EIM+EIM_nstep;
            Smelt=Smelt+Smelt_nstep; % combined ice and snow melt
            acc=acc+acc_nstep;
            rain=rain+rain_nstep;
        end
    end

    % Refreezing model: Simple dynamic parameterization or the refreezing process - Huybrechts & DeWolde 1999
    W_r = ESM + rain + par.f_EIM_ref .* EIM; % Available water mass - maximum amount of water that can be retained provided enough energy is available
    P_r = max(-(min(Ts,0))*par.d_ice*par.K/(par.kdif*par.rho*par.Latent),0); % Potential retention mass - maximum amount of liquid water that can be retained provided enough mass is available
    E_r = min(P_r, W_r); % Effective Refreezing - amount of water that is effectively retained at the end of the melting season
    runoff = (1 - par.f_EIM_ref) .* EIM + W_r - E_r;  % ice melt + rain and snowmelt liquid water beyond the snowpack saturation
 
end


