function Melt=ISMIP6OceanSlope(sina,par,ctr,ice_draft,icemask_shelves, ...
    TF_draft,basinNumber,deltaT_basin)

% Kori-ULB
% ISMIP6 ocean parameterization equation (1) + adapted to
% take the local slope into account
 
    cste = (par.rhow*par.cp0/(par.rho*par.Latent)).^2;  % in K^(-2)

    Nbasin=max(max(basinNumber))+1;

    Melt=zeros(ctr.imax,ctr.jmax);
    mean_TF=zeros(Nbasin,1);

    basinNumber=basinNumber+1;
    for i=1:Nbasin
        mean_TF(i)=mean(TF_draft(basinNumber==i & icemask_shelves>5.e-1 ...
            & ice_draft<=0.));
    end

    % Step3 calculation of melting rate
    % melt rate in m/yr (meters of pure water per year):
    % [*rhofw_SI/rhoi_SI to get it in meters of ice per year]

    for i=1:Nbasin
        Melt(basinNumber==i)=ctr.gammaT.*cste.*(TF_draft(basinNumber==i) ...
            +deltaT_basin(basinNumber==i)).*abs(mean_TF(i)+ ...
            deltaT_basin(basinNumber==i)).*sina((basinNumber==i));
    end
    Melt(TF_draft<=-5)=0;
    Melt=Melt*par.rhof/par.rho;
end


