function [Melt,butfac,H]=SubShelfMelt(ctr,fc,par,Tf,To,So,TF,butfac, ...
    HB,glMASK,H,B,ShelfN,numsh,shMASK,MASK,MASKlk,uxssa,uyssa,arcocn, ...
    MeltInv,cnt)

% Kori-ULB
% Options for sub-shelf melt calculation
% Quadratic, linear, PICO, PCIO, Plume, ISMIP6, ABUMIP, ...

    if islogical(TF)==1
        Tfm=max(To,Tf)-Tf;
    else
        Tfm=TF; % Use imposed thermal forcing
    end
    switch ctr.meltfunc
        case 1
            % Linear melt forcing
            Melt=ctr.gammaT*par.LatentMelt*Tfm*par.secperyear;
        case 2
            % Quadratic local melt forcing - Favier et al. (2019)
            Melt=ctr.gammaT*par.LatentMelt^2*abs(Tfm).*(Tfm)*par.secperyear;
        case 21
            % Quadratic semilocal melt forcing - Favier et al. (2019)
            mean_TF=zeros(ctr.imax,ctr.jmax);
            for i=1:numsh
                mean_TF(ShelfN==i)=mean(Tfm(ShelfN==i));
            end
            Melt=ctr.gammaT*par.LatentMelt^2*abs(mean_TF).*(Tfm)*par.secperyear;
        case 22
            % Quadratic semilocal melt forcing with local slope - Favier et
            % al. (2019)
            sina=ComputeAlphaLocal(HB, shMASK, ctr);
            mean_TF=zeros(ctr.imax,ctr.jmax);
            for i=1:numsh
                mean_TF(ShelfN==i)=mean(Tfm(ShelfN==i));
            end
            Melt=ctr.gammaT*par.LatentMelt^2*abs(mean_TF).*(Tfm).*sina* ...
                par.secperyear;
        case 23 % Quadratic local melt forcing with mean antarctic slope - Burgard et al. (2022)
            sin_theta = ComputeAlphaLocal(HB, shMASK, ctr);
            sin_theta=mean(sin_theta(shMASK==1)); % mean Antarctic slope - Burgard et al. (2022) value: 2.9e-3
            gammaT_B = (par.cp0/par.Latent) * par.BetaS .* So * par.g/(2*abs(par.f_coriolis)) * sin_theta;
            Melt = ctr.gammaT .* gammaT_B .* (par.rhow.*par.cp0./(par.rho.*par.Latent)) .* abs(Tfm) .* Tfm .*par.secperyear;
        case 24 % Quadratic local melt forcing with local slope - Burgard et al. (2022)
            sin_theta = ComputeAlphaLocal(HB, shMASK, ctr);
            gammaT_B = (par.cp0./par.Latent) .* par.BetaS .* So .* par.g./(2.*abs(par.f_coriolis)) .* sin_theta;
            Melt = ctr.gammaT .* gammaT_B .* (par.rhow.*par.cp0./(par.rho.*par.Latent)) .* abs(Tfm) .* Tfm .*par.secperyear;
        case {3,4}
            % PICO/PICOP
            [Bk,Ak,Bmax]=PICOsetup(glMASK,H,ctr,par,ShelfN,numsh,shMASK);
            [T0o,S0o,zb]=OceanVarBox(numsh,To,So,ShelfN,HB);
            T0o=max(T0o,par.lambda1*S0o+par.lambda2+par.lambda3*zb);
            [T,S]=BoxModel(S0o,T0o,numsh,par,Ak,Bk,ShelfN,HB, ...
                1/par.LatentMelt,ctr);
            if ctr.meltfunc==3
                % PICO
                Melt=PICOMelt(HB,shMASK,ShelfN,Ak,Bk,S0o,T0o,S,T,numsh, ...
                    Bmax,ctr,1/par.LatentMelt,par);
            else
                % PICOP
                [MASKpicop,MASKb]=PICOPmasks(MASK,MASKlk,H,ctr,par);
                [Zgl]=AdvecGL(HB,B,MASK,MASKpicop,MASKb,uxssa,uyssa,ctr);
                sina=ComputeAlphaLocal(HB, shMASK, ctr);
                Melt=PICOPmelt(HB,sina,Zgl,Bk,ShelfN,T,S,par,ctr,Bmax,numsh);
            end
        case 5
            % Plume model
            % Mean ocean properties by shelf cavity (Burgard et al. 2022)
            To_mean_cav=zeros(ctr.imax,ctr.jmax);
            So_mean_cav=zeros(ctr.imax,ctr.jmax);
            for b=1:numsh
                avg_To_int=mean(To(ShelfN==b));
                avg_So_int=mean(So(ShelfN==b));
                To_mean_cav(ShelfN==b)=avg_To_int;
                So_mean_cav(ShelfN==b)=avg_So_int;
            end
            [MASKpicop,MASKb]=PICOPmasks(MASK,MASKlk,H,ctr,par);
            [Zgl]=AdvecGL(HB,B,MASK,MASKpicop,MASKb,uxssa,uyssa,ctr);
            sina=ComputeAlphaLocal(HB, shMASK, ctr);
            Melt=PlumeMelt2019(To_mean_cav,So_mean_cav,HB,Zgl,sina,par,ctr);
        case 6
            % Cornford 2016 melt (function of H)
            Melt=max(min((H-100)*4/7,400),0);
        case 7
            % Uniform melting when butfac=1 - value of meltfac
            GeoBut=ones(ctr.imax,ctr.jmax)*butfac;
            butfac=1;
            Melt=ctr.meltfac*GeoBut;
        case 8
            % Ice shelf removal (needs butfac=[0])
            Melt=zeros(ctr.imax,ctr.jmax);
            if cnt==1
                H(MASK==0)=par.SeaIceThickness;
                Melt(MASK==0)=ctr.meltfac;
            else % can be modified for regrowth of ice shelves
                H(MASK==0)=par.SeaIceThickness;
                Melt(MASK==0)=ctr.meltfac;
            end
        case 9
            % ISMIP6 non-local melt rate parameterization
            Melt=ISMIP6ocean(par,ctr,HB,shMASK,Tfm,fc.basinNumber, ...
                fc.deltaT_basin);
        case 91
            % ISMIP6 non-local melt rate parameterization + dependy on local slope
            sina=ComputeAlphaLocal(HB, shMASK, ctr); 
            Melt=ISMIP6OceanSlope(sina,par,ctr,HB,shMASK,Tfm, ...
                fc.basinNumber,fc.deltaT_basin);
        case 92
            % non-local melt rate parameterization + mean ant slope
            sina=zeros(ctr.imax,ctr.jmax)+2.9e-3; % mean ant slope -- Burgard 2022
            deltaT_basin=zeros(ctr.imax,ctr.jmax); % no local correction
            Melt=ISMIP6OceanSlope(sina,par,ctr,HB,shMASK,Tfm, ...
                fc.basinNumber,deltaT_basin);
        case 10
            % MISMIP+ melt function
            Melt=0.2*tanh((HB-B)/75).*max(-100-HB,0)*ctr.meltfac;
            Melt(MASK==1)=0;
        case 11
            % run model with optimized sub-shelf melt
            Melt=MeltInv;
            Melt(MASK==1)=0;
    end
    switch ctr.meltfunc
        case {1,2,21,22,23,24,3,4,5,9,91,92}
            if par.ArcOcean==1
                Melt(Melt>0)=Melt(Melt>0).*max(0,min(1,(arcocn(Melt>0)-20)./20));
            end
            if ctr.meltfac>0
                Melt=Melt*ctr.meltfac; % Add meltfac to increase melt
            end
    end
    Melt(shMASK==0)=0; % Melt only for real shelves, Melt=0 for 'Lakes'
    Melt(glMASK==6)=0;

end


