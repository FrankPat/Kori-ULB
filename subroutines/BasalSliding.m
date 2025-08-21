function [Asf,Asfx,Asfy,Asfd,Neff,Wtil,r,expflw]=BasalSliding(ctr,par, ...
    A,As,Tbc,H,B,MASK,ub,Wd,Wtil,Bmelt,flw,Neff,expflw,kappa,updateHydro, ...
    bMASK,bMASKm,bMASKx,bMASKy)

% Kori-ULB
% Basal sliding module that determines Asf (prefactor in the sliding law,
% based on subglacial temperature and effective pressure
% extended with subglacial model from Kazmierczak et al. (2024)

    if ctr.SlidAdjust==1 && ctr.Tcalc>0
        if ctr.stdBexist==1
            Tr=min(par.TrTemp,-0.02*stdB); % Pollard (2015)
        else
            Tr=par.TrTemp;
        end
        r=max(0,min(1,(Tbc-Tr)./(-Tr))); % r on H-grid
        r(MASK==0)=1;
        if ctr.schoof>=1 || ctr.shelf==1
            r(glMASK==2)=1; % grounding line
        end
    else
        r=ones(ctr.imax,ctr.jmax);
    end

    if (updateHydro)
        expflw = ones(size(flw));
        Neff   = par.NeffScale * ones(size(flw));
        Po     = max(par.rho*par.g*H,1e5);

        % Parametrizations in Kazmierczak et al. (2022)
        if(ctr.subwaterflow <= 3)
            switch(ctr.subwaterflow)
                case 0
                    % HAB - Tsai et al. (2015)
                    pwv          = -par.PoreFrac*par.rhow*par.g*B;
                    pwv(MASK==0) =  par.PoreFrac*par.rho *par.g*H(MASK==0);
                    pwv(B>=0)    = 0;
                    Neff         = max(Po-r.*pwv,par.sigmat*Po);
                case 1
                    % SWD - Bueler and Brown (2009)
                    pwv  = par.PoreFrac*par.rho*par.g*H.*(Wd/par.Wdmax);
                    %pwv(MASK==0) =  par.PoreFrac*par.rho *par.g*H(MASK==0);
                    Neff = max(Po-r.*pwv,par.sigmat*Po);
                case 2
                    % TIL - Bueler and van Pelt (2015)
                    [Neff,Wtil] = TillWater(MASK,Wtil,Bmelt,ctr,Po,par);
                case 3
                    % SWF - Goeller et al. (2013)
                    expflw = exp(flw/par.flw0);
            end
        % HARD/SOFT/HYBRID
        else
            % Gradient of phi0
            phi0 = par.rho*par.g*H + par.rhow*par.g*B;
            Nconv = par.convWdwPhi;
            phi0smooth = conv2(phi0,(1/Nconv^2)*ones(Nconv),'same');
            dphi0dy = (phi0smooth-circshift(phi0smooth,[-1 0]))/ctr.delta;
            dphi0dx = (phi0smooth-circshift(phi0smooth,[0 -1]))/ctr.delta;
            dphi0ds  = sqrt(dphi0dx.^2 + dphi0dy.^2);
            
            % Far-field effective pressure
            switch(ctr.subwaterflow)
                case 4  % SCF - Schoof (2010)
                    NInf = SchoofWaterFarField(flw,dphi0ds,Po,A,ub,par);
                case 5  % HARD
                    NInf = HybridWater(flw,dphi0ds,Po,A,ub,par,0);
                case 6  % SOFT
                    NInf = HybridWater(flw,dphi0ds,Po,A,ub,par,1);
                case 7  % HYBRID
                    NInf = HybridWater(flw,dphi0ds,Po,A,ub,par,kappa);
            end

            % GL correction
            Neff = NInf.*(1 - erfc(sqrt(pi)/2*phi0./NInf));
            Neff(MASK == 0) = par.sigmat*Po(MASK == 0);
        end
    end
 
    effectHydro = max(par.effectHydroLimit, ((Neff/par.NeffScale).^(ctr.p))./expflw);
    Asf=par.AsScale*((1-r).*par.AsFroz+As.*r)./effectHydro;
    Asfx=0.5*(Asf+circshift(Asf,[0 -1]));
    Asfy=0.5*(Asf+circshift(Asf,[-1 0]));
    Asfd=h2d(Asf);
    if ctr.basin==1 %VL: necessary for SIA runs for basins
        Asf(bMASK==1)=par.AsScale*par.As0;
        Asfx(bMASKx==1)=par.AsScale*par.As0;
        Asfy(bMASKy==1)=par.AsScale*par.As0;
        Asfd(bMASKm==1)=par.AsScale*par.As0;
    end
end


