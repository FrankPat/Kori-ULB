function [Asf,Asfx,Asfy,Asfd,Tbc,Neff,pwv,Wtil,r]=BasalSliding(ctr,par, ...
    As,Tb,H,B,MASK,Wd,Wtil,Bmelt,flw,bMASK,bMASKm,bMASKx,bMASKy)

% Kori-ULB
% Basal sliding module that determines Asf (prefactor in the sliding law,
% based on subglacial temperature and effective pressure

    if ctr.Tcalc>0
        Tbc=Tb+par.pmp*H;
    else
        Tbc=false;
    end
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
    if ctr.subwaterflow==0
        % Tsai et al (2015)
        pwv=-par.PoreFrac*par.rhow*par.g*B;
        pwv(MASK==0)=par.PoreFrac*par.rho*par.g.*H(MASK==0);
        pwv(B>=0)=0;
    else
        pwv=par.PoreFrac*par.rho*par.g*H.*(Wd/par.Wdmax); % Bueler and Brown (2009)  
    end
    Po=max(par.rho*par.g*H,1e5);
    Neff=max(Po-r.*pwv,par.sigmat*Po);
    if ctr.subwaterflow==2
        [Neff,Wtil]=TillWater(MASK,Wtil,Bmelt,ctr,Po,par);
    end
    if ctr.subwaterflow==3
        Neff=ones(ctr.imax,ctr.jmax)*par.NeffScale^ctr.p;
        expflw=exp(flw/par.flw0);
    else
        expflw=1;
    end
    Asf=par.AsScale*((1-r).*par.AsFroz+As.*r).*expflw./((Neff.^ctr.p)/ ...
        (par.NeffScale^ctr.p));
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


