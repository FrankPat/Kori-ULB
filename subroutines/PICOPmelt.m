function Melt=PICOPmelt(HB,sina,Zgl,Bk,ShelfN,T,S,par,ctr,Bmax,numsh)

% Kori-ULB
% Sub-shelf melt according to the PICOP model (combination of PICO and
% Plume model. Note: the original Plume model is used here and not the
% revised version PICO2019

    % Define Ta, Sa for all ice shelves based on Box model
    Ta=zeros(ctr.imax,ctr.jmax)+par.Toi;
    Sa=zeros(ctr.imax,ctr.jmax)+par.Soi;
    for i=1:numsh
        for k=1:Bmax(i)
            Ta(Bk==k & ShelfN==i)=T(i,k);
            Sa(Bk==k & ShelfN==i)=S(i,k);
        end
    end
    % Limit on Ta
    Ta=max(Ta,par.lambda1*Sa+par.lambda2+par.lambda3*HB);
    
    Tfgl=par.lambda1*Sa+par.lambda2+par.lambda3*Zgl;
    CdGamTS=par.CdGamT*(par.gamma1+(par.gamma2*(Ta-Tfgl)/par.lambda3).*(par.Eo*sina./ ...
        (par.CdGamTS0+par.Eo*sina)));
    galfa=((sina./(par.Cd+par.Eo*sina)).^(0.5)).*((CdGamTS./(CdGamTS+par.Eo*sina)).^(0.5)).* ...
        (par.Eo*sina./(CdGamTS+par.Eo*sina));
    len=(Ta-Tfgl)/par.lambda3.*(par.x0*CdGamTS+par.Eo*sina)./(par.x0*(CdGamTS+par.Eo*sina));
    Xhat=(HB-Zgl)./len;
    Xhat=min(Xhat,1);
    Melt=ctr.M0*MscaledPoly(Xhat,par.pcof).*galfa.*(Ta-Tfgl).^2;
    
end

    
