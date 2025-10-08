function [Pr,Evp,runoff]=MbFunc(ctr,par,Prf,Evpf,runoff,runofff,Tsf,sn,S0,X,Y,Lj,Li,DeltaT,DeltaSL)

% Kori-ULB
% Surface mass balance parametrizations and corrections (when input files
% are read)

    switch ctr.MbType
        case 0
            % no correction for elevation changes
            Pr=Prf;
            Evp=Evpf;
            if ctr.PDDcalc==0
                runoff=runofff;
            end
        case 1
            % Correction for elevation changes - OPTION 1
            % correction of precipitation for changes in background temperature
            % Clausius-Clapeyron relationship
            Pr=Prf.*(1+par.Prfac).^DeltaT;
            if ctr.PDDcalc==0
                if ctr.runoffcorr==1 % option to correct externally-forced runoff field for elevation change
                    runoff=runofff+1.805*(exp(0.5745*(Tsf+par.Tlapse*(max(sn,DeltaSL)-S0)))-exp(0.5745*Tsf)); % inferred from MAR monthly data: runoff=1.805*exp(0.5745*t2m)
                else
                    runoff=runofff;
                end
            end
            Evp=Evpf; % Evp not corrected for elevation change
        case 2
            % Correction for elevation changes - OPTION 1 
            % Pollard et al., Garbe2020 - correction of precipitation for
            % elevation changes and change in background temperature
            Pr=Prf.*exp(0.05*(par.Tlapse*(max(sn,DeltaSL)-S0)+DeltaT));
            if ctr.PDDcalc==0
                if ctr.runoffcorr==1 % option to correct externally-forced runoff field for elevation change
                    runoff=runofff+1.805*(exp(0.5745*(Tsf+par.Tlapse*(max(sn,DeltaSL)-S0)))-exp(0.5745*Tsf)); % inferred from MAR monthly data: runoff=1.805*exp(0.5745*t2m)
                else
                    runoff=runofff;
                end
            end
            Evp=Evpf; % Evp not corrected for elevation change
        case 3
            % EISMINT fixed margin with forcing
            Pr=Prf+0.02*DeltaT;
            Evp=Evpf;
            runoff=runofff;
        case {4,5}
            % EISMINT moving margin
            s=1.e-2;
            Rel=450.;
            if ctr.MbType==5
                Rel=Rel+10*DeltaT;
            end
            dist=sqrt((X-Lj/2.).^2.+(Y-Li/2.).^2.);
            Pr=min(0.5,s*(Rel-dist));
            Evp=Evpf;
            runoff=runofff;
        case 6
            % McCall Glacier
            Pr=ctr.mbgrad*(sn-ctr.ELA); %0.0017
            Pr(sn>ctr.ELA)=ctr.mbgrad1*(sn(sn>ctr.ELA)-ctr.ELA);
            Evp=Evpf;
            runoff=runofff;
    end
end


