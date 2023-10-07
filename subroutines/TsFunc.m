function [Ts]=TsFunc(ctr,par,Tsf,S0,sn,DeltaSL,DeltaT)

% Kori-ULB
% Surface temperature parametrizations
% 0: no ice-elevation feedback
% 1: ice-elevation feedback   

    switch ctr.TsType
        case 0
            % no correction for elevation changes
            Ts=Tsf;
        case 1
            % correction for elevation changes
            Ts=Tsf+par.Tlapse*(max(sn,DeltaSL)-S0);
        case 2
            % EISMINT moving margin
            Tsf=270.-0.01*sn-par.T0+DeltaT;
            Ts=Tsf;
        case 3
            % McCall glacier
            Ts=-2.08-0.00335*sn;
            Ts(sn>ctr.ELA)=-6.3;
    end
end


