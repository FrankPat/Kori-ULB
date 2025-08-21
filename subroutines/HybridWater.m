function Neff=HybridWater(flw,dphi0,Po,A,ub,par,kappa)

% Kori-ULB
% Generalized effective pressure in hard/soft bed based on equilibrium 
% between melt, sliding, and creep.

    % Convert to SI units (per second instead of per year)
    flw = flw/par.secperyear;
    ub  = ub /par.secperyear;
    A   = A  /par.secperyear;
    
    % Parameters & coefficients
    alpha = par.alpha;         % Darcy-Weisbach exponent
    beta  = par.beta;          % Darcy-Weisbach exponent
    f     = par.f;             % Friction coefficient
    hb    = par.hb;            % Size of bed protrusions
    r     = par.distChannels;  % Typical spacing between channels
    hc    = par.hc;            % Thickness of canals
    Qc    = par.Qc;            % Critical flux for transition to canals

    c1    = 1/(par.rho*par.Latent);
    c2    = 2*A*(par.n^(-par.n));
    c3    = (2/pi)^(0.25)*sqrt((pi+2)/(par.rhow*f));
    
    % Flow in a single channel
    Q = flw*r;
    
    % Surface of a channel
    S = c3^(-1/alpha).*Q.^(1/alpha).*dphi0.^(-(beta-1)/alpha);

    % Equivalent thickness for soft bed
    hSoft = hc + (S.^(1/2)/par.FactDefTill - hc).*exp(-Q/Qc);

    % Equivalent thickness for hard bed
    hHard = S.^(1/2);

    % Equivalent thickness for hybrid hard/soft bed
    h = kappa.*hSoft + (1-kappa).*hHard;

    % Effective pressure
    Neff = (c1*Q.*dphi0 + ub*hb)./(c2.*(S./h).^2);
    Neff = Neff.^(1/par.n);

    % Limits
    Neff = max(min(Po,Neff),par.sigmat*Po);

end


