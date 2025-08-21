function Neff=SchoofWaterFarField(flw,dphi0,Po,A,ub,par)

% Kori-ULB
% Effective pressure in hard bed based on equilibrium between melt,
% sliding, and creep.

    % Convert to SI units (per second instead of per year)
    flw = flw/par.secperyear;
    ub  = ub /par.secperyear;
    A   = A  /par.secperyear;
    
    % Parameters & coefficients
    alpha = par.alpha;         % Darcy–Weisbach exponent
    beta  = par.beta;          % Darcy–Weisbach exponent
    f     = par.f;             % Friction coefficient
    hb    = par.hb;            % Size of bed protrusions
    r     = par.distChannels;  % Typical spacing between channels

    c1    = 1/(par.rho*par.Latent);
    c2    = 2*A*(par.n^(-par.n));
    c3    = (2/pi)^(0.25)*sqrt((pi+2)/(par.rhow*f));
    
    % Flow in a single channel
    Q = flw*r;
    
    % Effective pressure
    Neff = (c1*Q.*dphi0 + ub*hb)./(c2.*c3^(-1/alpha).*Q.^(1/alpha).*dphi0.^((1-beta)/alpha));
    Neff = Neff.^(1/par.n);

    % Limits
    Neff = max(min(Po,Neff),par.sigmat*Po);
end


