function default=InitDefault

% Kori-ULB
% Default values (different from zero) for control parameters in struct ctr

    default.plotGL=1; % plot grounding line position
    default.upstream=1; % upstream differences in solution of ice thickness equation
    default.ItSolv=1; %VL: new parameter
    default.BetaIter=5; % number of times that beta is updated within velocity nonlinear iteration
    default.snapshot=100;
    default.Ao=1.0e-16;
    default.m=1;
    default.kmax=11;
    default.Hinv=500;
    default.Tinv=200;
    default.stopoptim=0.1;
    default.HinvMelt=100;   %VL: new parameter
    default.TinvMelt=50;   %VL: new parameter
    default.radnorm=50;     %VL: new parameter
    default.shelftune=0.5;
    default.gammaTpico=3e-5;
    default.Cpico=1e6;
    default.gammaTquad=50e-5;
    default.gammaTlin=2e-5;
    default.gamma0_quad=1.447733676e4;
    default.gamma0_quad_slope=2060000;
    default.gammaTplume=5.9e-4;
    default.E0=3.6e-2;
    default.M0picop=10;
    default.meltfactor=0.3;
    default.G0=0.042;
    default.u0=1e12;
    
end


