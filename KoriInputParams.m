%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Kori-ULB (The ULB Ice Flow Model) is a 2.5-dimensional finite         %
% difference numerical ice sheet model of intermediate complexity.      %
%                                                                       %
% MIT License                                                           %
%                                                                       %
% Copyright (c) 2017-2023 Frank Pattyn                                  %
%                                                                       %
% Permission is hereby granted, free of charge, to any person obtaining %
% a copy of this software and associated documentation files (the       %
% "Software"), to deal in the Software without restriction, including   %
% without limitation the rights to use, copy, modify, merge, publish,   %
% distribute, sublicense, and/or sell copies of the Software, and to    %
% permit persons to whom the Software is furnished to do so, subject to %
% the following conditions:                                             %
%                                                                       %
% The above copyright notice and this permission notice shall be        %
% included in all copies or substantial portions of the Software.       %
%                                                                       %
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       %
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    %
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.%
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  %
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  %
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     %
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [par]=KoriInputParams(m,basin)

%-----------------------------------
% General control parameters
%-----------------------------------

par.dcolor='broc'; % Crameri colorscale (requires CrameriColourMaps7.0.mat) - vik
par.color='imola'; % Crameri colorscale (requires CrameriColourMaps7.0.mat) - -roma

%-----------------------------------
% Numerical control parameters
%-----------------------------------

par.maxspeed=40e3; % maximum ice speed limit (m/a)
par.omega=2.5; % Crank-Nicolson scale factor (0=explicit; 1=implicit; >1 over-implicit)
par.secperyear=31556926;
% 2d variables to be saved when timeslice=1
par.varlist={'MASK','H','B','ux','uy','flux','Tbc','SLR','Neff','Melt'};


%-----------------------------------
% Subglacial characteristics
%-----------------------------------

par.PoreFrac=0.96; % Fraction of water pressure to balance ice pressure (0.96)
par.longcoupwater=5; % distance in number of ice thicknesses over which
                 % hydraulic gradient coupling takes place
par.dirpp_war=[9 8 7 6 5 4 3 2 1];
par.waterviscosity=1.8e-3/par.secperyear;
par.NeffScale=5e6; % scale factor for Effective Pressure
par.Wdmin=1e-8; % minimum value for Wd and Wtil
par.Wdmax=0.015; % maximum value for Wd
par.Wmax=2; % maximum value for Wtil (2 m)
par.flw0=1e5; % maximum value for subglacial water flux
par.Cdr=1e-3; % background till drainage rate
par.Cc=0.12; % till compressibility (Tulaczyk et al., 2000a)
par.e0=0.69; % reference void ratio at N0(Tulaczyk et al., 2000a)
par.N0=1e3; % reference effective pressure (Tulaczyk et al., 2000a)
par.sigmat=0.02; % Ntil lower bound, as fraction of overburden pressure


%-----------------------------------
% Ice dynamics
%-----------------------------------

par.ShelfPinning=1; % sub-shelf pinning of ice shelves based on bedrock variability
par.g=9.81; % gravitational acceleration
par.rho=917.; % ice density
par.rhow=1027.; % sea water density
par.rhom=3370.; % mantle density
par.n=3; % flow law exponent
par.visciter=50; % Maximum number of iterations on the nonlinear part of the SSA equation (50)
par.visctol=5e-1; % Tolerance for calculation of the nonlinear part of the SSA equation (0.5)
par.veliter=50;  % Maximum number of iterations for the iterative SSA velocity solver (50)
par.veltol=1e-4;  % Tolerance for the iterative SSA velocity solver (1e-4)
if basin==1
    par.veltol=par.veltol/10;
end
par.Hiter=20; % max iteration number for iterative thickness solver (20)
par.Htol=1e-6;  % tolerance for iterative thickness solver (1e-6)
par.Z=2*(par.g*par.rho)^par.n; % SIA isothermal pre-term
par.dlim=0.3; % Limit on local crevasses depth (% of H)
par.damlim=0.7; % limit on total damage (% of H)

%-----------------------------------
% Ice-ocean interactions
%-----------------------------------

par.Latent=3.35e5; % Latent heat of freezing
par.cp0=3974.; % Heat capacity of ocean water
par.Soi=34.5; % ocean salinity for initialization
par.Toi=-1.7; % Ocean temperature for initialization
par.SeaIceThickness=0.1;
par.ArcOcean=0; % include Ocean Arc to control Melt and calving
par.LatentMelt=(par.rhow*par.cp0)/(par.rho*par.Latent);

%-----------------------------------
% PICO and plume model parameters
%-----------------------------------

par.nbox=5; % max number of ocean boxes
par.alphao=7.5e-5;
par.betao=7.7e-4;
par.rhoref=1033;

par.lambda1=-5.73e-2;
par.lambda2=8.32e-2;
par.lambda3=7.61e-4;
% par.gamma0=1.447733676e4; % seems not used
par.gamma1=0.545;
par.gamma2=3.5e-5;
par.CdGamT=1.1e-3;
par.CdGamTS0=6e-4;
par.Eo=3.6e-2;
par.Cd=2.5e-3;
par.x0=0.56;

% Melt function coefficients (Lazeroms)
par.pcof=[0.1371330075095435
    5.527656234709359e1
    -8.951812433987858e2
    8.927093637594877e3
    -5.563863123811898e4
    2.218596970948727e5
    -5.820015295669482e5
    1.015475347943186e6
    -1.166290429178556e6
    8.466870335320488e5
    -3.520598035764990e5
    6.387953795485420e4];

% Plume V2 constants (Lazeroms 2019)
par.C_eps_lazero = 0.6; % Slope correction parameter
par.alpha_coeff_lazero = 3.87e-5; % degC-1 Thermal expansion coefficient
par.beta_coeff_lazero = 7.86e-4; % psu-1 Haline contraction coefficient


%-----------------------------------
% Isostasy
%-----------------------------------

par.FlexRigid=1e25; % Flexural rigidity
par.bedrelax=3000.; % relaxation time astenosphere
par.nuB=0.25; % Poisson ratio in flexural rigidity

%-----------------------------------
% Local sea level (fingerprints)
%-----------------------------------

par.Re=6.3781e6; % Radius Earth
par.Me=5.972e24; % Mass Earth
par.Aoc=3.618e14; % ocean surface
par.geoidist=5000e3; % size of convolution filter
par.rhof=1000; % fresh water density
par.SLref=0; % reference sea level


%-----------------------------------
% Model initialization
%-----------------------------------

par.stdDevRegul=3.5; % standard deviation of the Gaussian filter for the regularization in grid cells (3.5)
if basin==1
    par.stdDevRegul=par.stdDevRegul+1;
end
par.invmin=1.e-10; % values valid for m=2; scaled with 150kPa for other m
par.invmax=1.e-3;
par.invmaxncor=1.e-5;
par.AsFroz=1e-11;
par.AsScale=1e5^(2-m);


%-----------------------------------
% Thermodynamics
%-----------------------------------

par.T0=273.15; % absolute temperature
par.K=2.1; % thermal conductivity
par.kdif=1.1487e-6; % kdif=K/(rho*cp), cp=2009; value of EISMINT
par.pmp=8.66e-4; % Clausius Clapeyron (Payne, 2000)
par.atune=1; % tuning factor in Arrhenius (1 for n=3; 1e-5 for n=4)
par.R=8.314; % gas constant
par.udfrac=0.25;
par.intT=10; % number of iterations for which tmp is calculated
par.TrTemp=-10; % Basal temperature for which ice is frozen to bed
par.Q1=78.2e3; % Arrhenius parameters from Ritz (1992)
par.Q2=95.45e3;
par.a1=1.66e-16;
par.a2=2e-16;

%-----------------------------------
% PDD model parameters
%-----------------------------------

par.PDDth=0; % PDD threshold temperature (0°C)
par.Train=2;
par.Tsnow=0;
par.snowfac=3/par.rho; % PDD factor for snow melt
par.icefac=8/par.rho; % PDD factor for ice melt
par.d_ice=5; % Maximum depth of refreezing of percolating meltwater (m)
par.Tlapse=-0.008; % Lapse rate for temperature correction with height
par.Tsigma=4; % standard deviation of mean T for PDD calculation
par.Psigma=3.5; % standard deviation of mean T for rain factor calculation
par.PDDsteps=48;

%-----------------------------------
% Basin model parameters
%-----------------------------------

par.As0=1e-20; % sliding coefficient outside basin
par.A0=1e-20; % Ice fluidity outside basin

end