%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Kori-ULB: The ULB ice flow model                      %
%                                                                       %
%                     Version 0.92 August 2025                          %
%                                                                       %
%                           Frank Pattyn                                %
%                    Laboratoire de Glaciologie                         %
%                  Universite libre de Bruxelles                        %
%                       Frank.Pattyn@ulb.be                             %
%                                                                       %
%                       co-developpers team                             %
%                        Javier Blasco Navarro                          %
%                         Kevin Bulthuis                                %
%                      Violaine Coulon                                  %
%                         Elise Kazmierczak                             %
%                        Daniel Moreno Parada                           %
%                        Thomas Gregov                                  %
%                        Olivia Raspoet                                 %
%                        Sainan Sun                                     %
%                          Lars Zipf                                    %
%                                                                       %
%                                                                       %
% Kori-ULB (The ULB Ice Flow Model) is a 2.5-dimensional finite         %
% difference numerical ice sheet model of intermediate complexity.      %
%                                                                       %
% MIT License                                                           %
%                                                                       %
% Copyright (c) 2017-2025 Frank Pattyn                                  %
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
% Other software packages used in Kori-ULB:                             %
%                                                                       %
%   - crameri: Fabio Crameri's scientific colormaps, version 8.0.       %
%              http://www.fabiocrameri.ch/colourmaps.php                %
%              C. Greene (UTIG, Texas) http://www.chadagreene.com       %
%   - convnfft, conv2fft:  Bruno Luong <brunoluong@yahoo.com>           %
%   - imagescn: C. Greene (UTIG, Texas) http://www.chadagreene.com      %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Today's scientists have substituted mathematics for experiments,
% and they wander off through equation after equation, and
% eventually build a structure which has no relation to reality.
%                                              Nikola Tesla (1934)
%
%-----------------------------------------------------------------
% Model Features
%-----------------------------------------------------------------
% -2.5D Finite difference ice sheet/ice shelf model
% -SSA-SIA hybrid velocity calculation (on Arakawa C grids)
% -SIA diffusive calculation (on Arakawa B-grid)
% -3D temperature field
% -Full thermomechanical coupling
% -Local and non-local isostatic adjustment (ELRA model) with spatially
%     varying flexural rigidity and asthenosphere viscosity
% -General slip law (viscous - power law - regularized Coulomb)
% -Grounding line parameterization with buttressing (optional)
% -Nudging method to determine spatially-varying basal slip coefficients
% -Nudging method to optimize sub-shelf mass balance for steady-state
% -PICO/PICOP/Plume ocean model for sub-shelf melt calculation
% -Calving, hydrofracturing and damage
% -Subglacial hydrology and till deformation
% -PDD model for surface melt
% -Colorblind-friendly output figures
%
%-----------------------------------------------------------------
% Model call
%-----------------------------------------------------------------
%
%   KoriModel(infile,outfile,ctr)
%
%               or
%
%   KoriModel(infile,outfile,ctr,fc)
%
%-----------------------------------------------------------------
% Model input
%-----------------------------------------------------------------
%
% Model input: input filename (infile); all other input is optional and
%   may contain ice geometry (H, B, MASK) or intial climate (Ts, Mb).
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   VERSION history                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% v0.92 (08/2025)
%
%------------------------------------------------------------------------


function varargout=KoriModel(infile,outfile,ctr,fc)

%-------------------
% model version
%-------------------

ctr.model='Kori-ULB';
ctr.version='v0.92';
fprintf('---%s %s---\n  [%s Frank.Pattyn@ulb.be]\n',ctr.model, ...
    ctr.version,char(169));

%------------------------------------------
% Run standard model test without arguments
%------------------------------------------
if nargin==0
    [ctr,outfile]=RunTest(ctr);
end

%---------------------------------------------------------------
% Determine whether and what forcing is used
%---------------------------------------------------------------
fc.forcingATM=0;
fc.forcingOCEAN=0;
if nargin==4
    if any(ismember(fields(fc),'atm_Ts_fname')) ||  ...
        any(ismember(fields(fc),'atm_Mb_fname')) || ...
        any(ismember(fields(fc),'atm_Pr_fname')) || ...
        any(ismember(fields(fc),'atm_Evp_fname')) || ...
        any(ismember(fields(fc),'atm_runoff_fname'))
        fc.forcingATM=1;
    end
    if any(ismember(fields(fc),'ocn_TF_fname')) || ...
        any(ismember(fields(fc),'ocn_To_fname')) || ...
        any(ismember(fields(fc),'ocn_So_fname'))
        fc.forcingOCEAN=1;
    end
end

%---------------------------------------------------------------
% Default values of control parameters (if different from zero)
%---------------------------------------------------------------

default=InitDefault;

%---------------------
% Input parameters
%---------------------

% Initialization of undefined control parameters (default values)
[ctr,fc]=InitCtr(ctr,fc,default);

% Impossible combination of parameters
if ctr.SSA==0 && ctr.shelf==1
    fprintf('If shelf=1, then SSA>0\n');
    return;
end

% Read model parameters
[par]=KoriInputParams(ctr.m,ctr.basin);

% Central differences and direct solver for SIA model
if ctr.SSA==0
    ctr.upstream=0;
    ctr.ItSolv=0;
end

% initialize the time slice counter
slicecount=0;

%--------------------------------------------------------------------
% Parameters that will be checked on whether they exist or not and
% initially set to 'false'
%--------------------------------------------------------------------

[Asor,stdB,v,vx,vy,tmp,Db,To,So,Tb,uxssa,uyssa,deltaZ,arcocn,arcocn0, ...
    E,wat,Epmp,Pr,Evp,runoff,MeltInv,lat,acc,Smelt,rain,TF,HAF,Hinit,ZB, ...
    flagHu,frb,kei,Ll,damage,CR,FMR,fluxmx,fluxmy]=deal(false);

%---------------------
% Initialization
%---------------------

[ctr.snapshot,plotst,cnt_atm,snp_atm,cnt_ocn,snp_ocn, ...
    Mb_update,Li,Lj,dtdx,dtdx2,X,Y,x,y,MASK,H,Ho,B,Bo,glMASK0, ...
    MASKo,Mb,Ts,As,G,u,ub,VAF,VA0,POV,SLC,Ag,Af,Btau,IVg,IVf,glflux, ...
    cfflux,dHdt,time,mbcomp,InvVol,ncor,dSLR,SLR,Wd,Wtil,Bmelt,NumStab, ...
    Dbw,CMB,FMB,flw,p,px,py,pxy,nodeu,nodev,nodes,node,VM,Tof,Sof, ...
    TFf,Tsf,Mbf,Prf,Evpf,runofff,Melt,shelftune,Melt_mean, ... 
    Bmelt_mean,Ts_mean,Mb_mean,To_mean,So_mean,TF_mean,CR_mean,FMR_mean, ...
    fluxmx_mean,fluxmy_mean,Smelt_mean,runoff_mean,rain_mean,acc_mean, ...
    Neff,expflw,kappa,etaD,beta2]=InitMatrices(ctr,par,default,fc);

%----------------------------------------------------------------------
% Read inputdata
% MASK=1: grounded ice sheet
% MASK=0: either ice shelf or open ocean
% MASK=3: ice shelf (Bedmap2/Bedmachine data Antarctica)
% Read optionally B, H, MASK, G, Mb, Ts, stdB, tmp, v, vx, vy, ...
%   Input initialized by zeros(B,H,Mb,Ts,stdB) and ones(MASK)
%----------------------------------------------------------------------

if nargin>0 && exist([infile,'.mat'],'file')
    load(infile);
end
[As,Mb0,Ts0,MASK0,MASK,bMASK,H,B,H0,B0]= ...
    InitInputData(ctr,par,As,Mb,Ts,MASK,MASKo,H,B);

%--------------------------------------------------
% Define basin-related variables (if basins exist)
%--------------------------------------------------

if islogical(ZB)==0
    mb_basin=zeros(ctr.nsteps,21,max(ZB(:)));
    SLC_basin=zeros(ctr.nsteps,max(ZB(:)));
    VAF_basin=zeros(ctr.nsteps,max(ZB(:)));
    IVg_basin=zeros(ctr.nsteps,max(ZB(:)));
    Ag_basin=zeros(ctr.nsteps,max(ZB(:)));
    if ctr.shelf==1
        IVf_basin=zeros(ctr.nsteps,max(ZB(:)));
        Af_basin=zeros(ctr.nsteps,max(ZB(:)));
    end
end

%------------------------------------------
% Test on existance of certain parameters
%------------------------------------------

[ctr,invmax2D,Asor,ncor,To,So,Pr0,Evp0,runoff0,Evp,Hinit]= ...
    ExistParams(ctr,par,ncor,Asor,stdB,v,uxssa,To,So,Db,B,MASK,As, ...
    Pr,Evp,runoff,Mb0,Hinit,Ho,damage);

%-------------------------------------
% Define grounded/floating ice sheet
%-------------------------------------

if (ctr.schoof>=1 || ctr.shelf==1) && ctr.inverse~=1
    [HAF,MASK,~,sn]=Floatation(par,B,SLR,H,MASK);
else % keep grounding line fixed
    [sn,~]=FixedFloatation(par,B,SLR,H,MASK);
end

%----------------------------------------------
% Copy initial datasets in separate variables
%----------------------------------------------

oldMASK=MASK; % used for interpolation between PICO calls
S0=sn; % initial surface elevation
HBo=Bo; % Original Bedmap/Bedmachine data to keep in each model run (Ho, Bo), if exist
HBo(MASK==0)=SLR(MASK==0)-par.rho/par.rhow*Ho(MASK==0);
HBo(HBo<Bo)=Bo(HBo<Bo);
sn0=HBo+Ho;
HBinit=Bo; 
HBinit(MASK==0)=SLR(MASK==0)-par.rho/par.rhow*Hinit(MASK==0);
HBinit(HBinit<Bo)=Bo(HBinit<Bo);
sninit=HBinit+Hinit; % Geometry from initialisation -- used for correction for elevation changes
Hn=H;
Bn=B;
To0=To;
So0=So;

%----------------------------------------
% Temperature parameter initialization
%----------------------------------------

cntT=0;
zeta=CalculateZeta(ctr.kmax,0.015);
if ctr.Tcalc>=1
    [tmp,Tb,dzc,dzp,dzm,E,Epmp]=InitTempParams(ctr,par,zeta,tmp,Ts,H,E);
    if ctr.Enthalpy==1
        CTSm=zeros(size(tmp));
        CTSp=zeros(size(tmp));
        Ht=zeros(size(Tb));
        Hw=zeros(size(Tb));
        Dfw=zeros(size(E));
        if ctr.Tinit==0
            Hw=max(0,min((Bmelt-par.Cdr)*ctr.dt*par.intT,par.Wmax));
            wat=max(0,(E-Epmp)/par.Latent);
            [CTSm,CTSp,Ht]=CalculateCTS(ctr,E,Epmp,MASK,H,zeta);
        end
    end
end

%--------------------------------------
% Initial Volume Above Floatation
% Initialization of geoid calculation
%--------------------------------------

SLR0=SLR;
VAF0=max(0,H0+min(B0-SLR,0)*(par.rhow/par.rho));
POV0=max(0,par.SLref-B0); % Potential ocean volume
if ctr.GeoidCalc==1
    [Pg0,geoide,frg]=InitGeoid(par,ctr,MASK,H0,B0,B,SLR);
end

%-------------------------------
% Initialization of bed loads
%-------------------------------

if ctr.BedAdj==1
    [B0,load0,frb,kei,Ll]=BedrockInit(ctr,par,H,B,SLR,MASK);
end
if ctr.BedAdj==2
    load0=zeros(ctr.imax,ctr.jmax);
    load0(MASK==1)=par.rho*H(MASK==1)/par.rhom;
    load0(MASK==0)=(-par.rhow)*(B(MASK==0)-SLR(MASK==0))/par.rhom;
    B0=B;
end

%--------------------------------------------------
% Initializing level set function for calving
%--------------------------------------------------

if exist('LSF','var')==0
    LSF=ones(ctr.imax,ctr.jmax);
    LSF(MASK==0 & H<=par.SeaIceThickness)=-1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                 %
%               ITERATION IN TIME                 %
%                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
cpu=cputime;
if ctr.runmode<2
    scrsz = get(groot,'ScreenSize');
    figure('Position',[30 scrsz(4)/2.5 scrsz(3)/1.4 scrsz(4)/2.2], ...
        'Name',outfile,'NumberTitle','off');
end
fprintf('... Running ...\n');

%-------------------------------------------------------------
% Restart procedure from previous 'toto' file when timeslice=1
%-------------------------------------------------------------

if ctr.restart==1 && (ctr.runmode==1 || ctr.runmode==3)
    outputname=[outfile,'_toto'];
    load(outputname); % read dumpfile
    ctr.restart=0; % reset restart
    cnt0=cnt+1; % advance one time step since dump
else
    cnt0=1;
end

for cnt=cnt0:ctr.nsteps

%-------------------------------------------------------------
% Print model progression on screen
%-------------------------------------------------------------

    cntT=cntT+1;
    if cnt>1
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    end
    fprintf('t = %12.2f', time(cnt));
    if cnt==ctr.nsteps
        fprintf('\n... End model run ...\n');
    end    

%------------------------------------------------------
% Copying Hn and Bn on H and B values
%------------------------------------------------------

    if cnt>1
        H=Hn; % assigning new ice thickness to previous ones
        if ctr.BedAdj>0
            B=Bn; % assigning new bed elevation to previous ones
        end
    end

%------------------------------------------------------
% Adjusting MASK based on floating condition
%------------------------------------------------------

    if ctr.GeoidCalc==0
        SLR=zeros(ctr.imax,ctr.jmax)+fc.DeltaSL(cnt);
    else
        SLR=dSLR+fc.DeltaSL(cnt);
    end
    if (ctr.schoof>=1 || ctr.shelf==1) && ctr.inverse~=1
        [HAF,MASK,HB,sn]=Floatation(par,B,SLR,H,MASK);
    else
        [sn,HB]=FixedFloatation(par,B,SLR,H,MASK);
    end
    
    
%------------------------------------------------------
% Surface mass balance and surface temperature
% parametrization and forcing
%------------------------------------------------------
    % Update atmospheric forcing based on external forcing data
    [Tsf,Mbf,Prf,Evpf,runofff,cnt_atm,snp_atm,Mb_update]= ...
        ATMupdate(fc,ctr,time,cnt,Ts0,Mb0,Pr0,Evp0,runoff0,cnt_atm, ...
        snp_atm,Mb_update,Tsf,Mbf,Prf,Evpf,runofff);

    % Correction/calculation of Ts 
    [Ts]=TsFunc(ctr,par,Tsf,sninit,sn,fc.DeltaSL(cnt),fc.DeltaT(cnt));
    % Correction/calculation of SMB components
    if Mb_update==0
        [Pr,Evp,runoff]=MbFunc(ctr,par,Prf,Evpf,runoff,runofff, ...
            Tsf,sn,sninit,X,Y,Lj,Li,fc.DeltaT(cnt),fc.DeltaSL(cnt));
    else
        [Mb,~,~]=MbFunc(ctr,par,Mbf,Evpf,runoff,runofff, ...
            Tsf,sn,sninit,X,Y,Lj,Li,fc.DeltaT(cnt),fc.DeltaSL(cnt));
    end
    if ctr.PDDcalc==1 % PDD MODEL
        if ctr.monthly==0
            [runoff,acc,rain,Smelt]=PDDyearly(ctr,par,Ts,Pr,lat,MASK,sn);
        else
            % Intra-annual (monthly or submonthly) variations are known
            % Evaluates mean yearly runoff at the beginning of the year
            if rem(time(cnt),1)==0
                [Ts_yc,Pr_yc]=ExtractAnnualCycle(fc,ctr,par,Pr, ...
                    Ts,sninit,sn,MASK,lat,fc.DeltaSL(cnt),fc.DeltaT(cnt),snp_atm);
                [runoff,acc,rain,Smelt]=PDDmonthly(Ts_yc,Pr_yc,par,ctr);
            end
            Pr=mean(Pr_yc,3);
        end
    end
    if Mb_update==0
        Mb=Pr-Evp-runoff; 
    end
    if ctr.PDDcalc==1 && ctr.PDD_anomaly==1
        aMb=Mb-fc.Mbref;
        Mb=Mb0+aMb;
    end

%------------------------------------------------------
% Geometry on staggered grids and grounding-line MASK
% Optionally define shelf mask and ice shelf numbering
%------------------------------------------------------

    [gradm,gradmx,gradmy,gradxy,gradsx,gradsy,gradHx,gradHy, ...
        Hm,Hmx,Hmy,Bmx,Bmy,signx,signy]=StaggeredGrid(sn,H,B,ctr);
    [bMASKm,bMASKx,bMASKy]=StaggeredBMASK(ctr,bMASK);

    % Searching GL on H-grid (define MASK & glMASK)
    [glMASK,MASK]=GroundingLineMask(MASK,bMASK,H,ctr,par);
    if ctr.schoof>=1 || ctr.shelf==1
        % sea ice thickness to control 'ocean' viscosity (0.1 m)
        H(MASK==0 & H<par.SeaIceThickness)=par.SeaIceThickness; 
    end
    if ctr.glMASKexist==1
        [ShelfN,numsh,shMASK,MASKlk]=ShelfSetup(MASK,glMASK,H,ctr,par);
    end

%--------------------------------------------------------------------------
% Thermomechanical coupling with method of Lliboutry (1979) and Ritz (1992)
%--------------------------------------------------------------------------

    [A,Ax,Ay,Ad,A3d,Tbc]=ThermoCoupling(ctr,par,Tb,H,bMASK,bMASKm, ...
        bMASKx,bMASKy,tmp,zeta,wat);

%--------------------------------------------------------------------------
% Basal sliding coefficients for different slip laws and basal hydrological
% models. All determined on h-grid
%--------------------------------------------------------------------------
    
    updateHydro=and(rem(cnt-1,ctr.FreqHydro)==0,ctr.FreqHydro<Inf);
    taud=par.rho*par.g*Hm.*sqrt(gradm); % taud on d-grid
    taudxy=par.rho*par.g*H.*sqrt(gradxy); % taud on h-grid
    [Asf,Asfx,Asfy,Asfd,Neff,Wtil,r,expflw]=BasalSliding(ctr,par, ...
        A,As,Tbc,H,B,MASK,ub,Wd,Wtil,Bmelt,flw,Neff,expflw,kappa,updateHydro, ...
        bMASK,bMASKm,bMASKx,bMASKy);

%-------------------------------------------------
% SIA velocity and diffusivities
% p=n for an isotherm ice sheet; p>=n otherwise
%-------------------------------------------------

    [d,udx,udy,ud,ubx,uby,ub,uxsia,uysia,p,pxy]= ...
        SIAvelocity(ctr,par,A,Ad,Ax,Ay,Asfd,Asfx,Asfy,taud,G,Tb, ...
        H,Hm,Hmx,Hmy,gradm,gradmx,gradmy,gradxy,signx,signy,MASK,p,px,py,pxy);
    
%---------------------------------------------
% Numerical solution of temperature field
%---------------------------------------------
    if ctr.Tcalc>=1 % on d-grid
        if cntT==1
            if (ctr.Tinit==1 && cnt==1) || ctr.Tinit==2
                [tmp,E,Epmp,wat]=InitTemp3d(G,taudxy,ub,ud,par,H,Mb,zeta,ctr,Ts, ...
                    MASK,fc.DeltaT(cnt));
            end
            if ctr.Tinit<2
                if ctr.SSA>=1 && (ctr.uSSAexist==1 || cnt>1)
                    if ctr.Enthalpy==1
                        [E,Epmp,wat,CTSm,CTSp,Bmelt,Dbw,Dfw,Hw,Ht,tmp]= ...
                            Enthalpy3d(par,ctr,E,Mb,Ts,G,A,H,pxy, ...
                            ctr.dt*par.intT,gradsx,gradsy,gradHx,gradHy, ...
                            gradxy,taudxy,udx,udy,vec2h(uxssa,uyssa), ...
                            uxssa,uyssa,zeta,dzc,dzm,dzp,fc.DeltaT,MASK, ...
                            Bmelt,Epmp,CTSm,CTSp,Hw,Ht,Dbw,Dfw,etaD,beta2,cnt);
                    else
                        [tmp,ctr]=Temperature3d(tmp,Mb,Ts,pxy,par, ...
                            ctr,ctr.dt*par.intT,gradsx,gradsy,gradHx,gradHy, ...
                            udx,udy,vec2h(uxssa,uyssa),uxssa,uyssa,zeta, ...
                            gradxy,H,dzc,dzm,dzp,G,taudxy, ...
                            A,fc.DeltaT,MASK,Bmelt,etaD,beta2,cnt);
                        Bmelt=BasalMelting(ctr,par,G,taudxy, ...
                            vec2h(uxssa,uyssa),H,tmp,dzm,MASK);
                    end
                else
                    if ctr.Enthalpy==1
                        [E,Epmp,wat,CTSm,CTSp,Bmelt,Dbw,Dfw,Hw,Ht,tmp]= ...
                            Enthalpy3d(par,ctr,E,Mb,Ts,G,A,H,pxy, ...
                            ctr.dt*par.intT,gradsx,gradsy,gradHx,gradHy, ...
                            gradxy,taudxy,udx,udy,ub,ubx,uby,zeta,dzc, ...
                            dzm,dzp,fc.DeltaT,MASK,Bmelt,Epmp,CTSm,CTSp, ...
                            Hw,Ht,Dbw,Dfw,etaD,beta2,cnt);
                    else
                        [tmp,ctr]=Temperature3d(tmp,Mb,Ts,pxy,par, ...
                            ctr,ctr.dt*par.intT,gradsx,gradsy,gradHx,gradHy, ...
                            udx,udy,ub,ubx,uby,zeta,gradxy,H,dzc,dzm,dzp, ...
                            G,taudxy,A,fc.DeltaT,MASK,Bmelt,etaD,beta2,cnt);
                        Bmelt=BasalMelting(ctr,par,G,taudxy,ub,H,tmp,dzm,MASK);
                    end
                end
            end
        end
        Tb=tmp(:,:,ctr.kmax)-par.T0;
    else
        if ctr.subwaterflow>0
            Bmelt=zeros(ctr.imax,ctr.jmax)+1e-8;
            Bmelt(MASK==0)=0;
        end
    end

%-------------------------------------
% Subglacial water flow
%-------------------------------------
    if not(or(ctr.subwaterflow==0, ctr.subwaterflow==2))
        if cntT==1
            [flw,Wd]=SubWaterFlux(ctr,par,H,HB,MASK,max(1e-8,Bmelt+Dbw));
            Wd0=Wd;
        else
            if ctr.inverse==0
                Wd=ExtrapolateWaterFlux(MASK,oldMASK,Wd,Wd0,ctr,par);
            end
        end
    end

%----------------------------------------------------------
% MASK for ice shelf/SSA (zero if SSA is not calculated)
%----------------------------------------------------------

    [MASKmx,MASKmy,HAFmx,HAFmy]=SSAmasks(ctr,par,SLR,Bmx,Bmy,Hmx, ...
        Hmy,bMASKx,bMASKy);

%---------------------------
% SSA ice shelf/ice stream
%---------------------------

    if ctr.SSA>=1
        if cnt==1
            su=zeros(ctr.imax*ctr.jmax*2,1);
            if ctr.uSSAexist==1
                su(nodeu)=uxssa;
                su(nodev)=uyssa;
            else
                uxssa=zeros(ctr.imax,ctr.jmax);
                uyssa=zeros(ctr.imax,ctr.jmax);
            end
        end
        [uxssa,uyssa,beta2,eta,etaD,dudx,dudy,dvdx,dvdy,su,ubx,uby,ux,uy, ...
            damage,NumStabVel]= ...
            SSAvelocity(ctr,par,su,Hmx,Hmy,gradmx,gradmy,signx,signy,zeta, ...
            uxssa,uyssa,etaD,A3d,H,HB,B,stdB,Asf,A,MASK,glMASK,HAF,HAFmx,HAFmy,cnt, ...
            nodeu,nodev,MASKmx,MASKmy,bMASK,uxsia,uysia,udx,udy,ubx,uby,node,nodes, ...
            Mb,Melt,dtdx,dtdx2,VM,damage,shelftune);
        if ctr.NumCheck==1
            NumStab(cnt,1:5)=NumStabVel;
        end
    else
        ux=uxsia;
        uy=uysia;
    end
    u=vec2h(ux,uy);
    ub=vec2h(ubx,uby);
    ub(MASK==0)=0;

%-----------------------------------------------
% Parameterization of grounding-line migration
%-----------------------------------------------

    if ctr.schoof>=1
       [uxsch,uysch]=GroundingLines(ctr,par,glMASK,HAF,SLR,Ax,Ay, ...
           fc.butfac(cnt),H,Hmx,Hmy,B,Bmx,Bmy,Asfx,Asfy,ux,uy, ...
           dudx,dvdy,dudy,dvdx,eta);
    else
        uxsch=ux;
        uysch=uy;
    end
    [uxsch,uysch,d]=DiffusiveCorrection(ctr,par,uxsch,uysch,udx,udy, ...
        d,Ad,p,Hm,gradm,bMASK);

%------------------------------------------------------------------
% Define distance to open ocean for calving and sub-shelf melting
%------------------------------------------------------------------

    if ctr.glMASKexist==1 && par.ArcOcean==1 && ctr.basin==0
        if cntT==1
            [arcocn,~,~]=OceanArc(MASK,H,MASKlk,ctr,par);
            arcocn0=arcocn;
        else
            arcocn=ExtrapolateArc(MASK,oldMASK,arcocn,arcocn0,ctr);
        end
    end

%------------------------------------------------------
% Sub-ice shelf melting
%------------------------------------------------------

	if ctr.glMASKexist==1 && ctr.inverse==0
		% Update ocean forcing based on external forcing data
		[Tof,Sof,TFf,cnt_ocn,snp_ocn]=OCEANupdate(fc,ctr,time,cnt, ...
            So0,To0,Tof,Sof,TFf,cnt_ocn,snp_ocn);

		% Extrapolate To and Tf to the depth of interest according 
        % to the chosen melt parameterization
		[To,So,TF]=OCEANdepthOfInt(fc,ctr,par,Tof,Sof,TFf,H, ...
            HB,B,MASK,glMASK,ShelfN,numsh);
	end

    if ctr.meltfunc>=1 && ctr.glMASKexist==1 && ctr.inverse==0
        Tf=par.lambda1*So+par.lambda2+par.lambda3*HB; % Ocean freezing point
        Tf(MASK==1)=0;
        [Melt,fc.butfac(cnt),H]=SubShelfMelt(ctr,fc,par,Tf,To,So,TF, ...
            fc.butfac(cnt),HB,glMASK,H,B,ShelfN,numsh,shMASK,MASK,MASKlk, ...
            uxssa,uyssa,arcocn,MeltInv,cnt);
    else
        Melt=zeros(ctr.imax,ctr.jmax);
        Melt(MASK==0)=ctr.meltfac;
    end
    if ctr.inverse==2 && cnt>1
        Melt=MeltInv;
    end

%---------------------------------------------------------------
% Calving and hydrofracturing (after Pollard et al., 2015)
% Melting at vertical face of calving front
% Use of LSF function (implemented by Vio - 2024)
%---------------------------------------------------------------

    if ctr.calving>=1 && ctr.shelf==1
        [he,fi]=DefineEdgeThickness(ctr,par,glMASK,H); % Pollard 2015
        [FMB,FMR]=VerticalFaceMelt(ctr,par,SLR,B,Melt,MASK,glMASK,he);
        [CMB,LSF,CR]=CalvingAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx, ...
            glMASK,H,A,uxssa,uyssa,arcocn,B,runoff,MASK,MASKo,Ho, ...
            bMASK,LSF,node,nodes,VM,cnt,ux,uy,Melt,he,fi,FMR);
    end

%---------------------------------------------------------
% Construct the finite-difference stiffness matrix 
% and solve continuity equation
%---------------------------------------------------------

    if ctr.diagnostic==0
        % sum of mass balance components for continuity equation
        Massb=Mb-Bmelt-Melt;
        if ctr.basin==1
            Massb(bMASK==1)=0; % only for ice thickness evolution
        end
        [Hn,flagH,relresH,iterH]=SparseSolverIceThickness(node,nodes, ...
            Massb,H,B,SLR,glMASK,dtdx,dtdx2, ...
            d,uxsch,uysch,ctr,cnt,bMASK,VM,par);
        if ctr.NumCheck==1
            NumStab(cnt,6:8)=[relresH,iterH,flagH];
        end
        if ctr.calving>=1 % remove icebergs
            Hn(LSF<0)=par.SeaIceThickness;
        end
        Hn(Hn<0)=0; % limit on minimal ice thickness
        dHdt(cnt)=mean(abs(Hn(:)-H(:)))/ctr.dt; % ice-sheet imbalance
    end
    % Ensure continuity of ice shelves during calving front advance/retreat
    % Grid point indices that have now become calving front and used
    % to be open sea. (DMP)
    if ctr.glMASKexist==1
        [row,col]=find((glMASK==5) & (glMASK0==6));
        % Ensure continuity if for such points.
        if ~isempty(row)
            [H,Hn]=IceShelfContinuity(ctr,row,col,H,Hn,glMASK);
        end
    end

%----------------------
% Bedrock adjustment
%----------------------

    if cntT==1 && ctr.BedAdj>0
        bload=BedrockAdjustment(ctr,par,load0,MASK,Hn,B,SLR,Db,VM, ...
            node,nodes,frb,kei,Ll);
        Bn=(B-B0+bload)*ctr.dt*par.intT./(-Btau)+B;
    end

%-----------------------------------------------------------------
% Inversion procedure to optimize basal sliding coefficients
% no optimization within the last 'stopoptim' part of nsteps
%-----------------------------------------------------------------

    if ctr.inverse>=1 && rem(cnt*ctr.dt,ctr.Tinv)==0 && cnt>1 && ...
            cnt<(1-ctr.stopoptim)*ctr.nsteps
        [As,deltaZ,Asor]=OptimizeIceSheet(ctr,par,cnt,Asor,MASK, ...
            bMASK,deltaZ,sn,sn0,r,ncor,B,stdB,vx,vy,ux,uy,invmax2D);
    end
    % optimization of basal melt rates based on Bernales (2017)
    if ctr.inverse==2 && rem(cnt*ctr.dt,ctr.TinvMelt)==0
        [MeltInv]=OptimizeIceShelf(ctr,MASK,glMASK,H,Ho,Melt,bMASK);
    end
    if ctr.inverse>=1
        % calculates misfit for ice grid cells (and only within drainage
        % basin if ctr.basin=1)
        InvVol(cnt,1)=sum(abs(sn(MASK==1 & bMASK==0)-sn0(MASK==1 & bMASK==0)));
        InvVol(cnt,2)=sum(sn(MASK==1 & bMASK==0)-sn0(MASK==1 & bMASK==0));
        if ctr.shelf==1
            InvVol(cnt,3)=mean(H(shMASK==1)-Ho(shMASK==1),'omitnan');
        end
    end

    
%------------------------------------------------------
% Volume above floatation, sea level contribution
% and geoid calculation
%------------------------------------------------------

    [VAF(cnt),POV(cnt),SLC(cnt),VA0(cnt)]=SeaLevel(ctr,par,SLR,B,H,H0,VAF0,POV0);
    if islogical(ZB)==0
        [VAF_basin(cnt,:),SLC_basin(cnt,:)]=SeaLevelBasin(ctr,par,SLR,B,H,H0,VAF0,POV0,ZB);
    end
    if ctr.GeoidCalc==1 && cntT==1
        [SLR,dSLR]=CalculateGeoid(par,ctr,fc.DeltaSL(cnt),SLC(cnt),SLR,SLR0, ...
            Bn,Hn,Pg0,frg,geoide);
    elseif ctr.GeoidCalc==2
        dSLR=zeros(ctr.imax,ctr.jmax)+SLC(cnt); % SLC added uniformly
        SLR=SLR0+dSLR+fc.DeltaSL(cnt);
    end
    
%------------------------------------------------------
% Time-dependent mass balance components
%------------------------------------------------------
    
    if islogical(ZB)==0
        mb_basin(cnt,:,:)=BasinFlux(ctr,par,acc,Smelt,runoff,rain,Mb,Pr, ...
            H,Hn,Bmelt,Melt,CMB,FMB,MASK,bMASK,squeeze(mb_basin(cnt,:,:)),B,Bn,SLR,ZB);
        mbcomp(cnt,:)=sum(squeeze(mb_basin(cnt,:,:)),2);
    else
        mbcomp(cnt,:)=MBcomponents(ctr,par,acc,Smelt,runoff,rain,Mb,Pr, ...
            H,Hn,Bmelt,Melt,CMB,FMB,MASK,bMASK,mbcomp(cnt,:),B,Bn,SLR);
    end
    IVg(cnt)=sum(H(MASK==1))*ctr.delta^2;
    Ag(cnt)=sum(MASK(H>0)==1)*ctr.delta^2;
    if ctr.glMASKexist==1
        Af(cnt)=sum(MASK(H>par.SeaIceThickness)==0)*ctr.delta^2;
        IVf(cnt)=sum(H(MASK==0 & H>par.SeaIceThickness))*ctr.delta^2;
        fluxmx=Hmx.*ux*ctr.delta;
        fluxmy=Hmy.*uy*ctr.delta;
        glflux(cnt)=sum(fluxmx(glMASK==2 & circshift(glMASK,[0 -1])>2))+ ...
            sum(fluxmy(glMASK==2 & circshift(glMASK,[-1 0])>2))+ ...
            sum(-fluxmx(glMASK>2 & circshift(glMASK,[0 -1])==2))+ ...
            sum(-fluxmy(glMASK>2 & circshift(glMASK,[-1 0])==2));
        cfflux(cnt)=sum(fluxmx(glMASK==5&circshift(glMASK,[0 -1])>5))+ ...
            sum(fluxmy(glMASK==5&circshift(glMASK,[-1 0])>5))+ ...
            sum(-fluxmx(glMASK>5&circshift(glMASK,[0 -1])==5))+ ...
            sum(-fluxmy(glMASK>5&circshift(glMASK,[-1 0])==5));
    end
    if islogical(ZB)==0
        for i=1:max(ZB(:))
            IVg_basin(cnt,i)=sum(H(MASK==1 & H>0 & ZB==i))*ctr.delta^2;
            Ag_basin(cnt,i)=sum(sum(MASK==1 & H>0 & ZB==i))*ctr.delta^2;
            if ctr.shelf==1
                Af_basin(cnt,i)=(sum(sum(MASK==0 & H>par.SeaIceThickness & ZB==i)))*ctr.delta^2;
                IVf_basin(cnt,i)=(sum(H(MASK==0 & H>par.SeaIceThickness & ZB==i)))*ctr.delta^2;
            end
        end
    end
    % Update glMASK to keep track of moving calving front (DMP)
    glMASK0=glMASK;

%------------------------------------
% NaN check
%------------------------------------   
    flagHu=NaNcheck(flagHu,Hn,ux,uy);
    if flagHu
        outputname=[outfile,'_toto'];
        save(outputname);
        fprintf('\n Run stopped at t = %12.2f\n\n', time(cnt));
        if nargout==1
            varargout{1}=flagHu;
        end
        return;
    end

%------------------------------------
% Save time-dependent matrices
%------------------------------------
    if ctr.YearlyMeans==1 && cnt==1
        [Melt_mean,Bmelt_mean,Ts_mean,Mb_mean,To_mean,So_mean,TF_mean, ...
            CR_mean,FMR_mean,fluxmx_mean,fluxmy_mean,Smelt_mean, ...
            runoff_mean,rain_mean,acc_mean]=InitYearlyMeans(Melt, ...
            Bmelt,Ts,Mb,To,So,TF,CR,FMR,fluxmx,fluxmy,Smelt,runoff,rain,acc);
    end
    if ctr.timeslice==1 
        if cnt==1 || (ctr.SnapList==1 && fc.snap_year(slicecount)==time(cnt)) ...
                || (ctr.SnapList==0 && rem(cnt-1,plotst)==0)
            slicecount=slicecount+1;
            fname=strcat(outfile,'_',num2str(slicecount-1,'%03i'));
            save(fname,par.varlist{1,1});
            for i=2:length(par.varlist)
                save(fname,par.varlist{1,i},'-append');
            end
            if (cnt>1 && ctr.SnapList==1 && ...
                    fc.snap_year(slicecount-1)==fc.snap_year(end))
                slicecount=length(fc.snap_year);
            end
        end
    end
    if ctr.YearlyMeans==1
        [Melt_mean,Bmelt_mean,Ts_mean,Mb_mean,To_mean,So_mean,TF_mean, ...
            CR_mean,FMR_mean,fluxmx_mean,fluxmy_mean,Smelt_mean, ...
            runoff_mean,rain_mean,acc_mean]=YearlyMeans(Melt, ...
            Bmelt,Ts,Mb,To,So,TF,CR,FMR,fluxmx,fluxmy,Smelt,runoff, ...
            rain,acc,cnt,ctr,Melt_mean,Bmelt_mean,Ts_mean,Mb_mean, ...
            To_mean,So_mean,TF_mean,CR_mean,FMR_mean,fluxmx_mean, ...
            fluxmy_mean,Smelt_mean,runoff_mean,rain_mean,acc_mean);
    end
    
%------------------------------------
% Main figure plot during run
%------------------------------------

    if ctr.runmode<2 && rem(cnt-1,plotst)==0
        PlotMainFigure(ctr,par,x,y,sn,S0,H,u,B,MASK,MASKo,glMASK,LSF);
    end
    
%------------------------------------
% Save intermediate output
%------------------------------------

    cntT(cntT>=par.intT)=0;
    oldMASK=MASK;
    if ctr.runmode==1 || ctr.runmode==3 || ctr.runmode==5
        if rem(cnt-1,plotst)==0 && cnt>1
            outputname=[outfile,'_toto'];
            save(outputname);
        end
        if ctr.runmode==5
            MeltDown; break;
        end
    end
        
%------------------------------------
% End of time loop
%------------------------------------
    
end

toc;
fprintf('CPU time: %.2f seconds\n', cputime-cpu);

%------------------------------------
% Save results
%------------------------------------

Mbend=Mb; % final values of Mb
Tsend=Ts; % final values of Ts;
Mb=Mb0; % reset output Mb to original value, since Mb=Mb-Melt-CMB
Ts=Ts0; % reset output Ts
To=To0; % reset output To
So=So0; % reset output So


save(outfile,'H','B','Ho','Bo','MASK','MASKo','As','G', ...
    'Ts','Mb');
if ctr.Tcalc>=1
    save(outfile,'tmp','Bmelt','-append');
end
if ctr.Enthalpy==1
    save(outfile,'E','-append'); %OR
end
if ctr.inverse>0
    save(outfile,'Asor','-append');
end
if ctr.stdBexist==1
    save(outfile,'stdB','-append');
end
if ctr.vexist==1
    save(outfile,'vx','vy','v','-append');
end
if ctr.Dbexist==1
    save(outfile,'Db','-append');
end
if exist('To','var')==1
    save(outfile,'To','So','-append');
end
if islogical(lat)==0
    save(outfile,'lat','lon','-append');
end
if ctr.SSA>=1
    save(outfile,'uxssa','uyssa','-append');
end
if ctr.damage>=1
    save(outfile,'damage','-append');
end
if islogical(Btau)==0
    save(outfile,'Btau','-append');
end
if exist('Bor','var')==1
    save(outfile,'Bor','-append');
end
if exist('ZB','var')==1
    save(outfile,'ZB','-append');
end
if islogical(Wd)==0
    save(outfile,'Wd','-append');
end
if islogical(Wtil)==0
    save(outfile,'Wtil','-append');
end
if islogical(MeltInv)==0
    save(outfile,'MeltInv','-append');
end
if islogical(Pr)==0
    Prend=Pr;
    Pr=Pr0;
    save(outfile,'Pr','-append');
    Evpend=Evp;
    Evp=Evp0;
    save(outfile,'Evp','-append');
end
if ctr.inverse>0
    Hinit=H;
    save(outfile,'Hinit','-append');
end

Ts=Tsend;
Mb=Mbend;

ctr=orderfields(ctr); % put in alphabetical order
outputname=[outfile,'_toto'];
save(outputname);
if nargout==1
    varargout{1}=flagHu;
end

%--------------------------
% Plots
%--------------------------
if ctr.runmode<2
    FinalPlots(ctr,outfile,time,InvVol,SLC,u,v,MASK,bMASK,glMASK);
end

end


