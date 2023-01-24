
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Kori-ULB: The ULB ice flow model                      %
%                                                                       %
%                      Version 0.9 January 2023                         %
%                                                                       %
%                           Frank Pattyn                                %
%                    Laboratoire de Glaciologie                         %
%                  Universite libre de Bruxelles                        %
%                       Frank.Pattyn@ulb.be                             %
%                                                                       %
%                       co-developpers team                             %
%                            Kevin Bulthuis                             %
%                         Violaine Coulon                               %
%                           Sainan Sun                                  %
%                             Lars Zipf                                 %
%                                                                       %
%                                                                       %
% Kori-ULB (The ULB Ice Flow Model) is a 2.5-dimensional finite         %
% difference numerical ice sheet model of intermediate complexity.      %
%                                                                       %
% MIT License                                                           %
%                                                                       %
% Copyright (c) 2023 Frank Pattyn                                       %
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
%   - crameri: Fabio Crameri's scientific colormaps, version 4.0.       %
%              http://www.fabiocrameri.ch/colourmaps.php                %
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
% Main Matlab files: KoriModel.m, KoriInputParams.m
%
% Model input: input filename (infile); all other input is optional and
%   may contain ice geometry (H, B, MASK) or intial climate (Ts, Mb).
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   VERSION history                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% v0.9 (12/2022)
%
%------------------------------------------------------------------------


function KoriModel(infile,outfile,ctr,fc)

%-------------------
% model version
%-------------------

ctr.model='Kori-ULB';
ctr.version='v0.9';
fprintf('---%s %s---\n  [%s Frank.Pattyn@ulb.be]\n',ctr.model, ...
    ctr.version,char(169));

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
    Pr,Evp,runoff,MeltInv,lat,acc,Smelt,rain,TF,HAF,Hinit]=deal(false);

%---------------------
% Initialization
%---------------------

[ctr.snapshot,plotst,cnt_atm,snp_atm,cnt_ocn,snp_ocn, ...
    Mb_update,Li,Lj,dtdx,dtdx2,X,Y,x,y,MASK,H,Ho,B,Bo, ...
    MASKo,Mb,Ts,As,G,u,VAF,VA0,POV,SLC,Ag,Af,Btau,IVg,IVf,glflux, ...
    dHdt,time,mbcomp,InvVol,ncor,dSLR,SLR,Wd,Wtil,Bmelt, ...
    CMB,FMB,flw,p,px,py,pxy,nodeu,nodev,nodes,node,VM,Tof,Sof, ...
    TFf,Tsf,Mbf,Prf,Evpf,runofff,Melt,damage]= ...
    InitMatrices(ctr,par,default,fc);

%----------------------------------------------------------------------
% Read inputdata
% MASK=1: grounded ice sheet
% MASK=0: either ice shelf or open ocean
% MASK=3: ice shelf (Bedmap2/Bedmachine data Antarctica)
% Read optionally B, H, MASK, G, Mb, Ts, stdB, tmp, v, vx, vy, ...
%   Input initialized by zeros(B,H,Mb,Ts,stdB) and ones(MASK)
%----------------------------------------------------------------------

if exist([infile,'.mat'],'file')
    load(infile);
end
[As,Mb0,Ts0,MASK0,MASK,bMASK,H,B,H0,B0]= ...
    InitInputData(ctr,par,As,Mb,Ts,MASK,MASKo,H,B);

%------------------------------------------
% Test on existance of certain parameters
%------------------------------------------

[ctr,invmax2D,Asor,ncor,To,So,Pr0,Evp0,runoff0,Evp,Hinit]= ...
    ExistParams(ctr,par,ncor,Asor,stdB,v,uxssa,To,So,Db,B,MASK,As, ...
    Pr,Evp,runoff,Mb0,Hinit,Ho);

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
if ctr.Tcalc>=1
    [tmp,Tb,zeta,dzc,dzp,dzm]=InitTempParams(ctr,par,tmp,Ts);
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
% Copying H and B on old values
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
    [Ts,Tsf]=TsFunc(ctr,par,Tsf,sninit,sn,fc.DeltaSL(cnt),fc.DeltaT(cnt));
    % Correction/calculation of SMB components
    if Mb_update==0
        [Pr,Evp,runoff]=MbFunc(ctr,Prf,Evpf,runoff,runofff, ...
            Tsf,Ts,sn,X,Y,Lj,Li,fc.DeltaT(cnt));
    else
        [Mb,~,~]=MbFunc(ctr,Mbf,Evpf,runoff,runofff, ...
            Tsf,Ts,sn,X,Y,Lj,Li,fc.DeltaT(cnt));
    end
    if ctr.PDDcalc==1 % PDD MODEL
        if ctr.monthly==0
            [runoff,acc,rain,Smelt]=PDDyearly(ctr,par,Ts,Pr,lat,MASK,sn);
        else
            % Intra-annual (monthly or submonthly) variations are known
            % Evaluates mean yearly runoff at the beginning of the year
            if rem(time(cnt),1)==0
                [Ts_yc,Pr_yc]=ExtractAnnualCycle(fc,ctr,par,Pr, ...
                    Ts,sninit,sn,MASK,lat,fc.DeltaT(cnt),snp_atm);
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
    [glMASK,MASK]=GroundingLineMask(MASK,bMASK,H,ctr);
    if ctr.schoof>=1 || ctr.shelf==1
        % sea ice thickness to control 'ocean' viscosity (0.1 m)
        H(MASK==0 & H<par.SeaIceThickness)=par.SeaIceThickness; 
    end
    if ctr.glMASKexist==1
        [ShelfN,numsh,shMASK,MASKlk]=ShelfSetup(MASK,glMASK,H,ctr,par);
    end

%--------------------------------------------------------------------------
% Basal sliding coefficients for different slip laws and basal hydrological
% models. All determined on h-grid
%--------------------------------------------------------------------------
    
    taud=par.rho*par.g*Hm.*sqrt(gradm); % taud on d-grid
    taudxy=par.rho*par.g*H.*sqrt(gradxy); % taud on h-grid
    [Asf,Asfx,Asfy,Asfd,Tbc,Neff,pwv,Wtil,r]=BasalSliding(ctr,par, ...
        As,Tb,H,B,MASK,Wd,Wtil,Bmelt,flw,bMASK,bMASKm,bMASKx,bMASKy);

%--------------------------------------------------------------------------
% Thermomechanical coupling with method of Lliboutry (1979) and Ritz (1992)
%--------------------------------------------------------------------------

    [A,Ax,Ay,Ad]=ThermoCoupling(ctr,par,Tb,Tbc,H,bMASK,bMASKm,bMASKx,bMASKy);

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
                tmp=InitTemp3d(G,taudxy,ub,ud,par,H,Mb,zeta,ctr,Ts, ...
                    MASK,fc.DeltaT(cnt));
            end
            if ctr.Tinit<2
                if ctr.shelf==1 && ctr.SSA>=1 && cnt>1
                    [tmp,ctr]=Temperature3d(tmp,Mb,Ts,pxy,par, ...
                        ctr,ctr.dt*par.intT,gradsx,gradsy,gradHx,gradHy, ...
                        udx,udy,vec2h(uxssa,uyssa),uxssa,uyssa,zeta, ...
                        gradxy,H,dzc,dzm,dzp,G,taudxy, ...
                        A,fc.DeltaT,MASK,Bmelt,cnt);
                    Bmelt=BasalMelting(ctr,par,G,taudxy, ...
                        vec2h(uxssa,uyssa),H,tmp,dzm,MASK);
                else
                    [tmp,ctr]=Temperature3d(tmp,Mb,Ts,pxy,par, ...
                        ctr,ctr.dt*par.intT,gradsx,gradsy,gradHx,gradHy, ...
                        udx,udy,ub,ubx,uby,zeta,gradxy,H,dzc,dzm,dzp, ...
                        G,taudxy,A,fc.DeltaT,MASK,Bmelt,cnt);
                    Bmelt=BasalMelting(ctr,par,G,taudxy,ub,H,tmp,dzm,MASK);
                end
            end
        end
        Tb=tmp(:,:,ctr.kmax)-par.T0;
    else
        Bmelt=zeros(ctr.imax,ctr.jmax)+5e-3; % mean value based on Pattyn (2010)
    end

%-------------------------------------
% Subglacial water flow
%-------------------------------------
    if ctr.subwaterflow==1 || ctr.subwaterflow==3
        if cntT==1
            [flw,Wd]=SubWaterFlux(ctr,par,H,HB,MASK,Bmelt);
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
        [uxssa,uyssa,beta2,eta,dudx,dudy,dvdx,dvdy,su,ubx,uby,ux,uy,damage]= ...
            SSAvelocity(ctr,par,su,Hmx,Hmy,gradmx,gradmy,signx,signy, ...
            uxssa,uyssa,H,HB,B,stdB,Asf,A,MASK,glMASK,HAF,HAFmx,HAFmy,cnt, ...
            nodeu,nodev,MASKmx,MASKmy,bMASK,uxsia,uysia,udx,udy,node,nodes, ...
            Mb,Melt,dtdx,dtdx2,VM,damage);
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

    if ctr.glMASKexist==1 && (par.ArcOcean==1 || ctr.calving==1 || ...
            ctr.calving==2) && ctr.basin==0
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
%---------------------------------------------------------------

    if ctr.calving>=1 && ctr.shelf==1
        [CMB,LSF,he]=CalvingAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,glMASK,H,A, ...
            uxssa,uyssa,arcocn,B,runoff,MASK,MASKo,Ho,bMASK,LSF,node,nodes,VM);
        if ctr.calving<5
            [FMB]=VerticalFaceMelt(ctr,par,SLR,B,Melt,MASK,glMASK,he);
        end
    end

%---------------------------------------------------------
% Construct the finite-difference stiffness matrix 
% and solve continuity equation
%---------------------------------------------------------

    if ctr.diagnostic==0
        % sum of mass balance components for continuity equation
        Massb=Mb-Bmelt-Melt-CMB-FMB;
        if ctr.basin==1
            Massb(bMASK==1)=0; % only for ice thickness evolution
        end
        Hn=SparseSolverIceThickness(node,nodes,Massb,H,B,SLR,glMASK,dtdx,dtdx2, ...
            d,uxsch,uysch,ctr,cnt,bMASK,VM,par);
        if ctr.calving==5
            % remove icebergs
            Hn(LSF<0)=par.SeaIceThickness;
        end
        Hn(Hn<0)=0; % limit on minimal ice thickness
        dHdt(cnt)=mean(abs(Hn(:)-H(:)))/ctr.dt; % ice-sheet imbalance
    end

%----------------------
% Bedrock adjustment
%----------------------

    if cntT==1
        if ctr.BedAdj==1
            bload=BedrockAdjustment(ctr,par,load0,MASK,Hn,B,SLR,Db,VM, ...
                node,nodes,frb,kei,Ll);
        end
        if ctr.BedAdj==2
            bload=zeros(ctr.imax,ctr.jmax);
            bload(MASK==1)=par.rho*Hn(MASK==1)/par.rhom-load0(MASK==1);
            bload(MASK==0)=(-par.rhow)*(B(MASK==0)-SLR(MASK==0))/ ...
                par.rhom-load0(MASK==0);
        end
        if ctr.BedAdj>0
            Bn=(B-B0+bload)*ctr.dt*par.intT./(-Btau)+B;
        end
    end

%-----------------------------------------------------------------
% Inversion procedure to optimize basal sliding coefficients
% no optimization within the last 'stopoptim' part of nsteps
%-----------------------------------------------------------------

    if ctr.inverse>=1 && rem(cnt*ctr.dt,ctr.Tinv)==0 && cnt>1 && ...
            cnt<(1-ctr.stopoptim)*ctr.nsteps
        [As,deltaZ,Asor]=OptimizeIceSheet(ctr,par,cnt,Asor,MASK, ...
            MASKo,bMASK,deltaZ,sn,sn0,r,ncor,B,stdB,vx,vy,ux,uy,invmax2D);
    end
    % optimization of basal melt rates based on Bernales (2017)
    if ctr.inverse==2 && rem(cnt*ctr.dt,ctr.TinvMelt)==0
        if ctr.GroundedMelt==1
            [MeltInv]=OptimizeIceShelf(ctr,MASKo,glMASK,H,Ho,Melt,bMASK);
        else
            [MeltInv]=OptimizeIceShelf(ctr,MASK,glMASK,H,Ho,Melt,bMASK);
        end
    end
    if ctr.inverse>=1
        InvVol(cnt,1)=sum(abs(sn(MASK==1)-sn0(MASK==1)));
        InvVol(cnt,2)=sum(sn(MASK==1)-sn0(MASK==1));
        if ctr.shelf==1
            InvVol(cnt,3)=mean(H(shMASK==1)-Ho(shMASK==1),'omitnan');
        end
    end
    
%------------------------------------------------------
% Volume above floatation, sea level contribution
% and geoid calculation
%------------------------------------------------------

    [VAF(cnt),POV(cnt),SLC(cnt),VA0(cnt)]=SeaLevel(ctr,par,SLR,B,H,H0,VAF0,POV0);
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
    
    mbcomp(cnt,:)=MBcomponents(ctr,par,acc,Smelt,runoff,rain,Mb,Pr, ...
        H,Hn,Bmelt,Melt,CMB,FMB,MASK,bMASK,mbcomp(cnt,:),B,Bn,SLR);
    flux=TotalFlux(H,ux,uy,ctr.delta);
    IVg(cnt)=sum(H(MASK==1))*ctr.delta^2;
    Ag(cnt)=sum(sum(MASK==1))*ctr.delta^2;
    if ctr.glMASKexist==1
        Af(cnt)=(sum(sum(MASK==0))-sum(sum(glMASK==6)))*ctr.delta^2;
        IVf(cnt)=(sum(H(MASK==0))-sum(H(glMASK==6)))*ctr.delta^2;
        glflux(cnt)=sum(flux(glMASK==2));
    end
    
%------------------------------------
% Save time-dependent matrices
%------------------------------------
    if ctr.timeslice==1 && rem(cnt-1,plotst)==0
        slicecount=slicecount+1;
        fname=strcat(outfile,'_',num2str(slicecount-1,'%03i'));
        save(fname,par.varlist{1,1});
        for i=2:length(par.varlist)
            save(fname,par.varlist{1,i},'-append');
        end
    end
    
%------------------------------------
% Main figure plot during run
%------------------------------------

    if ctr.runmode<2 && rem(cnt-1,plotst)==0
        PlotMainFigure(ctr,par,x,y,sn,S0,H,u,B,MASK,glMASK,LSF);
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
            if ctr.runmode==5
                MeltDown; break;
            end
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

save(outfile,'H','B','Ho','Bo','MASK','MASKo','As','G', ...
    'Ts','Mb');
if ctr.Tcalc>=1
    save(outfile,'tmp','Bmelt','-append');
end
if ctr.inverse>0    %VL: add Asor to init outputs
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

%--------------------------
% Plots
%--------------------------
if ctr.runmode<2
    FinalPlots(ctr,outfile,time,InvVol,SLC,u,v,MASK,bMASK,glMASK);
end

end


function [Tsf,Mbf,Prf,Evpf,runofff,cnt_atm,snp_atm,Mb_update]= ...
    ATMupdate(fc,ctr,time,cnt,Ts0,Mb0,Pr0,Evp0,runoff0,cnt_atm, ...
    snp_atm,Mb_update,Tsf,Mbf,Prf,Evpf,runofff)

% Kori-ULB
% Atmospheric forcing (precipitation, runoff, evaporation) based on
% sequential input files from (regional) climate models

    if fc.forcingATM==1
        if time(cnt)==fc.atm_Tinit % Initialise at first snapshot year
            cnt_atm=1;
            snp_atm=1;
        end
        if cnt_atm==1 % snapshot
            if time(cnt)<=fc.atm_Tend
                if any(ismember(fields(fc),'atm_Ts_fname'))
                    load([fc.atm_Ts_fname,num2str(snp_atm,'%03i')]);
                    Tsf=double(Ts); % make sure matrix is double
                else
                    Tsf=Ts0+fc.DeltaT(cnt); % simplified forcing otherwise
                end
                if any(ismember(fields(fc),'atm_Mb_fname'))
                    load([fc.atm_Mb_fname,num2str(snp_atm,'%03i')]);
                    Mbf=double(Mb); % make sure matrix is double
                    Mb_update=1;
                else % If Mb forcing does not exist, check for Mb component forcing
                    Mb_update=0;
                    Mbf=Mb0;
                    if any(ismember(fields(fc),'atm_Pr_fname'))
                        load([fc.atm_Pr_fname,num2str(snp_atm,'%03i')]);
                        Prf=double(Pr); % make sure matrix is double
                    else
                        Prf=Pr0;
                    end
                    if any(ismember(fields(fc),'atm_Evp_fname'))
                        load([fc.atm_Evp_fname,num2str(snp_atm,'%03i')]);
                        Evpf=double(Evp); % make sure matrix is double
                    else
                        Evpf=Evp0;
                    end
                    if any(ismember(fields(fc),'atm_runoff_fname')) && ctr.PDDcalc==0
                        load([fc.atm_runoff_fname,num2str(snp_atm,'%03i')]);
                        runofff=double(runoff); % make sure matrix is double
                    else
                        runofff=runoff0;
                    end
                end
                snp_atm=snp_atm+1;
            else % If forcing shorter than simulation time
                if snp_atm==fc.atm_snapshots+1
                    snp_atm=fc.atm_snapshots+1-fc.atm_nrep; % Re-initialise snapshot
                end
                if any(ismember(fields(fc),'atm_Ts_fname'))
                    load([fc.atm_Ts_fname,num2str(snp_atm,'%03i')]);
                    Tsf=double(Ts); % make sure matrix is double
                else
                    Tsf=Ts0+fc.DeltaT(cnt); % simplified forcing otherwise
                end
                if any(ismember(fields(fc),'atm_Mb_fname'))
                    load([fc.atm_Mb_fname,num2str(snp_atm,'%03i')]);
                    Mbf=double(Mb); % make sure matrix is double
                    Mb_update=1;
                else % If direct Mb forcing does not exist, check for Mb components forcing
                    Mb_update=0;
                    Mbf=Mb0;
                    if any(ismember(fields(fc),'atm_Pr_fname'))
                        load([fc.atm_Pr_fname,num2str(snp_atm,'%03i')]);
                        Prf=double(Pr); % make sure matrix is double
                    else
                        Prf=Pr0;
                    end
                    if any(ismember(fields(fc),'atm_Evp_fname'))
                        load([fc.atm_Evp_fname,num2str(snp_atm,'%03i')]);
                        Evpf=double(Evp); % make sure matrix is double
                    else
                        Evpf=Evp0;
                    end
                    if any(ismember(fields(fc),'atm_runoff_fname')) && ctr.PDDcalc==0
                        load([fc.atm_runoff_fname,num2str(snp_atm,'%03i')]);
                        runofff=double(runoff); % make sure matrix is double
                    else
                        runofff=runoff0;
                    end
                end
                snp_atm=snp_atm+1;
            end
        end
        cnt_atm=cnt_atm+1;
        cnt_atm(cnt_atm>fc.atm_cnt)=1;
    else
        Tsf=Ts0+fc.DeltaT(cnt); % simplified forcing otherwise
        Mbf=Mb0;
        Prf=Pr0;
        Evpf=Evp0;
        runofff=runoff0;
    end

end


function [Zgl]=AdvecGL(HB,B,MASK,MASKpicop,MASKb,vx,vy,ctr)

% Kori-ULB
% PICOP model Pelle et al. (2018)
% Advecting the grounding line thickness Zgl 
% Employs an upstream difference scheme

    epsilon=1e-5; % diffusion coefficient for stabilizing advection scheme

    VM=reshape(MASKpicop,ctr.imax*ctr.jmax,1);
    nodes=sum(sum(MASKpicop));
    node=zeros(ctr.imax,ctr.jmax);
    node(MASKpicop>0)=linspace(1,nodes,nodes);

    Zgl=HB;
    Zgl=min(Zgl,0);
    Zgl(MASK==0)=0;

    dtdx2=epsilon/(ctr.delta*ctr.delta);
    u1=circshift(vx,[0 1]); % vx(i,j-1)
    v1=circshift(vy,[1 0]); % vy(i-1,j)

    % Velocity sign masks
    MV=ones(ctr.imax,ctr.jmax); % u>=1, v>=1
    MV(vx>=0 & vy<0)=2;
    MV(vx<0 & vy>=0)=3;
    MV(vx<0 & vy<0)=4;
    % conditions for (u,v)>=0 (applied to whole matrix)
    V0=4*dtdx2+u1/ctr.delta+v1/ctr.delta; % i,j
    V1=zeros(ctr.imax,ctr.jmax)-dtdx2; % i,j+1
    V2=-u1/ctr.delta-dtdx2; % i,j-1
    V3=zeros(ctr.imax,ctr.jmax)-dtdx2; % i+1,j
    V4=-v1/ctr.delta-dtdx2; % i-1,j
    % conditions for u>=0 and v<0
    V0(MV==2)=4*dtdx2+u1(MV==2)/ctr.delta-vy(MV==2)/ctr.delta; % i,j
    V3(MV==2)=vy(MV==2)/ctr.delta-dtdx2; % i+1,j
    V4(MV==2)=-dtdx2; % i-1,j
    % conditions for u<0 and v>=0
    V0(MV==3)=4*dtdx2-vx(MV==3)/ctr.delta+v1(MV==3)/ctr.delta; % i,j
    V1(MV==3)=vx(MV==3)/ctr.delta-dtdx2; % i,j+1
    V2(MV==3)=-dtdx2; % i,j-1
    % conditions for u<0 and v<0
    V0(MV==4)=4*dtdx2-vx(MV==4)/ctr.delta-vy(MV==4)/ctr.delta;
    V1(MV==4)=vx(MV==4)/ctr.delta-dtdx2; % i,j+1
    V2(MV==4)=-dtdx2; % i,j-1
    V3(MV==4)=vy(MV==4)/ctr.delta-dtdx2; % i+1,j
    V4(MV==4)=-dtdx2; % i-1,j

    V0(MASKpicop==0)=0;
    V1(MASKpicop==0)=0;
    V2(MASKpicop==0)=0;
    V3(MASKpicop==0)=0;
    V4(MASKpicop==0)=0;

    R0=zeros(ctr.imax,ctr.jmax);

    V0(MASKb==1)=1;
    V1(MASKb==1)=0;
    V2(MASKb==1)=0;
    V3(MASKb==1)=0;
    V4(MASKb==1)=0;
    R0(MASKb==1)=Zgl(MASKb==1);

    V=[reshape(V0(VM==1),nodes,1)
        V1(V1~=0)
        V2(V2~=0)
        V3(V3~=0)
        V4(V4~=0)];

    row=[reshape(node(VM==1),nodes,1)
        node(V1~=0)
        node(V2~=0)
        node(V3~=0)
        node(V4~=0)];

    nodeV1=circshift(node,[0 -1]); % i,j+1
    nodeV2=circshift(node,[0 1]); % i,j-1
    nodeV3=circshift(node,[-1 0]); % i+1,j
    nodeV4=circshift(node,[1 0]); % i-1,j

    col=[reshape(node(VM==1),nodes,1)
    nodeV1(V1~=0)
    nodeV2(V2~=0)
    nodeV3(V3~=0)
    nodeV4(V4~=0)];

    R=reshape(R0(VM==1),nodes,1);

    % construct sparse matrix
    A=sparse(row,col,V);
    % Cholesky factor and solve
    s=A\R;
    Zgl(node>0)=s(node(node>0));
    Zgl(Zgl>HB)=HB(Zgl>HB);
    Zgl(Zgl<B)=B(Zgl<B);
    Zgl=min(Zgl,0);
end


function Bmelt=BasalMelting(ctr,par,G,taudxy,ub,H,tmp,dzm,MASK)

% Kori-ULB
% Basal melting underneath the grounded ice sheet based on geothermal heat
% flux and basal sliding
    
    Bgrad=(G+taudxy.*ub/par.secperyear).*H/par.K;
    Tgrad=(tmp(:,:,ctr.kmax)-tmp(:,:,ctr.kmax-1))./dzm(:,:,ctr.kmax);
    Bmelt=min(1,max(1e-8,(Bgrad-Tgrad)*par.K*par.secperyear./ ...
        (par.rho*par.Latent*H)));
    Bmelt(MASK==0)=0;
    
end


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


function stdB2=BedVar(B,ctr)

% Kori-ULB
% Calculate large scale bed variability for ncor

    sigma=2;
    AllB=zeros(ctr.imax,ctr.jmax,(2*sigma+1)^2);
    k=1;
    for j=-sigma:sigma
        for i=-sigma:sigma
            zdist = sqrt(j^2 + i^2); %Distance in grid cells
            if zdist<=sigma
                AllB(:,:,k)=circshift(B,[-i -j]);
            end
            k=k+1;
        end
    end
    stdB2=std(AllB,1,3);

end


function bload=BedrockAdjustment(ctr,par,load0,MASK,Hn,B,SLR,Db,VM,node, ...
    nodes,frb,kei,Ll)

% Kori-ULB
% Bedrock loads for isostatic adjustment for constant and variable flexural
% rigidity

    if ctr.Dbexist==1
        % variable flexural rigidity (Kevin)
        loadB=zeros(ctr.imax,ctr.jmax);
        loadB(MASK==1)=par.rho*par.g*Hn(MASK==1)-load0(MASK==1);
        loadB(MASK==0)=(-par.rhow)*par.g*(B(MASK==0)-SLR(MASK==0))-load0(MASK==0);
        % Takes into account geoid changes: !!! If GeoidCalc=1 & DeltaSL~=0, 
        % need to get rid of DeltaSL signal --> TO IMPLEMENT 
        bload=SparseSolverBedrock(Db,par,ctr,VM,node,nodes,loadB);                
    else
        % constant flexural rigidity
        P=zeros(ctr.imax,ctr.jmax);
        P(MASK==1)=par.rho*par.g*Hn(MASK==1)*ctr.delta^2.;
        P(MASK==0)=(-par.rhow)*par.g*(B(MASK==0)-SLR(MASK==0))*ctr.delta^2.;
        P0=NaN(ctr.imax+2*frb,ctr.jmax+2*frb);
        P0(frb+1:frb+ctr.imax,frb+1:frb+ctr.jmax)=P;
        P0(isnan(P0))=P(ctr.imax,ctr.jmax);
        bload=-xcorr2(P0,kei)*Ll^2./(2*pi*par.FlexRigid);
        bload=bload(2*frb+1:ctr.imax+2*frb,2*frb+1:ctr.jmax+2*frb)-load0;
    end
end


function [B0,load0,frb,kei,Ll]=BedrockInit(ctr,par,H,B,SLR,MASK)

% Kori-ULB
% Initialization of bedrock loads (considering initial bedrock in isostatic
% equilibrium) for constant and variable flexural rigidity

    if ctr.Dbexist==0
        % Constant flexural rigidity
        Ll=(par.FlexRigid/par.rhom/par.g)^0.25;
        frb=round(6*Ll/ctr.delta);
        P=zeros(ctr.imax,ctr.jmax);
        P(MASK==1)=par.rho*par.g*H(MASK==1)*ctr.delta^2.;
        P(MASK==0)=(-par.rhow)*par.g*(B(MASK==0)-SLR(MASK==0))*ctr.delta^2.;
        P0=NaN(ctr.imax+2*frb,ctr.jmax+2*frb);
        P0(frb+1:frb+ctr.imax,frb+1:frb+ctr.jmax)=P;
        P0(isnan(P0))=P(ctr.imax,ctr.jmax);
        [kerx,kery]=meshgrid(-frb*ctr.delta/Ll:ctr.delta/Ll:frb*ctr.delta/Ll, ...
            -frb*ctr.delta/Ll:ctr.delta/Ll:frb*ctr.delta/Ll);
        ker=sqrt((kerx.^2+kery.^2));
        ker(ker<1e-8)=1e-8;
        kei=imag(besselk(0,ker*(1+1i)/sqrt(2)));
        load0=-xcorr2(P0,kei)*Ll^2./(2*pi*par.FlexRigid);
        load0=load0(2*frb+1:ctr.imax+2*frb,2*frb+1:ctr.jmax+2*frb);
        B0=B;
    else
        frb=false;
        kei=false;
        Ll=false;
        % Varying flexural rigidity (Kevin)
        load0=zeros(ctr.imax,ctr.jmax);
        load0(MASK==1)=par.rho*par.g*H(MASK==1);
        load0(MASK==0)=(-par.rhow*par.g)*(B(MASK==0)-SLR(MASK==0));
        B0=B;
    end
end


function [T,S]=BoxModel(S0o,T0o,numsh,par,Ak,Bk,ShelfN,HB,nu,ctr)

% Kori-ULB
% Determination of temperature and salinity in different boxes of the PICO
% sub-shelf ocean circulation model (Reese et al.)

    % Initialization of boxes and ice shelves
    % numsh = num ice shelves; nbox = num boxes
    HBk=zeros(numsh,par.nbox); % mean ice thickness
    Ab=zeros(numsh,1); % surface of each box (all boxes within a shelf have equal size)
    T=zeros(numsh,par.nbox); % mean box temperature
    S=zeros(numsh,par.nbox); % mean box salinity
    for i=1:numsh
        Ab(i)=mean(Ak(ShelfN==i));
        for k=1:par.nbox
            HBk(i,k)=mean(HB(ShelfN==i & Bk==k));
        end
    end
    HBk(isnan(HBk))=0; % some shelves don't have all boxes - set HBk=0

    % Box I
    Tref=min(par.lambda1*S0o+par.lambda2+par.lambda3*HBk(:,1)-T0o,0); % make sure Tref<=0
    g1=Ab*ctr.gammaT;
    dif=g1./(ctr.C*par.rhoref*(par.betao*S0o/nu-par.alphao));
    x=(-0.5)*dif+0.5*sqrt(dif.^2-4*dif.*Tref);
    T(:,1)=T0o-x; % T1 should be lower than T0 in case of melting at GL
    S(:,1)=S0o.*(1-x/nu);

    % Box II to n
    q=ctr.C*par.rhoref*(par.betao*(S0o-S(:,1))-par.alphao*(T0o-T(:,1)));
    for i=2:par.nbox
        Tref=par.lambda1*S(:,i-1)+par.lambda2+par.lambda3*HBk(:,i)-T(:,i-1);
        x=-g1.*Tref./(q+g1.*(1-par.lambda1*S(:,i-1)/nu));
        T(:,i)=T(:,i-1)-x;
        S(:,i)=S(:,i-1).*(1-x/nu);
    end
end


function [arcocn,distocn,distgl]=CalcOceanArc(io,jo,MASK,imax,jmax, ...
    delta,narc,tanarc,ifi,idir)

% Kori-ULB
% Calculation of the angle of ice shelves to open ocean (Pollard and
% DeConto)
    
    arcocn=0;
    distocn=1e6;
    distgl=1e7;
    for mm=1:narc
        idirmm=idir(mm);
        tanarcmm=tanarc(mm);
        calcdist=0;
        if ifi(mm)==1
            ii=io;
            if idirmm==-1
                ttmax=io-1;
            else
                ttmax=imax-io-1;
            end
            for tt=1:ttmax
                ii=ii+idirmm;
                jj=round((ii-io)*tanarcmm)+jo;
                if jj<1 || jj>jmax
                    arcocn=arcocn+360./narc;
                    break
                end
                if MASK(ii,jj)==1
                    zdist=delta*sqrt((ii-io)^2+(jj-jo)^2);
                    distgl=min(distgl,zdist);
                    break
                elseif MASK(ii,jj)==0 && calcdist==0
                    zdist=delta*sqrt((ii-io)^2+(jj-jo)^2);
                    distocn=min(distocn,zdist);
                    calcdist=1;
                end
            end
        else
            jj=jo;
            if idirmm==-1
                ttmax=jo-1;
            else
                ttmax=jmax-jo-1;
            end
            for tt=1:ttmax
                jj=jj+idirmm;
                ii=round((jj-jo)/tanarcmm)+io;
                if ii<1 || ii>imax
                    arcocn=arcocn+360./narc;
                    break
                end
                if MASK(ii,jj)==1
                    zdist=delta*sqrt((ii-io)^2+(jj-jo)^2);
                    distgl=min(distgl,zdist);
                    break
                elseif MASK(ii,jj)==0 && calcdist==0
                    zdist=delta*sqrt((ii-io)^2+(jj-jo)^2);
                    distocn=min(distocn,zdist);
                    calcdist=1;
                end
            end
        end
        if tt==ttmax
            arcocn=arcocn+360./narc;
        end
    end
end


function [SLR,dSLR]=CalculateGeoid(par,ctr,DeltaSL,SLC,SLR,SLR0,Bn,Hn,Pg0,frg,geoide)

% Kori-ULB
% Calculation of geoid change due to mass changes of the ice sheet. The
% effect is stored in SLR (spatial variability of global sea level)

    [~,MASKn,~,~]=Floatation(par,Bn,SLR,Hn,zeros(ctr.imax,ctr.jmax));
    Pg1=zeros(ctr.imax,ctr.jmax);
    Pg1(MASKn==1)=par.rho*Hn(MASKn==1)*ctr.delta^2.;
    Pg1(MASKn==0)=(-par.rhow)*(Bn(MASKn==0)-SLR(MASKn==0)+DeltaSL)* ...
        ctr.delta^2.; % Takes into account previous geoid change but not
                      % DeltaSL (external SL forcing) as ocean mass
    Pg1=Pg1+par.rhom*Bn*ctr.delta^2.; % Takes into account land deformation
    DPg=Pg1-Pg0;
    % Fingerprint
    MaP=zeros(ctr.imax+2*frg,ctr.jmax+2*frg);
    MaP(frg+1:frg+ctr.imax,frg+1:frg+ctr.jmax)=MASKn;
    loadH=conv2fft(DPg,geoide); % FFT convolution
    % Determine local sea-level change
    Vd=sum(loadH(MaP==0))*ctr.delta^2; % Volume change in domain without mass addition
    loadG=loadH(frg+1:ctr.imax+frg,frg+1:ctr.jmax+frg);
    dSLR=loadG+SLC-Vd/par.Aoc; % previously VAF instead of SLC
    SLR=SLR0+dSLR+DeltaSL; % SLR = Initial SLR + Geoid perturbation + External SL forcing
        
end


function zeta=CalculateZeta(kmax,min_layer)

% Kori-ULB
% Calculation of the vertical scaled coordinate system for temperature
% calculation. The vertical layers are unevenly distributed with smaller
% layers near the bottom of the ice sheet.

    if (kmax-1)*min_layer > 1. % if number of layers is too high for min_layer
        zeta=0:(1/(kmax-1)):1;
    else
        xz=[1 kmax-1 kmax];
        yz=[0 1.-min_layer 1];
        p=polyfit(xz,yz,2);
        nxz=1:1:kmax;
        zeta=p(1)*nxz.^2+p(2)*nxz+p(3);
    end
    zeta(1)=0;
    zeta(kmax)=1;
    
end


function [CMB,LSF,he]=CalvingAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,glMASK,H,A, ...
    uxssa,uyssa,arcocn,B,runoff,MASK,MASKo,Ho,bMASK,LSF,node,nodes,VM)

% Kori-ULB
% Calving functions (UPDATE NEEDED)

    div=dudx+dvdy; % flow divergence
    [he,fi]=DefineEdgeThickness(ctr,par,glMASK,H); % Pollard 2015   %VL: add par
    Hshelf=H;
    Hshelf(fi<1)=he(fi<1);
    
    if ctr.calving==1 || ctr.calving==2
        dive=A.*(0.25*par.rho*par.g*Hshelf*(1.-par.rho/par.rhow)).^par.n;
        div(fi<1)=dive(fi<1); % divergence for shelf edge grid points
        div(div>dive)=dive(div>dive); % divergence not larger than maximum value
        d_s=2./(par.rho*par.g)*(max(0,div)./A).^(1/par.n); % dry surface crevasse depth
        d_b=2.*par.rho/((par.rhow-par.rho)*(par.rho*par.g))* ...
            (max(0,div)./A).^(1/par.n); % basal crevasse depth
        uxs1=circshift(uxssa,[0 1]); % ux(i+1,j)
        uys1=circshift(uyssa,[1 0]); % uy(i,j+1)
        usH=sqrt((0.5*(uxssa+uxs1)).^2+(0.5*(uyssa+uys1)).^2);
        d_a=Hshelf.*min(1,max(0,log(usH/1600))/log(1.2)); % additional crevasse deepening
        hc=150*max(0,min(1,(arcocn-70)./20)); % PD16: hc=200*max(0,min(1,(arcocn-40)./20));
        d_t=Hshelf.*max(0,min(1,(hc-Hshelf)/50)); % thin floating ice
        if ctr.HydroFrac==1
            RPDD=runoff; % annual surface melt + rainfall available after refreezing
            RPDD(runoff<1.5)=0; % R calculated as in DP16
            RPDD(runoff>3)=runoff(runoff>3).^2;
            RPDD(runoff<3 & runoff>1.5)=4*1.5*(runoff(runoff<3 & runoff>1.5)-1.5);
            d_w=100*RPDD;
        else
            d_w=zeros(ctr.imax,ctr.jmax);
        end
        % Calving
        ratio=(d_s+d_b+d_a+d_t+d_w)./he; 
        ratio_c=0.75; % Pollard(2015)
        CMB=3000*max(0,min(1,(ratio-ratio_c)./(1-ratio_c))).*Hshelf/ctr.delta;
        CMB(glMASK<4)=0;
    end
    if ctr.calving==2 % combination of 1 and 4 (making sure ice shelves don't get bigger)
        CMB(MASK==0 & MASKo==0 & Ho<par.SeaIceThickness)=50; % higher value than initial (LZ)
    end
    if ctr.calving==3
        wc=min(1,he/200.); % Pollard (2012)
        CMB=(1.-wc)*30+wc*3e5.*max(div,0).*he/ctr.delta; % Pollard (2012)
        if par.ArcOcean==1
            CMB=CMB.*max(0,min(1,(arcocn-70)./20)); % Vio: arcocean applied as in PD12
        end
        CMB(glMASK<5)=0;
    end
    if ctr.calving==4 % calving front kept at observed position (if ice is floating)
        CMB=zeros(ctr.imax,ctr.jmax);
        CMB(MASK==0 & MASKo==0 & Ho<par.SeaIceThickness)=50; % higher value than initial (LZ)
    end
    if ctr.calving>2 % Find holes in shelves + contour (where CMB applies)
        [MASKHole,HoleCT]=ShelfHole(MASK,H,par);
        CMB(MASKHole==1 | HoleCT==1)=0; % No calving in enclosed shelves holes
    end
    CMB(glMASK>=5 & B<-2700)=50; % Calving rule to prevent unrealistic areas of
                                 % thin ice extending seaward beyond continental
                                 % shelf (Vio)
    %VL: ice thinner than 5m can't persist --> equivalent to glMASK
    CMB(glMASK>=5 & H<5)=50;

    if ctr.calving==5 % LSF function calving (not operational)
        CMB=0;
        % determine eigenvalues of strain tensor
        StrTen=zeros(2,2);
        StrEig1=zeros(ctr.imax,ctr.jmax);
        StrEig2=zeros(ctr.imax,ctr.jmax);
        for i=1:ctr.imax
            for j=1:ctr.jmax
                StrTen(1,1)=dudx(i,j);
                StrTen(2,2)=dvdy(i,j);
                StrTen(1,2)=0.5*(dudy(i,j)+dvdx(i,j));
                StrTen(2,1)=StrTen(1,2);
                StrVal=eig(StrTen);
                StrEig1(i,j)=StrVal(1);
                StrEig2(i,j)=StrVal(2);
            end
        end
        EffStr=0.5*(max(0,StrEig1).^2+max(0,StrEig2).^2);
        sigmalim=0.8e5;
        q=sqrt(3)*A.^(-1/par.n).*EffStr.^(0.5/par.n)/sigmalim;
%             EigCalv=1e5*StrEig1.*StrEig2;
        FrontMelt=1000;
        wx=uxssa.*(1-q)-sign(uxssa).*FrontMelt; % w = u - c
        wy=uyssa.*(1-q)-sign(uyssa).*FrontMelt;
        if ctr.basin==1
            wx(bMASK==1)=0;
            wy(bMASK==1)=0;
        end
        wx=zeros(ctr.imax,ctr.jmax);
        wy=zeros(ctr.imax,ctr.jmax);
        LSF=LSFfunction(LSF,ctr,wx,wy,node,nodes,VM,MASK);
    end
end


function [M_hat]=ComputeMhat(X_hat)
 
% Kori-ULB
% Compute M_hat for plume parameterisation.
% This function computes the M_hat (or M0) for the plume parameterisation.
% This is Equation 26 in Lazeroms et al. 2019.
%
% Parameters
%     X_hat : scalar or array
%         Coordinate describing distance from plume origin.
% Returns
%     M_hat : scalar or array
%     Dimensionless melt rate, to be multiplied with the Mterm in Eq 28a.
% here M_hat = M0(X_tilde), Eq 26 in Lazeroms et al. 2019

    M_hat=1./(2*sqrt(2)).*(3*(1-X_hat).^(4/3)-1).*sqrt(1-(1-X_hat).^(4/3));
    
end


function [Mterm] = ComputeMterm(T_in,S_in,Tf,c_rho_1,c_tau,gamma,E0, ...
    thermal_forcing_factor,sina,par)

% Kori-ULB
% Compute M-term for plume parameterisation.
% This function computes the M-term for the plume parameterisation.
% This is the beginning of Equation 28(a) in Lazeroms et al. 2019.
%
% Parameters
%     T_in : scalar (or array?)
%         Ambient temperature in degrees C.
%     S_in : scalar (or array?)
%         Ambient salinity in psu.
%     Tf : scalar (or array?)
%         Freezing temperature
%     c_rho_1 : scalar
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     c_tau : scalar (or array?)
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     gamma: scalar
%         Effective thermal Stanton number. Can be modulated for tuning.
%     E0: scalar
%         Entrainment coefficient. Can be modulated for tuning.
%     alpha: scalar or array
%         Slope angle in rad (must be positive).
%     thermal_forcing_factor : scalar (or array?)
%         Factor to be multiplied to T0-Tf in the end of Eq. 28a - 
%         either thermal forcing or thermal forcing average.
% Returns
%     Mterm : scalar or array
%         Term to be multiplied with M_hat in Eq 28a.

    Mterm=sqrt((par.beta_coeff_lazero.*S_in.*par.g)./(par.lambda3.* ...
        (par.Latent./par.cp0).^3)).*sqrt(max(0,(1-c_rho_1.*gamma)./ ...
        (par.Cd + E0.*sina))).*((gamma.*E0.*sina)./(gamma+c_tau+ ...
        E0.*sina)).^(3/2).*(T_in-Tf).*thermal_forcing_factor;

end


function [sina]=ComputeSinSlope(HB,vx,vy,ctr)

% Kori-ULB
% Calculate local basal slop of ice shelves in Plume model

    s1=circshift(HB,[-1 0]); % i+1,j
    s2=circshift(HB,[1 0]); % i-1,j
    s3=circshift(HB,[0 -1]); % i,j+1
    s4=circshift(HB,[0 1]); % i,j-1
    MV=ones(ctr.imax,ctr.jmax); % u>=1, v>=1
    MV(vx>=0 & vy<0)=2;
    MV(vx<0 & vy>=0)=3;
    MV(vx<0 & vy<0)=4;
    sina=sin(sqrt(((s4-HB)/ctr.delta).^2+((s2-HB)/ctr.delta).^2));
    sina(MV==2)=sin(sqrt(((s4(MV==2)-HB(MV==2))/ctr.delta).^2+ ...
        ((s1(MV==2)-HB(MV==2))/ctr.delta).^2));
    sina(MV==3)=sin(sqrt(((s3(MV==3)-HB(MV==3))/ctr.delta).^2+ ...
        ((s2(MV==3)-HB(MV==3))/ctr.delta).^2));
    sina(MV==4)=sin(sqrt(((s3(MV==4)-HB(MV==4))/ctr.delta).^2+ ...
        ((s1(MV==4)-HB(MV==4))/ctr.delta).^2));
    sina=min(sina,0.05); % limit on slope of ice shelf
    
end


function [x_hat]=ComputeXhat(ice_draft_depth,zGL,T_in,Tf,E0,c_tau, ...
    gamma,sina,par)

% Kori-ULB
% Compute x_hat (or x_tilda) for plume parameterisation.
% This function computes x_hat (or x_tilda) for the plume parameterisation.
% It is a dimensionless coordinate describing distance from plume origin.
% This is Equation 28(b) in Lazeroms et al. 2019.
%
%  Parameters
%     ice_draft_depth : scalar (or array?)
%         Depth of the ice draft in m (depth is negative!).
%     zGL: scalar (or array?)
%         Depth of the grounding line where the source of the plume
%         is in m (depth is negative!).
%     T_in : scalar (or array?)
%         Ambient temperature in degrees C.
%     Tf : scalar (or array?)
%         Freezing temperature
%     E0: scalar
%         Entrainment coefficient. Can be modulated for tuning.
%     c_tau : scalar (or array?)
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     alpha: scalar or array
%         Slope angle in rad (must be positive).
%     gamma: scalar
%         Effective thermal Stanton number. Can be modulated for tuning.
% Returns
%     x_hat : scalar or array
%         Dimensionless coordinate describing distance from plume origin.
%         Has to be between 0 and 1
%     x_tilda in Eq 28b in Lazeroms et al. 2019

    x_hat=par.lambda3*(ice_draft_depth-zGL)./((T_in-Tf).*(1+ ...
        par.C_eps_lazero.*((E0*sina)./(gamma+c_tau+E0*sina)).^(3/4)));
    % all of this derivation is only valid for x in [0,1]
    x_hat=max(min(x_hat,1),0);
    x_hat(zGL>=0)=0;
    
end


function [c_rho_1,c_rho_2,c_tau]=Compute_c_rho_tau(gamma,S_in,par)

% Kori-ULB
% Compute the c-constants for plume parameterisation.
% This function computes c_rho_1, c_rho_2 and c_tau for the plume 
% parameterisation. They are constants given in Table 1 of 
% Lazeroms et al. 2019.
% 
% Parameters
%     gamma: scalar
%         Effective thermal Stanton number. Can be modulated for tuning.
%     S_in : scalar (or array?)
%         Ambient salinity in psu.
% Returns
%     c_rho_1 : scalar (or array?)
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     c_rho_2 : scalar
%         Constant given in Table 1 in Lazeroms et al. 2019.
%     c_tau : scalar (or array?)
%         Constant given in Table 1 in Lazeroms et al. 2019.
    
    c_rho_1=(par.Latent*par.alpha_coeff_lazero)./(par.cp0*gamma* ...
        par.beta_coeff_lazero.* S_in);
    c_rho_2=-par.lambda1*par.alpha_coeff_lazero./par.beta_coeff_lazero;
    c_tau=c_rho_2./c_rho_1;

end


function [contshelfMASK]=ContinentalShelfMASK(MASK,H,B,par)

% Kori-ULB
% Definition of contshelfMASK

    [~,contshelf]=bwboundaries(MASK==0 & B>-2000,'noholes',4);
    MASK(MASK==0 & H>par.SeaIceThickness)=3; % consider ice shelves (MASK==3)
    contshelfMASK=B*0;% Initializing contshelfMASK
    MASK1=circshift(MASK,[0 -1]); % MASK(i,j+1)
    MASK2=circshift(MASK,[0 1]); % MASK(i,j-1)
    MASK3=circshift(MASK,[-1 0]); % MASK(i+1,j)
    MASK4=circshift(MASK,[1 0]); % MASK(i-1,j)

    % Iterate through contshelf items and check if they are somehow connected
    % to either shelves or grounded ice sheet. If so, set contshelfMASK=1 for
    % that item
    for i=1:max(max(contshelf))
        if sum(MASK1(contshelf==i)~=0)+sum(MASK2(contshelf==i)~=0)+ ...
                sum(MASK3(contshelf==i)~=0)+sum(MASK4(contshelf==i)~=0)~=0
            contshelfMASK(contshelf==i)=1;
        end
    end
    
end


function [he,fi]=DefineEdgeThickness(ctr,par,glMASK,H)  %VL: add par
    
% Define ice thickness at the shelf edges used for calving.
% Based on Pollard et al (2015) and Pollard and DeConto (2012)
% Use of maximum thickness of surrounding points and not mean H

    IC=ones(ctr.imax,ctr.jmax);
    IC(glMASK>=5)=0; % only cells grounded/floated and not adjacent to ocean
    IC1=circshift(IC,[-1 0]);
    IC2=circshift(IC,[1 0]);
    IC3=circshift(IC,[0 1]);
    IC4=circshift(IC,[0 -1]);
    H1=circshift(H,[-1 0]);
    H2=circshift(H,[1 0]);
    H3=circshift(H,[0 1]);
    H4=circshift(H,[0 -1]);
    w1=1-min(1,H./(H1*exp(-ctr.delta/1e5)))*(1-exp(-ctr.delta/1e5));
    w2=1-min(1,H./(H2*exp(-ctr.delta/1e5)))*(1-exp(-ctr.delta/1e5));
    w3=1-min(1,H./(H3*exp(-ctr.delta/1e5)))*(1-exp(-ctr.delta/1e5));
    w4=1-min(1,H./(H4*exp(-ctr.delta/1e5)))*(1-exp(-ctr.delta/1e5));
    he=(H1.*w1.*IC1+H2.*w2.*IC2+H3.*w3.*IC3+H4.*w4.*IC4)./ ...
        (IC1+IC2+IC3+IC4);
    he(isnan(he))=H(isnan(he)); % no neighbouring ice-cells (iceberg)
    he(he<par.SeaIceThickness)=par.SeaIceThickness; %VL: make sure he can't become 0
    fi=min(1,H./he);
    fi(glMASK~=5)=1;
    fi(glMASK==6)=0;
    
end


function [uxsch,uysch,d]=DiffusiveCorrection(ctr,par,uxsch,uysch,udx,udy, ...
    d,Ad,p,Hm,gradm,bMASK)

% Kori-ULB
% Correction of diffusion coefficients in ice thickness equation for
% grounding line flux condition.
% This may become obsolete, as diffusion coefficients are not used for SSA

    if ctr.SSA>0
        uxsch=min(max(uxsch,-par.maxspeed),par.maxspeed);   %LZ2021
        uysch=min(max(uysch,-par.maxspeed),par.maxspeed);   %LZ2021
        if ctr.SSAdiffus>0 && ctr.SSA>1
            uxsch=uxsch-udx;
            uysch=uysch-udy;
            d=par.Z*Ad./(p+2).*Hm.^(par.n+2.).*gradm.^((par.n-1.)/2.);
        else
            if ctr.basin==1
                d(bMASK==0)=0;
            else
                d=zeros(ctr.imax,ctr.jmax);
            end
        end
    else % put Schoof/Tsai in d and set ux,uy zero for grounded part
        if ctr.schoof>=1
            dusch=velocity_on_grid(uxsch,uysch);
            d(gradm>1e-20)=dusch(gradm>1e-20).*Hm(gradm>1e-20)./ ...
                sqrt(gradm(gradm>1e-20));
        end
        uxsch=zeros(ctr.imax,ctr.jmax);
        uysch=zeros(ctr.imax,ctr.jmax);
    end
    
end
    

function flux=DpareaWarGds(ewr,nsr,MASK,Bmelt,kegdsx,kegdsy, ...
    par,ctr,flux,gdsmag,funcnt)

% Kori-ULB
% Subglacial meltwater routing function from LeBrocq et al.

    if MASK(ewr,nsr)>=1 && flux(ewr,nsr)<0.
        flux(ewr,nsr)=Bmelt(ewr,nsr)*ctr.delta^2;
        count=1;
        for j=nsr-1:nsr+1
            for i=ewr-1:ewr+1
                pp=zeros(9,1);
                if count~=5
                    if i>=1 && i<ctr.imax && j>=1 && j<ctr.jmax
                        if kegdsy(i,j)>0.
                            pp(6)=kegdsy(i,j)/gdsmag(i,j);
                        elseif kegdsy(i,j)<0.
                            pp(4)=abs(kegdsy(i,j))/gdsmag(i,j);
                        end
                        if kegdsx(i,j)<0.
                            pp(2)=abs(kegdsx(i,j))/gdsmag(i,j);
                        elseif kegdsx(i,j)>0.
                            pp(8)=kegdsx(i,j)/gdsmag(i,j);
                        end
                        if pp(par.dirpp_war(count))>0 && funcnt<=5e4
                            funcnt=funcnt+1;
                            flux=DpareaWarGds(i,j,MASK,Bmelt,kegdsx,kegdsy, ...
                                par,ctr,flux,gdsmag,funcnt);
                            flux(ewr,nsr)=flux(ewr,nsr)+ ...
                                pp(par.dirpp_war(count))*flux(i,j);
                        end
                    end
                end
                count=count+1;
            end
        end
    end
    
end


function [eta,dudx,dvdy,dudy,dvdx]=EffVisc(A,uxssa,uyssa,H,par,eta1,MASK,ctr)

% Kori-ULB
% Effective viscosity of the SSA solution. On the borders of the domain, a
% fixed value is determined based on a fixed ice thickness (Hshelf) that
% determines eta1

    dudx=(uxssa-circshift(uxssa,[0 1]))/ctr.delta;
    dvdy=(uyssa-circshift(uyssa,[1 0]))/ctr.delta;
    dudy=(circshift(uxssa,[-1 1])+circshift(uxssa,[-1 0])- ...
        circshift(uxssa,[1 1])-circshift(uxssa,[1 0]))/(4*ctr.delta);
    dvdx=(circshift(uyssa,[0 -1])+circshift(uyssa,[1 -1])- ...
        circshift(uyssa,[0 1])-circshift(uyssa,[1 1]))/(4*ctr.delta);
    EffStr=dudx.^2+dvdy.^2+dudx.*dvdy+0.25*(dudy+dvdx).^2;
    eta=0.5*H.*A.^(-1./par.n).*EffStr.^((1-par.n)/(2*par.n));
    eta(eta<=0)=NaN;
    eta(MASK==0)=eta(MASK==0)./ctr.shelftune(MASK==0);  %VL: 2D shelftune
    eta(isnan(eta))=1e7;
    
    % NEED TO CHECK THE VISCOSITY ON EDGES OF ICE SHELF/SEA ICE!!!
    % NOW TAKEN AS A CONSTANT VALUE FOR STABILITY AS MOSTLY SEA ICE
    
    if ctr.mismip>=1
        MASKb=ones(ctr.imax,ctr.jmax); % use constant eta on edges of ice shelf
        eta(:,1)=eta(:,3); % ice divide - symmetric
        eta(1,:)=eta(3,:); % periodic BC
        if ctr.mismip==1
            eta(ctr.imax,:)=eta(ctr.imax-2,:); % periodic BC
            MASKb(1:ctr.imax,1:ctr.jmax-2)=0;
        else
            MASKb(1:ctr.imax-1,1:ctr.jmax-2)=0;
        end
        eta(MASKb==1 & MASK==0)=eta1(MASKb==1 & MASK==0);
        % boundary conditions on strain components
        dudx(:,1)=-dudx(:,2);
        dudx(:,ctr.jmax)=dudx(:,ctr.jmax-1);
        dvdx(:,1)=-dvdx(:,2);
        dvdx(:,ctr.jmax)=dvdx(:,ctr.jmax-1);
        dvdy(1,:)=dvdy(3,:);
        dudy(1,:)=dudy(3,:);
        if ctr.mismip==1
            dvdy(ctr.imax,:)=dvdy(ctr.imax-3,:);
            dudy(ctr.imax,:)=dudy(ctr.imax-3,:);
        else
            dvdy(ctr.imax,:)=dvdy(ctr.imax-1,:);
            dudy(ctr.imax,:)=dudy(ctr.imax-1,:);
        end
    else
        MASKb=ones(ctr.imax,ctr.jmax); % use constant eta on edges of ice shelf
        MASKb(3:ctr.imax-2,3:ctr.jmax-2)=0;
        eta(MASKb==1 & MASK==0)=eta1(MASKb==1 & MASK==0);
    end
    eta=min(max(eta,1e5),1e15);
    
end


function [ctr,invmax2D,Asor,ncor,To,So,Pr0,Evp0,runoff0,Evp,Hinit]= ...
    ExistParams(ctr,par,ncor,Asor,stdB,v,uxssa,To,So,Db,B,MASK,As,Pr, ...
    Evp,runoff,Mb0,Hinit,Ho)

% Kori-ULB
% Test what parameters and matrices exist and initializes them accordingly

    % check whether model domain has grounding lines
    if ctr.shelf==1 || ctr.schoof>=1
        ctr.glMASKexist=1;
    else
        ctr.glMASKexist=0;
    end

    % check whether bedrock variability is defined
    if islogical(stdB)==1
        ctr.stdBexist=0; 
    else
        ctr.stdBexist=1;
    end

    % check whether observed velocity field exists
    if islogical(v)==1
        ctr.vexist=0;
    else
        ctr.vexist=1;
    end
    
    % check whether varying lithosphere thickness exists
    if islogical(Db)==1
        ctr.Dbexist=0;
    else
        ctr.Dbexist=1;
    end

    % check whether modelled SSA velocity field exists
    if islogical(uxssa)==1
        ctr.uSSAexist=0;
    else
        ctr.uSSAexist=1;
    end
    
    if islogical(To)==1
        To=zeros(ctr.imax,ctr.jmax)+par.Toi;
        So=zeros(ctr.imax,ctr.jmax)+par.Soi;
    else
        To(isnan(To))=par.Toi;
        So(isnan(So))=par.Soi;
    end
    if ctr.inverse>0
        if islogical(Asor)==1 %VL: check if Asor already exists
            Asor=As;
        end
        if ctr.stdBexist==1
            [invmax2D,ncor]=InitOptimization(ctr,par,ncor,MASK,B,stdB);
        else
            invmax2D=zeros(ctr.imax,ctr.jmax)+par.invmax;
        end
    else
        invmax2D=false;
    end
    if islogical(Pr)==1
        Pr0=Mb0;
    else
        Pr0=Pr;
    end
    if islogical(Evp)==1
        Evp0=Pr0-Mb0;
        Evp=Evp0;
    else
        Evp0=Evp; 
    end
    if islogical(runoff)==1
        runoff0=Pr0-Mb0-Evp0;
    else
        runoff0=runoff;
    end
    if islogical(Hinit)==1
        Hinit=Ho;
    end
end


function [Ts_yc,Pr_yc]=ExtractAnnualCycle(fc,ctr,par,Pr0,Ts0,S0,sn, ...
    MASK,lat,DeltaSL,snp_atm)

% Kori-ULB
% Extract annual cycles from input forcing data (atmospheric data)

    if any(ismember(fields(fc),'atm_Ts_fname'))
        Ts_yc=zeros(ctr.imax,ctr.jmax,fc.atm_yrstep);
        n=1;
        for i=snp_atm-1:snp_atm-1+fc.atm_yrstep-1
            load([fc.atm_Ts_fname,num2str(i,'%03i')]);
            Ts_yc(:,:,n)=Ts;
            n=n+1;
        end
        if ctr.TsType==1
            Ts_yc=Ts_yc+par.Tlapse*(max(sn,DeltaSL)-S0);
        end
    else
        fprintf('If ctr.monthly=1, intra-annual (monthly or submonthly) variations of the air temperature should be provided\n');
        % If intra-annual Ts not provided, then yearly cycle is 
        % approximated -- correction for elevation changes already applied
        
        % Yearly amplitude of the signal (applies to the Antarctic ice sheet)
        if islogical(lat)==1
            Ta=min(20+max(sn,0)/300.,30); % seasonal peak-to-peak amplitude in Ts
        else
            Ta=zeros(ctr.imax,ctr.jmax);
            Ta(MASK==1)=-11.06+0.003204*sn(MASK==1)-0.3584*lat(MASK==1);
            Ta(MASK==0)=-45.06-0.04598*sn(MASK==0)-0.969*lat(MASK==0);
            Ta(Ta<15)=15;
        end
    
        % Define idealized yearly cycle Tm:
        if any(ismember(fields(fc),'atm_yrstep')) 
            par.PDDsteps=fc.atm_yrstep;
        end
        t=repmat(reshape(linspace(1,365.25,par.PDDsteps),1,1,par.PDDsteps), ...
            [ctr.imax,ctr.jmax,1]);
        Ts_yc=repmat(Ts0,[1,1,par.PDDsteps])-0.5* ...
            repmat(Ta,[1,1,par.PDDsteps]).*sin(2*pi*t/365.25);
    end
    if any(ismember(fields(fc),'atm_Pr_fname'))
        Pr_yc=zeros(ctr.imax,ctr.jmax,fc.atm_yrstep);
        n=1;
        for i=snp_atm-1:snp_atm-1+fc.atm_yrstep-1
            load([fc.atm_Pr_fname,num2str(i,'%03i')]);
            Pr_yc(:,:,n)=Pr;
            n=n+1;
        end
        if ctr.MbType==1
            Pr_yc=Pr_yc.*exp(0.05*par.Tlapse*(max(sn,DeltaSL)-S0));
        elseif ctr.MbType==2
            Pr_yc=Pr_yc.*(1+0.053*(par.Tlapse*(max(sn,DeltaSL)-S0)));
        end
    else
        if any(ismember(fields(fc),'atm_yrstep'))
            par.PDDsteps=fc.atm_yrstep;
        end
        Pr_yc=zeros(ctr.imax,ctr.jmax,par.PDDsteps);
        Pr_yc=Pr_yc+Pr0; % If intra-annual Pr not provided, then considered 
                         % constant -- correction for elevation changes 
                         % already applied
    end
end


function Arc=ExtrapolateArc(MASK,oldMASK,Arc,Arc0,ctr)

% Extrapolate Arc between consecutive calls (to reduce
% calculation speed)

    dMASK=MASK-oldMASK; % difference between two consecutive
                        % masks (grounded - floated)
    Me=zeros(ctr.imax,ctr.jmax,8);
    OldArc=Arc;
    OldArc(MASK==1)=NaN;
    Me(:,:,1)=circshift(OldArc,[-1 1]);
    Me(:,:,2)=circshift(OldArc,[0 1]);
    Me(:,:,3)=circshift(OldArc,[1 1]);
    Me(:,:,4)=circshift(OldArc,[1 0]);
    Me(:,:,5)=circshift(OldArc,[-1 0]);
    Me(:,:,6)=circshift(OldArc,[-1 -1]);
    Me(:,:,7)=circshift(OldArc,[0 -1]);
    Me(:,:,8)=circshift(OldArc,[1 -1]);
    NewArc=mean(Me,3,'omitnan'); % Take mean of neighbours for which Arc exists
    NewArc(isnan(NewArc))=0;
    Arc(dMASK==-1)=NewArc(dMASK==-1); % Apply Arc for points becoming floated
    % Apply Arc for points that were previously determined
    Arc(dMASK==-1 & Arc0~=0)=Arc0(dMASK==-1 & Arc0~=0);
    Arc(dMASK==1)=0; % remove melt for points becoming grounded
end


function Wd=ExtrapolateWaterFlux(MASK,oldMASK,Wd,Wd0,ctr,par)

% Interpolate Wd between consecutive calls (to improve
% calculation speed)

    dMASK=MASK-oldMASK; % difference between two consecutive
                        % masks (grounded - floated)
    WD=zeros(ctr.imax,ctr.jmax,8);
    OldWD=Wd;
    OldWD(MASK==0)=NaN;
    WD(:,:,1)=circshift(OldWD,[-1 1]);
    WD(:,:,2)=circshift(OldWD,[0 1]);
    WD(:,:,3)=circshift(OldWD,[1 1]);
    WD(:,:,4)=circshift(OldWD,[1 0]);
    WD(:,:,5)=circshift(OldWD,[-1 0]);
    WD(:,:,6)=circshift(OldWD,[-1 -1]);
    WD(:,:,7)=circshift(OldWD,[0 -1]);
    WD(:,:,8)=circshift(OldWD,[1 -1]);
    NewWD=mean(WD,3,'omitnan'); % Take mean of neighbours for which Wd exists
    NewWD(isnan(NewWD))=par.Wdmin;
    Wd(dMASK==1)=NewWD(dMASK==1); % Apply depth for points becoming grounded
    % Apply flw depths for points that were previously determined
    %   with the subglacial water flow model
    Wd(dMASK==1 & Wd0~=par.Wdmin)=Wd0(dMASK==1 & Wd0~=par.Wdmin);
    Wd(dMASK==-1)=par.Wdmin; % set minimum Wd value for points becoming floated

end


function FinalPlots(ctr,outfile,time,InvVol,SLC,u,v,MASK,bMASK,glMASK)

% Kori-ULB
% Produce basic plots at the end of the model run
% When inverse>0, plots of the quality of the inversion
% When velocity data is availbale, plot of observed versus modelled
% velocity.

    if ctr.inverse>=1
        figure('Name',['InvVol: ' outfile],'NumberTitle','off');
        plot(time,InvVol(:,1)*ctr.delta^2.*1e-9);
        grid on; hold on;
        plot(time,InvVol(:,2)*ctr.delta^2.*1e-9,'r');
        plot(time,InvVol(:,3)*ctr.delta^2.*1e-9,'k');
        xlabel('Time (year)');
        ylabel('|H-H_0| (km^3)');
        legend('Abs. Misfit','Misfit');
    end

    figure('Name',['SLC: ' outfile],'NumberTitle','off');
    plot(time,SLC); grid on; 
    xlabel('Time (year)');
    ylabel('SLR contribution (m)');
    
    if ctr.vexist==1
        figure('Name',['Velocity: ' outfile],'NumberTitle','off');
        if ctr.basin==1
            loglog(v(MASK==1 & bMASK==0),u(MASK==1 & bMASK==0),'.');
        else
            loglog(v(MASK==1),u(MASK==1),'.');
        end
        hold on;
        if ctr.glMASKexist==1
            loglog(v(glMASK==4),u(glMASK==4),'r.');
        end
        plot([1e-1 4e3],[1e-1 4e3],'k-','linewidth',2);
        axis([1e-1 4e3 1e-1 4e3]);
        grid on;
        xlabel('Observed velocity (m a^{-1})');
        ylabel('Modelled velocity (m a^{-1})');
    end
end


function [sn,HB]=FixedFloatation(par,B,SLR,H,MASK)

% Kori-ULB
% height of the bottom of ice shelves based on floatation for a fixed MASK

    HB=B;
    HB(MASK==0)=SLR(MASK==0)-par.rho/par.rhow*H(MASK==0);
    HB(HB<B)=B(HB<B);
    sn=H+HB;
    
end


function [HAF,MASK,HB,sn]=Floatation(par,B,SLR,H,MASK)

% Kori-ULB
% Determination of MASK and height of the bottom of ice shelves based on
% floatation. Height above buoyancy is also returned.

    H(MASK==0 & H<=par.SeaIceThickness)=0;
    HAF=B-SLR+H*par.rho/par.rhow;
    MASK(HAF<0)=0;
    MASK(HAF>=0)=1;
    HB=max(SLR-par.rho/par.rhow*H,B);
    sn=HB+H;
    
end


function [uxsch,uysch]=GroundingLineFlux(ctr,glMASK,HAF,B,SLR,par, ...
    Ax,Ay,Tf,Txx,Tyy,Txy,butfac,ncorx,ncory,ux,uy,Hmx,Hmy,Asfx,Asfy)

% Kori-ULB
% Grounding line flux calculation according to Schoof following
% the method of Pollard and DeConto (2012; 2020)

% LZ2021: Implementation of Schoof following Pollard & DeConto (2012)
% + GL definition on H-grid + weighting factor [Pollard & DeConto (2020)]
% + correction GL definition (glMASK>=3 for adjacent floating grid cell
% instead glMASK==3)
% + correction for first floating grid cell:
%   --> advancing case (wx==0, wy==0): use minimum of velocity of GL grid 
%       cell and schoof flux
% 	--> retreating case (wx>0, wy>0): use schoof flux
% + stricter definition of GL grid cells: 
%   --> SSA velocity at GL has to have correct sign 
%   --> upstream H-grid cell has to be grounded
%   --> upstream velocity has to have correct sign
%   --> downstream H-grid cell has to be floating (if schoof vel. should be
%   applied at downstream grid cell)
% + apply Schoof downstream of GL only if shelf exists
% + orientation of GL (applied for buttressing and schoof velocity)
%       [Pollard & DeConto (2020)]

    uxsch=ux;
    uysch=uy;
    [jGL,iGL] = meshgrid(1:ctr.jmax,1:ctr.imax); %VL: grid for GL orientation

    % GL flux in x-direction
    q0=(Ax*(par.rho*par.g)^(par.n+1)*(1-par.rho/par.rhow)^par.n.* ...
        Asfx.^(1/ctr.m)/4^par.n);
    qs=q0.^(ctr.m/(ctr.m+1));
    qe0=ctr.m*(par.n+2)/(ctr.m+1);

    % conditions for GL in x-direction and ux>0
    qgx=zeros(ctr.imax,ctr.jmax);
    Mgl=zeros(ctr.imax,ctr.jmax);
    angnorm=zeros(ctr.imax,ctr.jmax);   %VL
    nx=ones(ctr.imax,ctr.jmax); %VL
    ny=zeros(ctr.imax,ctr.jmax); %VL
    Theta=ones(ctr.imax,ctr.jmax); %VL
    glMASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
    glMASK0=circshift(glMASK,[0 1]); % glMASK(i,j-1), upstream
    ux1=circshift(ux,[0 1]); % ux(i,j-1), upstream velocity
    HAF1=circshift(HAF,[0 -1]); % HAF(i,j+1)
    B1=circshift(B,[0 -1]); % B(i,j+1)
    Tf1=circshift(Tf,[0 -1]); %VL:  Tf(i,j+1)
    Txx1=circshift(Txx,[0 -1]); %VL:  Txx(i,j+1)
    Tyy1=circshift(Tyy,[0 -1]); %VL:  Tyy(i,j+1)
    Txy1=circshift(Txy,[0 -1]); %VL:  Txy(i,j+1)
    Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & ux>0 & ux1>0)=1; % u-grid(i,j) 
    fracgx=HAF./(HAF-HAF1);
    hbgx=(1.-fracgx).*B+fracgx.*B1;
    hgx=(SLR-hbgx)*par.rhow/par.rho;
    Mgl(hgx<0)=0;
    % GL in x-direction, u>0 --> jg=jGL, jh=jg+1, ih=ig=iGL
    [angnorm(Mgl==1)]=arrayfun(@(iGL,jGL) NormCalc(iGL,jGL,iGL,jGL+1, ...
        glMASK,ctr),iGL(Mgl==1),jGL(Mgl==1));
    nx(Mgl==1)=cos(angnorm(Mgl==1));
    ny(Mgl==1)=sin(angnorm(Mgl==1));
    Theta(Mgl==1)=max(min((butfac*(Txx1(Mgl==1).*nx(Mgl==1).^2+Tyy1(Mgl==1).* ...
        ny(Mgl==1).^2+Txy1(Mgl==1).*nx(Mgl==1).*ny(Mgl==1))+(1.-butfac)* ...
        Tf1(Mgl==1))./Tf1(Mgl==1),1),0);
    if ctr.schoof==1 % Schoof
        Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
    else % Pattyn
        Theta=Theta.^(par.n);
    end
    if ctr.schoof==1 % Schoof
        ugx=qs.*hgx.^qe0.*Theta;
    else % Pattyn
        ugx=q0.*hgx.^(par.n+3).*Theta;
    end
    ugx(Mgl==1)=ugx(Mgl==1).*abs(nx(Mgl==1));
    qgx(Mgl==1)=ugx(Mgl==1).*hgx(Mgl==1);
    wx=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
    wx(Mgl==1)=max(0,min(1,(ugx(Mgl==1)-ux(Mgl==1)).*hgx(Mgl==1)./10^5));
    % GL grid cell
    uxsch(Mgl==1)=ux(Mgl==1)+ncorx(Mgl==1).*wx(Mgl==1).*(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01)-ux(Mgl==1));
    % ice shelf grid cell downstream of GL
    Mgl=circshift(Mgl,[0 1]); % condition for i,j+1
    Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
    Mgl(Hmx==0)=0; % only apply if shelf exists
    qgx=circshift(qgx,[0 1]);
    ncx=circshift(ncorx,[0 1]);
    wx=circshift(wx,[0 1]);
    uxsch(Mgl==1)=ux(Mgl==1)+ncx(Mgl==1).*( ...
        sign(wx(Mgl==1)).*(1-wx(Mgl==1)).*(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01)-ux(Mgl==1)) ... % wx>0 (retreat)
        +(1-sign(wx(Mgl==1))).*(min(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01),ux1(Mgl==1))-ux(Mgl==1)));   % wx=0 (advance)
    
    % conditions for GL in x-direction and ux<0
    qgx=zeros(ctr.imax,ctr.jmax);
    Mgl=zeros(ctr.imax,ctr.jmax);
    angnorm=zeros(ctr.imax,ctr.jmax)+pi;   %VL
    nx=zeros(ctr.imax,ctr.jmax)-1; %VL
    ny=zeros(ctr.imax,ctr.jmax); %VL
    Theta=ones(ctr.imax,ctr.jmax); %VL
    glMASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
    glMASK0=circshift(glMASK,[0 -2]); % glMASK(i,j+2), upstream
    ux1=circshift(ux,[0 -1]);   % ux(i,j+1), upstream velocity
    HAF1=circshift(HAF,[0 -1]); % HAF(i,j+1)
    B1=circshift(B,[0 -1]); % B(i,j+1)
    Mgl(glMASK>=3 & glMASK1==2 & glMASK0<=2 & ux<0 & ux1<0)=1; % u-grid(i,j)
    fracgx=HAF1./(HAF1-HAF);
    hbgx=(1.-fracgx).*B1+fracgx.*B;
    hgx=(SLR-hbgx)*par.rhow/par.rho;
    Mgl(hgx<0)=0;
    %VL: GL in x-direction, u<0 --> jh=jGL, jg=jh+1, ih=ig=iGL
    [angnorm(Mgl==1)]=arrayfun(@(iGL,jGL) NormCalc(iGL,jGL+1,iGL,jGL, ...
        glMASK,ctr),iGL(Mgl==1),jGL(Mgl==1));  %VL
    nx(Mgl==1)=cos(angnorm(Mgl==1));
    ny(Mgl==1)=sin(angnorm(Mgl==1));
    Theta(Mgl==1)=max(min((butfac*(Txx(Mgl==1).*nx(Mgl==1).^2+Tyy(Mgl==1).* ...
        ny(Mgl==1).^2+Txy(Mgl==1).*nx(Mgl==1).*ny(Mgl==1))+(1.-butfac)* ...
        Tf(Mgl==1))./Tf(Mgl==1),1),0);
    if ctr.schoof==1 % Schoof
        Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
    else % Pattyn
        Theta=Theta.^(par.n);
    end
    if ctr.schoof==1 % Schoof
        ugx=qs.*hgx.^qe0.*Theta;
    else % Pattyn
        ugx=q0.*hgx.^(par.n+3).*Theta;
    end
    ugx(Mgl==1)=ugx(Mgl==1).*abs(nx(Mgl==1));
    qgx(Mgl==1)=-ugx(Mgl==1).*hgx(Mgl==1);
    wx=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
    wx(Mgl==1)=max(0,min(1,(ugx(Mgl==1)-abs(ux(Mgl==1))).*hgx(Mgl==1)./10^5));
    % GL grid cell
    uxsch(Mgl==1)=ux(Mgl==1)+ncorx(Mgl==1).*wx(Mgl==1).*(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01)-ux(Mgl==1));
    % ice shelf grid cell downstream of GL
    Mgl=circshift(Mgl,[0 -1]); % condition for i,j-1
    Mgl(glMASK<=2)=0; % only apply if downstream grid cell is floating
    Mgl(Hmx==0)=0; % only apply if shelf exists
    qgx=circshift(qgx,[0 -1]);
    ncx=circshift(ncorx,[0 -1]);
    wx=circshift(wx,[0 -1]);
    uxsch(Mgl==1)=ux(Mgl==1)+ncx(Mgl==1).*( ...
        sign(wx(Mgl==1)).*(1-wx(Mgl==1)).*(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01)-ux(Mgl==1)) ... % wx>0 (retreat)
        +(1-sign(wx(Mgl==1))).*(max(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01),ux1(Mgl==1))-ux(Mgl==1)));   % wx=0 (advance)
    
    
    % GL flux in y-direction
    q0=(Ay*(par.rho*par.g)^(par.n+1)*(1-par.rho/par.rhow)^par.n.* ...
        Asfy.^(1/ctr.m)/4^par.n);
    qs=q0.^(ctr.m/(ctr.m+1));

    % conditions for GL in y-direction and uy>0
    qgy=zeros(ctr.imax,ctr.jmax);
    Mgl=zeros(ctr.imax,ctr.jmax);
    angnorm=zeros(ctr.imax,ctr.jmax)+pi/2;   %VL
    nx=zeros(ctr.imax,ctr.jmax); %VL
    ny=ones(ctr.imax,ctr.jmax); %VL
    Theta=ones(ctr.imax,ctr.jmax); %VL
    glMASK1=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
    glMASK0=circshift(glMASK,[1 0]); % glMASK(i-1,j), upstream
    uy1=circshift(uy,[1 0]);    % uy(i-1,j), upstream velocity
    HAF1=circshift(HAF,[-1 0]); % HAF(i+1,j)
    B1=circshift(B,[-1 0]); % B(i+1,j)
    Tf1=circshift(Tf,[-1 0]); %VL:  Tf(i+1,j)
    Txx1=circshift(Txx,[-1 0]); %VL:  Txx(i+1,j)
    Tyy1=circshift(Tyy,[-1 0]); %VL:  Tyy(i+1,j)
    Txy1=circshift(Txy,[-1 0]); %VL:  Txy(i+1,j)
    Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & uy>0 & uy1>0)=1; % v-grid(i,j)
    fracgy=HAF./(HAF-HAF1);
    hbgy=(1.-fracgy).*B+fracgy.*B1;
    hgy=(SLR-hbgy)*par.rhow/par.rho;
    Mgl(hgy<0)=0;
    %VL: GL in y-direction, v>0 --> ig=iGL, ih=ig+1, jh=jg=jGL
    [angnorm(Mgl==1)]=arrayfun(@(iGL,jGL) NormCalc(iGL,jGL,iGL+1,jGL, ...
        glMASK,ctr),iGL(Mgl==1),jGL(Mgl==1));  %VL
    nx(Mgl==1)=cos(angnorm(Mgl==1));
    ny(Mgl==1)=sin(angnorm(Mgl==1));
    Theta(Mgl==1)=max(min((butfac*(Txx1(Mgl==1).*nx(Mgl==1).^2+ ...
        Tyy1(Mgl==1).*ny(Mgl==1).^2+Txy1(Mgl==1).*nx(Mgl==1).*ny(Mgl==1))+ ...
        (1.-butfac)*Tf1(Mgl==1))./Tf1(Mgl==1),1),0);
    if ctr.schoof==1 % Schoof
        Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
    else % Pattyn
        Theta=Theta.^(par.n);
    end
    if ctr.schoof==1 % Schoof
        ugy=qs.*hgy.^qe0.*Theta;
    else % Pattyn
        ugy=q0.*hgy.^(par.n+3).*Theta;
    end
    ugy(Mgl==1)=ugy(Mgl==1).*abs(ny(Mgl==1));
    qgy(Mgl==1)=ugy(Mgl==1).*hgy(Mgl==1);
    wy=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
    wy(Mgl==1)=max(0,min(1,(ugy(Mgl==1)-uy(Mgl==1)).*hgy(Mgl==1)./10^5));
    % GL grid cell
    uysch(Mgl==1)=uy(Mgl==1)+ncory(Mgl==1).*wy(Mgl==1).*(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01)-uy(Mgl==1));
    % ice shelf grid cell downstream of GL
    Mgl=circshift(Mgl,[1 0]); % condition for i+1,j
    Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
    Mgl(Hmy==0)=0; % only apply if shelf exists
    qgy=circshift(qgy,[1 0]);
    ncy=circshift(ncory,[1 0]);
    wy=circshift(wy,[1 0]);
    uysch(Mgl==1)=uy(Mgl==1)+ncy(Mgl==1).*( ...
        sign(wy(Mgl==1)).*(1-wy(Mgl==1)).*(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01)-uy(Mgl==1)) ...   % wy>0 (retreat)
        +(1-sign(wy(Mgl==1))).*(min(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01),uy1(Mgl==1))-uy(Mgl==1)));   % wy=0 (advance)
    
    % conditions for GL in y-direction and uy<0
    qgy=zeros(ctr.imax,ctr.jmax);
    Mgl=zeros(ctr.imax,ctr.jmax);
    angnorm=zeros(ctr.imax,ctr.jmax)-pi/2;   %VL
    nx=zeros(ctr.imax,ctr.jmax); %VL
    ny=zeros(ctr.imax,ctr.jmax)-1; %VL
    Theta=ones(ctr.imax,ctr.jmax); %VL
    glMASK1=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
    glMASK0=circshift(glMASK,[-2 0]); % glMASK(i+2,j), upstream
    uy1=circshift(uy,[-1 0]);   % uy(i+1,j), upstream velocity
    HAF1=circshift(HAF,[-1 0]); % HAF(i+1,j)
    B1=circshift(B,[-1 0]); % B(i+1,j)
    Mgl(glMASK>=3 & glMASK1==2 & glMASK0<=2 & uy<0 & uy1<0)=1; % v-grid(i,j)
    fracgy=HAF1./(HAF1-HAF);
    hbgy=(1.-fracgy).*B1+fracgy.*B;
    hgy=(SLR-hbgy)*par.rhow/par.rho;
    Mgl(hgy<0)=0;
    %VL: GL in y-direction, v<0 --> ih=iGL, ig=ih+1, jh=jg=jGL
    [angnorm(Mgl==1)]=arrayfun(@(iGL,jGL) NormCalc(iGL+1,jGL,iGL,jGL, ...
        glMASK,ctr),iGL(Mgl==1),jGL(Mgl==1));  %VL
    nx(Mgl==1)=cos(angnorm(Mgl==1));
    ny(Mgl==1)=sin(angnorm(Mgl==1));
    Theta(Mgl==1)=max(min((butfac*(Txx(Mgl==1).*nx(Mgl==1).^2+Tyy(Mgl==1).* ...
        ny(Mgl==1).^2+Txy(Mgl==1).*nx(Mgl==1).*ny(Mgl==1))+(1.-butfac)* ...
        Tf(Mgl==1))./Tf(Mgl==1),1),0);
    if ctr.schoof==1 % Schoof
        Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
    else % Pattyn
        Theta=Theta.^(par.n);
    end
    if ctr.schoof==1 % Schoof
        ugy=qs.*hgy.^qe0.*Theta;
    else % Pattyn
        ugy=q0.*hgy.^(par.n+3).*Theta;
    end
    ugy(Mgl==1)=ugy(Mgl==1).*abs(ny(Mgl==1));
    qgy(Mgl==1)=-ugy(Mgl==1).*hgy(Mgl==1);
    wy=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
    wy(Mgl==1)=max(0,min(1,(ugy(Mgl==1)-abs(uy(Mgl==1))).*hgy(Mgl==1)./10^5));
    % GL grid cell
    uysch(Mgl==1)=uy(Mgl==1)+ncory(Mgl==1).*wy(Mgl==1).*(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01)-uy(Mgl==1));
    % ice shelf grid cell downstream of GL
    Mgl=circshift(Mgl,[-1 0]); % condition for i-1,j
    Mgl(glMASK<=2)=0; % only apply if downstream grid cell is floating
    Mgl(Hmy==0)=0; % only apply if shelf exists
    qgy=circshift(qgy,[-1 0]);
    ncy=circshift(ncory,[-1 0]);
    wy=circshift(wy,[-1 0]);
    uysch(Mgl==1)=uy(Mgl==1)+ncy(Mgl==1).*( ...
        sign(wy(Mgl==1)).*(1-wy(Mgl==1)).*(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01)-uy(Mgl==1)) ...   % wy>0 (retreat)
        +(1-sign(wy(Mgl==1))).*(max(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01),uy1(Mgl==1))-uy(Mgl==1)));   % wy=0 (advance)
    
end


function [glMASK,MASK]=GroundingLineMask(MASK,bMASK,H,ctr)

% Kori-ULB
% Definition of glMASK when grounding lines appear in the domain
%   glMASK=1: grounded
%   glMASK=2: grounding line
%   glMASK=3: first floating grid point
%   glMASK=4: ice shelf
%   glMASK=5: calving front
%   glMASK=§: open ocean

    glMASK=MASK;
    MASK1=circshift(MASK,[0 -1]); % MASK(i,j+1)
    MASK2=circshift(MASK,[0 1]); % MASK(i,j-1)
    MASK3=circshift(MASK,[-1 0]); % MASK(i+1,j)
    MASK4=circshift(MASK,[1 0]); % MASK(i-1,j)

    glMASK(MASK==1 & (MASK1==0 | MASK2==0 | MASK3==0 | MASK4==0))=2;

    MASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
    MASK2=circshift(glMASK,[0 1]); % glMASK(i,j-1)
    MASK3=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
    MASK4=circshift(glMASK,[1 0]); % glMASK(i-1,j)

    % adjacent floating grid cell (=3)
    glMASK(MASK==0 & (MASK1==2 | MASK2==2 | MASK3==2 | MASK4==2))=3;

    if ctr.shelf==1
        glMASK(glMASK==0)=4; % ice shelf (=4)
        if ctr.calving>=1
            glMASK(glMASK==4 & H<5)=6;  % sea ice/ocean (=6)
            MASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
            MASK2=circshift(glMASK,[0 1]); % glMASK(i,j-1)
            MASK3=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
            MASK4=circshift(glMASK,[1 0]); % glMASK(i-1,j)
            glMASK((glMASK==3 | glMASK==4) & (MASK1==6 | MASK2==6 | ...
                MASK3==6 | MASK4==6))=5; % Calving front (=5)
        end
    end
    if ctr.mismip>=1
        glMASK(:,1)=glMASK(:,3);
        glMASK(1,:)=glMASK(3,:);
        if ctr.mismip==1
            glMASK(ctr.imax,:)=glMASK(ctr.imax-2,:);
        else
            glMASK(:,ctr.jmax)=glMASK(:,ctr.jmax-1);
            glMASK(ctr.imax,:)=glMASK(ctr.imax-1,:);
        end
    end
    if ctr.basin==1
        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(1,:)=1;
        MASKb(ctr.imax,:)=1;
        MASKb(:,1)=1;
        MASKb(:,ctr.jmax)=1;
        glMASK(MASKb==1 & MASK==1)=1;
        glMASK(MASKb==1 & MASK==0)=6;
        glMASK(bMASK==1)=1;
    end

end


function [uxsch,uysch]=GroundingLines(ctr,par,glMASK,HAF,SLR,Ax,Ay,butfac, ...
    H,Hmx,Hmy,B,Bmx,Bmy,Asfx,Asfy,ux,uy,dudx,dvdy,dudy,dvdx,eta)

% Kori-ULB
% Caluclation of buttress factors and parameterized flux condition at the
% grounding line for large-scale model simulations.

    if ctr.shelf==1
        Txx=2*eta.*(2.*dudx+dvdy)./max(par.SeaIceThickness,H);
        Tyy=2*eta.*(2.*dvdy+dudx)./max(par.SeaIceThickness,H);
        Txy=2*eta.*(dudy+dvdx)./max(par.SeaIceThickness,H); %VL: adaption to Pollard & deConto [2020]
        Tf=.5*par.rho*par.g*max(par.SeaIceThickness,H)*(1.-par.rho/par.rhow); % on h-grid
    else
        %VL: also provide Tf, Txx, Tyy, Txy for GL flux
        Tf=.5*par.rho*par.g*max(par.SeaIceThickness,H)*(1.-par.rho/par.rhow); % on h-grid
        Txx=Tf;
        Tyy=Tf;
        Txy=Tf;
    end
    % Applying Schoof BC and heuristic from Pollard
    % Weighting factor for Schoof correction as function bed
    ncorx=max(min(1-1/50*(Bmx+50),1),0);
    ncory=max(min(1-1/50*(Bmy+50),1),0);
    [uxsch,uysch]=GroundingLineFlux(ctr,glMASK,HAF,B,SLR,par,Ax,Ay, ...
        Tf,Txx,Tyy,Txy,butfac,ncorx,ncory,ux,uy,Hmx,Hmy,Asfx,Asfy); %VL: new function + new call!

end


function Melt=ISMIP6OceanSlope(sina,par,ctr,ice_draft,icemask_shelves, ...
    TF_draft,basinNumber,deltaT_basin)

% Kori-ULB
% ISMIP6 ocean parameterization equation (1) + adapted to
% take the local slope into account
 
    cste = (par.rhow*par.cp0/(par.rho*par.Latent)).^2;  % in K^(-2)

    Nbasin=max(max(basinNumber))+1;

    Melt=zeros(ctr.imax,ctr.jmax);
    mean_TF=zeros(Nbasin,1);

    basinNumber=basinNumber+1;
    for i=1:Nbasin
        mean_TF(i)=mean(TF_draft(basinNumber==i & icemask_shelves>5.e-1 ...
            & ice_draft<=0.));
    end

    % Step3 calculation of melting rate
    % melt rate in m/yr (meters of pure water per year):
    % [*rhofw_SI/rhoi_SI to get it in meters of ice per year]

    for i=1:Nbasin
        Melt(basinNumber==i)=ctr.gammaT.*cste.*(TF_draft(basinNumber==i) ...
            +deltaT_basin(basinNumber==i)).*abs(mean_TF(i)+ ...
            deltaT_basin(basinNumber==i)).*sina((basinNumber==i));
    end
    Melt(TF_draft<=-5)=0;
    Melt=Melt*par.rhof/par.rho;
end


function Melt=ISMIP6ocean(par,ctr,ice_draft,icemask_shelves,TF_draft, ...
    basinNumber,deltaT_basin)

% Kori-ULB
% ISMIP6 ocean parameterization equation (1)
 
    cste=(par.rhow*par.cp0/(par.rho*par.Latent)).^2;  % in K^(-2)
    Nbasin=max(max(basinNumber))+1;

    Melt=zeros(ctr.imax,ctr.jmax);
    mean_TF=zeros(Nbasin,1);

    basinNumber=basinNumber+1;
    for i=1:Nbasin
        mean_TF(i)=mean(TF_draft(basinNumber==i & icemask_shelves>5.e-1 & ...
            ice_draft<=0.));
    end

    % Step3 calculation of melting rate
    % melt rate in m/yr (meters of pure water per year):
    % [*rhofw_SI/rhoi_SI to get it in meters of ice per year]

    for i=1:Nbasin
        Melt(basinNumber==i)=ctr.gammaT.*cste.*(TF_draft(basinNumber==i) ...
            +deltaT_basin(basinNumber==i)).*abs(mean_TF(i)+ ...
            deltaT_basin(basinNumber==i));
    end
    Melt(TF_draft<=-5)=0;
    Melt=Melt*par.rhof/par.rho;
end

function [ctr,fc]=InitCtr(ctr,fc,default)

% Kori-ULB
% Initialization of ctr structure
% When not defined a priori, default values are assumed

    ctr.plotH(any(ismember(fields(ctr),'plotH'))==0)=0;
    ctr.SSAdiffus(any(ismember(fields(ctr),'SSAdiffus'))==0)=0;
        % (1) Use diffusion for deformational velocity in SSA=2
        % (2) idem, but hybrid model is addition of SSA (sliding) and deformational velocity
    ctr.runmode(any(ismember(fields(ctr),'runmode'))==0)=0;
    ctr.restart(any(ismember(fields(ctr),'restart'))==0)=0;
    ctr.diagnostic(any(ismember(fields(ctr),'diagnostic'))==0)=0;
    ctr.inverse(any(ismember(fields(ctr),'inverse'))==0)=0;
    ctr.timeslice(any(ismember(fields(ctr),'timeslice'))==0)=0;
    ctr.meltfac(any(ismember(fields(ctr),'meltfac'))==0)=0;
    ctr.meltfunc(any(ismember(fields(ctr),'meltfunc'))==0)=0;
    ctr.slidfac(any(ismember(fields(ctr),'slidfac'))==0)=1; % default no sliding perturbation
    ctr.schoof(any(ismember(fields(ctr),'schoof'))==0)=0;
    ctr.shelf(any(ismember(fields(ctr),'shelf'))==0)=0;
    ctr.SSA(any(ismember(fields(ctr),'SSA'))==0)=0;
    ctr.MbType(any(ismember(fields(ctr),'MbType'))==0)=0;
    ctr.TsType(any(ismember(fields(ctr),'TsType'))==0)=0;
    ctr.Tcalc(any(ismember(fields(ctr),'Tcalc'))==0)=0;
    ctr.Tinit(any(ismember(fields(ctr),'Tinit'))==0)=0;
    ctr.BedAdj(any(ismember(fields(ctr),'BedAdj'))==0)=0;
    ctr.m(any(ismember(fields(ctr),'m'))==0)=default.m; % default linear sliding
    ctr.p(any(ismember(fields(ctr),'p'))==0)=0;
    ctr.kmax(any(ismember(fields(ctr),'kmax'))==0)=default.kmax; % default number of z-levels    
    ctr.subwaterflow(any(ismember(fields(ctr),'subwaterflow'))==0)=0;
    ctr.SlidAdjust(any(ismember(fields(ctr),'SlidAdjust'))==0)=0;
    ctr.calving(any(ismember(fields(ctr),'calving'))==0)=0;
    ctr.HydroFrac(any(ismember(fields(ctr),'HydroFrac'))==0)=0;
    ctr.GeoidCalc(any(ismember(fields(ctr),'GeoidCalc'))==0)=0;
    ctr.starttime(any(ismember(fields(ctr),'starttime'))==0)=0;
    ctr.mismip(any(ismember(fields(ctr),'mismip'))==0)=0;
    ctr.basin(any(ismember(fields(ctr),'basin'))==0)=0;
    ctr.damage(any(ismember(fields(ctr),'damage'))==0)=0;
    ctr.GroundedMelt(any(ismember(fields(ctr),'GroundedMelt'))==0)=0;
    ctr.PDDcalc(any(ismember(fields(ctr),'PDDcalc'))==0)=0;
    ctr.monthly(any(ismember(fields(ctr),'monthly'))==0)=0;
    ctr.runoffcorr(any(ismember(fields(ctr),'runoffcorr'))==0)=0;
    ctr.PDD_anomaly(any(ismember(fields(ctr),'PDD_anomaly'))==0)=0;
    ctr.Hinv(any(ismember(fields(ctr),'Hinv'))==0)=default.Hinv;
    ctr.Tinv(any(ismember(fields(ctr),'Tinv'))==0)=default.Tinv;
    ctr.stopoptim(any(ismember(fields(ctr),'stopoptim'))==0)=default.stopoptim;
    ctr.HinvMelt(any(ismember(fields(ctr),'HinvMelt'))==0)=default.HinvMelt;
    ctr.TinvMelt(any(ismember(fields(ctr),'TinvMelt'))==0)=default.TinvMelt;
    ctr.radnorm(any(ismember(fields(ctr),'radnorm'))==0)=default.radnorm;
    ctr.snapshot(any(ismember(fields(ctr),'snapshot'))==0)=default.snapshot;
    ctr.BetaIter(any(ismember(fields(ctr),'BetaIter'))==0)=default.BetaIter;
    if any(ismember(fields(ctr),'shelftune'))==0
        ctr.shelftune=zeros(ctr.imax,ctr.jmax)+default.shelftune;
    elseif numel(ctr.shelftune)==1
        ctr.shelftune=zeros(ctr.imax,ctr.jmax)+ctr.shelftune;
    end
    ctr.meltfactor(any(ismember(fields(ctr),'meltfactor'))==0)=default.meltfactor;
    ctr.Ao(any(ismember(fields(ctr),'Ao'))==0)=default.Ao;
    ctr.u0(any(ismember(fields(ctr),'u0'))==0)=default.u0;
    ctr.plotGL(any(ismember(fields(ctr),'plotGL'))==0)=default.plotGL;
    ctr.upstream(any(ismember(fields(ctr),'upstream'))==0)=default.upstream;
    ctr.ItSolv(any(ismember(fields(ctr),'ItSolv'))==0)=default.ItSolv;
    if any(ismember(fields(ctr),'Asin'))==0
        ctr.Asin=zeros(ctr.imax,ctr.jmax)+1e-10;
    end
    if any(ismember(fields(ctr),'gammaT'))==0
        if ctr.meltfunc==1
            ctr.gammaT=default.gammaTlin;
        elseif ctr.meltfunc==2
            ctr.gammaT=default.gammaTquad;
        elseif ctr.meltfunc==9
            ctr.gammaT=default.gamma0_quad;
        elseif ctr.meltfunc==91
            ctr.gammaT=default.gamma0_quad_slope;
        else
            ctr.gammaT=default.gammaTpico;
        end
    end
    ctr.C(any(ismember(fields(ctr),'C'))==0)=default.Cpico;
    ctr.gammaTplume(any(ismember(fields(ctr),'gammaTplume'))==0)=default.gammaTplume;
    ctr.M0(any(ismember(fields(ctr),'M0'))==0)=default.M0picop;
    
    if any(ismember(fields(fc),'DeltaT'))==0
        fc.DeltaT=zeros(ctr.nsteps,1);
    end
    if any(ismember(fields(fc),'DeltaSL'))==0
        fc.DeltaSL=zeros(ctr.nsteps,1);
    end
    if any(ismember(fields(fc),'butfac'))==0
        fc.butfac=ones(ctr.nsteps,1);
    end
end


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


function [Pg0,geoide,frg]=InitGeoid(par,ctr,MASK,H0,B0,B,SLR)

% Kori-ULB
% Initialization of geoid fields

    Pg0=zeros(ctr.imax,ctr.jmax);
    Pg0(MASK==1)=par.rho*H0(MASK==1)*ctr.delta^2.;
    Pg0(MASK==0)=(-par.rhow)*(B(MASK==0)-SLR(MASK==0))*ctr.delta^2.;
    Pg0=Pg0+par.rhom*B0*ctr.delta^2.; % addition of bed change (if applicable)
    [gkx,gky]=meshgrid(-par.geoidist:ctr.delta:par.geoidist,-par.geoidist:ctr.delta:par.geoidist);
    frg=round(par.geoidist/ctr.delta-0.5);
    dist=sqrt(gkx.^2+gky.^2);
    geoide=par.Re./(2*par.Me*sin(dist/(2*par.Re)));
    geoidmax=par.Re./(2*par.Me*sin(10/(2*par.Re)));
    geoide(geoide>geoidmax)=0;
    geoide(geoide<1.5e-18)=0;
    
end


function [As,Mb0,Ts0,MASK0,MASK,bMASK,H,B,H0,B0]= ...
    InitInputData(ctr,par,As,Mb,Ts,MASK,MASKo,H,B)

% Kori-ULB
% Initialization of input matrices after reading the input file

    As=As*ctr.slidfac;
    Mb0=Mb; % initial mass balance if read from inputfile
    Ts0=Ts; % intial surface temperature if read from inputfile
    MASK0=MASK; % MASK at start of run (MASKo is origiginal BedMap mask, if exists)
    MASK(MASK==3)=0; % ice shelves in BedMap2/Bedmachine identified with MASK=3

    bMASK=zeros(ctr.imax,ctr.jmax);
    if ctr.basin==1
        bMASK(MASKo==-1)=1;
        MASK(bMASK==1)=1;
        As(bMASK==1)=par.As0;
        par.intT=1; % calculate temperature at each time step
        H(1,:)=H(2,:);
        B(1,:)=B(2,:);
        H(ctr.imax,:)=H(ctr.imax-1,:);
        B(ctr.imax,:)=B(ctr.imax-1,:);
        H(:,1)=H(:,2);
        B(:,1)=B(:,2);
        H(:,ctr.jmax)=H(:,ctr.jmax-1);
        B(:,ctr.jmax)=B(:,ctr.jmax-1);
    end

    H0=H; % initial ice thickness
    B0=B; % initial bed elevation
end


function [snapshot,plotst,cnt_atm,snp_atm,cnt_ocn,snp_ocn,Mb_update,Li,Lj, ...
    dtdx,dtdx2,X,Y,x,y,MASK,H,Ho,B,Bo,MASKo,Mb,Ts,As,G,u,VAF,VA0,POV, ...
    SLC,Ag,Af,Btau,IVg,IVf,glflux,dHdt,time,mbcomp,InvVol,ncor,dSLR,SLR, ...
    Wd,Wtil,Bmelt,CMB,FMB,flw,p,px,py,pxy,nodeu,nodev,nodes,node,VM,Tof, ...
    Sof,TFf,Tsf,Mbf,Prf,Evpf,runofff,Melt,damage]= ...
    InitMatrices(ctr,par,default,fc)
    
% Kori-ULB
% Initialization of main matrices used in the model

    plotst=floor(ctr.nsteps/ctr.snapshot);
    snapshot=floor(ctr.nsteps/plotst)+1;

    cnt_atm=0;
    snp_atm=0;
    cnt_ocn=0;
    snp_ocn=0;
    Mb_update=0;
    
    %-------------------------------------------
    % Numerical parameters and grid definition
    %-------------------------------------------
    
    Li=(ctr.imax-1)*ctr.delta/1e3; % length of domain in y
    Lj=(ctr.jmax-1)*ctr.delta/1e3; % length of domain in x
    dtdx=ctr.dt/(2.*ctr.delta*ctr.delta);
    dtdx2=ctr.dt/(2.*ctr.delta);
    [X,Y] = meshgrid(0:ctr.delta/1e3:Lj,0:ctr.delta/1e3:Li);
    x=0:ctr.delta/1e3:Lj;
    y=0:ctr.delta/1e3:Li;

    %------------------------------
    % Definition of matrix sizes
    %------------------------------
    
    [H,Ho,B,Bo,Mb,Ts,u,dSLR,CMB,FMB,Tof,Sof,TFf,Tsf,Mbf,Prf,Evpf, ...
        runofff,Melt,damage]=deal(zeros(ctr.imax,ctr.jmax));
    [MASK,MASKo,ncor]=deal(ones(ctr.imax,ctr.jmax));
    [VAF,VA0,POV,SLC,Ag,Af,IVg,IVf,glflux,dHdt]=deal(zeros(ctr.nsteps,1));
    G=zeros(ctr.imax,ctr.jmax)+default.G0;
    
    if ctr.BedAdj>0
        Btau=zeros(ctr.imax,ctr.jmax)+par.bedrelax;
    else
        Btau=false;
    end
    SLR=dSLR+fc.DeltaSL(1); % total sea level (background + fingerprint)
    time=(ctr.starttime:ctr.dt:ctr.starttime+(ctr.nsteps-1)*ctr.dt)';
    mbcomp=zeros(ctr.nsteps,21);
    if ctr.inverse>=1
        InvVol=zeros(ctr.nsteps,3);
    else
        InvVol=false;
    end
    p=zeros(ctr.imax,ctr.jmax)+par.n;
    As=ctr.Asin;
    if ctr.subwaterflow>0
        Wd=zeros(ctr.imax,ctr.jmax)+1e-8;
    else
        Wd=false;
    end
    [Wtil,Bmelt,flw]=deal(Wd);
    [px,py,pxy]=deal(p);
    
    % initialization of sparse matrix system for velocities
    nodeu=linspace(1,ctr.imax*ctr.jmax*2-1,ctr.imax*ctr.jmax)';
    nodeu=reshape(nodeu,[ctr.imax ctr.jmax]);
    nodev=linspace(2,ctr.imax*ctr.jmax*2,ctr.imax*ctr.jmax)';
    nodev=reshape(nodev,[ctr.imax ctr.jmax]);

    % Node count for sparse matrix for thickness solver
    nodes=ctr.imax*ctr.jmax;
    node=linspace(1,nodes,nodes)';
    node=reshape(node,[ctr.imax ctr.jmax]);
    VM=ones(nodes,1);

end


function [invmax2D,ncor]=init_optimization(ctr,par,ncor,MASK,B,stdB)

% Kori-ULB
% Initialization of control matrices for the optimization (initialization)
% of the ice sheet

    stdB2=BedVar(B,ctr); %VL: large scale bed variability for new ncor definition
    %VL: definition of new ncor
    ncor((stdB2>=prctile(stdB2(MASK==1),95) | ...
        stdB>=prctile(stdB(MASK==1),99)) & MASK==1 & B>0)=0;
    % Definition of invmax2D to decrease invmax in areas with high
    % stdB -- max variation of 3 orders of magnitude
    invmax2D=zeros(ctr.imax,ctr.jmax)+log10(par.invmax);
    invmax2D(ncor==0)=log10(par.invmaxncor)-3.*((stdB(ncor==0)- ...
        prctile(stdB(ncor==0),10))/(prctile(stdB(ncor==0),90)- ...
        prctile(stdB(ncor==0),10)));
    invmax2D=10.^(invmax2D);
    invmax2D(ncor==0 & invmax2D>par.invmaxncor)=par.invmaxncor;
    invmax2D(invmax2D<1e-8)=1e-8; % limit to 3 orders of magnitude

end


function tmp=InitTemp3d(G,taudxy,ub,ud,par,H,Mb,zeta,ctr,Ts,MASK,DeltaT)

% Kori-ULB
% Temperature field initialization with analytical solution for temperature
% profiles; different analytical solution for ice shelves

    uvel=ub+par.udfrac*ud;
    Tgrad=-(G+taudxy.*uvel/par.secperyear)/par.K;
    l=sqrt(2*par.kdif*(H+1e-8)./max(Mb,1e-8)*par.secperyear);
    repl=repmat(l,[1,1,ctr.kmax]);
    repH=repmat(H+1e-8,[1,1,ctr.kmax]);
    repTs=repmat(Ts,[1,1,ctr.kmax]);
    repTgrad=repmat(Tgrad,[1,1,ctr.kmax]);
    repz=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    tmp=repTs+sqrt(pi)*0.5*repl.*repTgrad.* ...
        (erf((1-repz).*repH./repl)-erf(repH./repl))+par.T0;
    Tp=par.pmp*repH.*repz;

    % Ice shelf and open ocean
    repmask=repmat(MASK,[1,1,ctr.kmax]);
    repM=repmat(max(Mb,1e-5),[1,1,ctr.kmax]);
    TBshelf=repmat(min(par.Toi+ctr.meltfactor*DeltaT- ...
        0.12e-3*par.rho*H/par.rhow,0),[1,1,ctr.kmax]);
    c2=exp(repM.*repH/(par.kdif*par.secperyear));
    c1=exp(repM.*repz.*repH/(par.kdif*par.secperyear));
    shelftmp=((repTs-TBshelf).*c1+TBshelf-repTs.*c2)./(1-c2)+par.T0;
    shelftmp(repH<=1)=TBshelf(repH<=1)+par.T0;
    tmp(repmask==0)=shelftmp(repmask==0);

    % Correction for pmp
    tmp(tmp>par.T0-Tp)=par.T0-Tp(tmp>par.T0-Tp);
    
end


function [tmp,Tb,zeta,dzc,dzp,dzm]=InitTempParams(ctr,par,tmp,Ts)

% Kori-ULB
% Initialization of vertical derivatives for temperature profiles in the 3d
% temperature calculation

    if islogical(tmp)==1
        tmp=repmat(Ts+par.T0,[1,1,ctr.kmax]);
        Tb=Ts;
    else
        Tb=tmp(:,:,ctr.kmax)-par.T0;
    end
    zeta=CalculateZeta(ctr.kmax,0.015);
    dzc=circshift(zeta,[0 -1])-circshift(zeta,[0 1]);
    dzp=circshift(zeta,[0 -1])-zeta;
    dzm=zeta-circshift(zeta,[0 1]);
    dzc=repmat(reshape(dzc,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    dzp=repmat(reshape(dzp,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    dzm=repmat(reshape(dzm,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);

end


function [var_int]=InterpToDepthOfInt(var,depth_levels,depth_of_int,ctr)

% Kori-ULB
% Interpolation of ocean variables to ice shelf depth

    var_int=zeros(ctr.imax,ctr.jmax);    
    var_up=var(:,:,1);
    var_int(depth_of_int>=depth_levels(1))=var_up(depth_of_int>=depth_levels(1));
    var_down=var(:,:,end);
    var_int(depth_of_int<=depth_levels(end))=var_down(depth_of_int<=depth_levels(end));
    for k=length(depth_levels):-1:2
        var_up=var(:,:,k);
        var_down=var(:,:,k-1);
        var_int(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1))=...
            ((depth_levels(k)-depth_of_int(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1))).*var_down(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1))...
            +(depth_of_int(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1))-depth_levels(k-1)).*var_up(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1)))...
            ./(depth_levels(k)-depth_levels(k-1));
    end
%     var_int(icemask_shelves<5.e-1 | depth_of_int>0.)=-9999.9;

end


function LSF=LSFfunction(LSF,ctr,u,v,node,nodes,VM,MASK)

% Kori-ULB
% Calculate the Level Set Function (LSF) for following the calving front.
% Used in the calving algorihms
% Still under development

    epsilon=1e-2;
    
    dtdx2=epsilon*ctr.dt/(ctr.delta*ctr.delta);
    dtdx=ctr.dt/ctr.delta;
    
    % velocity on LSF grid;
    cx=0.5*(u+circshift(u,[0 1]));
    cy=0.5*(v+circshift(v,[1 0]));

    MASKb=zeros(ctr.imax,ctr.jmax);
    MASKb(2:ctr.imax-1,2:ctr.jmax-1)=1;
    
    % conditions for diffusion scheme (init)
    V0=zeros(ctr.imax,ctr.jmax)+4*dtdx2; % i,j
    V1=zeros(ctr.imax,ctr.jmax)-dtdx2; % i,j+1
    V2=zeros(ctr.imax,ctr.jmax)-dtdx2; % i,j-1
    V3=zeros(ctr.imax,ctr.jmax)-dtdx2; % i+1,j
    V4=zeros(ctr.imax,ctr.jmax)-dtdx2; % i-1,j

    % Velocity sign masks
    MU=ones(ctr.imax,ctr.jmax);
    MU(cx<=0)=2;
    MV=ones(ctr.imax,ctr.jmax);
    MV(cy<=0)=2;
    MU(MASKb==0)=0;
    MV(MASKb==0)=0;

    V0a=zeros(ctr.imax,ctr.jmax);
    V1a=zeros(ctr.imax,ctr.jmax);
    V2a=zeros(ctr.imax,ctr.jmax);
    V3a=zeros(ctr.imax,ctr.jmax);
    V4a=zeros(ctr.imax,ctr.jmax);

    % conditions for MU=1 (grad(u)>0)
    V0a(MU==1)=cx(MU==1); % i,j
    V2a(MU==1)=-cx(MU==1); % i,j-1

    % conditions for MU=2 and (grad(u)<0)
    V0a(MU==2)=-cx(MU==2); % i,j
    V1a(MU==2)=cx(MU==2); % i,j+1

    % conditions for MV=1 (grad(v)>0)
    V0a(MV==1)=V0a(MV==1)+cy(MV==1); % i,j
    V4a(MV==1)=-cy(MV==1); % i-1,j

    % conditions for MV=2 (grad(v)<0)
    V0a(MV==2)=V0a(MV==2)-cy(MV==2); % i,j
    V3a(MV==2)=cy(MV==2); % i+1,j

    % Filling V-matrix
    V0=V0+V0a*dtdx;%.*(1.-alfa);
    V1=V1+V1a*dtdx;%.*(1.-alfa);
    V2=V2+V2a*dtdx;%.*(1.-alfa);
    V3=V3+V3a*dtdx;%.*(1.-alfa);
    V4=V4+V4a*dtdx;%.*(1.-alfa);

    R0=LSF;
    
    % boundaries
    V10=zeros(ctr.imax,ctr.jmax); % i=1 periodic boundary or ocean
    V11=zeros(ctr.imax,ctr.jmax); % i=imax periodic boundary or ocean

    wholemask=ctr.imax*ctr.jmax-sum(MASK(:));
    if wholemask~=0 && ctr.mismip>=1 % only when domain is not MASK=1 everywhere
        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(1,2:ctr.jmax-1)=1; % periodic BC at i=1 or ice divide
        V0(MASKb==1)=0;
        V10(MASKb==1)=-1;
        R0(MASKb==1)=0;
        
        if ctr.mismip==1
            MASKb=zeros(ctr.imax,ctr.jmax);
            MASKb(ctr.imax,2:ctr.jmax-1)=1; % periodic BC at i=imax or ocean
            V0(MASKb==1)=0;
            V11(MASKb==1)=-1;
            R0(MASKb==1)=0;
        else
            MASKb=zeros(ctr.imax,ctr.jmax);
            MASKb(2:ctr.imax-1,1)=1; % periodic BC at i=1 
            V0(MASKb==1)=0;
            V11(MASKb==1)=-1;
            R0(MASKb==1)=0;
        end
    end

    V=[reshape(V0(VM==1)+1,nodes,1)
        V1(V1~=0)
        V2(V2~=0)
        V3(V3~=0)
        V4(V4~=0)
        V10(V10~=0)
        V11(V11~=0)];

    row=[reshape(node(VM==1),nodes,1)
        node(V1~=0)
        node(V2~=0)
        node(V3~=0)
        node(V4~=0)
        node(V10~=0)
        node(V11~=0)];

    nodeV1=circshift(node,[0 -1]); % i,j+1
    nodeV2=circshift(node,[0 1]); % i,j-1
    nodeV3=circshift(node,[-1 0]); % i+1,j
    nodeV4=circshift(node,[1 0]); % i-1,j

    if ctr.mismip>=1
        nodeV10=circshift(node,[-2 0]); % 3,j - symmetry at i=1
        if ctr.mismip==1
            nodeV11=circshift(node,[2 0]); % n-2,j - PBC at i=imax
        elseif ctr.mismip==2
            nodeV11=circshift(node,[0 -2]); % n-2,j - PBC at i=imax
        else
            nodeV11=zeros(ctr.imax,ctr.jmax);
        end
    else
        nodeV10=zeros(ctr.imax,ctr.jmax);
        nodeV11=zeros(ctr.imax,ctr.jmax);
    end

    col=[reshape(node(VM==1),nodes,1)
        nodeV1(V1~=0)
        nodeV2(V2~=0)
        nodeV3(V3~=0)
        nodeV4(V4~=0)
        nodeV10(V10~=0)
        nodeV11(V11~=0)];

    R=reshape(R0(VM==1),nodes,1);

    % construct sparse matrix
    A=sparse(row,col,V);
    % Cholesky factor and solve
    tol=1e-3; % 1e-5
    maxit=40; % 15
    D=diag(diag(A));
    C1=tril(A);
    C2=D\triu(A);
    [s,flag]=pcg(A,R,tol,maxit,C1,C2);
    if flag>0
        s=A\R;
    end

    LSF(node>0)=s(node(node>0));
end


function mbcomp=MBcomponents(ctr,par,acc,Smelt,runoff,rain,Mb,Pr, ...
    H,Hn,Bmelt,Melt,CMB,FMB,MASK,bMASK,mbcomp,B,Bn,SLR)

% Kori-ULB
% Mass balance components of the ice sheet and ice shelf

% Global ice sheet components (ice sheet + ice shelf)
% (1) surface mass balance
% (2) dhdt
% (3) accumulation
% (4) surface melt
% (5) runoff
% (6) rain
% (7) basal melt grounded ice
% (8) sub-shelf melt
% (9) calving flux
% (10) ice shelf frontal melt
% (11) dynamical component

    if ctr.basin==1
        MASK=MASK+bMASK; % MASK=2 for outside basin domain and not accounted for
    end
    mbcomp(1)=VariabFlux(max(-H,Mb),H,MASK,par.SeaIceThickness); % smb
    mbcomp(2)=(sum(Hn(Hn>par.SeaIceThickness))-sum(H(H>par.SeaIceThickness)))/ctr.dt;
    if ctr.PDDcalc==1
        mbcomp(3)=VariabFlux(acc,H,MASK,par.SeaIceThickness); % acc
        mbcomp(4)=VariabFlux(min(H+acc,Smelt),H,MASK,par.SeaIceThickness); % Smelt
        mbcomp(5)=VariabFlux(min(H+acc,runoff),H,MASK,par.SeaIceThickness); % runoff
        mbcomp(6)=VariabFlux(rain,H,MASK,par.SeaIceThickness); % rain
    else
        mbcomp(3)=VariabFlux(Pr,H,MASK,par.SeaIceThickness); % acc
        mbcomp(5)=VariabFlux(min(H+Pr,runoff),H,MASK,par.SeaIceThickness); % runoff
    end
    if ctr.Tcalc>=1
        mbcomp(7)=VariabFlux(min(H+max(-H,Mb),Bmelt),H,MASK,par.SeaIceThickness); % Bmelt
    end
    if ctr.meltfunc>=1 && ctr.glMASKexist==1
        mbcomp(8)=VariabFlux(min(H+max(-H,Mb),Melt),H,MASK,par.SeaIceThickness);
    end
    if ctr.calving>=1 && ctr.shelf==1
        mbcomp(9)=VariabFlux(min(H+max(-H,Mb),CMB),H,MASK,par.SeaIceThickness); % CMB
        mbcomp(10)=VariabFlux(min(H+max(-H,Mb),FMB),H,MASK,par.SeaIceThickness); % FMB
    end
    mbcomp(11)=mbcomp(1)-mbcomp(2)-mbcomp(7)-mbcomp(8);
    
% Grounded ice sheet components
% (12) surface mass balance
% (13) dhdt
% (14) accumulation
% (15) surface melt
% (16) runoff
% (17) rain
% (18) basal melt grounded ice
% (19) grounding line flux
% (20) Rate of VAF
% (21) dynamical VAF component
    
    MASK(MASK==0)=NaN;
    mbcomp(12)=VariabFlux(max(-H,Mb),H,MASK,par.SeaIceThickness); % smb
    mbcomp(13)=(sum(Hn(MASK==1))-sum(H(MASK==1)))/ctr.dt;
    if ctr.PDDcalc==1
        mbcomp(14)=VariabFlux(acc,H,MASK,par.SeaIceThickness); % acc
        mbcomp(15)=VariabFlux(min(H+acc,Smelt),H,MASK,par.SeaIceThickness); % Smelt
        mbcomp(16)=VariabFlux(min(H+acc,runoff),H,MASK,par.SeaIceThickness); % runoff
        mbcomp(17)=VariabFlux(rain,H,MASK,par.SeaIceThickness); % rain
    else
        mbcomp(14)=VariabFlux(Pr,H,MASK,par.SeaIceThickness); % acc
        mbcomp(16)=VariabFlux(min(H+Pr,runoff),H,MASK,par.SeaIceThickness); % runoff
    end
    if ctr.Tcalc>=1
        mbcomp(18)=VariabFlux(min(H+max(-H,Mb),Bmelt),H,MASK,par.SeaIceThickness); % Bmelt
    end
    mbcomp(19)=mbcomp(12)-mbcomp(13)-mbcomp(18);
    VAF=max(0,H+min(B-SLR,0)*par.rhow/par.rho);
    VAFn=max(0,Hn+min(Bn-SLR,0)*par.rhow/par.rho);
    mbcomp(20)=(sum(sum(VAFn-VAF)))/ctr.dt;
    mbcomp(21)=mbcomp(12)-mbcomp(20)-mbcomp(18);
    % convert all components to Gt/a (water equivalent)
    mbcomp=mbcomp*ctr.delta^2*par.rho/1e12;
    
end


function [Pr,Evp,runoff]=MbFunc(ctr,Prf,Evpf,runoff,runofff,Ts0,Ts,sn,X,Y,Lj,Li,DeltaT)

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
            % Pollard et al., Garbe2020
            Pr=Prf.*exp(0.05*(Ts-Tsf));
            if ctr.PDDcalc==0
                if ctr.runoffcorr==1
                    runoff=runofff+1.805*(exp(0.5745*Ts)-exp(0.5745*Tsf)); % inferred from MAR monthly data: runoff=1.805*exp(0.5745*t2m)
                else
                    runoff=runofff;
                end
            end
            Evp=Evpf; % Evp not corrected for elevation change
        case 2
            % Correction for elevation changes - OPTION 2
            % Golledge 2015 (future runs, see Frieler)
            Pr=Prf.*(1+0.053*(Ts-Tsf));
            if ctr.PDDcalc==0
                if ctr.runoffcorr==1
                    runoff=runofff+1.805*(exp(0.5745*Ts)-exp(0.5745*Tsf)); % inferred from MAR monthly data: runoff=1.805*exp(0.5745*t2m)
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


function MeltDown()

% Kori-ULB
% No Comment

    figure;
    for i=1:20
        k=rand-0.5;
        hold on; clf;
        text(1+5*k,6+5*k,'All is lost','FontSize',24,'color',rand(1,3));
        text(1+5*k,4+5*k,'Complete meltdown','FontSize',24,'color',rand(1,3));
        axis([0 10 0 10]);
        axis off;
        hold off;
        pause(1);
    end
    
end


function [m]=MscaledPoly(x,p)

% Kori-ULB
% Polynomial function for Plume melt parameterization underneath ice
% shelves

    m=zeros(size(x));
    for i=1:12
        m=m+p(i)*x.^(i-1);
    end
    
end


function [angnorm]=NormCalc(ig,jg,ih,jh,glMASK,ctr)

% Kori-ULB
% Function normcalc (to calculate the direction of the grounding line)

nwin = round(ctr.radnorm*1e3/ctr.delta);  %radnorm in [km]
nwin = nwin + 2; %make sure to have all grid cells within radius radnorm

zavx = 0;
zavy = 0;

for j=jh-nwin:jh+nwin
    jj=max(1,min(j,ctr.jmax));
    for i=ih-nwin:ih+nwin
        ii=max(1,min(i,ctr.imax));
        if glMASK(ii,jj)==0 || glMASK(ii,jj)>2
            zdx = (j - 0.5*(jg+jh));
            zdy = (i - 0.5*(ig+ih));
            zdist = sqrt(zdx^2 + zdy^2)*ctr.delta/1e3; 
            if zdist>0 && zdist<=ctr.radnorm
                zavx = zavx + zdx; 
                zavy = zavy + zdy;
            end
        end
    end
end
angnorm = atan2(zavy,zavx);

end


function [damage]=NyeDamage(par,ctr,dudx,dvdy,dudy,dvdx,eta,H,HAF)

% Kori-ULB
% compute total crevasse depth d, either 0, the surface crevasse depth or
% the total, whichever is the larger, given the vertically integrated
% 'damaged' stretching stress t
%
% The crevasse depths are given by the Nye zero-stress rule, used recently
% in e.g Nick et al., 2010 JoG (eq.4 and 6)-our t is their Rxx
%
% After Sun et al. (2017)

    eps=1e-8;
    dw=zeros(ctr.imax,ctr.jmax);
    [lambda0,~]=PrincipalStrain(dudx,dvdy,dudy,dvdx);
    lambda0=2*lambda0.*eta; % convert strain to stress
    % basal crevasses
    db=(lambda0./(par.rho*par.g*(H+eps))-max(HAF,0))*par.rho/(par.rhow-par.rho);
    % surface crevasses
    ds=lambda0./(par.rho*par.g*(H+eps))+par.rhow*dw/par.rho;
%     damage=max(0,min(max(ds,ds+db),H*par.dlim));
    damage=max(0,min(max(ds,ds+db),H*par.dlim));
    
end


function [To,So,TF]=OCEANdepthOfInt(fc,ctr,par,Tof,Sof,TFf,H,HB,B,MASK,glMASK,ShelfN,numsh)

% Kori-ULB
% Interpolate ocean variables to depth of interest (bottom of ice shelf)

    if any(ismember(fields(fc),'z')) % Depth-profile of To and So data or TF
        if ctr.meltfunc==3 || ctr.meltfunc==4 
            % compute average depth of continental shelf for each ice shelf
            % front (Burgard et al. 2021 -- PROTECT)
            front_bot_depth_avg=zeros(ctr.imax,ctr.jmax)-500;
            for b=1:numsh
                avg_B=mean(B(ShelfN==b & glMASK==5));
                if isnan(avg_B)==1
                    avg_B=mean(B(ShelfN==b)); % some small shelves have front 
                                              % not well define, then average 
                                              % over shelf area
                end
                front_bot_depth_avg(ShelfN==b)=avg_B;
            end
            depth_of_int=front_bot_depth_avg;
            To=InterpToDepthOfInt(Tof,fc.z,depth_of_int,ctr);
            So=InterpToDepthOfInt(Sof,fc.z,depth_of_int,ctr);
            if islogical(TFf)==0
                TF=InterpToDepthOfInt(TFf,fc.z,depth_of_int,ctr);
            else
                TF=TFf;
            end
        else
            % compute deepest entrance depth of continental shelf for each 
            % ice shelf (Burgard et al. 2021 -- PROTECT)
            front_bot_depth_max=zeros(ctr.imax,ctr.jmax)-2500;
            for b=1:numsh
                max_B=min(B(ShelfN==b & glMASK==5));
                if isempty(max_B)==1
                    max_B=min(B(ShelfN==b)); % some small shelves have front 
                                             % not well define, then average 
                                             % over shelf area
                end
                front_bot_depth_max(ShelfN==b)=max_B;
            end
            depth_of_int=HB; % interpolate at ice draft
            depth_of_int(HB<front_bot_depth_max)= ...
                front_bot_depth_max(HB<front_bot_depth_max);
                % Limit depth of int to the deepest entrance depth 
                % of continental shelf
            if islogical(TFf)==0
                TF=InterpToDepthOfInt(TFf,fc.z,depth_of_int,ctr);
                To=Tof;
                So=Sof;
            else
                TF=TFf;
                To=InterpToDepthOfInt(Tof,fc.z,depth_of_int,ctr);
                So=InterpToDepthOfInt(Sof,fc.z,depth_of_int,ctr);
            end
        end
    else
        if (ctr.meltfunc==3 || ctr.meltfunc==4) && ...
                any(ismember(fields(fc),'PICO_basins')) 
            % PICO method of updating values by basin (Kreutzer et al. 2021)
            [contshelfMASK]=ContinentalShelfMASK(MASK,H,B,par);
                % update continental shelf mask
            SO=MASK*0;
            TO=MASK*0;
            for b=1:19
                % compute average of edge grid cells within a specific
                % PICO basin and assign that value to So or To field
                avg_SO=mean(mean(Sof(fc.PICO_basins==b & contshelfMASK==1)));
                avg_TO=mean(mean(Tof(fc.PICO_basins==b & contshelfMASK==1)));
                SO(fc.PICO_basins==b)=avg_SO;
                TO(fc.PICO_basins==b)=avg_TO;
            end
            To=TO;
            So=SO;
            TF=TFf;
        else
            TF=TFf;
            To=Tof;
            So=Sof;
        end
    end

end


function [Tof,Sof,TFf,cnt_ocn,snp_ocn]=OCEANupdate(fc,ctr,time,cnt,So0, ...
    To0,Tof,Sof,TFf,cnt_ocn,snp_ocn)

% Kori-ULB
% Update ocean forcing based on external forcing data

    if fc.forcingOCEAN==1
        if time(cnt)==fc.ocn_Tinit % Initialise at first snapshot year
            cnt_ocn=1;
            snp_ocn=1;
        end
        if cnt_ocn==1 % snapshot
            if time(cnt)<=fc.ocn_Tend
                if any(ismember(fields(fc),'ocn_TF_fname'))
                    load([fc.ocn_TF_fname,num2str(snp_ocn,'%03i')]);
                    TFf=double(TF); % make sure matrix is double
                    Tof=To0+ctr.meltfactor*fc.DeltaT(cnt); 
                    % simplified forcing otherwise
                    Sof=So0;
                else % If TF forcing does not exist, check for To or So forcing
                    TFf=false;
                    if any(ismember(fields(fc),'ocn_To_fname'))
                        load([fc.ocn_To_fname,num2str(snp_ocn,'%03i')]);
                        Tof=double(To); % make sure matrix is double
                    else
                        Tof=To0+ctr.meltfactor*fc.DeltaT(cnt); 
                        % simplified forcing otherwise
                    end
                    if any(ismember(fields(fc),'ocn_So_fname'))
                        load([fc.ocn_So_fname,num2str(snp_ocn,'%03i')]);
                        Sof=double(So); % make sure matrix is double
                    else
                        Sof=So0;
                    end
                end
                snp_ocn=snp_ocn+1;
            else % If forcing shorter than simulation time
                if snp_ocn==fc.ocn_snapshots+1
                    snp_ocn=fc.ocn_snapshots+1-fc.ocn_nrep; 
                    % Re-initialise snapshot
                end
                if any(ismember(fields(fc),'ocn_TF_fname'))
                    load([fc.ocn_TF_fname,num2str(snp_ocn,'%03i')]);
                    TFf=double(TF); % make sure matrix is double
                    Tof=To0+ctr.meltfactor*fc.DeltaT(cnt); 
                    % simplified forcing otherwise
                    Sof=So0;
                else % If TF forcing does not exist, check for To or So forcing
                    TFf=false;
                    if any(ismember(fields(fc),'ocn_To_fname'))
                        load([fc.ocn_To_fname,num2str(snp_ocn,'%03i')]);
                        Tof=double(To); % make sure matrix is double
                    else
                        Tof=To0+ctr.meltfactor*fc.DeltaT(cnt); 
                        % simplified forcing otherwise
                    end
                    if any(ismember(fields(fc),'ocn_So_fname'))
                        load([fc.ocn_So_fname,num2str(snp_ocn,'%03i')]);
                        Sof=double(So); % make sure matrix is double
                    else
                        Sof=So0;
                    end
                end
                snp_ocn=snp_ocn+1;
            end
        end
        cnt_ocn=cnt_ocn+1;
        cnt_ocn(cnt_ocn>fc.ocn_cnt)=1;
    else
        TFf=false;
        Tof=To0+ctr.meltfactor*fc.DeltaT(cnt); % simplified forcing otherwise
        Sof=So0;
    end

end


function [arcocn,distocn,distgl]=OceanArc(MASK,H,MASKlk,ctr,par)

% Kori-ULB
% Calculation of the angle of ice shelves to open ocean (Pollard and
% DeConto)

    arcocn=zeros(ctr.imax,ctr.jmax);
    distocn=zeros(ctr.imax,ctr.jmax)+1000e3;
    distgl=zeros(ctr.imax,ctr.jmax)+10000e3;
    narc=72; % number of directions (every 5 degree)
    angarc=zeros(1,narc);
    tanarc=zeros(1,narc);
    ifi=zeros(1,narc);
    idir=zeros(1,narc);
    MASK(MASK==0&H>par.SeaIceThickness)=3;  %mark shelves in MASK
    MASK1=circshift(MASK,[0 -1]); % MASK(i,j+1)
    MASK2=circshift(MASK,[0 1]); % MASK(i,j-1)
    MASK3=circshift(MASK,[-1 0]); % MASK(i+1,j)
    MASK4=circshift(MASK,[1 0]); % MASK(i-1,j)
    MASK5=circshift(MASK,[-1 -1]); % MASK(i+1,j+1)
    MASK6=circshift(MASK,[-1 1]); % MASK(i+1,j-1)
    MASK7=circshift(MASK,[1 -1]); % MASK(i-1,j+1)
    MASK8=circshift(MASK,[1 1]); % MASK(i-1,j-1)
    ifdo=zeros(ctr.imax,ctr.jmax);
    ifdo(MASK~=1&MASKlk==0&(MASK1~=0|MASK2~=0|MASK3~=0|MASK4~=0|MASK5~=0|MASK6~=0|MASK7~=0|MASK8~=0))=1;
    [io,jo] = meshgrid(1:ctr.imax,1:ctr.jmax);
    distocn(MASK==0)=0;

    for mm=1:narc
        angarc(mm)=-pi+(mm-0.5)*(2*pi)/narc;
        tanarc(mm)=tan(angarc(mm));
        if abs(angarc(mm))>=0.75*pi
            ifi(mm)=1;
            idir(mm)=-1;
        elseif angarc(mm)>=-0.75*pi && angarc(mm)<=-0.25*pi
            ifi(mm)=0;
            idir(mm)=-1;
        elseif abs(angarc(mm))<=0.25*pi
            ifi(mm)=1;
            idir(mm)=1;
        else
            ifi(mm)=0;
            idir(mm)=1;
        end
    end

    [arcocn(ifdo==1),distocn(ifdo==1),distgl(ifdo==1)]=arrayfun(@(io,jo) CalcOceanArc(io,jo,MASK,ctr.imax,ctr.jmax,ctr.delta,narc,tanarc,ifi,idir),io(ifdo==1),jo(ifdo==1));

    arcocn(MASK~=1&ifdo==0&MASKlk==0)=360;
    distocn(MASK==0)=0;
    distgl(MASK~=1&ifdo==0)=0;
end


function [T0o,S0o,zb]=OceanVarBox(numsh,To,So,ShelfN,HB)

% Kori-ULB
% Temperature and salinity in boxes of the PICO model

    T0o=zeros(numsh,1); % mean box T0
    S0o=zeros(numsh,1); % mean box S0
    zb=zeros(numsh,1); % mean box HB

    for k=1:numsh
        T0o(k)=mean(To(ShelfN==k));
        S0o(k)=mean(So(ShelfN==k));
        zb(k)=mean(HB(ShelfN==k));
    end
end


function [As,deltaZ,Asor]=OptimizeIceSheet(ctr,par,cnt, ...
    Asor,MASK,MASKo,bMASK,deltaZ,sn,sn0,r,ncor,B,stdB,vx,vy,ux,uy,invmax2D)

% Kori-ULB
% Optimization of basal sliding coefficients underneath grounded ice sheet
% using the nudging algorithm of Pollard and DeConto (2012)

    % introduce deltaZold to compare deltaZ between consecutive iterations
    deltaZold=zeros(ctr.imax,ctr.jmax);     %LZ
    % optimize only for grounded AND originally grounded grid cells
    if cnt*ctr.dt>ctr.Tinv
        deltaZold(MASK==1 & MASKo==1)=deltaZ(MASK==1 & MASKo==1);
    end
    deltaZ=zeros(ctr.imax,ctr.jmax);
    deltaZ(MASK==1 & MASKo==1)=max(-1.5,min(1.5,(sn(MASK==1 & MASKo==1) ...
        -sn0(MASK==1 & MASKo==1))/ctr.Hinv));
    if ctr.SlidAdjust==1 && ctr.Tcalc>0
        Asor(r>0 & abs(deltaZ)>=abs(deltaZold))=Asor(r>0 & ...
            abs(deltaZ)>=abs(deltaZold)).*10.^deltaZ(r>0 & ...
            abs(deltaZ)>=abs(deltaZold));
    else
        Asor(abs(deltaZ)>=abs(deltaZold))=Asor(abs(deltaZ)>= ...
            abs(deltaZold)).*10.^deltaZ(abs(deltaZ)>=abs(deltaZold));
    end
    Asor(Asor<par.invmin/10)=par.invmin/10;
    Asor(Asor>par.invmax*1000)=par.invmax*1000;
    Asor(Asor>invmax2D & ncor==0)=invmax2D(Asor>invmax2D & ncor==0);
    if ctr.vexist==1
        As = RegularizationNew(Asor,B,stdB,vx,vy,ctr,ncor,par.stdDevRegul);
    else
        As = RegularizationNew(Asor,B,stdB,ux,uy,ctr,ncor,par.stdDevRegul);
    end
    As(As<par.invmin)=par.invmin; % first limits
    As(As>par.invmax)=par.invmax;
    if ctr.basin==1 %VL: As on h-grid
        As(bMASK==1)=par.As0;
    end
end


function [Melt]=OptimizeIceShelf(ctr,MASKo,glMASK,H,Ho,Melt,bMASK)

% Kori-ULB
% Optimization of basal melt and accretion underneath ice shelves according
% to the method by Bernales

    Ftan=1.725;
    Melt=Melt+Ftan.*tan(max(-1.3,min(1.3,(H-Ho)./ctr.HinvMelt)));
    % smaller limit as Bernales (changes of Melt were too big for limit 
    % 1.5 --> created holes in Ross)
    Melt(MASKo==1)=0; % No Melt for grounded grid cells
    Melt=min(100,max(-100,Melt)); % limit melt
    Melt(glMASK==6)=0;
    if ctr.basin==1
        Melt(bMASK==1)=0;
    end
    
end


function [runoff,acc,rain,Smelt]=PDDmonthly(Ts_yc,Pr_yc,par,ctr)

% Kori-ULB
% Scheme for monthly PDD calculation

    Tm=Ts_yc; % Ts yearly cycle
    Ts=mean(Ts_yc,3); % Yearly mean
    Pr=Pr_yc; % Pr yearly cycle
    par.PDDsteps=size(Tm,3);

    % PDD calculation from Calov and Greve (2005)
    nT=(Tm-par.PDDth)/(sqrt(2)*par.Tsigma);
    T=(par.Tsigma/sqrt(2*pi)*exp(-(nT.^2))+0.5*(Tm-par.PDDth).*erfc(-nT)) ...
        +par.PDDth;

    % Rain fraction calculation
    nP=(Tm-par.Tsnow)/(sqrt(2)*par.Psigma);
    P=(par.Psigma/sqrt(2*pi)*exp(-(nP.^2))+0.5*(Tm-par.Tsnow).*erfc(-nP)) ...
        +par.Tsnow;
    PDD=T*365.25/par.PDDsteps;
    R=(1-max(0,min((par.Train-P)/(par.Train-par.Tsnow),1)));
    
    % Melt model
    snow_depth=zeros(ctr.imax,ctr.jmax);
    nyear=2; % Sufficient to initialize snow depth

    % Loop over timesteps through one year
    for y=1:nyear % run through several years to equilibrium (Tsai 2020)
        Smelt=zeros(ctr.imax,ctr.jmax); % Annual ice  and snow melt
        acc=zeros(ctr.imax,ctr.jmax);
        rain=zeros(ctr.imax,ctr.jmax);
        ESM=zeros(ctr.imax,ctr.jmax);
        EIM=zeros(ctr.imax,ctr.jmax);
        for nstep=1:par.PDDsteps
            PR=Pr(:,:,nstep)/par.PDDsteps; % Total Precipitation during nstep
            acc_nstep=max(0,PR.*(1-R(:,:,nstep))); 
            % Total Snow Accumulation during nstep. R is the fraction 
            % of precipitation that falls as rain
            rain_nstep=max(0,PR.*R(:,:,nstep)); % Total Rain during nstep
            snow_depth=max(0,snow_depth+acc_nstep);
            ESM_nstep=min(max(snow_depth,0),par.snowfac.*PDD(:,:,nstep)); 
            % Effective snow melt (cannot exceed the amount of snow) 
            % in m ice equivalent
            PDDr=max(PDD(:,:,nstep)-ESM_nstep./par.snowfac,0); 
            % remaining PDDs after snow melt
            snow_depth=max(0,snow_depth-ESM_nstep);
            EIM_nstep=par.icefac.*PDDr; % Effective ice melt
            Smelt_nstep=ESM_nstep+EIM_nstep;

            % Increment annual quantities
            ESM=ESM+ESM_nstep;
            EIM=EIM+EIM_nstep;
            Smelt=Smelt+Smelt_nstep; % combined ice and snow melt
            acc=acc+acc_nstep;
            rain=rain+rain_nstep;
        end
    end
    Pr=mean(Pr,3);

    % Refreezing model
    W_r=ESM+rain; % Available water mass - maximum amount of water that 
                  % can be retained provided enough energy is available
    P_r=max(-(min(Ts,0))*par.d_ice*par.K/(par.kdif*par.rho*par.Latent),0); 
    % Potential retention mass - maximum amount of liquid water that 
    % can be retained provided enough mass is available
    E_r=min(Pr,min(P_r,W_r)); 
    % Effective Refreezing - amount of water that is effectively 
    % retained at the end of the melting season
    % Er is limited by the annual precipitation
    runoff=EIM+W_r-E_r; % ice melt + rain and snowmelt liquid water 
                        % beyond the snowpack saturation
 
end


function [runoff,acc,rain,Smelt]=PDDyearly(ctr,par,Ts,Pr,lat,MASK,sn)

% Kori-ULB
% Scheme for yearly PDD model caculation

    % Yearly amplitude of the signal (applies to the Antarctic ice sheet)
    if islogical(lat)==1
        Ta=min(20+max(sn,0)/300.,30); % seasonal peak-to-peak amplitude in Ts
    else
        Ta=zeros(ctr.imax,ctr.jmax);
        Ta(MASK==1)=-11.06+0.003204*sn(MASK==1)-0.3584*lat(MASK==1);
        Ta(MASK==0)=-45.06-0.04598*sn(MASK==0)-0.969*lat(MASK==0);
        Ta(Ta<15)=15;
    end
    
    % Define idealized yearly cycle Tm:
    t=repmat(reshape(linspace(1,365.25,par.PDDsteps),1,1,par.PDDsteps), ...
        [ctr.imax,ctr.jmax,1]);
    Tm=repmat(Ts,[1,1,par.PDDsteps])-0.5*repmat(Ta,[1,1,par.PDDsteps]).* ...
        sin(2*pi*t/365.25);
    % PDD calculation from Calov and Greve (2005)
    nT=(Tm-par.PDDth)/(sqrt(2)*par.Tsigma);
    T=(par.Tsigma/sqrt(2*pi)*exp(-(nT.^2))+0.5*(Tm-par.PDDth).* ...
        erfc(-nT))+par.PDDth;

    % Rain fraction calculation
    nP=(Tm-par.Tsnow)/(sqrt(2)*par.Psigma);
    P=(par.Psigma/sqrt(2*pi)*exp(-(nP.^2))+0.5*(Tm-par.Tsnow).* ...
        erfc(-nP))+par.Tsnow;

    PDD=sum(T,3)*365.25/par.PDDsteps;
    R=sum(1-max(0,min((par.Train-P)/(par.Train-par.Tsnow),1)),3)/par.PDDsteps;
    
    % Melt model
    acc=max(0,Pr.*(1-R)); % Snow Accumulation Rate. R is the fraction of 
                          % precipitation that falls as rain
    rain=max(0,Pr.*R);
    ESM=min(max(acc,0),par.snowfac.*PDD); 
    % Effective snow melt (cannot exceed the amount of snow)
    PDDr=max(PDD-ESM./par.snowfac,0); % remaining PDDs after snow melt
    EIM=par.icefac.*PDDr; % Effective ice melt
    Smelt=ESM+EIM; % total surface melt

    % Refreezing model
    W_r=ESM+rain; % Available water mass - maximum amount of water that 
                  % can be retained provided enough energy is available
    P_r=max(-(min(Ts,0))*par.d_ice*par.K/(par.kdif*par.rho*par.Latent),0); 
    % Potential retention mass - maximum amount of liquid water that can 
    % be retained provided enough mass is available
    E_r=min(Pr,min(P_r,W_r)); 
    % Effective Refreezing - amount of water that is effectively 
    % retained at the end of the melting season
    % Er is limited by the annual precipitation
    runoff=EIM+W_r-E_r; % ice melt + rain and snowmelt liquid water beyond 
                        % the snowpack saturation

end


function [Melt]=PICOMelt(HB,shMASK,ShelfN,Ak,Bk,S0o,T0o,S,T, ...
    numsh,Bmax,ctr,nu,par)

% Kori-ULB
% Define sub-shelf melt underneath each ice shelf based on the results for each
% ocean box
    
    % Results for local shelf system
    Melt=zeros(ctr.imax,ctr.jmax);
    HB(shMASK==0)=0;
    
    % Box I
    g1=Ak*ctr.gammaT;
    for i=1:numsh
        Tref=min(par.lambda1*S0o(i)+par.lambda2+par.lambda3*HB-T0o(i),0); 
        % make sure Tref<=0
        dif=g1./(ctr.C*par.rhoref*(par.betao*S0o(i)/nu-par.alphao));
        x=-0.5*dif+0.5*sqrt(dif.^2-4*dif.*Tref);
        Tsh=T0o(i)-x; % T1 should be lower than T0
        Ssh=S0o(i)*(1-x/nu);
        m=-ctr.gammaT*(par.lambda1*Ssh+par.lambda2+par.lambda3*HB-Tsh)/ ...
            nu*par.secperyear;
        Melt(Bk==1 & ShelfN==i)=m(Bk==1 & ShelfN==i);
    end
    
    % Box II to nbox
    for i=1:numsh
        q=ctr.C*par.rhoref*(par.betao*(S0o(i)-S(i,1))-par.alphao*(T0o(i)-T(i,1)));
        for k=2:Bmax(i)
            Tref=par.lambda1*S(i,k-1)+par.lambda2+par.lambda3*HB-T(i,k-1);
            x=-g1.*Tref./(q+g1*(1-par.lambda1*S(i,k-1)/nu));
            Tsh=T(i,k-1)-x;
            Ssh=S(i,k-1)*(1-x/nu);
            m=-ctr.gammaT*(par.lambda1*Ssh+par.lambda2+par.lambda3*HB-Tsh)/ ...
                nu*par.secperyear;
            Melt(Bk==k & ShelfN==i)=m(Bk==k & ShelfN==i);
        end
    end
    
end


function [MASKpicop,MASKb]=PICOPmasks(MASK,MASKlk,H,ctr,par)

% Kori-ULB
% Defines the MASK for PICOP domain (ice shelf + grounding line + first
% ocean grid point) and the boundary (grounding line and first ocean
% grid point)

    MASKpicop=zeros(ctr.imax,ctr.jmax); % initialize
    MASKb=zeros(ctr.imax,ctr.jmax); % ice shelf boundary points
    MASK=MASK+MASKlk; % FP: make sure that lakes are removed
    % grounding line (=2)
    MASK1=circshift(MASK,[0 -1]);
    MASK2=circshift(MASK,[0 1]);
    MASK3=circshift(MASK,[-1 0]);
    MASK4=circshift(MASK,[1 0]);
    MASKpicop(MASK==1 & (MASK1==0 | MASK2==0 | MASK3==0 | MASK4==0))=1;
    MASKb(MASK==1 & (MASK1==0 | MASK2==0 | MASK3==0 | MASK4==0))=1;
    % ice shelves (=1)
    MASKpicop(MASK==0 & H>par.SeaIceThickness)=1;
    % adjacent floating grid cell
    MASK1=circshift(MASKpicop,[0 -1]);
    MASK2=circshift(MASKpicop,[0 1]);
    MASK3=circshift(MASKpicop,[-1 0]);
    MASK4=circshift(MASKpicop,[1 0]);
    MASKb(MASKpicop==0 & MASK==0 & (MASK1==1 | MASK2==1 | MASK3==1 | MASK4==1))=1;
    MASKpicop(MASKpicop==0 & MASK==0 & (MASK1==1 | MASK2==1 | MASK3==1 | MASK4==1))=1;
    % remove ice shelves near border of domain
%     MASKpicop(:,1:3)=0;
%     MASKpicop(:,ctr.jmax-2:ctr.jmax)=0;
%     MASKpicop(1:3,:)=0;
%     MASKpicop(ctr.imax-2:ctr.imax,:)=0;

end


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

    
function [Bk,Ak,Bmax]=PICOsetup(glMASK,H,ctr,par,ShelfN,numsh,shMASK)

% Kori-ULB
% Setup of the PICO model to determine ice shelf area
% (Ak) and the ocean boxes underneath each shelf (Bk)

    % Add first ocean point in glMASK
    MASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
    MASK2=circshift(glMASK,[0 1]); % glMASK(i,j-1)
    MASK3=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
    MASK4=circshift(glMASK,[1 0]); % glMASK(i-1,j)
    glMASK(glMASK==6 & (MASK1==5 | MASK2==5 | ...
        MASK3==5 | MASK4==5))=7; % First floating point (=7)
        
    % relative distance to GL: gMASK
    gMASK=ones(ctr.imax,ctr.jmax)-2;
    gMASK(glMASK==2)=0; % start at GL position
    gMASK(glMASK==3 | glMASK==4 | glMASK==5 | glMASK==7)=-2;
    for i=0:1000
        sg=sum(sum(gMASK==-2));
        if i>0 && (sg-sg0)==0
            break;
        end
        MASK1=circshift(gMASK,[0 -1]); % glMASK(i,j+1)
        MASK2=circshift(gMASK,[0 1]); % glMASK(i,j-1)
        MASK3=circshift(gMASK,[-1 0]); % glMASK(i+1,j)
        MASK4=circshift(gMASK,[1 0]); % glMASK(i-1,j)
        gMASK(gMASK==-2 & (MASK1==i | MASK2==i | MASK3==i | MASK4==i))=i+1;
        sg0=sg;
    end

    % relative distance from front: fMASK
    fMASK=ones(ctr.imax,ctr.jmax)-2;
    fMASK(glMASK==5)=0; % start at front position
    fMASK(glMASK==2 | glMASK==3 | glMASK==4)=-2;
    for i=0:1000
        sg=sum(sum(fMASK==-2));
        if i>0 && (sg-sg0)==0
            break;
        end
        MASK1=circshift(fMASK,[0 -1]); % glMASK(i,j+1)
        MASK2=circshift(fMASK,[0 1]); % glMASK(i,j-1)
        MASK3=circshift(fMASK,[-1 0]); % glMASK(i+1,j)
        MASK4=circshift(fMASK,[1 0]); % glMASK(i-1,j)
        fMASK(fMASK==-2 & (MASK1==i | MASK2==i | MASK3==i | MASK4==i))=i+1;
        sg0=sg;
    end

    gMASK(gMASK<0 | H<5*par.SeaIceThickness)=0;
    fMASK(fMASK<0 | H<5*par.SeaIceThickness)=0;
    gMASK=gMASK*ctr.delta;
    fMASK=fMASK*ctr.delta;

    % define boxes and relative distance r
    dmax=max(max(gMASK));
    rd=gMASK./(gMASK+fMASK); % relative distance between GL and front
    
    % define number of boxes for different ice shelves
    nB=1+round(sqrt(gMASK/dmax)*(par.nbox-1));
    nD=zeros(ctr.imax,ctr.jmax);
    for i=1:numsh
        nD(ShelfN==i)=max(nB(ShelfN==i));
    end

    Bk=zeros(ctr.imax,ctr.jmax);
    for k=1:par.nbox
        LL=1-sqrt(abs((nD-k+1)./nD));
        UL=1-sqrt(abs((nD-k)./nD));
        Bk(rd>=LL & rd<=UL)=k;
        Bk(Bk>nD)=nD(Bk>nD);
    end
    Bk(glMASK==3)=1; % make sure that boxes near GL are always Box 1
    Bmax=zeros(numsh,1);
    for i=1:numsh
        Bmax(i)=max(Bk(ShelfN==i));
    end
    % Calculate surface of each box (Ak)
    Ak=shMASK; % size of each box within particular ice shelf
    for i=1:numsh
        Ak(ShelfN==i)=Ak(ShelfN==i)*sum(sum(ShelfN==i));
    end
    Bk(Ak==1)=1;
    Ak(nD>0)=Ak(nD>0)*ctr.delta^2./nD(nD>0);
end


function PlotMainFigure(ctr,par,x,y,sn,S0,H,u,B,MASK,glMASK,LSF)

% Kori-ULB
% Plot of the main figure during the model run with changes in ice
% thickness and total velocity field of the ice sheet

    if ctr.plotH==1
        sn1=sn;
    elseif ctr.plotH==2
        sn1=H;
        sn1(H==0)=NaN;
    else
        sn1=sn-S0;
        plim=max(max(abs(sn-S0)));
        plim(plim==0)=1;
    end
    u1=u;
    if ctr.glMASKexist==1
        sn1(glMASK==6)=NaN;
        u1(glMASK==6)=NaN;
    else
        sn1(MASK==0)=NaN;
        u1(MASK==0)=NaN;
    end
    ax1=subplot(1,2,1);
    if ctr.plotH>=1
        imagescn(x,y,sn1);
        if ctr.plotH==2
            title('Ice thickness (m)');
        else
            title('Surface elevation (m)');
        end
        colormap(ax1,crameri(par.color));
        colorbar;
    else
        imagescn(x,y,sn1,[-plim plim]);
        title('Ice thickness change (m)');
        colormap(ax1,crameri(par.dcolor));
        colorbar;
    end
    if ctr.glMASKexist==1 && ctr.plotGL==1
        hold on;
        contour(x,y,MASK,1,'LineColor','k','LineWidth',0.5);
        if ctr.calving==5
            contour(x,y,LSF,[0 0],'LineColor','b','LineWidth',0.5);
        end
        hold off;
    end
    if ctr.mismip>=1
        hold on;
        contour(x,y,B,'LineColor','w','LineWidth',0.5);
        hold off;
    end
    axis xy;
    if ctr.mismip==0 || ctr.mismip==2
        axis equal;
    end
    axis tight;
    xlabel('x (km)');
    ylabel('y (km)');

    ax2=subplot(1,2,2);
    imagescn(x,y,log10(u1),[-1 3.5]);
    axis xy;
    if ctr.mismip==0 || ctr.mismip==2
        axis equal;
    end
    axis tight;
    colormap(ax2,crameri(par.color));
    colorbar;
    if ctr.glMASKexist==1 && ctr.plotGL==1
        hold on;
        contour(x,y,MASK,1,'LineColor','k','LineWidth',0.5);
        hold off;
    end
    if ctr.mismip>=1
        hold on;
        contour(x,y,B,'LineColor','w','LineWidth',0.5);
        hold off;
    end
    title('log_{10}(ice velocity) log_{10}(m a^{-1})');
    xlabel('x (km)');
    ylabel('y (km)');
    pause(0.0001);
end


function [Melt]=PlumeMelt2019(T_in, S_in, ice_draft_depth, zGL, sina, par, ctr)

% Kori-ULB
% Apply the plume parametrization.
% This function computes the basal melt based on a plume parametrization
% (see Lazeroms et al. 2018 and Lazeroms et al. 2019).
%
% Parameters
%   T_in : scalar (or array?)
%       Ambient temperature in degrees C.
%   S_in : scalar (or array?)
%       Ambient salinity in psu.
%   ice_draft_depth : scalar or array
%       Depth of the ice draft in m (depth is negative!).
%   zGL: scalar or array
%       Depth of the grounding line where the source of the plume is in m (depth is negative!).
%   alpha: scalar or array
%       Slope angle in rad (must be positive).
%   gamma: scalar
%       Effective thermal Stanton number. Can be modulated for tuning.
%   E0: scalar
%       Entrainment coefficient. Can be modulated for tuning.
%   picop: Boolean
%       Option defining which Mterm function to use.
% Returns
%   melt_rate : scalar or array
%       Melt rate in m ice per second.

    E0=par.Eo;
    gamma=ctr.gammaTplume;

    [c_rho_1,c_rho_2,c_tau]=Compute_c_rho_tau(gamma,S_in,par);

    % freezing temperature at the grounding line
    Tf=par.lambda1*S_in+par.lambda2+par.lambda3*zGL; % Ocean freezing point
    thermal_forcing=T_in-Tf;

    x_hat=ComputeXhat(ice_draft_depth,zGL,T_in,Tf,E0,c_tau,gamma,sina,par);
    M_hat=ComputeMhat(x_hat);
    Mterm=ComputeMterm(T_in,S_in,Tf,c_rho_1,c_tau,gamma,E0, ...
        thermal_forcing,sina,par);
    Melt=Mterm.*M_hat.*par.secperyear; % m ice per year

end


function [pot]=PotentialFilling(pot,ctr)

% Kori-ULB
% Hollow filling for potential gradient

    for k=1:10
        pool=zeros(ctr.imax,ctr.jmax);
        p1=circshift(pot,[0 -1]);
        p2=circshift(pot,[0 1]);
        p3=circshift(pot,[-1 0]);
        p4=circshift(pot,[1 0]);
        pool(pot<p1 & pot<p2 & pot<p3 & pot<p4)=1;
        pot(pool==1)=(p1(pool==1)+p2(pool==1)+p3(pool==1)+p4(pool==1))/4;
    end

end


function [lambda0,lambda1]=PrincipalStrain(dudx,dvdy,dudy,dvdx)

% Kori-ULB
% Calculation of principal strain based on Cauchy stress

    Exx=(2*dudx+dvdy);
    Eyy=(2*dvdy+dudx);
    Exy=0.5*(dudy+dvdx);
    b=0.5*(Exx+Eyy);
    d=sqrt((0.5*(Exx-Eyy)).^2+Exy.^2);
    lambda0=b+d;
    lambda1=b-d;

end


function As = RegularizationNew(Asor,B,stdB,vx,vy,ctr,ncor,sigma)

% Kori-ULB
% New regularization function based on Gaussian filter for the optimization
% of subglacial sliding coefficients of the grounded ice sheet

    DEM = log10(Asor);
    % sigma is the standard deviation of the Gaussian distribution
    % span defines the size of the Gaussian filter
    nwin=round(sigma*2);
    span=nwin+0.5;  

    vx1=0.5*(vx+circshift(vx,[0 1]));
    vy1=0.5*(vy+circshift(vy,[1 0]));
    v=sqrt(vx1.^2+vy1.^2);
    VelAng=atan2(vy1,vx1);

    TotWeightFac=zeros(ctr.imax,ctr.jmax);
    DEMnew=zeros(ctr.imax,ctr.jmax);

    for j=-nwin:nwin
        for i=-nwin:nwin
            zdist = sqrt(j^2 + i^2); %Distance in grid cells
            if zdist<=span
                zang=atan2(i,j);
                GaussFiltWeight=normpdf(zdist/sigma);
                VelMagWeight=min(1,max(0,log10(v/100)*zdist/sigma));
                VelWeight=1-abs(sin(zang-VelAng)).*VelMagWeight;
                BedWeight=1./2.^(0.5*(max(1,log10(abs(circshift(B,[-i -j]) ...
                    -B)))-1));
                if ctr.stdBexist==1
                    stdBweight=1./2.^(max(1,log10(circshift(stdB,[-i -j])))-1);
                else
                    stdBweight=1;
                end
                DEMnew=DEMnew+GaussFiltWeight*VelWeight.*BedWeight.* ...
                    stdBweight.*circshift(DEM,[-i -j]);
                TotWeightFac=TotWeightFac+GaussFiltWeight*VelWeight.* ...
                    BedWeight.*stdBweight; %for normalization --> sum of all Probs=1
            end
        end
    end
    As=Asor;
    As(ncor>0.5)=10.^(DEMnew(ncor>0.5)./TotWeightFac(ncor>0.5));

end


function [d,udx,udy,ud,ubx,uby,ub,uxsia,uysia,p,pxy]= ...
    SIAvelocity(ctr,par,A,Ad,Ax,Ay,Asfd,Asfx,Asfy,taud,G,Tb, ...
    H,Hm,Hmx,Hmy,gradm,gradmx,gradmy,gradxy,signx,signy,MASK,p,px,py,pxy)

% Kori-ULB
% Deformationald and basal velocity according to the shallow-ice
% approximation. Parameters for the thermomechanical coupling are also
% defined based on basal temperature gradients.

    if ctr.Tcalc==2
        p=par.n-1+par.Q1*G.*Hm./(par.K*par.R.*(h2d(Tb)+par.T0).^2);
        % on staggered d-grid
        pxy=par.n-1+par.Q1*G.*H./(par.K*par.R.*(Tb+par.T0).^2); % on h-grid     
        px=par.n-1+par.Q1*G.*Hmx./(par.K*par.R.*(0.5* ...
            (Tb+circshift(Tb,[0 -1]))+par.T0).^2); % on h-grid     
        py=par.n-1+par.Q1*G.*Hmy./(par.K*par.R.*(0.5* ...
            (Tb+circshift(Tb,[-1 0]))+par.T0).^2); % on h-grid     
    end
    if ctr.u0>1e10
        ds=1; % disable correction when u0 is default value
    else
        ds=max(0.01,min(1-Asfd./ctr.u0.*taud.^ctr.m,1)); 
        % Zoet correction factor %VL: use Asfd
    end
    d=par.Z*Ad./(p+2).*Hm.^(par.n+2.).*gradm.^((par.n-1.)/2.)+ ...
        Asfd*(par.rho*par.g)^ctr.m.*Hm.^(ctr.m+1.).* ...
        gradm.^((ctr.m-1.)/2.)./ds;
    udx=par.Z*Ax./(px+2).*signx.*Hmx.^(par.n+1.).*gradmx.^(par.n/2.);
    udy=par.Z*Ay./(py+2).*signy.*Hmy.^(par.n+1.).*gradmy.^(par.n/2.);
    ud=par.Z*A./(pxy+2).*H.^(par.n+1.).*gradxy.^(par.n/2.); 
    % only used for temperature (on h-grid)
    ubx=Asfx.*signx.*(par.rho*par.g*Hmx).^ctr.m.*gradmx.^(ctr.m/2.); 
    % only for temperature
    uby=Asfy.*signy.*(par.rho*par.g*Hmy).^ctr.m.*gradmy.^(ctr.m/2.);
    uxsia=min(max(-par.maxspeed,udx+ubx),par.maxspeed);
    uysia=min(max(-par.maxspeed,udy+uby),par.maxspeed);
    ub=vec2h(ubx,uby); % only for temperature (on h grid)
    ub(MASK==0)=0;
    
end


function [MASKmx,MASKmy,HAFmx,HAFmy]=SSAmasks(ctr,par,SLR,Bmx, ...
    Bmy,Hmx,Hmy,bMASKx,bMASKy)

% Kori-ULB
% Ice shelf masks on staggered grids for SSA computation

    HAFmx=Bmx-SLR+Hmx*par.rho/par.rhow;
    HAFmy=Bmy-SLR+Hmy*par.rho/par.rhow;
    MASKmx=zeros(ctr.imax,ctr.jmax);
    MASKmy=zeros(ctr.imax,ctr.jmax);
    if ctr.shelf==1 || ctr.SSA>=1
        if ctr.shelf==1
            MASKmx(HAFmx<0)=1; % floating
            MASKmy(HAFmy<0)=1;
        end
        if ctr.SSA>=1
            MASKmx(HAFmx>=0)=1; % grounded
            MASKmy(HAFmy>=0)=1;
        end
        if ctr.basin==1
            MASKmx(bMASKx==1)=0;
            MASKmy(bMASKy==1)=0;
        end
    end
end


function [uxssa,uyssa,beta2,eta,dudx,dudy,dvdx,dvdy,su,ubx,uby,ux,uy,damage]= ...
    SSAvelocity(ctr,par,su,Hmx,Hmy,gradmx,gradmy,signx,signy, ...
    uxssa,uyssa,H,HB,B,stdB,Asf,A,MASK,glMASK,HAF,HAFmx,HAFmy,cnt, ...
    nodeu,nodev,MASKmx,MASKmy,bMASK,uxsia,uysia,udx,udy,node,nodes, ...
    Mb,Melt,dtdx,dtdx2,VM,damage)

% Kori-ULB
% Iterative solution to the SSA velocity (both pure SSA and hybrid model

    eps=1e-8;
    taudx=par.rho*par.g*Hmx.*sqrt(gradmx).*signx;
    taudy=par.rho*par.g*Hmy.*sqrt(gradmy).*signy;
    ussa=vec2h(uxssa,uyssa);    %VL: ussa on h-grid
    if par.ShelfPinning==1 && ctr.stdBexist==1 && ctr.inverse==0 
        fg=max(0,1-(HB-B)./stdB); % subgrid pinning points in ice shelf
    else
        fg=1;
    end
    ussa=max(ussa,1e-3); % lower limit to prevent absence of sliding in ideal cases
    if ctr.u0>1e10
        beta2=fg.*(ussa.^(1/ctr.m-1)).*Asf.^(-1/ctr.m);
    else
        beta2=fg.*(ussa.^(1/ctr.m-1)).*((ussa+ctr.u0).* ...
            Asf/ctr.u0).^(-1/ctr.m);
    end
    beta2=min(beta2,1e8);
    beta2(MASK==0)=0;
    betax=0.5*(beta2+circshift(beta2,[0 -1]));
    betay=0.5*(beta2+circshift(beta2,[-1 0]));

    if ctr.mismip>=1
        betax(:,1)=betax(:,2); % symmetric divide
        betax(1,:)=betax(3,:); % symmetry axis
        betax(ctr.imax,:)=betax(ctr.imax-2,:); % periodic BC
        betax(:,ctr.jmax)=0; % ocean
        betay(:,1)=betay(:,3); % symmetric divide
        betay(1,:)=betay(2,:); % symmetry axis
        betay(ctr.imax,:)=betay(ctr.imax-1,:); % periodic BC
        if ctr.mismip==2 % Thule setup
            betax(ctr.imax,:)=0;
            betay(ctr.imax,:)=0;
        end
    end

    Hshelf=200; % Mean ice shelf thickness to control viscosity on domain edge
    Tf=.5*par.rho*par.g*Hshelf*(1.-par.rho/par.rhow); % on h-grid
    eta1=0.5*Hshelf.*Tf.^(1.-par.n)./A; 
    % Viscosity on edge (or when viscosity cannot be calculated)
    if ctr.SSA>1 % hydrid model
        udx(HAFmx<0)=0; % no deformation for floating cells
        udy(HAFmy<0)=0;
    else % pure SSA
        udx=zeros(ctr.imax,ctr.jmax);
        udy=zeros(ctr.imax,ctr.jmax);
    end
    if ctr.damage==1 && cnt>1
        dtr=TransportDamage(node,nodes,damage,Mb,Melt,H,glMASK,dtdx,dtdx2, ...
            uxssa,uyssa,ctr,cnt,bMASK,VM,par);
    end
    for ll=1:par.visciter % iteration over effective viscosity
        [eta,dudx,dvdy,dudy,dvdx]=EffVisc(A,uxssa,uyssa,H,par,eta1,MASK,ctr);
        if ctr.damage==1 && cnt>1
%             if ll==1
                damage=NyeDamage(par,ctr,dudx,dvdy,dudy,dvdx,eta,H,HAF);
                damage=min(par.damlim*H,max(damage,dtr));
                scale_eta=(H-min(damage,H-eps))./(H+eps);
%             end
        else
            scale_eta=1;
            damage=zeros(ctr.imax,ctr.jmax);
        end
        eta=eta.*scale_eta;
        if ctr.shelf==1 || ctr.schoof>0
            eta(glMASK==6)=1e7;
        end
        [uxs1,uys1,su]=SparseSolverSSA(nodeu,nodev,su,MASKmx,MASKmy,bMASK, ...
            H,eta,betax,betay,uxssa,uyssa,uxsia,uysia,udx,udy,taudx, ...
            taudy,ctr,par);           
        duxs=sqrt((uxs1-uxssa).^2+(uys1-uyssa).^2);
        duxs(isnan(duxs))=0;
        uxssa=uxs1;
        uyssa=uys1;
        limit=sum(duxs(:))/(ctr.imax*ctr.jmax);
        %---------iterative beta---------
        if cnt<=ctr.BetaIter
            ussa=vec2h(uxssa,uyssa); %VL: ussa on h-grid
            if ctr.u0>1e10
                beta2=fg.*(ussa.^(1/ctr.m-1)).*Asf.^(-1/ctr.m);
            else
                beta2=fg.*(ussa.^(1/ctr.m-1)).*((ussa+ctr.u0).*Asf ...
                    /ctr.u0).^(-1/ctr.m);
            end
            beta2=min(beta2,1e8);
            beta2(MASK==0)=0;
            betax=0.5*(beta2+circshift(beta2,[0 -1]));
            betay=0.5*(beta2+circshift(beta2,[-1 0]));
        end
        %--------------------------------
        if limit<par.visctol % Limit on convergence
            break;
        end
    end
    uxssa(Hmx==0 & HAFmx>0)=0; % only for grounded ice sheet with H=0
    uyssa(Hmy==0 & HAFmy>0)=0;
    uxssa=min(max(-par.maxspeed,uxssa),par.maxspeed);
    uyssa=min(max(-par.maxspeed,uyssa),par.maxspeed);
    ux=uxssa;   %LZ2021
    uy=uyssa;   %LZ2021
    if ctr.SSAdiffus==2 % SSA as basal velocity only
        ux=ux+udx;
        uy=uy+udy;
    end
    if ctr.SSA==1 % basal sliding velocity SSA
        ubx=ux;
        uby=uy;
    else
        ubx=ux-udx;
        uby=uy-udy;
    end
end


function [VAF,POV,SLC,VA0]=SeaLevel(ctr,par,SLR,B,H,H0,VAF0,POV0)

% Kori-ULB
% Calculation of local sea level

    VAFi=max(0,H+min(B-SLR,0)*(par.rhow/par.rho));
    % VAF variation (SL equivalent in ocean water)
    VAF=sum(sum((VAF0-VAFi)*ctr.delta^2.))*par.rho/(par.Aoc*par.rhow);
    VA0i=max(0,H+min(B-par.SLref,0)*(par.rhow/par.rho));
    % VAreferenceSL variation (SL equivalent in ocean water)
    VA0=sum(sum((VAF0-VA0i)*ctr.delta^2.))*par.rho/(par.Aoc*par.rhow);
    POVi=max(0,par.SLref-B);
    % Potential Ocean volume variation (SL equivalent in ocean water)
    POV=sum(sum((POV0-POVi)*ctr.delta^2.))/par.Aoc;
    % Density correction for transformation from ice to freshwater and 
    % not ocean water (SL equivalent)
    DENScorr=sum(sum((H0-H)*ctr.delta^2.))*((par.rho/par.rhow)- ...
        (par.rho/par.rhof))/par.Aoc;
    % Sea-level contribution from modelled ice-sheet
    SLC=VA0+POV-DENScorr;

end


function [MASKHole,HoleCT]=ShelfHole(MASK,H,par)

% Kori-ULB
% Checking for eventual holes in ice shelves and returning a MASK that
% identifies these

    MASK(MASK==0 & H>par.SeaIceThickness)=3;
    [~,MASKHole]=bwboundaries(MASK>0,4);
    MASKHole(MASK==1 | H>par.SeaIceThickness)=0;
    MASKHole(MASKHole~=0)=1;
    H1=circshift(MASKHole,[-1 0]); % i+1,j
    H2=circshift(MASKHole,[1 0]); % i-1,j
    H3=circshift(MASKHole,[0 -1]); % i,j+1
    H4=circshift(MASKHole,[0 1]); % i,j-1
    HoleCT=H1==1 | H2==1 | H3==1 | H4==1;
end


function [ShelfN,numsh,shMASK,MASKlk]=ShelfSetup(MASK,glMASK,H,ctr,par)

% Kori-ULB
% Delimit ice shelves and number them separately (ShelfN). Also returns a
% shelf MASK and identifies lakes (locally floating conditions within
% grounded ice sheet

    %LZ: find all 'lakes'
    [~,MASKlk]=bwboundaries(MASK,4);
    MASKlk(MASK==1)=0;
    MASKlk(MASKlk~=0)=1;

    shMASK=zeros(ctr.imax,ctr.jmax); % Mask of ice shelf
    shMASK(MASKlk==0 & (glMASK==3 | glMASK==4 | glMASK==5))=1;
    shMASK(H<5*par.SeaIceThickness)=0;
    
    %LZ: new function for calculation of ShelfN
    ShelfN=bwlabel(shMASK);
    numsh=max(ShelfN(:));

end


function bload=SparseSolverBedrock(Db,par,ctr,VM,node,nodes,loadB)

% Kori-ULB
% Solving thin-plate equation with spatially-varying flexural rigidity for
% isostatic adjustment

    Db1=circshift(Db,[0 -1]); % d(i,j+1)
    Db2=circshift(Db,[0  1]); % d(i,j-1)
    Db3=circshift(Db,[-1 0]); % d(i+1,j)
    Db4=circshift(Db,[1 0]); % d(i-1,j)
    Db5=circshift(Db,[-1 -1]); % d(i+1,j+1)
    Db6=circshift(Db,[1  -1]); % d(i-1,j+1)
    Db7=circshift(Db,[-1 1]); % d(i+1,j-1)
    Db8=circshift(Db,[1 1]); % d(i-1,j-1)

    MASKb=zeros(ctr.imax,ctr.jmax);
    MASKb(1:2,:)=1;
    MASKb(ctr.imax-1:ctr.imax,:)=1;
    MASKb(:,1:2)=1;
    MASKb(:,ctr.jmax-1:ctr.jmax)=1;

    nabla2Db = (-4*Db+Db1+Db2+Db3+Db4);
    dDbx = (Db1-Db2);
    dDby = (Db3-Db4);
    dDbx2 = (-2*Db+Db1+Db2);
    dDbxy = (Db5-Db6-Db7+Db8)/4;
    dDby2 = (-2*Db+Db3+Db4);

    V0 = (20*Db-4*nabla2Db-(1-par.nuB)*(-2*dDbx2-2*dDby2))/ctr.delta^4+ ...
        par.rhom*par.g;
    V1 = (-8*Db-dDbx-dDbx+nabla2Db-(1-par.nuB)*dDby2)/ctr.delta^4;
    V2 = (-8*Db+dDbx+dDbx+nabla2Db-(1-par.nuB)*dDby2)/ctr.delta^4;
    V3 = (-8*Db-dDby-dDby+nabla2Db-(1-par.nuB)*dDbx2)/ctr.delta^4; 
    V4 = (-8*Db+dDby+dDby+nabla2Db-(1-par.nuB)*dDbx2)/ctr.delta^4;         
    V5 = (2*Db+0.5*dDbx+0.5*dDby-0.5*(1-par.nuB)*(-dDbxy))/ctr.delta^4;
    V6 = (2*Db+0.5*dDbx-0.5*dDby-0.5*(1-par.nuB)*(dDbxy))/ctr.delta^4;
    V7 = (2*Db-0.5*dDbx+0.5*dDby-0.5*(1-par.nuB)*(dDbxy))/ctr.delta^4;
    V8 = (2*Db-0.5*dDbx-0.5*dDby-0.5*(1-par.nuB)*(-dDbxy))/ctr.delta^4;
    V9 = (Db+0.5*dDbx)/ctr.delta^4;
    V10 = (Db-0.5*dDbx)/ctr.delta^4;
    V11 = (Db+0.5*dDby)/ctr.delta^4;
    V12 = (Db-0.5*dDby)/ctr.delta^4;

    R0 = loadB;

    V0(MASKb==1)=1;
    V1(MASKb==1)=0;
    V2(MASKb==1)=0;
    V3(MASKb==1)=0;
    V4(MASKb==1)=0;
    V5(MASKb==1)=0;
    V6(MASKb==1)=0;
    V7(MASKb==1)=0;
    V8(MASKb==1)=0;
    V9(MASKb==1)=0;
    V9(MASKb==1)=0;
    V10(MASKb==1)=0;
    V11(MASKb==1)=0;
    V12(MASKb==1)=0;
    R0(MASKb==1)=0;

    V=[reshape(V0(VM==1),nodes,1)
        V1(V1~=0)
        V2(V2~=0)
        V3(V3~=0)
        V4(V4~=0)
        V5(V5~=0)
        V6(V6~=0)
        V7(V7~=0)
        V8(V8~=0)
        V9(V9~=0)
        V10(V10~=0)
        V11(V11~=0)
        V12(V12~=0)
        ];

    row=[reshape(node(VM==1),nodes,1)
        node(V1~=0)
        node(V2~=0)
        node(V3~=0)
        node(V4~=0)
        node(V5~=0)
        node(V6~=0)
        node(V7~=0)
        node(V8~=0)
        node(V9~=0)
        node(V10~=0)
        node(V11~=0)
        node(V12~=0)
        ];

    nodeV1=circshift(node,[0 -1]); % i,j+1
    nodeV2=circshift(node,[0 1]); % i,j-1
    nodeV3=circshift(node,[-1 0]); % i+1,j
    nodeV4=circshift(node,[1 0]); % i-1,j
    nodeV5=circshift(node,[-1 -1]); % i+1,j+1
    nodeV6=circshift(node,[1 -1]); % i-1,j+1
    nodeV7=circshift(node,[-1 1]); % i+1,j-1
    nodeV8=circshift(node,[1 1]); % i-1,j-1
    nodeV9=circshift(node,[0 -2]); % i,j+2
    nodeV10=circshift(node,[0 2]); % i,j-2
    nodeV11=circshift(node,[-2 0]); % i+2,j
    nodeV12=circshift(node,[2 0]); % i-2,j

    col=[reshape(node(VM==1),nodes,1)
        nodeV1(V1~=0)
        nodeV2(V2~=0)
        nodeV3(V3~=0)
        nodeV4(V4~=0)
        nodeV5(V5~=0)
        nodeV6(V6~=0)
        nodeV7(V7~=0)
        nodeV8(V8~=0)
        nodeV9(V9~=0)
        nodeV10(V10~=0)
        nodeV11(V11~=0)
        nodeV12(V12~=0)
        ];

    R=reshape(R0(VM==1),nodes,1);

    % construct sparse matrix
    A=sparse(row,col,V);

    % solve
    s=A\R;
    bload = zeros(ctr.imax,ctr.jmax);
    bload(node>0)=s(node(node>0));

end


function [H]=SparseSolverIceThickness(node,nodes,Mb,H,B,SLR,MASK,dtdx,dtdx2, ...
    d,u,v,ctr,cnt,bMASK,VM,par)

% Kori-ULB
% Sparse solver of the ice thickness equation

    d=min(d,1e20); % limit due to NaNs in initialization for high resolution

    d1=circshift(d,[0 1]); % d(i,j-1)
    d2=circshift(d,[1 0]); % d(i-1,j)
    d3=circshift(d,[1 1]); % d(i-1,j-1)
    dMASK=zeros(ctr.imax,ctr.jmax); %dMASK (floating=0, grounded=1)
    dMASK(MASK==1 | MASK==2)=1;
    if ctr.shelf==0
        MASK=dMASK;
    end
    um1=circshift(u,[0 1]); % u(i,j-1)
    vm1=circshift(v,[1 0]); % v(i-1,j)
    if ctr.upstream==1
        up1=circshift(u,[0 -1]); % u(i,j+1)
        vp1=circshift(v,[-1 0]); % v(i+1,j)
        um2=circshift(u,[0 2]); % u(i,j-2)
        vm2=circshift(v,[2 0]); % v(i-2,j)
    end

    dipx=dtdx*(d+d2);
    dimx=dtdx*(d1+d3);
    dipy=dtdx*(d+d1);
    dimy=dtdx*(d2+d3);
    MASKfac1=dMASK+(1-dMASK)*(1-par.rho/par.rhow);
    MASKfac2=B.*dMASK+SLR.*(1-dMASK);

    if ctr.upstream==1
        % conditions for diffusion scheme (init)
        V0=2*(d+d1+d2+d3).*MASKfac1*dtdx; % i,j
        V1=-dipx.*circshift(MASKfac1,[0 -1]); % i,j+1
        V2=-dimx.*circshift(MASKfac1,[0 1]); % i,j-1
        V3=-dipy.*circshift(MASKfac1,[-1 0]); % i+1,j
        V4=-dimy.*circshift(MASKfac1,[1 0]); % i-1,j
        V5=zeros(ctr.imax,ctr.jmax); % i,j-2
        V6=zeros(ctr.imax,ctr.jmax); % i-2,j
        V7=zeros(ctr.imax,ctr.jmax); % i,j+2
        V8=zeros(ctr.imax,ctr.jmax); % i+2,j

        % Velocity sign masks
        MU=zeros(ctr.imax,ctr.jmax);
        MU(u>=0 & um1>=0 & um2>=0)=1;
        MU(u<=0 & um1<=0 & up1<=0)=2;
        MV=zeros(ctr.imax,ctr.jmax);
        MV(v>=0 & vm1>=0 & vm2>=0)=1;
        MV(v<=0 & vm1<=0 & vp1<=0)=2;

        if ctr.basin==1
            MU(bMASK==1)=0;
            MV(bMASK==1)=0;
        end

        V0a=zeros(ctr.imax,ctr.jmax);
        V1a=zeros(ctr.imax,ctr.jmax);
        V2a=zeros(ctr.imax,ctr.jmax);
        V3a=zeros(ctr.imax,ctr.jmax);
        V4a=zeros(ctr.imax,ctr.jmax);

        % conditions for MU=0 (central difference)
        V0a(MU==0)=u(MU==0)-um1(MU==0); % i,j
        V1a(MU==0)=u(MU==0); % i,j+1
        V2a(MU==0)=-um1(MU==0); % i,j-1

        % conditions for MU=1 (grad(u)>0)
        V0a(MU==1)=2*u(MU==1)+um1(MU==1); % i,j
        V2a(MU==1)=-3*um1(MU==1)-um2(MU==1); % i,j-1
        V5(MU==1)=um2(MU==1); % i,j-2

        % conditions for MU=2 and (grad(u)<0)
        V0a(MU==2)=-u(MU==2)-2*um1(MU==2); % i,j
        V1a(MU==2)=3*u(MU==2)+up1(MU==2); % i,j+1
        V7(MU==2)=(-up1(MU==2)); % i,j+2

        % conditions for MV=0 (central difference)
        V0a(MV==0)=V0a(MV==0)+v(MV==0)-vm1(MV==0); % i,j
        V3a(MV==0)=v(MV==0); % i+1,j
        V4a(MV==0)=-vm1(MV==0); % i-1,j

        % conditions for MV=1 (grad(v)>0)
        V0a(MV==1)=V0a(MV==1)+2*v(MV==1)+vm1(MV==1); % i,j
        V4a(MV==1)=-3*vm1(MV==1)-vm2(MV==1); % i-1,j
        V6(MV==1)=vm2(MV==1); % i-2,j

        % conditions for MV=2 (grad(v)<0)
        V0a(MV==2)=V0a(MV==2)-v(MV==2)-2*vm1(MV==2); % i,j
        V3a(MV==2)=3*v(MV==2)+vp1(MV==2); % i+1,j
        V8(MV==2)=-vp1(MV==2); % i+2,j

        % Filling V-matrix
        V0=V0+V0a*dtdx2;%.*(1.-alfa);
        V1=V1+V1a*dtdx2;%.*(1.-alfa);
        V2=V2+V2a*dtdx2;%.*(1.-alfa);
        V3=V3+V3a*dtdx2;%.*(1.-alfa);
        V4=V4+V4a*dtdx2;%.*(1.-alfa);
        V5=V5*dtdx2;%.*(1.-alfa);
        V6=V6*dtdx2;%.*(1.-alfa);
        V7=V7*dtdx2;%.*(1.-alfa);
        V8=V8*dtdx2;%.*(1.-alfa);
    else
        V0=2*(d+d1+d2+d3).*MASKfac1*dtdx+dtdx2*(u-um1+v-vm1); % i,j
        V1=-dipx.*circshift(MASKfac1,[0 -1])+dtdx2*u; % i,j+1
        V2=-dimx.*circshift(MASKfac1,[0 1])-dtdx2*um1; % i,j-1
        V3=-dipy.*circshift(MASKfac1,[-1 0])+dtdx2*v; % i+1,j
        V4=-dimy.*circshift(MASKfac1,[1 0])-dtdx2*vm1; % i-1,j
    end

    R0=Mb*ctr.dt+H-H*(1-par.omega).*V0-circshift(H,[0 -1])*(1-par.omega) ...
        .*V1-circshift(H,[0 1])*(1-par.omega).*V2-circshift(H,[-1 0])* ...
        (1-par.omega).*V3-circshift(H,[1 0])*(1-par.omega).*V4- ...
        MASKfac2.*(d+d1+d2+d3)*2*dtdx+circshift(MASKfac2,[0 -1]).*dipx+ ...
        circshift(MASKfac2,[0 1]).*dimx+circshift(MASKfac2,[-1 0]).*dipy+ ...
        circshift(MASKfac2,[1 0]).*dimy;
    if ctr.upstream==1
        R0=R0-circshift(H,[0 2])*(1-par.omega).*V5-circshift(H,[2 0])* ...
            (1-par.omega).*V6-circshift(H,[0 -2])*(1-par.omega).*V7- ...
            circshift(H,[-2 0])*(1-par.omega).*V8;
    end

    V0(MASK==0)=0; % note that for shelf=1, MASK=glMASK in the call
    V1(MASK==0)=0;
    V2(MASK==0)=0;
    V3(MASK==0)=0;
    V4(MASK==0)=0;
    if ctr.upstream==1
        V5(MASK==0)=0;
        V6(MASK==0)=0;
        V7(MASK==0)=0;
        V8(MASK==0)=0;
    end
    R0(MASK==0)=H(MASK==0);

    % boundaries
    V9=zeros(ctr.imax,ctr.jmax); % ice divide or ocean
    V10=zeros(ctr.imax,ctr.jmax); % i=1 periodic boundary or ocean
    V11=zeros(ctr.imax,ctr.jmax); % i=imax periodic boundary or ocean
    V12=zeros(ctr.imax,ctr.jmax); % j=jmax ocean contact

    wholemask=ctr.imax*ctr.jmax-sum(MASK(:));
    if wholemask~=0 % only when domain is not MASK=1 everywhere
        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(:,1)=1; % symmetric divide or ocean
        V0(MASKb==1)=0;
        V9(MASKb==1)=-1;
        R0(MASKb==1)=0;

        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(1,2:ctr.jmax-1)=1; % periodic BC at i=1 or ocean
        V0(MASKb==1)=0;
        V10(MASKb==1)=-1;
        R0(MASKb==1)=0;

        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(ctr.imax,2:ctr.jmax-1)=1; % periodic BC at i=imax or ocean
        V0(MASKb==1)=0;
        V11(MASKb==1)=-1;
        R0(MASKb==1)=0;

        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(1:ctr.imax,ctr.jmax)=1; % ocean
        V0(MASKb==1)=0;
        V12(MASKb==1)=-1;
        R0(MASKb==1)=0;
    end

    MASKb=zeros(ctr.imax,ctr.jmax);
    MASKb(1,:)=1;
    MASKb(ctr.imax,:)=1;
    MASKb(:,1)=1;
    MASKb(:,ctr.jmax)=1;
    V1(MASKb==1)=0;
    V2(MASKb==1)=0;
    V3(MASKb==1)=0;
    V4(MASKb==1)=0;
    if ctr.upstream==1
        V5(MASKb==1)=0;
        V6(MASKb==1)=0;
        V7(MASKb==1)=0;
        V8(MASKb==1)=0;
    end
    if wholemask==0
        V0(MASKb==1)=0;
        R0(MASKb==1)=0;
    end

    if ctr.upstream==1
        V=[reshape(V0(VM==1)*par.omega+1,nodes,1)
            V1(V1~=0)*par.omega
            V2(V2~=0)*par.omega
            V3(V3~=0)*par.omega
            V4(V4~=0)*par.omega
            V5(V5~=0)*par.omega
            V6(V6~=0)*par.omega
            V7(V7~=0)*par.omega
            V8(V8~=0)*par.omega
            V9(V9~=0)
            V10(V10~=0)
            V11(V11~=0)
            V12(V12~=0)];

        row=[reshape(node(VM==1),nodes,1)
            node(V1~=0)
            node(V2~=0)
            node(V3~=0)
            node(V4~=0)
            node(V5~=0)
            node(V6~=0)
            node(V7~=0)
            node(V8~=0)
            node(V9~=0)
            node(V10~=0)
            node(V11~=0)
            node(V12~=0)];
    else
        V=[reshape(V0(VM==1)*par.omega+1,nodes,1)
            V1(V1~=0)*par.omega
            V2(V2~=0)*par.omega
            V3(V3~=0)*par.omega
            V4(V4~=0)*par.omega
            V9(V9~=0)
            V10(V10~=0)
            V11(V11~=0)
            V12(V12~=0)];

        row=[reshape(node(VM==1),nodes,1)
            node(V1~=0)
            node(V2~=0)
            node(V3~=0)
            node(V4~=0)
            node(V9~=0)
            node(V10~=0)
            node(V11~=0)
            node(V12~=0)];
    end

    nodeV1=circshift(node,[0 -1]); % i,j+1
    nodeV2=circshift(node,[0 1]); % i,j-1
    nodeV3=circshift(node,[-1 0]); % i+1,j
    nodeV4=circshift(node,[1 0]); % i-1,j
    if ctr.upstream==1
        nodeV5=circshift(node,[0 2]); % i,j-2
        nodeV6=circshift(node,[2 0]); % i-2,j
        nodeV7=circshift(node,[0 -2]); % i,j+2
        nodeV8=circshift(node,[-2 0]); % i+2,j
    end
    if ctr.mismip>=1
        nodeV9=circshift(node,[0 -2]); % i,j+2 - divide
        nodeV10=circshift(node,[-2 0]); % 3,j - symmetry at i=1
        if ctr.mismip==1
            nodeV11=circshift(node,[2 0]); % n-2,j - PBC at i=imax
        else
            nodeV11=circshift(node,[1 0]);
        end
        nodeV12=circshift(node,[0 1]); % i,jmax-1 - ocean
    else
        nodeV9=circshift(node,[0 -1]); % i,2 - ocean
        nodeV10=circshift(node,[-1 0]); % 2,j - ocean
        nodeV11=circshift(node,[1 0]); % imax-1,j - ocean
        nodeV12=circshift(node,[0 1]); % i,jmax-1 - ocean
    end

    if ctr.upstream==1
        col=[reshape(node(VM==1),nodes,1)
            nodeV1(V1~=0)
            nodeV2(V2~=0)
            nodeV3(V3~=0)
            nodeV4(V4~=0)
            nodeV5(V5~=0)
            nodeV6(V6~=0)
            nodeV7(V7~=0)
            nodeV8(V8~=0)
            nodeV9(V9~=0)
            nodeV10(V10~=0)
            nodeV11(V11~=0)
            nodeV12(V12~=0)];
    else
        col=[reshape(node(VM==1),nodes,1)
            nodeV1(V1~=0)
            nodeV2(V2~=0)
            nodeV3(V3~=0)
            nodeV4(V4~=0)
            nodeV9(V9~=0)
            nodeV10(V10~=0)
            nodeV11(V11~=0)
            nodeV12(V12~=0)];
    end

    R=reshape(R0(VM==1),nodes,1);

    % construct sparse matrix
    A=sparse(row,col,V);
    % Cholesky factor and solve
    if ctr.inverse==1 || ctr.ItSolv==0
        s=A\R;
    else
        D=diag(diag(A));
        C1=tril(A);
        C2=D\triu(A);
        [s,flag]=pcg(A,R,par.Htol,par.Hiter,C1,C2);
        if flag>0 || cnt==1
            s=A\R;
        end
    end

    H(node>0)=s(node(node>0));

end


function [u,v,s]=SparseSolverSSA(nodeu,nodev,s0,MASKmx,MASKmy,bMASK, ...
    H,eta,betax,betay,u,v,usia,vsia,udx,udy,taudx,taudy,ctr,par)

% Kori-ULB
% Solving the SSA equation (both pure and hybrid SSA)
% sparse matrix solution for two two-dimensional ice shelf velocity field
% with kinematic boundary conditions
% u-field on u-grid, viscosity on h-grid
% u-velocities: quadrant U and Uv of solution matrix
% V-velocities: quadrant V and Vu of solution matrix
% Subsequent interleaving of u and v velocities (Quiquet et al., 2018) to
% improve stability and speed of the algorithm

    limit=1e-5; % limit on effective viscosity gradients

    if ctr.SSAdiffus==2
        udx=zeros(ctr.imax,ctr.jmax);
        udy=zeros(ctr.imax,ctr.jmax);
    end

    R0=zeros(ctr.imax,ctr.jmax);
    U0=zeros(ctr.imax,ctr.jmax); % u(i,j)
    U1=zeros(ctr.imax,ctr.jmax); % u(i,j+1)
    U2=zeros(ctr.imax,ctr.jmax); % u(i,j-1)
    U3=zeros(ctr.imax,ctr.jmax); % u(i+1,j)
    U4=zeros(ctr.imax,ctr.jmax); % u(i-1,j)
    U5=zeros(ctr.imax,ctr.jmax); % u(i+1,j+1)
    U6=zeros(ctr.imax,ctr.jmax); % u(i+1,j-1)
    U7=zeros(ctr.imax,ctr.jmax); % u(i-1,j+1)
    U8=zeros(ctr.imax,ctr.jmax); % u(i-1,j-1)
    U9=zeros(ctr.imax,ctr.jmax); % periodic BC i=1
    U10=zeros(ctr.imax,ctr.jmax); % periodic BC i=imax
    Uv0=zeros(ctr.imax,ctr.jmax); % v(i,j)
    Uv1=zeros(ctr.imax,ctr.jmax); % v(i,j+1)
    Uv2=zeros(ctr.imax,ctr.jmax); % v(i-1,j)
    Uv3=zeros(ctr.imax,ctr.jmax); % v(i-1,j+1)

    eta1=circshift(eta,[0 -1]); % eta(i,j+1)
    H1=circshift(H,[0 -1]); % H(i,j+1)

    eta2=circshift(eta,[-1 0]); % eta(i+1,j)
    eta3=circshift(eta,[-1 -1]); % eta(i+1,j+1)
    eta4=circshift(eta,[1 0]); % eta(i-1,j)
    eta5=circshift(eta,[1 -1]); % eta(i-1,j+1)
    dmudx=(eta1-eta)/ctr.delta;
    dmudy=0.25*(eta2+eta3-eta4-eta5)/ctr.delta;

    dmudx=min(limit,max(dmudx,-limit));
    dmudy=min(limit,max(dmudy,-limit));

    MASKb=zeros(ctr.imax,ctr.jmax); % domain mask for SSA
    MASKb(2:ctr.imax-1,2:ctr.jmax-2)=1;
    MASKb(MASKmx==0)=0;

    U0(MASKb==1)=-5.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)- ...
        betax(MASKb==1);
    U1(MASKb==1)=2.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)+ ...
        2.*dmudx(MASKb==1)/ctr.delta;
    U2(MASKb==1)=2.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)- ...
        2.*dmudx(MASKb==1)/ctr.delta;
    U3(MASKb==1)=0.5*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)+ ...
        0.5*dmudy(MASKb==1)/ctr.delta;
    U4(MASKb==1)=0.5*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)- ...
        0.5*dmudy(MASKb==1)/ctr.delta;
    Uv0(MASKb==1)=-1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)+ ...
        (dmudx(MASKb==1)-0.5*dmudy(MASKb==1))/ctr.delta;
    Uv1(MASKb==1)=1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)+ ...
        (dmudx(MASKb==1)+0.5*dmudy(MASKb==1))/ctr.delta;
    Uv2(MASKb==1)=1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)- ...
        (dmudx(MASKb==1)+0.5*dmudy(MASKb==1))/ctr.delta;
    Uv3(MASKb==1)=-1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)- ...
        (dmudx(MASKb==1)-0.5*dmudy(MASKb==1))/ctr.delta;
    R0(MASKb==1)=-taudx(MASKb==1)-betax(MASKb==1).*udx(MASKb==1);

    MASKb=zeros(ctr.imax,ctr.jmax); % domain mask for non-SSA
    MASKb(2:ctr.imax-1,2:ctr.jmax-2)=1;
    MASKb(MASKmx>0)=0;

    U0(MASKb==1)=1;
    R0(MASKb==1)=usia(MASKb==1);

    % boundary conditions

    if ctr.shelf==1 && ctr.mismip==0
        % j=1; contact with ocean (upwinding in x)
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(2:ctr.imax-1,1)=1;
        if ctr.basin==1
            U0(MASKb==1)=1;
            R0(MASKb==1)=u(MASKb==1);
            MASKb(bMASK==1)=0;
        end
        U0(MASKb==1)=-4.*eta1(MASKb==1)/ctr.delta;
        U1(MASKb==1)=4.*eta1(MASKb==1)/ctr.delta;
        Uv1(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
        Uv3(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;
        R0(MASKb==1)=0.5*par.rho*par.g*H1(MASKb==1).^2.*(1.- ...
            par.rho/par.rhow);

        % j=jmax-1; contact with ocean
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(2:ctr.imax-1,ctr.jmax-1)=1;
        if ctr.basin==1
            U0(MASKb==1)=1;
            R0(MASKb==1)=u(MASKb==1);
            MASKb(bMASK==1)=0;
        end
        U0(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
        U2(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta;
        Uv0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
        Uv2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;
        R0(MASKb==1)=0.5*par.rho*par.g*H(MASKb==1).^2.*(1.- ...
            par.rho/par.rhow);

        % i=1; contact with ocean (v-direction)
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(1,2:ctr.jmax-2)=1;
        if ctr.basin==1
            U0(MASKb==1)=1;
            R0(MASKb==1)=u(MASKb==1);
            MASKb(bMASK==1)=0;
        end
        U0(MASKb==1)=-1/ctr.delta;
        U3(MASKb==1)=1./ctr.delta;
        Uv0(MASKb==1)=-1./ctr.delta;
        Uv1(MASKb==1)=1./ctr.delta;
        R0(MASKb==1)=0;

        % i=imax; contact with ocean (v-direction)
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(ctr.imax,2:ctr.jmax-2)=1;
        if ctr.basin==1
            U0(MASKb==1)=1;
            R0(MASKb==1)=u(MASKb==1);
            MASKb(bMASK==1)=0;
        end
        U0(MASKb==1)=1./ctr.delta;
        U4(MASKb==1)=-1./ctr.delta;
        Uv2(MASKb==1)=-1./ctr.delta;
        Uv3(MASKb==1)=1./ctr.delta;
        R0(MASKb==1)=0;

        % j=jmax;
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(:,ctr.jmax)=1;
        U0(MASKb==1)=1;
        R0(MASKb==1)=NaN;

        % Model corners
        U0(1,1)=1;
        U1(1,1)=-1;
        U3(1,1)=-1;
        U5(1,1)=1;
        R0(1,1)=0;
        U0(1,ctr.jmax-1)=1;
        U2(1,ctr.jmax-1)=-1;
        U3(1,ctr.jmax-1)=-1;
        U6(1,ctr.jmax-1)=1;
        R0(1,ctr.jmax-1)=0;
        U0(ctr.imax,1)=1;
        U1(ctr.imax,1)=-1;
        U4(ctr.imax,1)=-1;
        U7(ctr.imax,1)=1;
        R0(ctr.imax,1)=0;
        U0(ctr.imax,ctr.jmax-1)=1;
        U2(ctr.imax,ctr.jmax-1)=-1;
        U4(ctr.imax,ctr.jmax-1)=-1;
        U8(ctr.imax,ctr.jmax-1)=1;
        R0(ctr.imax,ctr.jmax-1)=0;

    elseif ctr.shelf==1 && ctr.mismip>=1

        % j=1: ice divide (symmetric)
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(1:ctr.imax,1)=1;
        U0(MASKb==1)=1;
        U1(MASKb==1)=1;
        R0(MASKb==1)=0;

        % j=jmax-1; contact with ocean
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(2:ctr.imax-1,ctr.jmax-1)=1;
        U0(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
        U2(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta;
        Uv0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
        Uv2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;
        R0(MASKb==1)=0.5*par.rho*par.g*H(MASKb==1).^2.*(1.-par.rho/par.rhow);

        % i=1: periodic boundary condition
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(1,2:ctr.jmax)=1;
        U0(MASKb==1)=1;
        U9(MASKb==1)=-1;
        R0(MASKb==1)=0;

        if ctr.mismip==1
            % i=imax: periodic boundary condition
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(ctr.imax,2:ctr.jmax)=1;
            U0(MASKb==1)=1;
            U10(MASKb==1)=-1;
            R0(MASKb==1)=0;
        else
            % i=imax; contact with ocean (v-direction)
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(ctr.imax,2:ctr.jmax-2)=1;
            U0(MASKb==1)=1./ctr.delta;
            U4(MASKb==1)=-1./ctr.delta;
            Uv2(MASKb==1)=-1./ctr.delta;
            Uv3(MASKb==1)=1./ctr.delta;
            R0(MASKb==1)=0;
            U0(ctr.imax,ctr.jmax-1)=1;
            U2(ctr.imax,ctr.jmax-1)=-1;
            U4(ctr.imax,ctr.jmax-1)=-1;
            U8(ctr.imax,ctr.jmax-1)=1;
            R0(ctr.imax,ctr.jmax-1)=0;
        end

        % j=jmax;
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(:,ctr.jmax)=1;
        U0(MASKb==1)=1;
        R0(MASKb==1)=NaN;

    else
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(1,1:ctr.jmax-1)=1;
        MASKb(ctr.imax,1:ctr.jmax-1)=1;
        MASKb(:,1)=1;
        MASKb(:,ctr.jmax-1)=1;
        U0(MASKb==1)=1;
        R0(MASKb==1)=0;

        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(:,ctr.jmax)=1;
        U0(MASKb==1)=1;
        R0(MASKb==1)=NaN;
    end

    % v-velocities: quadrant V and Vu of solution matrix

    S0=zeros(ctr.imax,ctr.jmax);
    V0=zeros(ctr.imax,ctr.jmax); % v(i,j)
    V1=zeros(ctr.imax,ctr.jmax); % v(i+1,j)
    V2=zeros(ctr.imax,ctr.jmax); % v(i-1,j)
    V3=zeros(ctr.imax,ctr.jmax); % v(i,j+1)
    V4=zeros(ctr.imax,ctr.jmax); % v(i,j-1)
    V5=zeros(ctr.imax,ctr.jmax); % v(i+1,j+1)
    V6=zeros(ctr.imax,ctr.jmax); % v(i-1,j+1)
    V7=zeros(ctr.imax,ctr.jmax); % v(i+1,j-1)
    V8=zeros(ctr.imax,ctr.jmax); % v(i-1,j-1)
    V9=zeros(ctr.imax,ctr.jmax); % Periodic BC on i=1
    V10=zeros(ctr.imax,ctr.jmax); % Periodic BC on i=imax-1
    V11=zeros(ctr.imax,ctr.jmax); % Symmetric ice divide v(i,j+3)
    Vu0=zeros(ctr.imax,ctr.jmax); % u(i,j)
    Vu1=zeros(ctr.imax,ctr.jmax); % u(i+1,j)
    Vu2=zeros(ctr.imax,ctr.jmax); % u(i,j-1)
    Vu3=zeros(ctr.imax,ctr.jmax); % u(i+1,j-1)

    eta1=circshift(eta,[-1 0]); % eta(i+1,j)
    H1=circshift(H,[-1 0]); % H(i+1,j)

    eta2=circshift(eta,[0 -1]); % eta(i,j+1)
    eta3=circshift(eta,[-1 -1]); % eta(i+1,j+1)
    eta4=circshift(eta,[0 1]); % eta(i,j-1)
    eta5=circshift(eta,[-1 1]); % eta(i+1,j-1)
    dmudy=(eta1-eta)/ctr.delta;
    dmudx=0.25*(eta2+eta3-eta4-eta5)/ctr.delta;

    dmudx=min(limit,max(dmudx,-limit));
    dmudy=min(limit,max(dmudy,-limit));

    MASKb=zeros(ctr.imax,ctr.jmax); % domain mask for SSA
    MASKb(2:ctr.imax-2,2:ctr.jmax-1)=1;
    MASKb(MASKmy==0)=0;

    V0(MASKb==1)=-5.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)- ...
        betay(MASKb==1);
    V1(MASKb==1)=2.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2.)+ ...
        2.*dmudy(MASKb==1)/ctr.delta;
    V2(MASKb==1)=2.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2.)- ...
        2.*dmudy(MASKb==1)/ctr.delta;
    V3(MASKb==1)=0.5*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2.)+ ...
        0.5*dmudx(MASKb==1)/ctr.delta;
    V4(MASKb==1)=0.5*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2.)- ...
        0.5*dmudx(MASKb==1)/ctr.delta;
    Vu0(MASKb==1)=-1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)+ ...
        (dmudy(MASKb==1)-0.5*dmudx(MASKb==1))/ctr.delta;
    Vu1(MASKb==1)=1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)+ ...
        (dmudy(MASKb==1)+0.5*dmudx(MASKb==1))/ctr.delta;
    Vu2(MASKb==1)=1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)- ...
        (dmudy(MASKb==1)+0.5*dmudx(MASKb==1))/ctr.delta;
    Vu3(MASKb==1)=-1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)- ...
        (dmudy(MASKb==1)-0.5*dmudx(MASKb==1))/ctr.delta;
    S0(MASKb==1)=-taudy(MASKb==1)-betay(MASKb==1).*udy(MASKb==1);

    MASKb=zeros(ctr.imax,ctr.jmax); % domain mask for non-SSA
    MASKb(2:ctr.imax-2,2:ctr.jmax-1)=1;
    MASKb(MASKmy>0)=0;
    V0(MASKb==1)=1;
    S0(MASKb==1)=vsia(MASKb==1);

    % boundary conditions

    if ctr.shelf==1 && ctr.mismip==0
        % i=1; contact with ocean
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(1,2:ctr.jmax-1)=1;
        if ctr.basin==1
            V0(MASKb==1)=1;
            S0(MASKb==1)=v(MASKb==1);
            MASKb(bMASK==1)=0;
        end
        V0(MASKb==1)=-4.*eta1(MASKb==1)/ctr.delta;
        V1(MASKb==1)=4.*eta1(MASKb==1)/ctr.delta;
        Vu1(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
        Vu3(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;
        S0(MASKb==1)=0.5*par.rho*par.g*H1(MASKb==1).^2.*(1.-par.rho/par.rhow);

        % i=imax-1; contact with ocean
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(ctr.imax-1,2:ctr.jmax-1)=1;
        if ctr.basin==1
            V0(MASKb==1)=1;
            S0(MASKb==1)=v(MASKb==1);
            MASKb(bMASK==1)=0;
        end
        V0(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
        V2(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta;
        Vu0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
        Vu2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;
        S0(MASKb==1)=0.5*par.rho*par.g*H(MASKb==1).^2.*(1.-par.rho/par.rhow);

        % j=1; contact with ocean (u-direction)
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(2:ctr.imax-2,1)=1;
        if ctr.basin==1
            V0(MASKb==1)=1;
            S0(MASKb==1)=v(MASKb==1);
            MASKb(bMASK==1)=0;
        end
        V0(MASKb==1)=-1/ctr.delta;
        V3(MASKb==1)=1./ctr.delta;
        Vu0(MASKb==1)=-1./ctr.delta;
        Vu1(MASKb==1)=1./ctr.delta;
        S0(MASKb==1)=0;

        % j=jmax; contact with ocean (u-direction)
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(2:ctr.imax-2,ctr.jmax)=1;
        if ctr.basin==1
            V0(MASKb==1)=1;
            S0(MASKb==1)=v(MASKb==1);
            MASKb(bMASK==1)=0;
        end
        V0(MASKb==1)=1./ctr.delta;
        V4(MASKb==1)=-1./ctr.delta;
        Vu2(MASKb==1)=-1./ctr.delta;
        Vu3(MASKb==1)=1./ctr.delta;
        S0(MASKb==1)=0;

        % Model corners
        V0(1,1)=1;
        V1(1,1)=-1;
        V3(1,1)=-1;
        V5(1,1)=1;
        S0(1,1)=0;
        V0(ctr.imax-1,1)=1;
        V2(ctr.imax-1,1)=-1;
        V3(ctr.imax-1,1)=-1;
        V6(ctr.imax-1,1)=1;
        S0(ctr.imax-1,1)=0;
        V0(1,ctr.jmax)=1;
        V1(1,ctr.jmax)=-1;
        V4(1,ctr.jmax)=-1;
        V7(1,ctr.jmax)=1;
        S0(1,ctr.jmax)=0;
        V0(ctr.imax-1,ctr.jmax)=1;
        V2(ctr.imax-1,ctr.jmax)=-1;
        V4(ctr.imax-1,ctr.jmax)=-1;
        V8(ctr.imax-1,ctr.jmax)=1;
        S0(ctr.imax-1,ctr.jmax)=0;

        % i=imax;
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(ctr.imax,:)=1;
        V0(MASKb==1)=1;
        S0(MASKb==1)=NaN;

    elseif ctr.shelf==1 && ctr.mismip>=1

        % j=1: ice divide (symmetric)
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(1:ctr.imax-1,1)=1;
        V0(MASKb==1)=1;
        V11(MASKb==1)=-1;
        S0(MASKb==1)=0;

        % j=jmax; contact with ocean (u-direction)
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(2:ctr.imax-2,ctr.jmax)=1;
        V0(MASKb==1)=1./ctr.delta;
        V4(MASKb==1)=-1./ctr.delta;
        Vu2(MASKb==1)=-1./ctr.delta;
        Vu3(MASKb==1)=1./ctr.delta;
        S0(MASKb==1)=0;

        % i=1: periodic boundary condition
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(1,2:ctr.jmax)=1;
        V0(MASKb==1)=1;
        V9(MASKb==1)=1;
        S0(MASKb==1)=0;

        if ctr.mismip==1
            % i=imax-1: periodic boundary condition
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(ctr.imax-1,2:ctr.jmax)=1;
            V0(MASKb==1)=1;
            V10(MASKb==1)=1;
            S0(MASKb==1)=0;
        else
            % i=imax-1; contact with ocean
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(ctr.imax-1,2:ctr.jmax-1)=1;
            V0(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
            V2(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta;
            Vu0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            Vu2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;
            S0(MASKb==1)=0.5*par.rho*par.g*H(MASKb==1).^2.*(1.- ...
                par.rho/par.rhow);
            V0(ctr.imax-1,ctr.jmax)=1;
            V2(ctr.imax-1,ctr.jmax)=-1;
            V4(ctr.imax-1,ctr.jmax)=-1;
            V8(ctr.imax-1,ctr.jmax)=1;
            S0(ctr.imax-1,ctr.jmax)=0;
        end

        % i=imax;
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(ctr.imax,:)=1;
        V0(MASKb==1)=1;
        S0(MASKb==1)=NaN;

    else
        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(1:ctr.imax-1,1)=1;
        MASKb(1:ctr.imax-1,ctr.jmax)=1;
        MASKb(1,:)=1;
        MASKb(ctr.imax-1,:)=1;
        V0(MASKb==1)=1;
        S0(MASKb==1)=0;

        MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
        MASKb(ctr.imax,:)=1;
        V0(MASKb==1)=1;
        S0(MASKb==1)=NaN;
    end

    nodes=ctr.imax*ctr.jmax;
    V=[reshape(U0,nodes,1)
        reshape(V0,nodes,1)
        U1(U1~=0)
        U2(U2~=0)
        U3(U3~=0)
        U4(U4~=0)
        U5(U5~=0)
        U6(U6~=0)
        U7(U7~=0)
        U8(U8~=0)
        U9(U9~=0)
        U10(U10~=0)
        Uv0(Uv0~=0)
        Uv1(Uv1~=0)
        Uv2(Uv2~=0)
        Uv3(Uv3~=0)
        V1(V1~=0)
        V2(V2~=0)
        V3(V3~=0)
        V4(V4~=0)
        V5(V5~=0)
        V6(V6~=0)
        V7(V7~=0)
        V8(V8~=0)
        V9(V9~=0)
        V10(V10~=0)
        V11(V11~=0)
        Vu0(Vu0~=0)
        Vu1(Vu1~=0)
        Vu2(Vu2~=0)
        Vu3(Vu3~=0)
      ];

    row=[reshape(nodeu,nodes,1)
        reshape(nodev,nodes,1)
        nodeu(U1~=0)
        nodeu(U2~=0)
        nodeu(U3~=0)
        nodeu(U4~=0)
        nodeu(U5~=0)
        nodeu(U6~=0)
        nodeu(U7~=0)
        nodeu(U8~=0)
        nodeu(U9~=0)
        nodeu(U10~=0)
        nodeu(Uv0~=0)
        nodeu(Uv1~=0)
        nodeu(Uv2~=0)
        nodeu(Uv3~=0)
        nodev(V1~=0)
        nodev(V2~=0)
        nodev(V3~=0)
        nodev(V4~=0)
        nodev(V5~=0)
        nodev(V6~=0)
        nodev(V7~=0)
        nodev(V8~=0)
        nodev(V9~=0)
        nodev(V10~=0)
        nodev(V11~=0)
        nodev(Vu0~=0)
        nodev(Vu1~=0)
        nodev(Vu2~=0)
        nodev(Vu3~=0)
      ];
    nodeU1=circshift(nodeu,[0 -1]); %i,j+1
    nodeU2=circshift(nodeu,[0 1]); %i,j-1
    nodeU3=circshift(nodeu,[-1 0]); %i+1,j
    nodeU4=circshift(nodeu,[1 0]); %i-1,j
    nodeU5=circshift(nodeu,[-1 -1]); %i+1,j+1
    nodeU6=circshift(nodeu,[-1 1]); %i+1,j-1
    nodeU7=circshift(nodeu,[1 -1]); %i-1,j+1
    nodeU8=circshift(nodeu,[1 1]); %i-1,j-1
    nodeU9=circshift(nodeu,[-2 0]); % periodic BC i=1
    nodeU10=circshift(nodeu,[2 0]); %periodic BC i=imax
    nodeUv1=circshift(nodev,[0 -1]); %i,j+1
    nodeUv2=circshift(nodev,[1 0]); %i-1,j
    nodeUv3=circshift(nodev,[1 -1]); %i-1,j+1

    nodeV1=circshift(nodev,[-1 0]); %i+1,j
    nodeV2=circshift(nodev,[1 0]); %i-1,j
    nodeV3=circshift(nodev,[0 -1]); %i,j+1
    nodeV4=circshift(nodev,[0 1]); %i,j-1
    nodeV5=circshift(nodev,[-1 -1]); %i+1,j+1
    nodeV6=circshift(nodev,[1 -1]); %i-1,j+1
    nodeV7=circshift(nodev,[-1 1]); %i+1,j-1
    nodeV8=circshift(nodev,[1 1]); %i-1,j-1
    nodeV9=circshift(nodev,[-1 0]); % periodic BC i=1
    nodeV10=circshift(nodev,[1 0]); %periodic BC i=imax
    nodeV11=circshift(nodev,[0 -2]); % ice divide
    nodeVu1=circshift(nodeu,[-1 0]); %i+1,j
    nodeVu2=circshift(nodeu,[0 1]); %i,j-1
    nodeVu3=circshift(nodeu,[-1 1]); %i+1,j-1

    col=[reshape(nodeu,nodes,1)
        reshape(nodev,nodes,1)
        nodeU1(U1~=0)
        nodeU2(U2~=0)
        nodeU3(U3~=0)
        nodeU4(U4~=0)
        nodeU5(U5~=0)
        nodeU6(U6~=0)
        nodeU7(U7~=0)
        nodeU8(U8~=0)
        nodeU9(U9~=0)
        nodeU10(U10~=0)
        nodev(Uv0~=0)
        nodeUv1(Uv1~=0)
        nodeUv2(Uv2~=0)
        nodeUv3(Uv3~=0)
        nodeV1(V1~=0)
        nodeV2(V2~=0)
        nodeV3(V3~=0)
        nodeV4(V4~=0)
        nodeV5(V5~=0)
        nodeV6(V6~=0)
        nodeV7(V7~=0)
        nodeV8(V8~=0)
        nodeV9(V9~=0)
        nodeV10(V10~=0)
        nodeV11(V11~=0)
        nodeu(Vu0~=0)
        nodeVu1(Vu1~=0)
        nodeVu2(Vu2~=0)
        nodeVu3(Vu3~=0)
      ];

    % R=[reshape(R0,nodes,1)
    %     reshape(S0,nodes,1)];
    R(nodeu)=R0;
    R(nodev)=S0;
    R=R';
    R(isnan(R))=0;

    % construct sparse matrix
    A=sparse(row,col,V);

    % Cholesky factor and solve
    if ctr.ItSolv==1
        D=diag(diag(A));
        C1=tril(A);
        C2=D\triu(A);
        [s,flag]=bicgstab(A,R,par.veltol,par.veliter,C1,C2,s0);
        if flag>0
            s=A\R;
        end
    else
        s=A\R;
    end
    u=s(nodeu);
    v=s(nodev);

end


function [Ax,Ay,Ad]=StaggeredA(A)

% Kori-ULB
% Flow parameter A on different staggered grids. Originally, A is determined
% on the H-grid

    A1=circshift(A,[0 -1]); % A(i,j+1)
    A2=circshift(A,[-1 0]); % A(i+1,j)
    
    Ax=(A+A1)/2.; % A on u-grid
    Ay=(A+A2)/2.; % A on v-grid
    Ad=h2d(A); % A on d-grid

end


function [bMASKm,bMASKx,bMASKy]=StaggeredBMASK(ctr,bMASK)

% Kori-ULB
% MASK for basin calculations on staggered grid

    if ctr.basin==1
        bMASKm=round((bMASK+circshift(bMASK,[-1 0])+circshift(bMASK,[0 -1])+ ...
            circshift(bMASK,[-1 -1]))/4); 
        bMASKx=round((bMASK+circshift(bMASK,[0 -1]))/2);
        bMASKy=round((bMASK+circshift(bMASK,[-1 0]))/2);
    else
        bMASKm=false;
        bMASKx=false;
        bMASKy=false;
    end
    
end


function [gradm,gradmx,gradmy,gradxy,gradsx,gradsy,gradHx,gradHy, ...
    Hm,Hmx,Hmy,Bmx,Bmy,signx,signy]=StaggeredGrid(sn,H,B,ctr)

% Kori-ULB
% Variables on staggered grids

    sn1=circshift(sn,[-1 0]); % sn(i+1,j)
    sn2=circshift(sn,[0 -1]); % sn(i,j+1)
    sn3=circshift(sn,[-1 -1]); % sn(i+1,j+1)
    sn4=circshift(sn,[0 1]); % sn(i,j-1)
    sn5=circshift(sn,[1 0]); % sn(i-1,j)

    Hm=(H+circshift(H,[-1 0])+circshift(H,[0 -1])+circshift(H,[-1 -1]))/4.; 
    gradm=((sn2+sn3-sn-sn1)/(2*ctr.delta)).^2+((sn1+sn3-sn-sn2)/(2*ctr.delta)).^2;
    gradxy=((sn2-sn4)/(2*ctr.delta)).^2+((sn1-sn5)/(2*ctr.delta)).^2;
    
    H1=circshift(H,[-1 0]); % sn(i+1,j)
    H2=circshift(H,[0 -1]); % sn(i,j+1)
    H4=circshift(H,[0 1]); % sn(i,j-1)
    H5=circshift(H,[1 0]); % sn(i-1,j)
    gradsx=(sn2-sn4)/(2*ctr.delta);
    gradsy=(sn1-sn5)/(2*ctr.delta);
    gradHx=(H2-H4)/(2*ctr.delta);
    gradHy=(H1-H5)/(2*ctr.delta);
    
    gradmx=(sn2-sn)/ctr.delta;
    Hmx=(H+circshift(H,[0 -1]))/2.;
    Bmx=(B+circshift(B,[0 -1]))/2.;
    gradmy=(sn1-sn)/ctr.delta;
    Hmy=(H+circshift(H,[-1 0]))/2.;
    Bmy=(B+circshift(B,[-1 0]))/2.;
    
    signx=sign(-gradmx);
    signy=sign(-gradmy);
    gradmx=gradmx.^2;
    gradmy=gradmy.^2;
    
    if ctr.mismip>0
        gradxy(:,1)=gradxy(:,3);
        gradxy(:,end)=gradxy(:,end-1);
        gradxy(1,:)=gradxy(3,:);
        gradxy(end,:)=gradxy(end-3,:);
    end

end


function [Melt,butfac,H]=SubShelfMelt(ctr,fc,par,Tf,To,So,TF,butfac, ...
    HB,glMASK,H,B,ShelfN,numsh,shMASK,MASK,MASKlk,uxssa,uyssa,arcocn, ...
    MeltInv,cnt)

% Kori-ULB
% Options for sub-shelf melt calculation
% Quadratic, linear, PICO, PCIO, Plume, ISMIP6, ABUMIP, ...

    if islogical(TF)==1
        Tfm=max(To,Tf)-Tf;
    else
        Tfm=TF; % Use imposed thermal forcing
    end
    switch ctr.meltfunc
        case 1
            % Linear melt forcing
            Melt=ctr.gammaT*par.LatentMelt*Tfm*par.secperyear;
        case 2
            % Quadratic local melt forcing with mean antarctic slope
            Melt=ctr.gammaT*par.LatentMelt^2*abs(Tfm).*(Tfm)*par.secperyear;
        case 21
            % Quadratic semilocal melt forcing with mean Antarctic slope
            mean_TF=zeros(ctr.imax,ctr.jmax);
            for i=1:numsh
                mean_TF(ShelfN==i)=mean(Tfm(ShelfN==i));
            end
            Melt=ctr.gammaT*par.LatentMelt^2*abs(mean_TF).*(Tfm)*par.secperyear;
        case 22
            % Quadratic semilocal melt forcing with local slope
            sina=ComputeSinSlope(HB,uxssa,uyssa,ctr); 
            mean_TF=zeros(ctr.imax,ctr.jmax);
            for i=1:numsh
                mean_TF(ShelfN==i)=mean(Tfm(ShelfN==i));
            end
            Melt=ctr.gammaT*par.LatentMelt^2*abs(mean_TF).*(Tfm).*sina* ...
                par.secperyear;
        case {3,4}
            % PICO/PICOP
            [Bk,Ak,Bmax]=PICOsetup(glMASK,H,ctr,par,ShelfN,numsh,shMASK);
            [T0o,S0o,zb]=OceanVarBox(numsh,To,So,ShelfN,HB);
            T0o=max(T0o,par.lambda1*S0o+par.lambda2+par.lambda3*zb);
            [T,S]=BoxModel(S0o,T0o,numsh,par,Ak,Bk,ShelfN,HB, ...
                1/par.LatentMelt,ctr);
            if ctr.meltfunc==3
                % PICO
                Melt=PICOMelt(HB,shMASK,ShelfN,Ak,Bk,S0o,T0o,S,T,numsh, ...
                    Bmax,ctr,1/par.LatentMelt,par);
            else
                % PICOP
                [MASKpicop,MASKb]=PICOPmasks(MASK,MASKlk,H,ctr,par);
                [Zgl]=AdvecGL(HB,B,MASK,MASKpicop,MASKb,uxssa,uyssa,ctr);
                sina=ComputeSinSlope(HB,uxssa,uyssa,ctr);
                Melt=PICOPmelt(HB,sina,Zgl,Bk,ShelfN,T,S,par,ctr,Bmax,numsh);
            end
        case 5
            % Plume model
            % Mean ocean properties by shelf cavity (Burgard et al. 2022)
            To_mean_cav=zeros(ctr.imax,ctr.jmax);
            So_mean_cav=zeros(ctr.imax,ctr.jmax);
            for b=1:numsh
                avg_To_int=mean(To(ShelfN==b));
                avg_So_int=mean(So(ShelfN==b));
                To_mean_cav(ShelfN==b)=avg_To_int;
                So_mean_cav(ShelfN==b)=avg_So_int;
            end
            [MASKpicop,MASKb]=PICOPmasks(MASK,MASKlk,H,ctr,par);
            [Zgl]=AdvecGL(HB,B,MASK,MASKpicop,MASKb,uxssa,uyssa,ctr);
            sina=ComputeSinSlope(HB,uxssa,uyssa,ctr);
            Melt=PlumeMelt2019(To_mean_cav,So_mean_cav,HB,Zgl,sina,par,ctr);
        case 6
            % Cornford 2016 melt (function of H)
            Melt=max(min((H-100)*4/7,400),0);
        case 7
            % Uniform melting when butfac=1 - value of meltfac
            GeoBut=ones(ctr.imax,ctr.jmax)*butfac;
            butfac=1;
            Melt=ctr.meltfac*GeoBut;
        case 8
            % Ice shelf removal (needs butfac=[0])
            Melt=zeros(ctr.imax,ctr.jmax);
            if cnt==1
                H(MASK==0)=par.SeaIceThickness;
                Melt(MASK==0)=ctr.meltfac;
            else % can be modified for regrowth of ice shelves
                H(MASK==0)=par.SeaIceThickness;
                Melt(MASK==0)=ctr.meltfac;
            end
        case 9
            % ISMIP6 non-local melt rate parameterization
            Melt=ISMIP6ocean(par,ctr,HB,shMASK,Tfm,fc.basinNumber, ...
                fc.deltaT_basin);
        case 91
            % ISMIP6 non-local melt rate parameterization + dependy on local slope
            sina=ComputeSinSlope(HB,uxssa,uyssa,ctr); 
            Melt=ISMIP6OceanSlope(sina,par,ctr,HB,shMASK,Tfm, ...
                fc.basinNumber,fc.deltaT_basin);
        case 92
            % non-local melt rate parameterization + mean ant slope
            sina=zeros(ctr.imax,ctr.jmax)+2.9e-3; % mean ant slope -- Burgard 2022
            deltaT_basin=zeros(ctr.imax,ctr.jmax); % no local correction
            Melt=ISMIP6OceanSlope(sina,par,ctr,HB,shMASK,Tfm, ...
                fc.basinNumber,deltaT_basin);
        case 10
            % MISMIP+ melt function
            Melt=0.2*tanh((HB-B)/75).*max(-100-HB,0)*ctr.meltfac;
            Melt(MASK==1)=0;
        case 11
            % run model with optimized sub-shelf melt
            Melt=MeltInv;
            Melt(MASK==1)=0;
    end
    switch ctr.meltfunc
        case {1,2,3,4,5,9}
            if par.ArcOcean==1
                Melt(Melt>0)=Melt(Melt>0).*max(0,min(1,(arcocn(Melt>0)-20)./20));
            end
            if ctr.meltfac>0
                Melt=Melt*ctr.meltfac; % Add meltfac to increase melt
            end
    end
    Melt(shMASK==0)=0; % Melt only for real shelves, Melt=0 for 'Lakes'
    Melt(glMASK==6)=0;

end


function [flw,Wd]=SubWaterFlux(ctr,par,H,HB,MASK,Bmelt)

% Kori-ULB
% Steady-state subglacial water flow model
% Method due to LeBrocq (GDS Warner model)
% Simplified implementation for ice thickness classes as a function of
% the distance over hich hydraulic gradient coupling takes place
% Water flux and depth calculated on h-grid

    flux=zeros(ctr.imax,ctr.jmax)-1;
    pot=par.rho*par.g*H+par.rhof*par.g*HB;
    pot(MASK==0)=0;
    [pot]=PotentialFilling(pot,ctr);
    % Correction of ice thickness due to filling of 'pot'
    H=(pot-par.rhof*par.g*HB)/(par.rho*par.g);
    % Hydraulic potential gradients on staggered d-grid
    gdsx=-par.rho*par.g*(circshift(H,[0 -1])-circshift(H,[0 1]))/ ...
        (2*ctr.delta)-par.rhof*par.g*(circshift(HB,[0 -1])- ...
        circshift(HB,[0 1]))/(2*ctr.delta);
    gdsy=-par.rho*par.g*(circshift(H,[-1 0])-circshift(H,[1 0]))/ ...
        (2*ctr.delta)-par.rhof*par.g*(circshift(HB,[-1 0])- ...
        circshift(HB,[1 0]))/(2*ctr.delta);

    % define multiplier for GDS
    Hlev=max(mean(H(MASK==1)),10);
    scale=Hlev*par.longcoupwater*2.; % scale of GDS window
    width=2.*scale;
    scale(width<=ctr.delta)=ctr.delta/2+1;
    maxlevel=2*round(width/ctr.delta-0.5)+1;
    frb=(maxlevel-1)/2;
    [ew,ns]=meshgrid(1:maxlevel,1:maxlevel); % size of GDS window

    kegdsx0=zeros(ctr.imax+2*frb,ctr.jmax+2*frb);
    kegdsy0=zeros(ctr.imax+2*frb,ctr.jmax+2*frb);
    kegdsx0(frb+1:frb+ctr.imax,frb+1:frb+ctr.jmax)=gdsx;
    kegdsy0(frb+1:frb+ctr.imax,frb+1:frb+ctr.jmax)=gdsy;

    % Define mulipliers as a function of ice thickness (scale) and convolute
    dist=sqrt((ctr.delta*(ew-frb-1)).^2+(ctr.delta*(ns-frb-1)).^2)/scale;
    multiplier=max(0.,1.-dist/2.);
    multiplier=multiplier/sum(sum((multiplier)));
    kegx=xcorr2(kegdsx0,multiplier);
    kegy=xcorr2(kegdsy0,multiplier);
    kegdsx=kegx(2*frb+1:ctr.imax+2*frb,2*frb+1:ctr.jmax+2*frb);
    kegdsy=kegy(2*frb+1:ctr.imax+2*frb,2*frb+1:ctr.jmax+2*frb);

    % Recursive algorithm
    gdsmag=abs(kegdsx)+abs(kegdsy);
    for i=1:ctr.imax-1
        for j=1:ctr.jmax-1
            if MASK(i,j)==1
                flux=DpareaWarGds(i,j,MASK,Bmelt, ...
                    kegdsx,kegdsy,par,ctr,flux,gdsmag,1);
            end
        end
    end

    corfac=zeros(ctr.imax,ctr.jmax)+1;
    denom=sqrt(kegdsx.^2+kegdsy.^2);
    corfac(gdsmag>0)=gdsmag(gdsmag>0)./denom(gdsmag>0);
    flw=min(max(flux./corfac/ctr.delta,0),1e5);
    flw(MASK==0)=0;
    Wd=min(par.Wdmax,max(par.Wdmin,(12*par.waterviscosity*flw/ ...
        mean(gdsmag(MASK==1))).^(1/3)));
end


function [tmp,ctr]=Temperature3d(tmp,Mb,Ts,pxy,par, ...
    ctr,dt,gradsx,gradsy,gradHx,gradHy,udx,udy,ub,ubx,uby,zeta,gradxy,H, ...
    dzc,dzm,dzp,G,taudxy,A,DeltaT,MASK,Bmelt,cnt)

% Kori-ULB
% 3d englacial temperature calculation in ice sheet and ice shelves

    tmp(:,:,1)=Ts+par.T0;
    
    % horizontal velocities on H grid
    uxdt=0.5*(udx+circshift(udx,[0 1]));
    uydt=0.5*(udy+circshift(udy,[1 0]));
    uxbt=0.5*(ubx+circshift(ubx,[0 1]));
    uybt=0.5*(uby+circshift(uby,[1 0]));
    pl=repmat(pxy,[1,1,ctr.kmax]);
    zl=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    ushllib=(pl+2)./(pl+1).*(1-zl.^(pl+1));
    ut=repmat(uxdt,[1,1,ctr.kmax]).*ushllib+repmat(uxbt,[1,1,ctr.kmax]);
    vt=repmat(uydt,[1,1,ctr.kmax]).*ushllib+repmat(uybt,[1,1,ctr.kmax]);
    
    % horizontal advection
    dTdxm=(circshift(tmp,[0 -1])-tmp)/ctr.delta;
    dTdxp=(tmp-circshift(tmp,[0 1]))/ctr.delta;
    advecx=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    advecx(ut>0)=ut(ut>0).*dTdxp(ut>0)*dt;
    advecx(ut<0)=ut(ut<0).*dTdxm(ut<0)*dt;
    dTdym=(circshift(tmp,[-1 0])-tmp)/ctr.delta;
    dTdyp=(tmp-circshift(tmp,[1 0]))/ctr.delta;
    advecy=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    advecy(vt>0)=vt(vt>0).*dTdyp(vt>0)*dt;
    advecy(vt<0)=vt(vt<0).*dTdym(vt<0)*dt;
    
    % adjusted vertical velocity (old f.ETISh version)
%     ws=-max(Mb,1e-5)-(ud.*ushllib(:,:,1)+ub).*sqrt(gradxy);
%     wshllib=1-(pl+2).*zl./(pl+1)+1./(pl+1).*zl.^(pl+2);
%     w=repmat(ws,[1,1,ctr.kmax]).*wshllib;
    
    % vertical velocity according to Pattyn (2010) with Lliboutry shape
    % function
    % FP: check whether +Bmelt (Huybrechts) or -Bmelt (Pattyn, 2010) - it
    % should be -Bmelt I think
    wshllib=1-(pl+2).*zl./(pl+1)+1./(pl+1).*zl.^(pl+2);
    w=repmat(-max(Mb,1e-5),[1,1,ctr.kmax]).*wshllib-Bmelt+ut.* ...
        (repmat(gradsx,[1,1,ctr.kmax])-zl.*repmat(gradHx,[1,1,ctr.kmax])) ...
        +vt.*(repmat(gradsy,[1,1,ctr.kmax])-zl.*repmat(gradHy,[1,1,ctr.kmax]));
    
%     % integrated verical velocity according to Pattyn (2003)
%     w=zeros(ctr.imax,ctr.jmax,ctr.kmax);
%     u=repmat(udx,[1,1,ctr.kmax]).*ushllib+repmat(ubx,[1,1,ctr.kmax]);
%     v=repmat(udy,[1,1,ctr.kmax]).*ushllib+repmat(uby,[1,1,ctr.kmax]);
%     dudx=(u-circshift(u,[0 1 0]))/ctr.delta;
%     dudx(u<0)=(circshift(u(u<0),[0 -1 0])-u(u<0))/ctr.delta;
%     dvdy=(v-circshift(v,[0 1 0]))/ctr.delta;
%     dvdy(v<0)=(circshift(v(v<0),[0 -1 0])-v(v<0))/ctr.delta;
%     w(:,:,ctr.kmax)=u(:,:,ctr.kmax).*(gradsx-gradHx)+v(:,:,ctr.kmax).* ...
%         (gradsy-gradHy)-Bmelt;
%     for k=ctr.kmax-1:-1:1
%         u1=0.5*H.*(dudx(:,:,k)+dudx(:,:,k+1));
%         v1=0.5*H.*(dvdy(:,:,k)+dvdy(:,:,k+1));
%         u2=(u(:,:,k+1)-u(:,:,k)).*(gradsx-0.5*(zeta(k)+zeta(k+1))*gradHx);
%         v2=(v(:,:,k+1)-v(:,:,k)).*(gradsy-0.5*(zeta(k)+zeta(k+1))*gradHy);
%         w(:,:,k)=(u1+u2+v1+v2)*(zeta(k)-zeta(k+1))+w(:,:,k+1);
%     end
    
    % Internal heating
    repz=repmat(reshape(zeta,1,1,ctr.kmax),[ctr.imax,ctr.jmax,1]);
    dudz=repmat(2*A.*taudxy.^par.n.*H.*(pxy+2)/(par.n+2),[1,1,ctr.kmax]).* ...
        repz.^pl;
    fric=par.rho*par.g*par.kdif*dt*repz.*dudz.*repmat(sqrt(gradxy), ...
        [1,1,ctr.kmax])/par.K;
    repmask=repmat(MASK,[1,1,ctr.kmax]);
    fric(repmask==0)=0; % no frictional heat in ice shelves
    extraterm=max(min(fric-advecx-advecy,5),-5); %10

    % Temperature solution
    repH=repmat(H+1e-8,[1,1,ctr.kmax]);
    Tp=par.pmp*repH.*repz;
    atp=(2*par.kdif*par.secperyear./(repH.*dzm)-w)*dt./(repH.*dzc);
    btp=1+2*par.kdif*par.secperyear*dt./((repH.^2).*dzp.*dzm);
    ctp=(2*par.kdif*par.secperyear./(repH.*dzp)+w)*dt./(repH.*dzc);
    ftp=ones(ctr.imax,ctr.jmax,ctr.kmax);
    gtp=zeros(ctr.imax,ctr.jmax,ctr.kmax);
    % basal BC with strain heating (correction with ud)
    gtp(:,:,ctr.kmax)=(G+taudxy.*ub/par.secperyear).*H.*dzm(:,:,ctr.kmax)/par.K;
    
    % ice shelves basal boundary condition
    mask=zeros(ctr.imax,ctr.jmax);
    mask(MASK==0)=1;
    repmask=repmat(mask,[1,1,ctr.kmax]);
    repmask(:,:,1:ctr.kmax-1)=0;
    TBshelf=par.T0+repmat(min(par.Toi+ctr.meltfactor*DeltaT(cnt)- ...
        0.12e-3*par.rho*H/par.rhow,0),[1,1,ctr.kmax]);
    ftp(repmask==1)=0;
    gtp(repmask==1)=TBshelf(repmask==1);
    
    % Temperature solution
    for k=ctr.kmax-1:-1:2
        ftp(:,:,k)=atp(:,:,k)./(btp(:,:,k)-ctp(:,:,k).*ftp(:,:,k+1));
        gtp(:,:,k)=(tmp(:,:,k)+extraterm(:,:,k)+ ...
            ctp(:,:,k).*gtp(:,:,k+1))./(btp(:,:,k)-ctp(:,:,k).*ftp(:,:,k+1));
    end
    for k=2:ctr.kmax
        tmp(:,:,k)=tmp(:,:,k-1).*ftp(:,:,k)+gtp(:,:,k);
    end
    c3=ceil(par.rhom/par.rho)*10+ceil(par.rhow/par.rho);
    ctr.runmode(DeltaT(cnt)==c3)=5;
    
    % Correction for unstable temperature profiles
    % find an anomaly where temperature decreases with depth
    [ipos,jpos]=find((tmp(:,:,ctr.kmax-1)-tmp(:,:,ctr.kmax))>3 | ...
        (tmp(:,:,ctr.kmax-2)-tmp(:,:,ctr.kmax-1))>3 | ...
        (tmp(:,:,ctr.kmax-3)-tmp(:,:,ctr.kmax-2))>3);
    nMASK=zeros(size(MASK));
    nMASK(sub2ind(size(nMASK),ipos,jpos))=1;
    % use analytical solution at the domain boundary when using basins
    if ctr.basin==1
        nMASK(1,:)=1;
        nMASK(ctr.imax,:)=1;
        nMASK(:,1)=1;
        nMASK(:,ctr.jmax)=1;
    end
    % apply linear temperature profile when MASK=0
    nMASK(nMASK==1 & MASK==0)=2;
    Tgrad=-G/par.K;
    l=sqrt(2*par.kdif*(H+1e-8)./max(Mb,1e-8)*par.secperyear);
    nMASKz=repmat(nMASK,[1,1,ctr.kmax]);
    repl=repmat(l,[1,1,ctr.kmax]);
    repTs=repmat(Ts,[1,1,ctr.kmax]);
    repTgrad=repmat(Tgrad,[1,1,ctr.kmax]);
    tmpb=repTs+sqrt(pi)*0.5*repl.*repTgrad.* ...
        (erf((1-repz).*repH./repl)-erf(repH./repl))+par.T0;
    tmp(nMASKz==1)=tmpb(nMASKz==1);
    Tshelf=repTs+par.T0-(repTs+par.T0-TBshelf).*repz;
    tmp(nMASKz==2)=Tshelf(nMASKz==2); 

    % Correction for pmp
    tmp(tmp>par.T0-Tp)=par.T0-Tp(tmp>par.T0-Tp);

    mintemp=min(min(Ts))+min(DeltaT)+par.T0-5;
    tmp(tmp<mintemp)=mintemp;
    
    if ctr.mismip>0
        tmp(:,1,:)=tmp(:,3,:);
        tmp(:,end,:)=tmp(:,end-1,:);
        tmp(1,:,:)=tmp(3,:,:);
        tmp(end,:,:)=tmp(end-2,:,:);
    end
    
end


function [A,Ax,Ay,Ad]=ThermoCoupling(ctr,par,Tb,Tbc,H,bMASK,bMASKm,bMASKx,bMASKy)

% Kori-ULB
% Thermomechanical coupling using Arrhenius relationship

    A=zeros(ctr.imax,ctr.jmax)+ctr.Ao;
    if ctr.Tcalc==2
        A=0.5*par.atune*par.a1*exp(par.Q1./par.R*(1./(par.T0+par.pmp*H)-1./ ...
            (Tbc+par.T0)));
        A(Tb>-6.5)=0.5*par.atune*par.a2*exp(par.Q2./par.R*(1./(par.T0+ ...
            par.pmp*H(Tb>-6.5))-1./(Tbc(Tb>-6.5)+par.T0)));
        [Ax,Ay,Ad]=StaggeredA(A);
    else
        Ax=A;
        Ay=A;
        Ad=A;
    end
    if ctr.basin==1
        A(bMASK==1)=par.A0; % A on h-grid
        Ax(bMASKx==1)=par.A0;
        Ay(bMASKy==1)=par.A0;
        Ad(bMASKm==1)=par.A0; % A on d-grid
    end

end


function [Ntil,Wtil]=TillWater(MASK,Wtil,Bmelt,ctr,Po,par)

% Kori-ULB
% Effective pressure in till based on subglacial water saturation

    Wtil=max(par.Wdmin,min(Wtil+(Bmelt-par.Cdr)*ctr.dt,par.Wmax));
    s=Wtil/par.Wmax;
    s(MASK==0)=1;
    Ntil=par.N0*((par.sigmat*Po/par.N0).^s).*10.^((par.e0/par.Cc).*(1.-s));
    Ntil=max(min(Po,Ntil),par.sigmat*Po);

end


function [flux]=TotalFlux(H,ux,uy,delta)

% Kori-ULB
% Calculation of ice sheet flux using staggered grids

    ux1=circshift(ux,[0 1]); % ux(i,j-1)
    uy1=circshift(uy,[1 0]); % uy(i-1,j)
    flux=H.*sqrt((0.5*(ux+ux1)).^2+(0.5*(uy+uy1)).^2)*delta;
    
end


function [dtr]=TransportDamage(node,nodes,dtr,Mb,Melt,H,MASK,dtdx,dtdx2, ...
    u,v,ctr,cnt,bMASK,VM,par)

% Kori-ULB
% Sparse solver of damage transport (based on ice thickness solver)

    epsilon=1e-10; % artificial diffusion
    MASK(MASK==6)=0;

    um1=circshift(u,[0 1]); % u(i,j-1)
    vm1=circshift(v,[1 0]); % v(i-1,j)
    if ctr.upstream==1
        up1=circshift(u,[0 -1]); % u(i,j+1)
        vp1=circshift(v,[-1 0]); % v(i+1,j)
        um2=circshift(u,[0 2]); % u(i,j-2)
        vm2=circshift(v,[2 0]); % v(i-2,j)
    end

    if ctr.upstream==1
        % conditions for diffusion scheme (init)
        V0=zeros(ctr.imax,ctr.jmax)+8*epsilon*dtdx; % i,j
        V1=zeros(ctr.imax,ctr.jmax)-2*epsilon*dtdx; % i,j+1
        V2=V1; % i,j-1
        V3=V1; % i+1,j
        V4=V1; % i-1,j
        V5=zeros(ctr.imax,ctr.jmax); % i,j-2
        V6=zeros(ctr.imax,ctr.jmax); % i-2,j
        V7=zeros(ctr.imax,ctr.jmax); % i,j+2
        V8=zeros(ctr.imax,ctr.jmax); % i+2,j

        % Velocity sign masks
        MU=zeros(ctr.imax,ctr.jmax);
        MU(u>=0 & um1>=0 & um2>=0)=1;
        MU(u<=0 & um1<=0 & up1<=0)=2;
        MV=zeros(ctr.imax,ctr.jmax);
        MV(v>=0 & vm1>=0 & vm2>=0)=1;
        MV(v<=0 & vm1<=0 & vp1<=0)=2;

        if ctr.basin==1
            MU(bMASK==1)=0;
            MV(bMASK==1)=0;
        end

        V0a=zeros(ctr.imax,ctr.jmax);
        V1a=zeros(ctr.imax,ctr.jmax);
        V2a=zeros(ctr.imax,ctr.jmax);
        V3a=zeros(ctr.imax,ctr.jmax);
        V4a=zeros(ctr.imax,ctr.jmax);

        % conditions for MU=0 (central difference)
        V0a(MU==0)=u(MU==0)-um1(MU==0); % i,j
        V1a(MU==0)=u(MU==0); % i,j+1
        V2a(MU==0)=-um1(MU==0); % i,j-1

        % conditions for MU=1 (grad(u)>0)
        V0a(MU==1)=2*u(MU==1)+um1(MU==1); % i,j
        V2a(MU==1)=-3*um1(MU==1)-um2(MU==1); % i,j-1
        V5(MU==1)=um2(MU==1); % i,j-2

        % conditions for MU=2 and (grad(u)<0)
        V0a(MU==2)=-u(MU==2)-2*um1(MU==2); % i,j
        V1a(MU==2)=3*u(MU==2)+up1(MU==2); % i,j+1
        V7(MU==2)=(-up1(MU==2)); % i,j+2

        % conditions for MV=0 (central difference)
        V0a(MV==0)=V0a(MV==0)+v(MV==0)-vm1(MV==0); % i,j
        V3a(MV==0)=v(MV==0); % i+1,j
        V4a(MV==0)=-vm1(MV==0); % i-1,j

        % conditions for MV=1 (grad(v)>0)
        V0a(MV==1)=V0a(MV==1)+2*v(MV==1)+vm1(MV==1); % i,j
        V4a(MV==1)=-3*vm1(MV==1)-vm2(MV==1); % i-1,j
        V6(MV==1)=vm2(MV==1); % i-2,j

        % conditions for MV=2 (grad(v)<0)
        V0a(MV==2)=V0a(MV==2)-v(MV==2)-2*vm1(MV==2); % i,j
        V3a(MV==2)=3*v(MV==2)+vp1(MV==2); % i+1,j
        V8(MV==2)=-vp1(MV==2); % i+2,j

        % Filling V-matrix
        V0=V0+V0a*dtdx2;%.*(1.-alfa);
        V1=V1+V1a*dtdx2;%.*(1.-alfa);
        V2=V2+V2a*dtdx2;%.*(1.-alfa);
        V3=V3+V3a*dtdx2;%.*(1.-alfa);
        V4=V4+V4a*dtdx2;%.*(1.-alfa);
        V5=V5*dtdx2;%.*(1.-alfa);
        V6=V6*dtdx2;%.*(1.-alfa);
        V7=V7*dtdx2;%.*(1.-alfa);
        V8=V8*dtdx2;%.*(1.-alfa);
    else
        V0=8*epsilon*dtdx+dtdx2*(u-um1+v-vm1); % i,j
        V1=-2*epsilon*dtdx+dtdx2*u; % i,j+1
        V2=-2*epsilon*dtdx-dtdx2*um1; % i,j-1
        V3=-2*epsilon*dtdx+dtdx2*v; % i+1,j
        V4=-2*epsilon*dtdx-dtdx2*vm1; % i-1,j
    end

    R0=dtr-dtr*(1-par.omega).*V0-circshift(dtr,[0 -1])*(1-par.omega) ...
        .*V1-circshift(dtr,[0 1])*(1-par.omega).*V2-circshift(dtr,[-1 0])* ...
        (1-par.omega).*V3-circshift(dtr,[1 0])*(1-par.omega).*V4 ...
        +(max(Mb,0)+max(Melt,0))*ctr.dt./max(H,1e-5);
    if ctr.upstream==1
        R0=R0-circshift(dtr,[0 2])*(1-par.omega).*V5-circshift(dtr,[2 0])* ...
            (1-par.omega).*V6-circshift(dtr,[0 -2])*(1-par.omega).*V7- ...
            circshift(dtr,[-2 0])*(1-par.omega).*V8;
    end

    V0(MASK==0)=0; % note that for shelf=1, MASK=glMASK in the call
    V1(MASK==0)=0;
    V2(MASK==0)=0;
    V3(MASK==0)=0;
    V4(MASK==0)=0;
    if ctr.upstream==1
        V5(MASK==0)=0;
        V6(MASK==0)=0;
        V7(MASK==0)=0;
        V8(MASK==0)=0;
    end
    R0(MASK==0)=dtr(MASK==0);

    % boundaries
    V9=zeros(ctr.imax,ctr.jmax); % ice divide or ocean
    V10=zeros(ctr.imax,ctr.jmax); % i=1 periodic boundary or ocean
    V11=zeros(ctr.imax,ctr.jmax); % i=imax periodic boundary or ocean
    V12=zeros(ctr.imax,ctr.jmax); % j=jmax ocean contact

    wholemask=ctr.imax*ctr.jmax-sum(MASK(:));
    if wholemask~=0 % only when domain is not MASK=1 everywhere
        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(:,1)=1; % symmetric divide or ocean
        V0(MASKb==1)=0;
        V9(MASKb==1)=-1;
        R0(MASKb==1)=0;

        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(1,2:ctr.jmax-1)=1; % periodic BC at i=1 or ocean
        V0(MASKb==1)=0;
        V10(MASKb==1)=-1;
        R0(MASKb==1)=0;

        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(ctr.imax,2:ctr.jmax-1)=1; % periodic BC at i=imax or ocean
        V0(MASKb==1)=0;
        V11(MASKb==1)=-1;
        R0(MASKb==1)=0;

        MASKb=zeros(ctr.imax,ctr.jmax);
        MASKb(1:ctr.imax,ctr.jmax)=1; % ocean
        V0(MASKb==1)=0;
        V12(MASKb==1)=-1;
        R0(MASKb==1)=0;
    end

    MASKb=zeros(ctr.imax,ctr.jmax);
    MASKb(1,:)=1;
    MASKb(ctr.imax,:)=1;
    MASKb(:,1)=1;
    MASKb(:,ctr.jmax)=1;
    V1(MASKb==1)=0;
    V2(MASKb==1)=0;
    V3(MASKb==1)=0;
    V4(MASKb==1)=0;
    if ctr.upstream==1
        V5(MASKb==1)=0;
        V6(MASKb==1)=0;
        V7(MASKb==1)=0;
        V8(MASKb==1)=0;
    end
    if wholemask==0
        V0(MASKb==1)=0;
        R0(MASKb==1)=0;
    end

    if ctr.upstream==1
        V=[reshape(V0(VM==1)*par.omega+1,nodes,1)
            V1(V1~=0)*par.omega
            V2(V2~=0)*par.omega
            V3(V3~=0)*par.omega
            V4(V4~=0)*par.omega
            V5(V5~=0)*par.omega
            V6(V6~=0)*par.omega
            V7(V7~=0)*par.omega
            V8(V8~=0)*par.omega
            V9(V9~=0)
            V10(V10~=0)
            V11(V11~=0)
            V12(V12~=0)];

        row=[reshape(node(VM==1),nodes,1)
            node(V1~=0)
            node(V2~=0)
            node(V3~=0)
            node(V4~=0)
            node(V5~=0)
            node(V6~=0)
            node(V7~=0)
            node(V8~=0)
            node(V9~=0)
            node(V10~=0)
            node(V11~=0)
            node(V12~=0)];
    else
        V=[reshape(V0(VM==1)*par.omega+1,nodes,1)
            V1(V1~=0)*par.omega
            V2(V2~=0)*par.omega
            V3(V3~=0)*par.omega
            V4(V4~=0)*par.omega
            V9(V9~=0)
            V10(V10~=0)
            V11(V11~=0)
            V12(V12~=0)];

        row=[reshape(node(VM==1),nodes,1)
            node(V1~=0)
            node(V2~=0)
            node(V3~=0)
            node(V4~=0)
            node(V9~=0)
            node(V10~=0)
            node(V11~=0)
            node(V12~=0)];
    end

    nodeV1=circshift(node,[0 -1]); % i,j+1
    nodeV2=circshift(node,[0 1]); % i,j-1
    nodeV3=circshift(node,[-1 0]); % i+1,j
    nodeV4=circshift(node,[1 0]); % i-1,j
    if ctr.upstream==1
        nodeV5=circshift(node,[0 2]); % i,j-2
        nodeV6=circshift(node,[2 0]); % i-2,j
        nodeV7=circshift(node,[0 -2]); % i,j+2
        nodeV8=circshift(node,[-2 0]); % i+2,j
    end
    if ctr.mismip>=1
        nodeV9=circshift(node,[0 -2]); % i,j+2 - divide
        nodeV10=circshift(node,[-2 0]); % 3,j - symmetry at i=1
        if ctr.mismip==1
            nodeV11=circshift(node,[2 0]); % n-2,j - PBC at i=imax
        else
            nodeV11=circshift(node,[1 0]);
        end
        nodeV12=circshift(node,[0 1]); % i,jmax-1 - ocean
    else
        nodeV9=circshift(node,[0 -1]); % i,2 - ocean
        nodeV10=circshift(node,[-1 0]); % 2,j - ocean
        nodeV11=circshift(node,[1 0]); % imax-1,j - ocean
        nodeV12=circshift(node,[0 1]); % i,jmax-1 - ocean
    end

    if ctr.upstream==1
        col=[reshape(node(VM==1),nodes,1)
            nodeV1(V1~=0)
            nodeV2(V2~=0)
            nodeV3(V3~=0)
            nodeV4(V4~=0)
            nodeV5(V5~=0)
            nodeV6(V6~=0)
            nodeV7(V7~=0)
            nodeV8(V8~=0)
            nodeV9(V9~=0)
            nodeV10(V10~=0)
            nodeV11(V11~=0)
            nodeV12(V12~=0)];
    else
        col=[reshape(node(VM==1),nodes,1)
            nodeV1(V1~=0)
            nodeV2(V2~=0)
            nodeV3(V3~=0)
            nodeV4(V4~=0)
            nodeV9(V9~=0)
            nodeV10(V10~=0)
            nodeV11(V11~=0)
            nodeV12(V12~=0)];
    end

    R=reshape(R0(VM==1),nodes,1);

    % construct sparse matrix
    A=sparse(row,col,V);
    % Cholesky factor and solve
    if ctr.inverse==1 || ctr.ItSolv==0
        s=A\R;
    else
        D=diag(diag(A));
        C1=tril(A);
        C2=D\triu(A);
        [s,flag]=pcg(A,R,par.Htol,par.Hiter,C1,C2);
        if flag>0 || cnt==1
            s=A\R;
        end
    end

    dtr(node>0)=s(node(node>0));
    dtr=max(0,min(H-eps,dtr)); % FP: put limit on maximum dtr as 1/20 of ice thickness

end


function [Ts,Tsf]=TsFunc(ctr,par,Tsf,S0,sn,DeltaSL,DeltaT)

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


function flux=VariabFlux(variab,H,MASK,seaice)

% Kori-ULB
% Calculate mass balance component fluxes

    flux=sum(variab(MASK==1 | (MASK==0 & H>seaice)));
    
end


function [FMB]=VerticalFaceMelt(ctr,par,SLR,B,Melt,MASK,glMASK,he)

% Kori-ULB
% Calculation of ice shelf frontal melt

    Dz=SLR-B;
    Dz(MASK==0)=par.rho*he(MASK==0)/par.rhow;
    FMB=max(Dz.*Melt/ctr.delta,0);
    FMB(glMASK~=5)=0;
    
end

        
function C = conv2fft(varargin)
% C = conv2fft(A, B)
% C = conv2fft(H1, H2, A)
%
%   C = CONV2FFT(A, B) performs the 2-D convolution of matrices A and B.
%   If [ma,na] = size(A), [mb,nb] = size(B), and [mc,nc] = size(C), then
%   mc = max([ma+mb-1,ma,mb]) and nc = max([na+nb-1,na,nb]).
%
%   C = CONV2FFT(H1, H2, A) first convolves each column of A with the vector
%   H1 and then convolves each row of the result with the vector H2.  If
%   n1 = length(H1), n2 = length(H2), and [mc,nc] = size(C) then
%   mc = max([ma+n1-1,ma,n1]) and nc = max([na+n2-1,na,n2]).
%   CONV2(H1, H2, A) is equivalent to CONV2FFT(H1(:)*H2(:).', A) up to
%   round-off.
%
%   C = CONV2FFT(..., SHAPE) returns a subsection of the 2-D
%   convolution with size specified by SHAPE:
%     'full'  - (default) returns the full 2-D convolution,
%     'same'  - returns the central part of the convolution
%               that is the same size as A.
%     'valid' - returns only those parts of the convolution
%               that are computed without the zero-padded edges.
%               size(C) = max([ma-max(0,mb-1),na-max(0,nb-1)],0).
%
%   See also CONV2, CONVN, CONVNFFT
%
%   Author: Bruno Luong <brunoluong@yahoo.com>
%   History:
%       Original: 21-April-2014
if length(varargin)>=3 && isnumeric(varargin{3})
    [H1, H2, A] = deal(varargin{1:3});
    varargin = varargin(4:end);
    C = convnfft(H1(:), A, varargin{:});
    C = convnfft(H2(:).', C, varargin{:});
else
    C = convnfft(varargin{:});
end
end % conv2fft


function A = convnfft(A, B, shape, dims, options)
% CONVNFFT  FFT-BASED N-dimensional convolution.
%   C = CONVNFFT(A, B) performs the N-dimensional convolution of
%   matrices A and B. If nak = size(A,k) and nbk = size(B,k), then
%   size(C,k) = max([nak+nbk-1,nak,nbk]);
% 
%   C = CONVNFFT(A, B, SHAPE) controls the size of the answer C:
%       'full'   - (default) returns the full N-D convolution
%       'same'   - returns the central part of the convolution that
%                  is the same size as A.
%       'valid'  - returns only the part of the result that can be
%                  computed without assuming zero-padded arrays.
%                  size(C,k) = max([nak-max(0,nbk-1)],0).
%
%   C = CONVNFFT(..., SHAPE, DIMS) with DIMS is vector of dimensions where
%       the convolution will be carried out. By default DIMS is
%       [1:max(ndims(A),ndims(B))] (all dimensions). A and B must have the
%       same lengths on other dimensions.
%   C = CONVNFFT(..., SHAPE, DIMS, GPU)
%       GPU is boolean flag, see next
%
%   C = CONVNFFT(..., SHAPE, DIMS, OPTIONS)
%
%   OPTIONS is structure with following optional fields
%       - 'GPU', boolean. If GPU is TRUE Jacket/GPU FFT engine will be used
%       By default GPU is FALSE.
%       - 'Power2Flag', boolean. If it is TRUE, use FFT with length rounded
%       to the next power-two. It is faster but requires more memory.
%       Default value is TRUE.
%
% Class support for inputs A,B:
% float: double, single
%
% METHOD: CONVNFFT uses Fourier transform (FT) convolution theorem, i.e.
%         FT of the convolution is equal to the product of the FTs of the
%         input functions.
%         In 1-D, the complexity is O((na+nb)*log(na+nb)), where na/nb are
%         respectively the lengths of A and B.
%
% Usage recommendation:
%         In 1D, this function is faster than CONV for nA, nB > 1000.
%         In 2D, this function is faster than CONV2 for nA, nB > 20.
%         In 3D, this function is faster than CONVN for nA, nB > 5.
% 
% See also conv, conv2, convn.
% 
%   Author: Bruno Luong <brunoluong@yahoo.com>
%   History:
%       Original: 21-Jun-2009
%       23-Jun-2009: correct bug when ndims(A)<ndims(B)
%       02-Sep-2009: GPU/JACKET option
%       04-Sep-2009: options structure
%       16-Sep-2009: inplace product
if nargin<3 || isempty(shape)
    shape = 'full';
end
if nargin<5 || isempty(options)
    options = struct();
elseif ~isstruct(options) % GPU options
    options = struct('GPU', options);
end
nd = max(ndims(A),ndims(B));
% work on all dimensions by default
if nargin<4 || isempty(dims)
    dims = 1:nd;
end
dims = reshape(dims, 1, []); % row (needed for for-loop index)
% GPU enable flag
GPU = getoption(options, 'GPU', false);
% Check if Jacket is installed
GPU = GPU && ~isempty(which('ginfo'));
% IFUN function will be used later to truncate the result
% M and N are respectively the length of A and B in some dimension
switch lower(shape)
    case 'full',
        ifun = @(m,n) 1:m+n-1;
    case 'same',
        ifun = @(m,n) ceil((n-1)/2)+(1:m);
    case 'valid',
        ifun = @(m,n) n:m;
    otherwise
        error('convnfft: unknown shape %s', shape);
end
classA = class(A);
classB = class(B);
ABreal = isreal(A) && isreal(B);
% Special case, empty convolution, try to follow MATLAB CONVN convention
if any(size(A)==0) || any(size(B)==0)
    szA = zeros(1,nd); szA(1:ndims(A))=size(A);
    szB = zeros(1,nd); szB(1:ndims(B))=size(B);
    % Matlab wants these:
    szA = max(szA,1); szB = max(szB,1);
    szC = szA;
    for dim=dims
        szC(dim) = length(ifun(szA(dim),szB(dim)));
    end
    A = zeros(szC,classA); % empty -> return zeros
    return
end
power2flag = getoption(options, 'Power2Flag', true);
if power2flag
    % faster FFT if the dimension is power of 2
    lfftfun = @(l) 2^nextpow2(l);
else
    % slower, but smaller temporary arrays
    lfftfun = @(l) l;
end
if GPU % GPU/Jacket FFT
    if strcmp(classA,'single')
        A = gsingle(A);
    else
        A = gdouble(A);
    end
    if strcmp(classB,'single')
        B = gsingle(B);
    else
        B = gdouble(B);
    end
    % Do the FFT
    subs(1:ndims(A)) = {':'};
    for dim=dims
        m = size(A,dim);
        n = size(B,dim);
        % compute the FFT length
        l = lfftfun(m+n-1);
        % We need to swap dimensions because GPU FFT works along the
        % first dimension
        if dim~=1 % do the work when only required
            swap = 1:nd;
            swap([1 dim]) = swap([dim 1]);
            A = permute(A, swap);
            B = permute(B, swap);
        end
        A = fft(A,l);
        B = fft(B,l);
        subs{dim} = ifun(m,n);
    end
else % Matlab FFT
    % Do the FFT
    subs(1:ndims(A)) = {':'};
    for dim=dims
        m = size(A,dim);
        n = size(B,dim);
        % compute the FFT length
        l = lfftfun(m+n-1);
        A = fft(A,l,dim);
        B = fft(B,l,dim);
        subs{dim} = ifun(m,n);
    end
end
 
if GPU
    A = A.*B;
    clear B
else
%     % inplace product to save 1/3 of the memory
%     inplaceprod(A,B);
    A = A.*B;
end
% Back to the non-Fourier space
if GPU % GPU/Jacket FFT
    for dim=dims(end:-1:1) % reverse loop
        A = ifft(A,[]);
        % Swap back the dimensions
        if dim~=1 % do the work when only required
            swap = 1:nd;
            swap([1 dim]) = swap([dim 1]);
            A = permute(A, swap);
        end        
    end   
else % Matlab IFFT  
    for dim=dims
        A = ifft(A,[],dim);
    end
end
% Truncate the results
if ABreal
    % Make sure the result is real
    A = real(A(subs{:}));
else
    A = A(subs{:});
end
% GPU/Jacket
if GPU
    % Cast the type back
    if strcmp(class(A),'gsingle')
        A = single(A);
    else
        A = double(A);
    end
end
end % convnfft


function cmap = crameri(ColormapName,varargin) 
% crameri returns perceptually-uniform scientific colormaps created
% by Fabio Crameri. 
% 
% Syntax 
% 
%  crameri 
%  cmap = crameri('ColormapName') 
%  cmap = crameri('-ColormapName') 
%  cmap = crameri(...,NLevels)
%  cmap = crameri(...,'pivot',PivotValue) 
%  crameri(...)
% 
% Description 
% 
% crameri without any inputs displays the options for colormaps. 
% 
% cmap = crameri('ColormapName') returns a 256x3 colormap.  For a visual
% depiction of valid colormap names, type |crameri|. 
%
% cmap = crameri('-ColormapName') a minus sign preceeding any ColormapName flips the
% order of the colormap. 
%
% cmap = crameri(...,NLevels) specifies a number of levels in the colormap.  Default
% value is 256. 
%
% cmap = crameri(...,'pivot',PivotValue) centers a diverging colormap such that white 
% corresponds to a given value and maximum extents are set using current caxis limits. 
% If no PivotValue is set, 0 is assumed. 
%
% crameri(...) without any outputs sets the current colormap to the current axes.  
% 
% Examples 
% For examples, type: 
% 
%  showdemo crameri_documentation
%
% Author Info 
% This function was written by Chad A. Greene of the University of Texas
% Institute for Geophysics (UTIG), August 2018, using Fabio Crameri's 
% scientific colormaps, version 4.0. http://www.fabiocrameri.ch/colourmaps.php
% 
% Citing this colormap: 
% Please acknowledge the free use of these colormaps by citing
% 
% Crameri, F. (2018). Scientific colour-maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
% 
% Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation and 
% StagLab 3.0, Geosci. Model Dev., 11, 2541-2562, doi:10.5194/gmd-11-2541-2018.
% 
% For more on choosing effective and accurate colormaps for science, be sure
% to enjoy this fine beach reading: 
% 
% Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True 
% colors of oceanography: Guidelines for effective and accurate colormap selection. 
% Oceanography 29(3):9-13, http://dx.doi.org/10.5670/oceanog.2016.66.
% 
% See also colormap and caxis.
%
% Display colormap options: 
if nargin==0
   figure('menubar','none','numbertitle','off','Name','crameri options:')
   
   if license('test','image_toolbox')
      imshow(imread('crameri7.0.png')); 
   else
      axes('pos',[0 0 1 1])
      image(imread('crameri7.0.png')); 
      axis image off
   end
   
   return
end
% Error checks: 
assert(isnumeric(ColormapName)==0,'Input error: ColormapName must be a string.') 
% Set defaults: 
NLevels = 256; 
autopivot = false; 
PivotValue = 0; 
InvertedColormap = false; 
% Parse inputs: 
% Does user want to flip the colormap direction? 
dash = strncmp(ColormapName,'-',1); 
if any(dash) 
   InvertedColormap = true; 
   ColormapName(dash) = []; 
end
% Standardize all colormap names to lowercase: 
ColormapName = lower(ColormapName); 
% Oleron's too hard for me to remember, so I'm gonna use dem or topo. 
if ismember(ColormapName,{'dem','topo'})
   ColormapName = 'oleron'; 
end
% Does the user want to center a diverging colormap on a specific value? 
% This parsing support original 'zero' syntax and current 'pivot' syntax. 
tmp = strncmpi(varargin,'pivot',3); 
if any(tmp) 
   autopivot = true; 
   try
      if isscalar(varargin{find(tmp)+1})
         PivotValue = varargin{find(tmp)+1}; 
         tmp(find(tmp)+1) = 1; 
      end
   end
   varargin = varargin(~tmp); 
end
% Has user requested a specific number of levels? 
tmp = isscalar(varargin); 
if any(tmp) 
   NLevels = varargin{tmp}; 
end
% Load RGB values and interpolate to NLevels: 
try
   S = load('CrameriColourMaps7.0.mat',ColormapName); 
   cmap = S.(ColormapName); 
catch
   error(['Unknown colormap name ''',ColormapName,'''. Try typing crameri with no inputs to check the options and try again.'])
end
% Interpolate if necessary: 
if NLevels~=size(cmap,1) 
   cmap = interp1(1:size(cmap,1), cmap, linspace(1,size(cmap,1),NLevels),'linear');
end
% Invert the colormap if requested by user: 
if InvertedColormap
   cmap = flipud(cmap); 
end
% Adjust values to current caxis limits? 
if autopivot
   clim = caxis; 
   maxval = max(abs(clim-PivotValue)); 
   cmap = interp1(linspace(-maxval,maxval,size(cmap,1))+PivotValue, cmap, linspace(clim(1),clim(2),size(cmap,1)),'linear');
end
% Clean up 
if nargout==0
   colormap(gca,cmap) 
   clear cmap  
end

% Code to collect Fabio's data into a single .mat file: 
% Unzip the latest folder, navigate to that filepath, and run this.
% Update the file list as needed. 
%
% clear all
% f = {'acton','bam','bamO','bamako','batlow','batlowK','batlowW','berlin','bilbao','broc','brocO','buda','bukavu','cork',...
%    'corkO','davos','devon','fes','grayC','hawaii','imola','lajolla','lapaz','lisbon',...
%    'nuuk','oleron','oslo','roma','romaO','tofino','tokyo','turku','vanimo','vik','vikO'}; 
% 
% for k = 1:length(f)
%    load([f{k},'/',f{k},'.mat'])
% end
% 
% clear f k 
% save('CrameriColourMaps7.0.mat')
end


function value = getoption(options, name, defaultvalue)

% Get default option in convnfft
% function value = getoption(options, name, defaultvalue)

    value = defaultvalue;
    fields = fieldnames(options);
    found = strcmpi(name,fields);
    if any(found)
        i = find(found,1,'first');
        if ~isempty(options.(fields{i}))
            value = options.(fields{i});
        end
    end
end


function y=h2d(x)

% Kori-ULB
% Interpolate values initialiiy on d-grid to h-grid

    y=0.25*(x+circshift(x,[0 -1])+circshift(x,[-1 -1])+circshift(x,[-1 0]));
    
end


function h = imagescn(varargin) 

% imagescn behaves just like imagesc, but makes NaNs transparent, sets
% axis to xy (aka ydirection normal) if xdata and ydata are included, and has a little more 
% error checking than imagesc. 
% 
% Syntax 
% 
%  imagescn(C) 
%  imagescn(x,y,C) 
%  imagescn(x,y,C,clims) 
%  imagescn('PropertyName',PropertyValue,...) 
%  h = imagescn(...) 
% 
% Description 
% 
% imagescn(C) displays the data in array C as an image that uses the full range of colors in the colormap. 
% Each element of C specifies the color for 1 pixel of the image. The resulting image is an m-by-n grid of 
% pixels where m is the number of columns and n is the number of rows in C. The row and column indices of 
% the elements determine the centers of the corresponding pixels. NaN values in C appear transparent. 
% 
% imagescn(x,y,C) specifies x and y locations of the centers of the pixels in C. If x and y are two-element
% vectors, the outside rows and columns of C are centered on the values in x and y. Mimicking imagesc, if 
% x or y are vectors with more than two elements, only the first and last elements of of the vectors are 
% considered, and spacing is automatically scaled as if you entered two-element arrays. The imagescn function
% takes this one step further, and allows you to enter x and y as 2D grids the same size as C. If x and y
% are included, the imagescn function automatically sets axes to cartesian xy rather than the (reverse) ij axes. 
% 
% imagescn(x,y,C,clims) specifies the data values that map to the first and last elements of the colormap. 
% Specify clims as a two-element vector of the form [cmin cmax], where values less than or equal to cmin 
% map to the first color in the colormap and values greater than or equal to cmax map to the last color in 
% the colormap.
% 
% imagescn('PropertyName',PropertyValue,...) specifies image properties as name-value pairs. 
% 
% h = imagescn(...) retrns a handle of the object created. 
% 
% Differences between imagesc, imagescn, and pcolor
% The imagescn function plots data with imagesc, but after plotting, sets NaN pixels to an 
% alpha value of 0. The imagesc function allows input coordinates x and y to be grids, which 
% are assumed to be evenly-spaced and monotonic as if created by meshgrid. If x and y data 
% are included when calling imagescn, y axis direction is changed from reverse to normal. 
% 
% The imagescn function is faster than pcolor. Pcolor (nonsensically) deletes an outside row 
% and column of data, and pcolor also refuses to plot data points closest to any NaN holes. 
% The imagescn function does not delete any data.  However, you may still sometimes wish to 
% use pcolor if x,y coordinates are not evenly spaced or if you want interpolated shading. 
% 
% Examples 
% For examples, type 
% 
%  cdt imagescn
% 
% Author Info 
% 
% This function was written by Chad A. Greene of the University of 
% Texas Institute for Geophysics (UTIG), January 2017. 
% http://www.chadagreene.com 
% 
% See also imagesc, image, and pcolor.
% The imagesc function does not have error checking regarding number of elements
% in xdata, ydata versus number of elements in the input image, so I'm gonna add
% some error checking: 

% Check inputs: 

xydata = false; 
if nargin>2
   if all([isnumeric(varargin{1}) isnumeric(varargin{2})]) 
      % This is an assumption that should typically be safe to make: 
      xydata = true; 
      
      % Determine if input coordinates are meshgrid type and if so, convert to vector: 
      if isequal(size(varargin{1}),size(varargin{2}),size(varargin{3}))
         X = varargin{1}; 
         Y = varargin{2}; 
         
         varargin{1} = [X(1,1) X(end,end)]; 
         varargin{2} = [Y(1,1) Y(end,end)]; 
      end
   end
end

% Plot

% Plot imagesc: 
h = imagesc(varargin{:}); 

% Make NaNs transparent: 
cd = get(h,'CData'); 
set(h,'alphadata',isfinite(cd)); 

if xydata
   axis xy
end

if nargout==0
   clear h
end

end


function u=vec2h(ux,uy) %VL: velocity on h-grid

    ux1=circshift(ux,[0 1]); % ux(i,j-1)
    uy1=circshift(uy,[1 0]); % uy(i-1,j)
    u=sqrt((0.5*(ux+ux1)).^2+(0.5*(uy+uy1)).^2);
    
end

