function [snapshot,plotst,cnt_atm,snp_atm,cnt_ocn,snp_ocn,Mb_update,Li,Lj, ...
    dtdx,dtdx2,X,Y,x,y,MASK,H,Ho,B,Bo,MASKo,Mb,Ts,As,G,u,VAF,VA0,POV, ...
    SLC,Ag,Af,Btau,IVg,IVf,glflux,cfflux,dHdt,time,mbcomp,InvVol,ncor,dSLR,SLR, ...
    Wd,Wtil,Bmelt,NumStab,CMB,FMB,flw,p,px,py,pxy,nodeu,nodev,nodes,node, ...
    VM,Tof,Sof,TFf,Tsf,Mbf,Prf,Evpf,runofff,Melt,damage,shelftune,Melt_mean, ... 
    Bmelt_mean,Ts_mean,Mb_mean,To_mean,So_mean,TF_mean,CMB_mean,FMB_mean, ...
    fluxmx_mean,fluxmy_mean]= ...
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
        runofff,Melt,damage,Melt_mean, ... 
    Bmelt_mean,Ts_mean,Mb_mean,To_mean,So_mean,TF_mean,CMB_mean,FMB_mean, ...
    fluxmx_mean,fluxmy_mean]=deal(zeros(ctr.imax,ctr.jmax));
    Mb=Mb+ctr.MbConst;
    Ts=Ts+ctr.TsConst;
    [MASK,MASKo,ncor]=deal(ones(ctr.imax,ctr.jmax));
    [VAF,VA0,POV,SLC,Ag,Af,IVg,IVf,glflux,cfflux,dHdt]=deal(zeros(ctr.nsteps,1));
    G=zeros(ctr.imax,ctr.jmax)+default.G0;
    
    if ctr.NumCheck==1
        NumStab=zeros(ctr.nsteps,8);
    else
        NumStab=false;
    end
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
    if length(ctr.Asin)>1
        As=ctr.Asin;
    else
        As=zeros(ctr.imax,ctr.jmax)+ctr.Asin;
    end
    if length(ctr.shelftune)>1
        shelftune=ctr.shelftune;
    else
        shelftune=zeros(ctr.imax,ctr.jmax)+ctr.shelftune;
    end
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


