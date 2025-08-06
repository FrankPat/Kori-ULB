function [u,v,s,flag,relres,iter]=SparseSolverSSA(nodeu,nodev,s0, ...
    MASKmx,MASKmy,bMASK,glMASK, ...
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
    
    magh=max(sqrt(u.^2+v.^2),1e-10); % velocity magnitude on h-grid
    angle=atan2(v./magh,u./magh);
    nx=cos(angle);
    ny=sin(angle);
    
    % Define masks for the ocean/ice boundary conditions
    
    BCMASK=zeros(ctr.imax,ctr.jmax);
    if ctr.shelfBC==1
        BCMASK(glMASK==4 & circshift(glMASK,[0 -1])==5 & nx>0.5)=1;
        BCMASK(glMASK==5 & circshift(glMASK,[0 1])==6 & nx<-0.5)=2;
    end
    if ctr.mismip>0
        BCMASK(1,:)=0;
        BCMASK(ctr.imax,:)=0;
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
        R0(MASKb==1)=1.5*0.5*par.rho*par.g*H1(MASKb==1).^2.*(1.- ...
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
        R0(MASKb==1)=1.5*0.5*par.rho*par.g*H(MASKb==1).^2.*(1.- ...
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
        R0(MASKb==1)=0;

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
        R0(MASKb==1)=1.5*0.5*par.rho*par.g*(1.-par.rho/par.rhow)*H(MASKb==1).^2;

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
        R0(MASKb==1)=0;

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
        R0(MASKb==1)=0;
    end
    
    % Boundaries for the calving front (implies overriding original values)
    R0(BCMASK>0)=0;
    U0(BCMASK>0)=0;
    U1(BCMASK>0)=0;
    U2(BCMASK>0)=0;
    U3(BCMASK>0)=0;
    U4(BCMASK>0)=0;
    U5(BCMASK>0)=0;
    U6(BCMASK>0)=0;
    U7(BCMASK>0)=0;
    U8(BCMASK>0)=0;
    U9(BCMASK>0)=0;
    U10(BCMASK>0)=0;
    Uv0(BCMASK>0)=0;
    Uv1(BCMASK>0)=0;
    Uv2(BCMASK>0)=0;
    Uv3(BCMASK>0)=0;
    
    % ice-ocean boundary (BCMASK=1)
    U0(BCMASK==1)=4.*eta(BCMASK==1)/ctr.delta;
    U2(BCMASK==1)=-4.*eta(BCMASK==1)/ctr.delta;
    Uv0(BCMASK==1)=2.*eta(BCMASK==1)/ctr.delta;
    Uv2(BCMASK==1)=-2.*eta(BCMASK==1)/ctr.delta;
    R0(BCMASK==1)=1.5*0.5*par.rho*par.g*(1.-par.rho/par.rhow)*H(BCMASK==1).^2;
    
    % ocean-ice boundary (BCMASK=2)
    U0(BCMASK==2)=-4.*eta1(BCMASK==2)/ctr.delta;
    U1(BCMASK==2)=4.*eta1(BCMASK==2)/ctr.delta;
    Uv1(BCMASK==2)=2.*eta(BCMASK==2)/ctr.delta;
    Uv3(BCMASK==2)=-2.*eta(BCMASK==2)/ctr.delta;
    R0(BCMASK==2)=1.5*0.5*par.rho*par.g*(1.-par.rho/par.rhow)*H1(BCMASK==2).^2;
    
    

    % v-velocities: quadrant V and Vu of solution matrix

    % Define masks for the ocean/ice boundary conditions   
    BCMASK=zeros(ctr.imax,ctr.jmax);
    if ctr.shelfBC==1
        BCMASK(glMASK==4 & circshift(glMASK,[-1 0])==5 & ny>0.5)=1;
        BCMASK(glMASK==5 & circshift(glMASK,[1 0])==6 & ny<-0.5)=2;
    end
    if ctr.mismip>0
        BCMASK(1,:)=0;
        BCMASK(ctr.imax,:)=0;
    end

    
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
        S0(MASKb==1)=1.5*0.5*par.rho*par.g*H1(MASKb==1).^2.*(1.-par.rho/par.rhow);

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
        S0(MASKb==1)=1.5*0.5*par.rho*par.g*H(MASKb==1).^2.*(1.-par.rho/par.rhow);

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
        S0(MASKb==1)=0;

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
            S0(MASKb==1)=1.5*0.5*par.rho*par.g*H(MASKb==1).^2.*(1.- ...
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
        S0(MASKb==1)=0;

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
        S0(MASKb==1)=0;
    end

    % Boundaries for the calving front (implies overriding original values)
    S0(BCMASK>0)=0;
    V0(BCMASK>0)=0;
    V1(BCMASK>0)=0;
    V2(BCMASK>0)=0;
    V3(BCMASK>0)=0;
    V4(BCMASK>0)=0;
    V5(BCMASK>0)=0;
    V6(BCMASK>0)=0;
    V7(BCMASK>0)=0;
    V8(BCMASK>0)=0;
    V9(BCMASK>0)=0;
    V10(BCMASK>0)=0;
    V11(BCMASK>0)=0;
    Vu0(BCMASK>0)=0;
    Vu1(BCMASK>0)=0;
    Vu2(BCMASK>0)=0;
    Vu3(BCMASK>0)=0;
    
    % ice-ocean boundary (BCMASK=1)
    V0(BCMASK==1)=4.*eta(BCMASK==1)/ctr.delta;
    V2(BCMASK==1)=-4.*eta(BCMASK==1)/ctr.delta;
    Vu0(BCMASK==1)=2.*eta(BCMASK==1)/ctr.delta;
    Vu2(BCMASK==1)=-2.*eta(BCMASK==1)/ctr.delta;
    S0(BCMASK==1)=1.5*0.5*par.rho*par.g*(1.-par.rho/par.rhow)*H(BCMASK==1).^2;
    
    % ocean-ice boundary (BCMASK=2)
    V0(BCMASK==2)=-4.*eta1(BCMASK==2)/ctr.delta;
    V1(BCMASK==2)=4.*eta1(BCMASK==2)/ctr.delta;
    Vu1(BCMASK==2)=2.*eta(BCMASK==2)/ctr.delta;
    Vu3(BCMASK==2)=-2.*eta(BCMASK==2)/ctr.delta;
    S0(BCMASK==2)=1.5*0.5*par.rho*par.g*(1.-par.rho/par.rhow)*H1(BCMASK==2).^2;

    
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
        [s,flag,relres,iter]=bicgstab(A,R,par.veltol,par.veliter,C1,C2,s0);
        if flag>0
            s=A\R;
        end
    else
        s=A\R;
        [flag,relres,iter]=deal(false);
    end
    u=s(nodeu);
    v=s(nodev);

end


