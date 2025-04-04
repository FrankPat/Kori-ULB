function [H,flag,relres,iter]=SparseSolverIceThickness(node,nodes,Mb, ...
    H,B,SLR,MASK,dtdx,dtdx2,d,u,v,ctr,cnt,bMASK,VM,par)

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
        [flag,relres,iter]=deal(false);
    else
        D=diag(diag(A));
        C1=tril(A);
        C2=D\triu(A);
        [s,flag,relres,iter]=pcg(A,R,par.Htol,par.Hiter,C1,C2);
        if flag>0 || cnt==1
            s=A\R;
        end
    end

    H(node>0)=s(node(node>0));

end


