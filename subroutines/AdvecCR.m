function [CRadvec]=AdvecCR(CR,H,glMASK,MASK,vx,vy,ctr,par)

% Kori-ULB
% Advection of calving rate into ocean 
% Employs an upstream difference scheme
% based on advection scheme from Pelle (for PICOP and Plume)

    epsilon=1e-5; % diffusion coefficient for stabilizing advection scheme

    % ocean + last ice point (calving front or last grounded)
    MASKadvec=glMASK*0;
    MASKadvec(glMASK==6)=1;
    % Ice areas
    MASKice=MASK;
    MASKice(MASK==0 & H>par.SeaIceThickness)=1;
    % Last ice point
    MASK1=circshift(MASKice,[0 -1]);
    MASK2=circshift(MASKice,[0 1]);
    MASK3=circshift(MASKice,[-1 0]);
    MASK4=circshift(MASKice,[1 0]);
    MASKadvec(MASKice==1 & (MASK1==0 | MASK2==0 | MASK3==0 | MASK4==0))=1;
    % domain edge + last ice point
    MASKb=MASK.*0;
    MASKb(MASKice==1 & (MASK1==0 | MASK2==0 | MASK3==0 | MASK4==0))=1;
    MASKb(1:end,1:2)=1;
    MASKb(1:2,1:end)=1;
    MASKb(1:end,end-1:end)=1;
    MASKb(end-1:end,1:end)=1;

    VM=reshape(MASKadvec,ctr.imax*ctr.jmax,1);
    nodes=sum(sum(MASKadvec));
    node=zeros(ctr.imax,ctr.jmax);
    node(MASKadvec>0)=linspace(1,nodes,nodes);

    CRadvec=CR;
%     CRadvec=min(CRadvec,0);
%     CRadvec(MASK==0)=0;

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

    V0(MASKadvec==0)=0;
    V1(MASKadvec==0)=0;
    V2(MASKadvec==0)=0;
    V3(MASKadvec==0)=0;
    V4(MASKadvec==0)=0;

    R0=zeros(ctr.imax,ctr.jmax);

    V0(MASKb==1)=1;
    V1(MASKb==1)=0;
    V2(MASKb==1)=0;
    V3(MASKb==1)=0;
    V4(MASKb==1)=0;
    R0(MASKb==1)=CRadvec(MASKb==1);

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
    CRadvec(node>0)=s(node(node>0));
    CRadvec=max(CRadvec,0);
end


