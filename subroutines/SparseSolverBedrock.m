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


