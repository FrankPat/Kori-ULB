function [dtr]=TransportDamage(node,nodes,dtr,Mb,Melt,ThinComp,H,MASK,dtdx,dtdx2, ...
                                     u,v,ctr,cnt,bMASK,VM,par)

    % Kori-ULB
    % Sparse solver of damage transport (based on ice thickness solver)
    % Damage transport based on Kachuck et al., (2022) -> Based on Bassis & Ma (2015)
    % DOI:0.1017/jog.2022.12

    epsilon=1e-10; % artificial diffusion
    MASK(MASK==6)=0;

    um1=circshift(u,[0 1]);        % u(i,j-1)
    vm1=circshift(v,[1 0]);        % v(i-1,j)

    % 
    V0=8*epsilon*dtdx+dtdx2*(u-um1+v-vm1)+(max(ThinComp,0)+max(Melt,0)+max(Mb,0))*ctr.dt./max(H,1e-5); % i,j
    V1=-2*epsilon*dtdx+dtdx2*u;   % i,j+1
    V2=-2*epsilon*dtdx-dtdx2*um1; % i,j-1
    V3=-2*epsilon*dtdx+dtdx2*v;   % i+1,j
    V4=-2*epsilon*dtdx-dtdx2*vm1; % i-1,j

    R0=dtr;

    V0(MASK==0)=0; % note that for shelf=1, MASK=glMASK in the call
    V1(MASK==0)=0;
    V2(MASK==0)=0;
    V3(MASK==0)=0;
    V4(MASK==0)=0;
    
    R0(MASK==0)=dtr(MASK==0);

    % boundaries
    V9=zeros(ctr.imax,ctr.jmax);  % ice divide or ocean
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
    if wholemask==0
        V0(MASKb==1)=0;
        R0(MASKb==1)=0;
    end

    V=[reshape(V0(VM==1)+1,nodes,1)
       V1(V1~=0)
       V2(V2~=0)
       V3(V3~=0)
       V4(V4~=0)
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

    nodeV1=circshift(node,[0 -1]);  % i,j+1
    nodeV2=circshift(node,[0 1]);   % i,j-1
    nodeV3=circshift(node,[-1 0]);  % i+1,j
    nodeV4=circshift(node,[1 0]);   % i-1,j
    nodeV9=circshift(node,[0 -1]);  % i,2 - ocean
    nodeV10=circshift(node,[-1 0]); % 2,j - ocean
    nodeV11=circshift(node,[1 0]);  % imax-1,j - ocean
    nodeV12=circshift(node,[0 1]);  % i,jmax-1 - ocean

    col=[reshape(node(VM==1),nodes,1)
         nodeV1(V1~=0)
         nodeV2(V2~=0)
         nodeV3(V3~=0)
         nodeV4(V4~=0)
         nodeV9(V9~=0)
         nodeV10(V10~=0)
         nodeV11(V11~=0)
         nodeV12(V12~=0)];

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
