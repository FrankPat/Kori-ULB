function LSF=LSFfunction(LSF,ctr,u,v,node,nodes,VM,MASK)

% Kori-ULB
% Calculate the Level Set Function (LSF) for following the calving front.
% Used in the calving algorihms
% Still under development

    epsilon=1e-2;
    
    dtdx2=epsilon*ctr.dt/(ctr.delta*ctr.delta);
    dtdx=ctr.dt/ctr.delta;
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
    MU(u<=0)=2;
    MV=ones(ctr.imax,ctr.jmax);
    MV(v<=0)=2;
    MU(MASKb==0)=0;
    MV(MASKb==0)=0;

    V0a=zeros(ctr.imax,ctr.jmax);
    V1a=zeros(ctr.imax,ctr.jmax);
    V2a=zeros(ctr.imax,ctr.jmax);
    V3a=zeros(ctr.imax,ctr.jmax);
    V4a=zeros(ctr.imax,ctr.jmax);

    % conditions for MU=1 (grad(u)>0)
    V0a(MU==1)=u(MU==1); % i,j
    V2a(MU==1)=-u(MU==1); % i,j-1

    % conditions for MU=2 and (grad(u)<0)
    V0a(MU==2)=-u(MU==2); % i,j
    V1a(MU==2)=u(MU==2); % i,j+1

    % conditions for MV=1 (grad(v)>0)
    V0a(MV==1)=V0a(MV==1)+v(MV==1); % i,j
    V4a(MV==1)=-v(MV==1); % i-1,j

    % conditions for MV=2 (grad(v)<0)
    V0a(MV==2)=V0a(MV==2)-v(MV==2); % i,j
    V3a(MV==2)=v(MV==2); % i+1,j

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
    
%     % Avoid numerical issues when calving front coincides with grounding line.
%     % Allow for a couple of grid cells of calving front between GL and open ocean.
%     M1 = circshift(MASK,[3 3]);
%     M2 = circshift(MASK,[3 -3]);
%     M3 = circshift(MASK,[-3 3]);
%     M4 = circshift(MASK,[-3 -3]);
%     a = (MASK==1)|(M1==1)|(M2==1)|(M3==1)|(M4==1);
%     M(a) = 1;
% 
%     % Calving front cannot retreat further than the GL by definition.
%     LSF(M==1) = R0(M==1);

end


