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


