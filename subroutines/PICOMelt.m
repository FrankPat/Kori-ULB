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


