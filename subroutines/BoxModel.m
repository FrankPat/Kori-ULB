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


