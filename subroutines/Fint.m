function [F2,F2x,F2y]=Fint(ctr,etaD,H,zeta)
    
% Kori-ULB
% Numerical integration of effective viscosity and F terms for the DIVA model

    % Preallocate arrays
    F2=zeros(ctr.imax,ctr.jmax);

    % Vertical integration
    for k=2:ctr.kmax
        etalayer=max(1e-8,(etaD(:,:,k)+etaD(:,:,k-1))/2);
        zetalayer=(zeta(k)+zeta(k-1))/2;
        F2=F2+zetalayer^2.*(zeta(k)-zeta(k-1))./etalayer;
    end
    F2(isnan(F2))=0;
    F2=F2.*H;
    % Stagger since it is used to compute velocity and beta.
    F2x=0.5*(F2+circshift(F2,[0 -1]));    % (i,j+1)
    F2y=0.5*(F2+circshift(F2,[-1 0]));    % (i+1,j)

end


