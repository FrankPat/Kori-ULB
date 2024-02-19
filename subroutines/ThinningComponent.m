function [ThinComp] = ThinningComponent(ctr,par,dudx,dvdy,dudy,dvdx,eta,H) 
    % Thinning component of basal crevasses.
    % Term comes from Bassis and Ma (2015) and is applied in Kachuck et al., (2022)
    % Applies to long wavelenghts: i.e. the widths of basal crevasses are large compared to the ice thickness

    eps=1e-8; % avoid zero values

    % Initialize to zeros
    ThinComp = zeros(ctr.imax,ctr.jmax); 
    % 1st/2nd principal strain
    [lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); 
    tau1=2*lambda1.*eta./(H+eps);
    % Thinning component
    alpha    = lambda2./lambda1;
    n_star   = (4*par.n*(1+alpha+alpha.*alpha))./(4*(1+alpha+alpha.*alpha)+3*(par.n-1)*alpha.*alpha);
    So       = (par.rho*(par.rhow-par.rho)*par.g)./(2*tau1.*par.rhow);
    ThinComp = max(H,1e-5).*n_star.*(1-So).*lambda1;

    for i=1:ctr.imax
            for j=1:ctr.jmax
                if tau1(i,j)<par.tauice
                    ThinComp(i,j)=0.0;
        end
    end

end
