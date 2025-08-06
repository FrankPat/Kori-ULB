function ThinComp=ThinningComponent(ctr,par,dudx,dvdy,dudy,dvdx,eta,H, ...
    damage,bMASK,glMASK)

% Kori-ULB
% Thinning component of basal crevasses.
% Term comes from Bassis and Ma (2015) and is applied in Kachuck et al., (2022)
% Applies to long wavelenghts: i.e. the widths of basal crevasses are large compared to the ice thickness

    eps=1e-8; % avoid zero values

    % 1st/2nd principal strain
    [lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); 
    tau1=2*lambda1.*eta./(H+eps);

    % Thinning component
    alpha=lambda2./lambda1;
    n_star=(4*par.n*(1+alpha+alpha.*alpha))./(4*(1+alpha+alpha.*alpha)+3*(par.n-1)*alpha.*alpha);
    So=(par.rho*(par.rhow-par.rho)*par.g.*H)./(2*tau1.*par.rhow);
    ThinComp=max(H,1e-5).*n_star.*(1-So).*lambda1;

    % no thinning if no local nor transported damage
    ThinComp(damage==0 | tau1<ctr.tauice)=0;
    if ctr.basin==1
        ThinComp(bMASK==1)=0;
    end
    ThinComp(glMASK==6)=0;
end


