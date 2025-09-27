function [eta,dudx,dvdy,dudy,dvdx]=EffVisc(A,uxssa,uyssa,H,H0,par,MASK, ...
    glMASK,shelftune,damage,ctr)

% Kori-ULB
% Effective viscosity of the SSA solution. On the borders of the domain, a
% fixed value is determined either based on the mean viscosity of the ice
% shelf or a specific value for the ocean.

    eps=1e-8;
    [dudx,dudy,dvdx,dvdy]=VelocityGradients(uxssa,uyssa,ctr);

    EffStr=dudx.^2+dvdy.^2+dudx.*dvdy+0.25*(dudy+dvdx).^2;
    EffStr=max(EffStr,1e-12);
    scale_eta=(H-min(min(H.*ctr.damlim,damage),H-eps))./(H+eps);
    if ctr.damexist==1 && ctr.damage==0
        scale_eta=(H0-min(min(H0.*ctr.damlim,damage),H0-eps))./(H0+eps);
    end
    Astar=A.*scale_eta.^(-par.n);
    eta=0.5*H.*Astar.^(-1./par.n).*EffStr.^((1-par.n)/(2*par.n));
    eta(MASK==0)=eta(MASK==0)./shelftune(MASK==0);  %VL: 2D shelftune
    
    if ctr.shelf==1 || ctr.schoof>0
        MASKb=ones(ctr.imax,ctr.jmax); % use constant eta on edges of ice shelf
        if ctr.mismip==0
            MASKb(3:ctr.imax-2,3:ctr.jmax-2)=0;
        elseif ctr.mismip==1
            MASKb(:,1:ctr.jmax-2)=0;
        elseif ctr.mismip==2
            MASKb(1:ctr.imax-2,1:ctr.jmax-2)=0;
        end
        MASKb(MASK==1)=0;
        eta(MASKb==1)=trimmean(eta(MASKb==1),10);
        
        % Instead of calculating effective viscosity on the sea ice,
        % keep constant viscosity. Need to further check how to deal
        % with this. May be quoted
        eta(glMASK==6)=ctr.OceanVisc; % Default: 1.0e7. Daniel: 8.0e9
        % 1e8 is approx for ice shelf of 200m thick
    end
    
end


