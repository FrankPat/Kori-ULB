function [d,udx,udy,ud,ubx,uby,ub,uxsia,uysia,p,pxy]= ...
    SIAvelocity(ctr,par,A,Ad,Ax,Ay,Asfd,Asfx,Asfy,taud,G,Tb, ...
    H,Hm,Hmx,Hmy,gradm,gradmx,gradmy,gradxy,signx,signy,MASK,p,px,py,pxy)

% Kori-ULB
% Deformational and basal velocity according to the shallow-ice
% approximation. Parameters for the thermomechanical coupling are also
% defined based on basal temperature gradients.

    if ctr.Tcalc==2
        p=par.n-1+par.Q1*G.*Hm./(par.R.*(h2d(Tb)+par.T0).^2);
        % on staggered d-grid
        pxy=par.n-1+par.Q1*G.*H./(par.R.*(Tb+par.T0).^2); % on h-grid     
        px=par.n-1+par.Q1*G.*Hmx./(par.R.*(0.5* ...
            (Tb+circshift(Tb,[0 -1]))+par.T0).^2); % on h-grid     
        py=par.n-1+par.Q1*G.*Hmy./(par.R.*(0.5* ...
            (Tb+circshift(Tb,[-1 0]))+par.T0).^2); % on h-grid
    end
    if ctr.u0>1e10
        ds=1; % disable correction when u0 is default value
    else
        ds=max(0.01,min(1-Asfd./ctr.u0.*taud.^ctr.m,1)); 
        % Zoet correction factor %VL: use Asfd
    end
    d=par.Z*Ad./(p+2).*Hm.^(par.n+2.).*gradm.^((par.n-1.)/2.)+ ...
        Asfd*(par.rho*par.g)^ctr.m.*Hm.^(ctr.m+1.).* ...
        gradm.^((ctr.m-1.)/2.)./ds;
    udx=par.Z*Ax./(px+2).*signx.*Hmx.^(par.n+1.).*gradmx.^(par.n/2.);
    udy=par.Z*Ay./(py+2).*signy.*Hmy.^(par.n+1.).*gradmy.^(par.n/2.);
    ud=par.Z*A./(pxy+2).*H.^(par.n+1.).*gradxy.^(par.n/2.);
    ubx=Asfx.*signx.*(par.rho*par.g*Hmx).^ctr.m.*gradmx.^(ctr.m/2.); 
    uby=Asfy.*signy.*(par.rho*par.g*Hmy).^ctr.m.*gradmy.^(ctr.m/2.);
    uxsia=min(max(-par.maxspeed,udx+ubx),par.maxspeed);
    uysia=min(max(-par.maxspeed,udy+uby),par.maxspeed);
    ub=vec2h(ubx,uby); % only for temperature (on h grid)
    ub(MASK==0)=0;
    
end


