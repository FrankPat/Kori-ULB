function [ub,ubx,uby]=BasalVelocitySSA(uxssa,uyssa,beta2,etaD,H,zeta,MASK,ctr)

% Kori-ULB
% Basal velocity according to the SSA, hybrid and DIVA models.

    if ctr.SSA==1 || ctr.SSA==2
        ubx=uxssa;
        uby=uyssa;
    elseif ctr.SSA==3
        betax=0.5*(beta2+circshift(beta2,[0 -1]));
        betay=0.5*(beta2+circshift(beta2,[-1 0]));
        [~,F2x,F2y]=Fint(ctr,etaD,H,zeta);
        ubx=uxssa./(1+betax.*F2x);
        uby=uyssa./(1+betay.*F2y);
    end
    ub=vec2h(ubx,uby); % only for temperature (on h grid)
    ub(MASK==0)=0;

end


