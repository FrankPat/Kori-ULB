function [Ntil,Wtil]=TillWater(MASK,Wtil,Bmelt,ctr,Po,par)

% Kori-ULB
% Effective pressure in till based on subglacial water saturation

    Wtil=max(par.Wdmin,min(Wtil+(Bmelt-par.Cdr)*ctr.dt,par.Wmax));
    s=Wtil/par.Wmax;
    s(MASK==0)=1;
    Ntil=par.N0*((par.sigmat*Po/par.N0).^s).*10.^((par.e0/par.Cc).*(1.-s));
    Ntil=max(min(Po,Ntil),par.sigmat*Po);

end


