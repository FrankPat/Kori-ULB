function [uxsch,uysch,d]=DiffusiveCorrection(ctr,par,uxsch,uysch,udx,udy, ...
    d,Ad,p,Hm,gradm,bMASK)

% Kori-ULB
% Correction of diffusion coefficients in ice thickness equation for
% grounding line flux condition.
% This may become obsolete, as diffusion coefficients are not used for SSA

    if ctr.SSA>0
        uxsch=min(max(uxsch,-par.maxspeed),par.maxspeed);   %LZ2021
        uysch=min(max(uysch,-par.maxspeed),par.maxspeed);   %LZ2021
        if ctr.SSAdiffus>0 && ctr.SSA>1
            uxsch=uxsch-udx;
            uysch=uysch-udy;
            d=par.Z*Ad./(p+2).*Hm.^(par.n+2.).*gradm.^((par.n-1.)/2.);
        else
            if ctr.basin==1
                d(bMASK==0)=0;
            else
                d=zeros(ctr.imax,ctr.jmax);
            end
        end
    else % put Schoof/Tsai in d and set ux,uy zero for grounded part
        if ctr.schoof>=1
            dusch=velocity_on_grid(uxsch,uysch);
            d(gradm>1e-20)=dusch(gradm>1e-20).*Hm(gradm>1e-20)./ ...
                sqrt(gradm(gradm>1e-20));
        end
        uxsch=zeros(ctr.imax,ctr.jmax);
        uysch=zeros(ctr.imax,ctr.jmax);
    end
    
end
    

