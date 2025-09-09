function [uxsch,uysch,d]=DiffusiveCorrection(ctr,par,uxsch,uysch,d,bMASK)

% Kori-ULB
% Correction of diffusion coefficients in ice thickness equation for
% grounding line flux condition.
% This may become obsolete, as diffusion coefficients are not used for SSA

    if ctr.SSA>0
        uxsch=min(max(uxsch,-par.maxspeed),par.maxspeed);   %LZ2021
        uysch=min(max(uysch,-par.maxspeed),par.maxspeed);   %LZ2021
        if ctr.basin==1
            % use SIA and diffusion scheme outside basin
            d(bMASK==0)=0;
        else
            d=zeros(ctr.imax,ctr.jmax);
        end
    else % set ux,uy zero for grounded part, only diffusion for SIA
        uxsch=zeros(ctr.imax,ctr.jmax);
        uysch=zeros(ctr.imax,ctr.jmax);
    end
    
end
    

