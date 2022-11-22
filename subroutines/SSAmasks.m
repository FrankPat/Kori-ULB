function [MASKmx,MASKmy,HAFmx,HAFmy]=SSAmasks(ctr,par,SLR,Bmx, ...
    Bmy,Hmx,Hmy,bMASKx,bMASKy)

% Kori-ULB
% Ice shelf masks on staggered grids for SSA computation

    HAFmx=Bmx-SLR+Hmx*par.rho/par.rhow;
    HAFmy=Bmy-SLR+Hmy*par.rho/par.rhow;
    MASKmx=zeros(ctr.imax,ctr.jmax);
    MASKmy=zeros(ctr.imax,ctr.jmax);
    if ctr.shelf==1 || ctr.SSA>=1
        if ctr.shelf==1
            MASKmx(HAFmx<0)=1; % floating
            MASKmy(HAFmy<0)=1;
        end
        if ctr.SSA>=1
            MASKmx(HAFmx>=0)=1; % grounded
            MASKmy(HAFmy>=0)=1;
        end
        if ctr.basin==1
            MASKmx(bMASKx==1)=0;
            MASKmy(bMASKy==1)=0;
        end
    end
end


