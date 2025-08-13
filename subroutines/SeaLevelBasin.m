function [VAF_basin,SLC_basin]=SeaLevelBasin(ctr,par,SLR,B,H,H0,VAF0,POV0,ZB)

% Kori-ULB
% Sea level contribution per basin

    for i=1:max(ZB(:))
        VAFi=max(0,H+min(B-SLR,0)*(par.rhow/par.rho));
        % VAF variation (SL equivalent in ocean water)
        VAF=sum(sum((VAF0(ZB==i)-VAFi(ZB==i))*ctr.delta^2.))*par.rho/(par.Aoc*par.rhow);
        VA0i=max(0,H+min(B-par.SLref,0)*(par.rhow/par.rho));
        % VAreferenceSL variation (SL equivalent in ocean water)
        VA0=sum(sum((VAF0(ZB==i)-VA0i(ZB==i))*ctr.delta^2.))*par.rho/(par.Aoc*par.rhow);
        POVi=max(0,par.SLref-B);
        % Potential Ocean volume variation (SL equivalent in ocean water)
        POV=sum(sum((POV0(ZB==i)-POVi(ZB==i))*ctr.delta^2.))/par.Aoc;
        % Density correction for transformation from ice to freshwater and not ocean water (SL equivalent)
        DENScorr=sum(sum((H0(ZB==i)-H(ZB==i))*ctr.delta^2.))*((par.rho/par.rhow)-(par.rho/par.rhof))/par.Aoc;
        % Sea-level contribution from modelled ice-sheet
        VAF_basin(i)=VAF;
        SLC_basin(i)=VA0+POV-DENScorr;
    end

end


