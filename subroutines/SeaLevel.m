function [VAF,POV,SLC,VA0]=SeaLevel(ctr,par,SLR,B,H,H0,VAF0,POV0)

% Kori-ULB
% Calculation of local sea level

    VAFi=max(0,H+min(B-SLR,0)*(par.rhow/par.rho));
    % VAF variation (SL equivalent in ocean water)
    VAF=sum(sum((VAF0-VAFi)*ctr.delta^2.))*par.rho/(par.Aoc*par.rhow);
    VA0i=max(0,H+min(B-par.SLref,0)*(par.rhow/par.rho));
    % VAreferenceSL variation (SL equivalent in ocean water)
    VA0=sum(sum((VAF0-VA0i)*ctr.delta^2.))*par.rho/(par.Aoc*par.rhow);
    POVi=max(0,par.SLref-B);
    % Potential Ocean volume variation (SL equivalent in ocean water)
    POV=sum(sum((POV0-POVi)*ctr.delta^2.))/par.Aoc;
    % Density correction for transformation from ice to freshwater and 
    % not ocean water (SL equivalent)
    DENScorr=sum(sum((H0-H)*ctr.delta^2.))*((par.rho/par.rhow)- ...
        (par.rho/par.rhof))/par.Aoc;
    % Sea-level contribution from modelled ice-sheet
    SLC=VA0+POV-DENScorr;

end


