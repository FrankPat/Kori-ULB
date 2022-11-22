function [uxsch,uysch]=GroundingLines(ctr,par,glMASK,HAF,SLR,Ax,Ay,butfac, ...
    H,Hmx,Hmy,B,Bmx,Bmy,Asfx,Asfy,ux,uy,dudx,dvdy,dudy,dvdx,eta)

% Kori-ULB
% Caluclation of buttress factors and parameterized flux condition at the
% grounding line for large-scale model simulations.

    if ctr.shelf==1
        Txx=2*eta.*(2.*dudx+dvdy)./max(par.SeaIceThickness,H);
        Tyy=2*eta.*(2.*dvdy+dudx)./max(par.SeaIceThickness,H);
        Txy=2*eta.*(dudy+dvdx)./max(par.SeaIceThickness,H); %VL: adaption to Pollard & deConto [2020]
        Tf=.5*par.rho*par.g*max(par.SeaIceThickness,H)*(1.-par.rho/par.rhow); % on h-grid
    else
        %VL: also provide Tf, Txx, Tyy, Txy for GL flux
        Tf=.5*par.rho*par.g*max(par.SeaIceThickness,H)*(1.-par.rho/par.rhow); % on h-grid
        Txx=Tf;
        Tyy=Tf;
        Txy=Tf;
    end
    % Applying Schoof BC and heuristic from Pollard
    % Weighting factor for Schoof correction as function bed
    ncorx=max(min(1-1/50*(Bmx+50),1),0);
    ncory=max(min(1-1/50*(Bmy+50),1),0);
    [uxsch,uysch]=GroundingLineFlux(ctr,glMASK,HAF,B,SLR,par,Ax,Ay, ...
        Tf,Txx,Tyy,Txy,butfac,ncorx,ncory,ux,uy,Hmx,Hmy,Asfx,Asfy); %VL: new function + new call!

end


