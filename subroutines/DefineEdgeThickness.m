function [he,fi]=DefineEdgeThickness(ctr,par,glMASK,H)  %VL: add par
    
% Define ice thickness at the shelf edges used for calving.
% Based on Pollard et al (2015) and Pollard and DeConto (2012)
% Use of maximum thickness of surrounding points and not mean H

    IC=ones(ctr.imax,ctr.jmax);
    IC(glMASK>=5)=0; % only cells grounded/floated and not adjacent to ocean
    IC1=circshift(IC,[-1 0]);
    IC2=circshift(IC,[1 0]);
    IC3=circshift(IC,[0 1]);
    IC4=circshift(IC,[0 -1]);
    H1=circshift(H,[-1 0]);
    H2=circshift(H,[1 0]);
    H3=circshift(H,[0 1]);
    H4=circshift(H,[0 -1]);
    w1=1-min(1,H./(H1*exp(-ctr.delta/1e5)))*(1-exp(-ctr.delta/1e5));
    w2=1-min(1,H./(H2*exp(-ctr.delta/1e5)))*(1-exp(-ctr.delta/1e5));
    w3=1-min(1,H./(H3*exp(-ctr.delta/1e5)))*(1-exp(-ctr.delta/1e5));
    w4=1-min(1,H./(H4*exp(-ctr.delta/1e5)))*(1-exp(-ctr.delta/1e5));
    he=(H1.*w1.*IC1+H2.*w2.*IC2+H3.*w3.*IC3+H4.*w4.*IC4)./ ...
        (IC1+IC2+IC3+IC4);
    he(glMASK==6)=NaN;
    he(isnan(he))=H(isnan(he)); % no neighbouring ice-cells (iceberg)
    he(he<par.SeaIceThickness)=par.SeaIceThickness; %VL: make sure he can't become 0
    fi=min(1,H./he);
    fi(glMASK~=5)=1;
    fi(glMASK==6)=0;
    
end


