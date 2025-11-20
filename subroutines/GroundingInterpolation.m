function [beta2]=GroundingInterpolation()

    B1=circshift(B,[0 1]); % i,j-1
    B2=circshift(B,[0 -1]); % i,j+1
    B3=circshift(B,[1 0]); % i-1,j
    B4=circshift(B,[1 -1]); % i-1,j+1
    B5=circshift(B,[1 1]); % i-1,j-1
    B6=circshift(B,[-1 0]); % i+1,j
    B7=circshift(B,[-1 -1]); % i+1,j+1
    B8=circshift(B,[-1 1]); % i+1,j-1
    
    BSW=0.25*(B+B1+B3+B5);
    BSE=0.25*(B+B2+B3+B4);
    BNW=0.25*(B+B1+B6+B8);
    BNE=0.25*(B+B2+B6+B7);
    
    H1=circshift(H,[0 1]); % i,j-1
    H2=circshift(H,[0 -1]); % i,j+1
    H3=circshift(H,[1 0]); % i-1,j
    H4=circshift(H,[1 -1]); % i-1,j+1
    H5=circshift(H,[1 1]); % i-1,j-1
    H6=circshift(H,[-1 0]); % i+1,j
    H7=circshift(H,[-1 -1]); % i+1,j+1
    H8=circshift(H,[-1 1]); % i+1,j-1
    
    HSW=0.25*(H+H1+H3+H5);
    HSE=0.25*(H+H2+H3+H4);
    HNW=0.25*(H+H1+H6+H8);
    HNE=0.25*(H+H2+H6+H7);
    
    fSW=-BSW-par.rho*HSW/par.rhow;
    fSE=-BSE-par.rho*HSE/par.rhow;
    fNW=-BNW-par.rho*HNW/par.rhow;
    fNE=-BNE-par.rho*HNE/par.rhow;
    
    a=fSW;
    b=fSE-fSW;
    c=fNW-fSW;
    d=fNE+fSW-fNW-fSE;
    
    phig=((b.*c-a.*d).*log(abs(1-a.*d./(b.*c+1e-8)))+a*d)./d.^2;



end

