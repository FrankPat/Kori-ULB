function [lambda0,lambda1]=PrincipalStrain(dudx,dvdy,dudy,dvdx)

% Kori-ULB
% Calculation of principal strain based on Cauchy stress

    Exx=(2*dudx+dvdy);
    Eyy=(2*dvdy+dudx);
    Exy=0.5*(dudy+dvdx);
    b=0.5*(Exx+Eyy);
    d=sqrt((0.5*(Exx-Eyy)).^2+Exy.^2);
    lambda0=b+d;
    lambda1=b-d;

end
