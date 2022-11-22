function [uxsch,uysch]=GroundingLineFlux(ctr,glMASK,HAF,B,SLR,par, ...
    Ax,Ay,Tf,Txx,Tyy,Txy,butfac,ncorx,ncory,ux,uy,Hmx,Hmy,Asfx,Asfy)

% Kori-ULB
% Grounding line flux calculation according to Schoof following
% the method of Pollard and DeConto (2012; 2020)

% LZ2021: Implementation of Schoof following Pollard & DeConto (2012)
% + GL definition on H-grid + weighting factor [Pollard & DeConto (2020)]
% + correction GL definition (glMASK>=3 for adjacent floating grid cell
% instead glMASK==3)
% + correction for first floating grid cell:
%   --> advancing case (wx==0, wy==0): use minimum of velocity of GL grid 
%       cell and schoof flux
% 	--> retreating case (wx>0, wy>0): use schoof flux
% + stricter definition of GL grid cells: 
%   --> SSA velocity at GL has to have correct sign 
%   --> upstream H-grid cell has to be grounded
%   --> upstream velocity has to have correct sign
%   --> downstream H-grid cell has to be floating (if schoof vel. should be
%   applied at downstream grid cell)
% + apply Schoof downstream of GL only if shelf exists
% + orientation of GL (applied for buttressing and schoof velocity)
%       [Pollard & DeConto (2020)]

    uxsch=ux;
    uysch=uy;
    [jGL,iGL] = meshgrid(1:ctr.jmax,1:ctr.imax); %VL: grid for GL orientation

    % GL flux in x-direction
    q0=(Ax*(par.rho*par.g)^(par.n+1)*(1-par.rho/par.rhow)^par.n.* ...
        Asfx.^(1/ctr.m)/4^par.n);
    qs=q0.^(ctr.m/(ctr.m+1));
    qe0=ctr.m*(par.n+2)/(ctr.m+1);

    % conditions for GL in x-direction and ux>0
    qgx=zeros(ctr.imax,ctr.jmax);
    Mgl=zeros(ctr.imax,ctr.jmax);
    angnorm=zeros(ctr.imax,ctr.jmax);   %VL
    nx=ones(ctr.imax,ctr.jmax); %VL
    ny=zeros(ctr.imax,ctr.jmax); %VL
    Theta=ones(ctr.imax,ctr.jmax); %VL
    glMASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
    glMASK0=circshift(glMASK,[0 1]); % glMASK(i,j-1), upstream
    ux1=circshift(ux,[0 1]); % ux(i,j-1), upstream velocity
    HAF1=circshift(HAF,[0 -1]); % HAF(i,j+1)
    B1=circshift(B,[0 -1]); % B(i,j+1)
    Tf1=circshift(Tf,[0 -1]); %VL:  Tf(i,j+1)
    Txx1=circshift(Txx,[0 -1]); %VL:  Txx(i,j+1)
    Tyy1=circshift(Tyy,[0 -1]); %VL:  Tyy(i,j+1)
    Txy1=circshift(Txy,[0 -1]); %VL:  Txy(i,j+1)
    Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & ux>0 & ux1>0)=1; % u-grid(i,j) 
    fracgx=HAF./(HAF-HAF1);
    hbgx=(1.-fracgx).*B+fracgx.*B1;
    hgx=(SLR-hbgx)*par.rhow/par.rho;
    Mgl(hgx<0)=0;
    % GL in x-direction, u>0 --> jg=jGL, jh=jg+1, ih=ig=iGL
    [angnorm(Mgl==1)]=arrayfun(@(iGL,jGL) NormCalc(iGL,jGL,iGL,jGL+1, ...
        glMASK,ctr),iGL(Mgl==1),jGL(Mgl==1));
    nx(Mgl==1)=cos(angnorm(Mgl==1));
    ny(Mgl==1)=sin(angnorm(Mgl==1));
    Theta(Mgl==1)=max(min((butfac*(Txx1(Mgl==1).*nx(Mgl==1).^2+Tyy1(Mgl==1).* ...
        ny(Mgl==1).^2+Txy1(Mgl==1).*nx(Mgl==1).*ny(Mgl==1))+(1.-butfac)* ...
        Tf1(Mgl==1))./Tf1(Mgl==1),1),0);
    if ctr.schoof==1 % Schoof
        Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
    else % Pattyn
        Theta=Theta.^(par.n);
    end
    if ctr.schoof==1 % Schoof
        ugx=qs.*hgx.^qe0.*Theta;
    else % Pattyn
        ugx=q0.*hgx.^(par.n+3).*Theta;
    end
    ugx(Mgl==1)=ugx(Mgl==1).*abs(nx(Mgl==1));
    qgx(Mgl==1)=ugx(Mgl==1).*hgx(Mgl==1);
    wx=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
    wx(Mgl==1)=max(0,min(1,(ugx(Mgl==1)-ux(Mgl==1)).*hgx(Mgl==1)./10^5));
    % GL grid cell
    uxsch(Mgl==1)=ux(Mgl==1)+ncorx(Mgl==1).*wx(Mgl==1).*(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01)-ux(Mgl==1));
    % ice shelf grid cell downstream of GL
    Mgl=circshift(Mgl,[0 1]); % condition for i,j+1
    Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
    Mgl(Hmx==0)=0; % only apply if shelf exists
    qgx=circshift(qgx,[0 1]);
    ncx=circshift(ncorx,[0 1]);
    wx=circshift(wx,[0 1]);
    uxsch(Mgl==1)=ux(Mgl==1)+ncx(Mgl==1).*( ...
        sign(wx(Mgl==1)).*(1-wx(Mgl==1)).*(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01)-ux(Mgl==1)) ... % wx>0 (retreat)
        +(1-sign(wx(Mgl==1))).*(min(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01),ux1(Mgl==1))-ux(Mgl==1)));   % wx=0 (advance)
    
    % conditions for GL in x-direction and ux<0
    qgx=zeros(ctr.imax,ctr.jmax);
    Mgl=zeros(ctr.imax,ctr.jmax);
    angnorm=zeros(ctr.imax,ctr.jmax)+pi;   %VL
    nx=zeros(ctr.imax,ctr.jmax)-1; %VL
    ny=zeros(ctr.imax,ctr.jmax); %VL
    Theta=ones(ctr.imax,ctr.jmax); %VL
    glMASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
    glMASK0=circshift(glMASK,[0 -2]); % glMASK(i,j+2), upstream
    ux1=circshift(ux,[0 -1]);   % ux(i,j+1), upstream velocity
    HAF1=circshift(HAF,[0 -1]); % HAF(i,j+1)
    B1=circshift(B,[0 -1]); % B(i,j+1)
    Mgl(glMASK>=3 & glMASK1==2 & glMASK0<=2 & ux<0 & ux1<0)=1; % u-grid(i,j)
    fracgx=HAF1./(HAF1-HAF);
    hbgx=(1.-fracgx).*B1+fracgx.*B;
    hgx=(SLR-hbgx)*par.rhow/par.rho;
    Mgl(hgx<0)=0;
    %VL: GL in x-direction, u<0 --> jh=jGL, jg=jh+1, ih=ig=iGL
    [angnorm(Mgl==1)]=arrayfun(@(iGL,jGL) NormCalc(iGL,jGL+1,iGL,jGL, ...
        glMASK,ctr),iGL(Mgl==1),jGL(Mgl==1));  %VL
    nx(Mgl==1)=cos(angnorm(Mgl==1));
    ny(Mgl==1)=sin(angnorm(Mgl==1));
    Theta(Mgl==1)=max(min((butfac*(Txx(Mgl==1).*nx(Mgl==1).^2+Tyy(Mgl==1).* ...
        ny(Mgl==1).^2+Txy(Mgl==1).*nx(Mgl==1).*ny(Mgl==1))+(1.-butfac)* ...
        Tf(Mgl==1))./Tf(Mgl==1),1),0);
    if ctr.schoof==1 % Schoof
        Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
    else % Pattyn
        Theta=Theta.^(par.n);
    end
    if ctr.schoof==1 % Schoof
        ugx=qs.*hgx.^qe0.*Theta;
    else % Pattyn
        ugx=q0.*hgx.^(par.n+3).*Theta;
    end
    ugx(Mgl==1)=ugx(Mgl==1).*abs(nx(Mgl==1));
    qgx(Mgl==1)=-ugx(Mgl==1).*hgx(Mgl==1);
    wx=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
    wx(Mgl==1)=max(0,min(1,(ugx(Mgl==1)-abs(ux(Mgl==1))).*hgx(Mgl==1)./10^5));
    % GL grid cell
    uxsch(Mgl==1)=ux(Mgl==1)+ncorx(Mgl==1).*wx(Mgl==1).*(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01)-ux(Mgl==1));
    % ice shelf grid cell downstream of GL
    Mgl=circshift(Mgl,[0 -1]); % condition for i,j-1
    Mgl(glMASK<=2)=0; % only apply if downstream grid cell is floating
    Mgl(Hmx==0)=0; % only apply if shelf exists
    qgx=circshift(qgx,[0 -1]);
    ncx=circshift(ncorx,[0 -1]);
    wx=circshift(wx,[0 -1]);
    uxsch(Mgl==1)=ux(Mgl==1)+ncx(Mgl==1).*( ...
        sign(wx(Mgl==1)).*(1-wx(Mgl==1)).*(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01)-ux(Mgl==1)) ... % wx>0 (retreat)
        +(1-sign(wx(Mgl==1))).*(max(qgx(Mgl==1)./ ...
        max(Hmx(Mgl==1),0.01),ux1(Mgl==1))-ux(Mgl==1)));   % wx=0 (advance)
    
    
    % GL flux in y-direction
    q0=(Ay*(par.rho*par.g)^(par.n+1)*(1-par.rho/par.rhow)^par.n.* ...
        Asfy.^(1/ctr.m)/4^par.n);
    qs=q0.^(ctr.m/(ctr.m+1));

    % conditions for GL in y-direction and uy>0
    qgy=zeros(ctr.imax,ctr.jmax);
    Mgl=zeros(ctr.imax,ctr.jmax);
    angnorm=zeros(ctr.imax,ctr.jmax)+pi/2;   %VL
    nx=zeros(ctr.imax,ctr.jmax); %VL
    ny=ones(ctr.imax,ctr.jmax); %VL
    Theta=ones(ctr.imax,ctr.jmax); %VL
    glMASK1=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
    glMASK0=circshift(glMASK,[1 0]); % glMASK(i-1,j), upstream
    uy1=circshift(uy,[1 0]);    % uy(i-1,j), upstream velocity
    HAF1=circshift(HAF,[-1 0]); % HAF(i+1,j)
    B1=circshift(B,[-1 0]); % B(i+1,j)
    Tf1=circshift(Tf,[-1 0]); %VL:  Tf(i+1,j)
    Txx1=circshift(Txx,[-1 0]); %VL:  Txx(i+1,j)
    Tyy1=circshift(Tyy,[-1 0]); %VL:  Tyy(i+1,j)
    Txy1=circshift(Txy,[-1 0]); %VL:  Txy(i+1,j)
    Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & uy>0 & uy1>0)=1; % v-grid(i,j)
    fracgy=HAF./(HAF-HAF1);
    hbgy=(1.-fracgy).*B+fracgy.*B1;
    hgy=(SLR-hbgy)*par.rhow/par.rho;
    Mgl(hgy<0)=0;
    %VL: GL in y-direction, v>0 --> ig=iGL, ih=ig+1, jh=jg=jGL
    [angnorm(Mgl==1)]=arrayfun(@(iGL,jGL) NormCalc(iGL,jGL,iGL+1,jGL, ...
        glMASK,ctr),iGL(Mgl==1),jGL(Mgl==1));  %VL
    nx(Mgl==1)=cos(angnorm(Mgl==1));
    ny(Mgl==1)=sin(angnorm(Mgl==1));
    Theta(Mgl==1)=max(min((butfac*(Txx1(Mgl==1).*nx(Mgl==1).^2+ ...
        Tyy1(Mgl==1).*ny(Mgl==1).^2+Txy1(Mgl==1).*nx(Mgl==1).*ny(Mgl==1))+ ...
        (1.-butfac)*Tf1(Mgl==1))./Tf1(Mgl==1),1),0);
    if ctr.schoof==1 % Schoof
        Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
    else % Pattyn
        Theta=Theta.^(par.n);
    end
    if ctr.schoof==1 % Schoof
        ugy=qs.*hgy.^qe0.*Theta;
    else % Pattyn
        ugy=q0.*hgy.^(par.n+3).*Theta;
    end
    ugy(Mgl==1)=ugy(Mgl==1).*abs(ny(Mgl==1));
    qgy(Mgl==1)=ugy(Mgl==1).*hgy(Mgl==1);
    wy=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
    wy(Mgl==1)=max(0,min(1,(ugy(Mgl==1)-uy(Mgl==1)).*hgy(Mgl==1)./10^5));
    % GL grid cell
    uysch(Mgl==1)=uy(Mgl==1)+ncory(Mgl==1).*wy(Mgl==1).*(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01)-uy(Mgl==1));
    % ice shelf grid cell downstream of GL
    Mgl=circshift(Mgl,[1 0]); % condition for i+1,j
    Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
    Mgl(Hmy==0)=0; % only apply if shelf exists
    qgy=circshift(qgy,[1 0]);
    ncy=circshift(ncory,[1 0]);
    wy=circshift(wy,[1 0]);
    uysch(Mgl==1)=uy(Mgl==1)+ncy(Mgl==1).*( ...
        sign(wy(Mgl==1)).*(1-wy(Mgl==1)).*(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01)-uy(Mgl==1)) ...   % wy>0 (retreat)
        +(1-sign(wy(Mgl==1))).*(min(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01),uy1(Mgl==1))-uy(Mgl==1)));   % wy=0 (advance)
    
    % conditions for GL in y-direction and uy<0
    qgy=zeros(ctr.imax,ctr.jmax);
    Mgl=zeros(ctr.imax,ctr.jmax);
    angnorm=zeros(ctr.imax,ctr.jmax)-pi/2;   %VL
    nx=zeros(ctr.imax,ctr.jmax); %VL
    ny=zeros(ctr.imax,ctr.jmax)-1; %VL
    Theta=ones(ctr.imax,ctr.jmax); %VL
    glMASK1=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
    glMASK0=circshift(glMASK,[-2 0]); % glMASK(i+2,j), upstream
    uy1=circshift(uy,[-1 0]);   % uy(i+1,j), upstream velocity
    HAF1=circshift(HAF,[-1 0]); % HAF(i+1,j)
    B1=circshift(B,[-1 0]); % B(i+1,j)
    Mgl(glMASK>=3 & glMASK1==2 & glMASK0<=2 & uy<0 & uy1<0)=1; % v-grid(i,j)
    fracgy=HAF1./(HAF1-HAF);
    hbgy=(1.-fracgy).*B1+fracgy.*B;
    hgy=(SLR-hbgy)*par.rhow/par.rho;
    Mgl(hgy<0)=0;
    %VL: GL in y-direction, v<0 --> ih=iGL, ig=ih+1, jh=jg=jGL
    [angnorm(Mgl==1)]=arrayfun(@(iGL,jGL) NormCalc(iGL+1,jGL,iGL,jGL, ...
        glMASK,ctr),iGL(Mgl==1),jGL(Mgl==1));  %VL
    nx(Mgl==1)=cos(angnorm(Mgl==1));
    ny(Mgl==1)=sin(angnorm(Mgl==1));
    Theta(Mgl==1)=max(min((butfac*(Txx(Mgl==1).*nx(Mgl==1).^2+Tyy(Mgl==1).* ...
        ny(Mgl==1).^2+Txy(Mgl==1).*nx(Mgl==1).*ny(Mgl==1))+(1.-butfac)* ...
        Tf(Mgl==1))./Tf(Mgl==1),1),0);
    if ctr.schoof==1 % Schoof
        Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
    else % Pattyn
        Theta=Theta.^(par.n);
    end
    if ctr.schoof==1 % Schoof
        ugy=qs.*hgy.^qe0.*Theta;
    else % Pattyn
        ugy=q0.*hgy.^(par.n+3).*Theta;
    end
    ugy(Mgl==1)=ugy(Mgl==1).*abs(ny(Mgl==1));
    qgy(Mgl==1)=-ugy(Mgl==1).*hgy(Mgl==1);
    wy=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
    wy(Mgl==1)=max(0,min(1,(ugy(Mgl==1)-abs(uy(Mgl==1))).*hgy(Mgl==1)./10^5));
    % GL grid cell
    uysch(Mgl==1)=uy(Mgl==1)+ncory(Mgl==1).*wy(Mgl==1).*(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01)-uy(Mgl==1));
    % ice shelf grid cell downstream of GL
    Mgl=circshift(Mgl,[-1 0]); % condition for i-1,j
    Mgl(glMASK<=2)=0; % only apply if downstream grid cell is floating
    Mgl(Hmy==0)=0; % only apply if shelf exists
    qgy=circshift(qgy,[-1 0]);
    ncy=circshift(ncory,[-1 0]);
    wy=circshift(wy,[-1 0]);
    uysch(Mgl==1)=uy(Mgl==1)+ncy(Mgl==1).*( ...
        sign(wy(Mgl==1)).*(1-wy(Mgl==1)).*(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01)-uy(Mgl==1)) ...   % wy>0 (retreat)
        +(1-sign(wy(Mgl==1))).*(max(qgy(Mgl==1)./ ...
        max(Hmy(Mgl==1),0.01),uy1(Mgl==1))-uy(Mgl==1)));   % wy=0 (advance)
    
end


