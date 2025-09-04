function [As,deltaZ,Asor]=OptimizeIceSheet(ctr,par,cnt, ...
    Asor,MASK,bMASK,deltaZ,sn,sn0,r,ncor,B,stdB,vx,vy,ux,uy,invmax2D)

% Kori-ULB
% Optimization of basal sliding coefficients underneath grounded ice sheet
% using the nudging algorithm of Pollard and DeConto (2012)

    % introduce deltaZold to compare deltaZ between consecutive iterations
    deltaZold=zeros(ctr.imax,ctr.jmax);     %LZ
    % optimize only for grounded grid cells
    if cnt*ctr.dt>ctr.Tinv
        deltaZold(MASK==1)=deltaZ(MASK==1);
    end
    deltaZ=zeros(ctr.imax,ctr.jmax);
    deltaZ(MASK==1)=max(-1.5,min(1.5,(sn(MASK==1) ...
        -sn0(MASK==1))/ctr.Hinv));
    if ctr.SlidAdjust==1 && ctr.Tcalc>0
        Asor(r>0 & abs(deltaZ)>=abs(deltaZold))=Asor(r>0 & ...
            abs(deltaZ)>=abs(deltaZold)).*10.^deltaZ(r>0 & ...
            abs(deltaZ)>=abs(deltaZold));
    else
        Asor(abs(deltaZ)>=abs(deltaZold))=Asor(abs(deltaZ)>= ...
            abs(deltaZold)).*10.^deltaZ(abs(deltaZ)>=abs(deltaZold));
    end
    Asor(Asor<par.invmin/10)=par.invmin/10;
    Asor(Asor>par.invmax*1000)=par.invmax*1000;
    Asor(Asor>invmax2D & ncor==0)=invmax2D(Asor>invmax2D & ncor==0);
    if ctr.vexist==1
        As = RegularizationNew(Asor,B,stdB,vx,vy,ctr,ncor,par.stdDevRegul);
    else
        As = RegularizationNew(Asor,B,stdB,ux,uy,ctr,ncor,par.stdDevRegul);
    end
    As(As<par.invmin)=par.invmin; % first limits
    As(As>par.invmax)=par.invmax;
    if ctr.basin==1 %VL: As on h-grid
        As(bMASK==1)=par.As0;
    end
end


