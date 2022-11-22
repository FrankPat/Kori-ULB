function As = RegularizationNew(Asor,B,stdB,vx,vy,ctr,ncor,sigma)

% Kori-ULB
% New regularization function based on Gaussian filter for the optimization
% of subglacial sliding coefficients of the grounded ice sheet

    DEM = log10(Asor);
    % sigma is the standard deviation of the Gaussian distribution
    % span defines the size of the Gaussian filter
    nwin=round(sigma*2);
    span=nwin+0.5;  

    vx1=0.5*(vx+circshift(vx,[0 1]));
    vy1=0.5*(vy+circshift(vy,[1 0]));
    v=sqrt(vx1.^2+vy1.^2);
    VelAng=atan2(vy1,vx1);

    TotWeightFac=zeros(ctr.imax,ctr.jmax);
    DEMnew=zeros(ctr.imax,ctr.jmax);

    for j=-nwin:nwin
        for i=-nwin:nwin
            zdist = sqrt(j^2 + i^2); %Distance in grid cells
            if zdist<=span
                zang=atan2(i,j);
                GaussFiltWeight=normpdf(zdist/sigma);
                VelMagWeight=min(1,max(0,log10(v/100)*zdist/sigma));
                VelWeight=1-abs(sin(zang-VelAng)).*VelMagWeight;
                BedWeight=1./2.^(0.5*(max(1,log10(abs(circshift(B,[-i -j]) ...
                    -B)))-1));
                if ctr.stdBexist==1
                    stdBweight=1./2.^(max(1,log10(circshift(stdB,[-i -j])))-1);
                else
                    stdBweight=1;
                end
                DEMnew=DEMnew+GaussFiltWeight*VelWeight.*BedWeight.* ...
                    stdBweight.*circshift(DEM,[-i -j]);
                TotWeightFac=TotWeightFac+GaussFiltWeight*VelWeight.* ...
                    BedWeight.*stdBweight; %for normalization --> sum of all Probs=1
            end
        end
    end
    As=Asor;
    As(ncor>0.5)=10.^(DEMnew(ncor>0.5)./TotWeightFac(ncor>0.5));

end


