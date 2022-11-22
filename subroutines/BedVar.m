function stdB2=BedVar(B,ctr)

% Kori-ULB
% Calculate large scale bed variability for ncor

    sigma=2;
    AllB=zeros(ctr.imax,ctr.jmax,(2*sigma+1)^2);
    k=1;
    for j=-sigma:sigma
        for i=-sigma:sigma
            zdist = sqrt(j^2 + i^2); %Distance in grid cells
            if zdist<=sigma
                AllB(:,:,k)=circshift(B,[-i -j]);
            end
            k=k+1;
        end
    end
    stdB2=std(AllB,1,3);

end


