function [invmax2D,ncor]=InitOptimization(ctr,par,ncor,MASK,B,stdB)

% Kori-ULB
% Initialization of control matrices for the optimization (initialization)
% of the ice sheet

    stdB2=BedVar(B,ctr); %VL: large scale bed variability for new ncor definition
    %VL: definition of new ncor
    ncor((stdB2>=prctile(stdB2(MASK==1),95) | ...
        stdB>=prctile(stdB(MASK==1),99)) & MASK==1 & B>0)=0;
    % Definition of invmax2D to decrease invmax in areas with high
    % stdB -- max variation of 3 orders of magnitude
    invmax2D=zeros(ctr.imax,ctr.jmax)+log10(par.invmax);
    invmax2D(ncor==0)=log10(par.invmaxncor)-3.*((stdB(ncor==0)- ...
        prctile(stdB(ncor==0),10))/(prctile(stdB(ncor==0),90)- ...
        prctile(stdB(ncor==0),10)));
    invmax2D=10.^(invmax2D);
    invmax2D(ncor==0 & invmax2D>par.invmaxncor)=par.invmaxncor;
    invmax2D(invmax2D<1e-8)=1e-8; % limit to 3 orders of magnitude

end


