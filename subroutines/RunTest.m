function [ctr,outfile]=RunTest(ctr)

    ctr.imax=31;
    ctr.jmax=31;
    ctr.delta=50.e3;
    ctr.nsteps=501;
    ctr.dt=50.;
    ctr.plotH=1;
    ctr.MbConst=0.3;
    outfile='KoriModelTest';
    
end
    