function [ctr,invmax2D,Asor,ncor,To,So,Pr0,Evp0,runoff0,Evp,Hinit]= ...
    ExistParams(ctr,par,ncor,Asor,stdB,v,uxssa,To,So,Db,B,MASK,As,Pr, ...
    Evp,runoff,Mb0,Hinit,Ho,damage)

% Kori-ULB
% Test what parameters and matrices exist and initializes them accordingly

    % check whether model domain has grounding lines
    if ctr.shelf==1 || ctr.schoof>=1
        ctr.glMASKexist=1;
    else
        ctr.glMASKexist=0;
    end

    % check whether bedrock variability is defined
    if islogical(stdB)==1
        ctr.stdBexist=0; 
    else
        ctr.stdBexist=1;
    end

    % check whether observed velocity field exists
    if islogical(v)==1
        ctr.vexist=0;
    else
        ctr.vexist=1;
    end
    
    % check whether varying lithosphere thickness exists
    if islogical(Db)==1
        ctr.Dbexist=0;
    else
        ctr.Dbexist=1;
    end

    % check whether modelled SSA velocity field exists
    if islogical(uxssa)==1
        ctr.uSSAexist=0;
    else
        ctr.uSSAexist=1;
    end
    
    % check whether damage field exists
    if islogical(damage)==1
        ctr.damexist=0;
    else
        ctr.damexist=1;
    end
    
    if islogical(To)==1
        To=zeros(ctr.imax,ctr.jmax)+par.Toi;
        So=zeros(ctr.imax,ctr.jmax)+par.Soi;
    else
        To(isnan(To))=par.Toi;
        So(isnan(So))=par.Soi;
    end
    if ctr.inverse>0
        if islogical(Asor)==1 %VL: check if Asor already exists
            Asor=As;
        end
        if ctr.stdBexist==1
            [invmax2D,ncor]=InitOptimization(ctr,par,ncor,MASK,B,stdB);
        else
            invmax2D=zeros(ctr.imax,ctr.jmax)+par.invmax;
        end
    else
        invmax2D=false;
    end
    if islogical(Pr)==1
        Pr0=Mb0;
    else
        Pr0=Pr;
    end
    if islogical(Evp)==1
        Evp0=Pr0-Mb0;
        Evp=Evp0;
    else
        Evp0=Evp; 
    end
    if islogical(runoff)==1
        runoff0=Pr0-Mb0-Evp0;
    else
        runoff0=runoff;
    end
    if islogical(Hinit)==1
        Hinit=Ho;
    end
end


