function [Ts_yc,Pr_yc]=ExtractAnnualCycle(fc,ctr,par,Pr0,Ts0,S0,sn, ...
    MASK,lat,DeltaSL,DeltaT,snp_atm)

% Kori-ULB
% Extract annual cycles from input forcing data (atmospheric data)

    if any(ismember(fields(fc),'atm_Ts_fname'))
        Ts_yc=zeros(ctr.imax,ctr.jmax,fc.atm_yrstep);
        n=1;
        for i=snp_atm-1:snp_atm-1+fc.atm_yrstep-1
            load([fc.atm_Ts_fname,num2str(i,'%03i')]);
            Ts_yc(:,:,n)=Ts;
            n=n+1;
        end
        if ctr.TsType==1
            Ts_yc=Ts_yc+par.Tlapse*(max(sn,DeltaSL)-S0);
        end
    else
        fprintf('If ctr.monthly=1, intra-annual (monthly or submonthly) variations of the air temperature should be provided\n');
        % If intra-annual Ts not provided, then yearly cycle is 
        % approximated -- correction for elevation changes already applied
        
        % Yearly amplitude of the signal (applies to the Antarctic ice sheet)
        if islogical(lat)==1
            Ta=min(20+max(sn,0)/300.,30); % seasonal peak-to-peak amplitude in Ts
        else
            Ta=zeros(ctr.imax,ctr.jmax);
            Ta(MASK==1)=-11.06+0.003204*sn(MASK==1)-0.3584*lat(MASK==1);
            Ta(MASK==0)=-45.06-0.04598*sn(MASK==0)-0.969*lat(MASK==0);
            Ta(Ta<15)=15;
        end
    
        % Define idealized yearly cycle Tm:
        if any(ismember(fields(fc),'atm_yrstep')) 
            par.PDDsteps=fc.atm_yrstep;
        end
        t=repmat(reshape(linspace(1,365.25,par.PDDsteps),1,1,par.PDDsteps), ...
            [ctr.imax,ctr.jmax,1]);
        Ts_yc=repmat(Ts0,[1,1,par.PDDsteps])-0.5* ...
            repmat(Ta,[1,1,par.PDDsteps]).*sin(2*pi*t/365.25);
    end
    if any(ismember(fields(fc),'atm_Pr_fname'))
        Pr_yc=zeros(ctr.imax,ctr.jmax,fc.atm_yrstep);
        n=1;
        for i=snp_atm-1:snp_atm-1+fc.atm_yrstep-1
            load([fc.atm_Pr_fname,num2str(i,'%03i')]);
            Pr_yc(:,:,n)=Pr;
            n=n+1;
        end
        if ctr.MbType==1
            Pr_yc=Pr_yc.*exp(0.05*(par.Tlapse*(max(sn,DeltaSL)-S0)));
        elseif ctr.MbType==2
            Pr_yc=Pr_yc.*(1+0.053*(par.Tlapse*(max(sn,DeltaSL)-S0)));
        end
    else
        if any(ismember(fields(fc),'atm_yrstep'))
            par.PDDsteps=fc.atm_yrstep;
        end
        Pr_yc=zeros(ctr.imax,ctr.jmax,par.PDDsteps);
        Pr_yc=Pr_yc+Pr0; % If intra-annual Pr not provided, then considered 
                         % constant -- correction for elevation changes 
                         % already applied
    end
end


