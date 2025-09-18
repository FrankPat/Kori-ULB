function [As,Mb0,Ts0,MASK0,MASK,bMASK,H,B,H0,B0]= ...
    InitInputData(ctr,par,As,Mb,Ts,MASK,MASKo,H,B)

% Kori-ULB
% Initialization of input matrices after reading the input file

    As=As*ctr.slidfac;
    Mb0=Mb; % initial mass balance if read from inputfile
    Ts0=Ts; % intial surface temperature if read from inputfile
    MASK0=MASK; % MASK at start of run (MASKo is origiginal BedMap mask, if exists)
    MASK(MASK==3)=0; % ice shelves in BedMap2/Bedmachine identified with MASK=3

    bMASK=zeros(ctr.imax,ctr.jmax);
    if ctr.basin==1
        bMASK(MASKo==-1)=1;
        MASK(bMASK==1)=1;
        As(bMASK==1)=par.As0;
        H(1,:)=H(2,:);
        B(1,:)=B(2,:);
        H(ctr.imax,:)=H(ctr.imax-1,:);
        B(ctr.imax,:)=B(ctr.imax-1,:);
        H(:,1)=H(:,2);
        B(:,1)=B(:,2);
        H(:,ctr.jmax)=H(:,ctr.jmax-1);
        B(:,ctr.jmax)=B(:,ctr.jmax-1);
    end

    H0=H; % initial ice thickness
    B0=B; % initial bed elevation
end


