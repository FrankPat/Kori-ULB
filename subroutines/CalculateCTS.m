function [CTSm,CTSp,Ht,E,idx]=CalculateCTS(ctr,E,Epmp,MASK,H,zeta)

% Kori-ULB
% Calculate Cold Temperate Transitions Surface (CTS)
    
    % Cold-temperate transition surface
    CTS=zeros(ctr.imax,ctr.jmax,ctr.kmax); CTSm=CTS; CTSp=CTS;
    
    % only keep the first value for the CTS where the condition apply
    % --> assumption that only one CTS exists for a given ice column
    first_one=false(size(CTS, 1), size(CTS, 2));
    second_one=false(size(CTS, 1), size(CTS, 2));

    for k=ctr.kmax-1:-1:2
        % first grid point where E>=Epmp
        CTSm(:,:,k)=E(:,:,k)>=Epmp(:,:,k)&E(:,:,k-1)<Epmp(:,:,k-1) & ~first_one;    
        first_one=first_one | CTSm(:,:,k);
        % last grid point where E<Epmp
        CTSp(:,:,k)=E(:,:,k)<Epmp(:,:,k)&E(:,:,k+1)>=Epmp(:,:,k+1) & ~second_one;
        second_one = second_one | CTSp(:,:,k);
    end

    % Find index of the CTS
    [idx_i,idx_j,idx_k]=ind2sub(size(CTSm), find(CTSm == 1));
    idx=zeros(ctr.imax,ctr.jmax);
    if ~isempty(idx_i)
        lin_idx = sub2ind(size(idx), idx_i, idx_j);
        idx(lin_idx) = idx_k;
    end

    % Basal conditions
    Cld=(E(:,:,ctr.kmax)<Epmp(:,:,ctr.kmax)); % cold base
    Abv=(E(:,:,ctr.kmax)>=Epmp(:,:,ctr.kmax)&E(:,:,ctr.kmax-1)<Epmp(:,:,ctr.kmax-1)); % temperate base
    % Blw=(E(:,:,ctr.kmax)>=Epmp(:,:,ctr.kmax)&E(:,:,ctr.kmax-1)>=Epmp(:,:,ctr.kmax-1)); % temperate layer

    % Temperate layer thickness
    zpad=[1, zeta]; % correction to take idx+1 and avoid the error of 0 indexing
    Ht=H.*(1-zpad(idx+1));
    Ht(Abv==1|Cld==1)=0;
    Ht(MASK==0)=0;
    % No CTS if cold layer of ice above the bed
    CTSp(repmat(Abv|Cld,[1 1 ctr.kmax]))=0;
    CTSm(repmat(Abv|Cld,[1 1 ctr.kmax]))=0;
    % Correction when CTS is at the bed-ice interface
    repCTSm=CTSm(:,:,ctr.kmax); repCTSm(Abv==1)=1; CTSm(:,:,ctr.kmax)=repCTSm;
    repE=E(:,:,ctr.kmax); repEpmp=Epmp(:,:,ctr.kmax);
    repE(Abv==1)=repEpmp(Abv==1); E(:,:,ctr.kmax)=repE;

    for k=2:ctr.kmax
    % Find if ice is cold below the CTS while it should be temperate
    mask1 = (k>idx) & (E(:,:,k)<Epmp(:,:,k)) & Ht>0;
    % Find if ice is temperate above the CTS while it should be cold
    mask2 = (k<idx) & (E(:,:,k)>Epmp(:,:,k)) & Ht>0;
    % Correct E if necessary
    repE=E(:,:,k); 
    repEpmp=Epmp(:,:,k);
    repE(mask1|mask2)=repEpmp(mask1|mask2);
    E(:,:,k)=repE;
    end
end


