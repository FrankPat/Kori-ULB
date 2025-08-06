function [CTSm,CTSp,Ht]=CalculateCTS(ctr,E,Epmp,MASK,H,zeta)

% Kori-ULB
% Calculate Cold Transitions Surface (CTS)

    % Temperate ice layer thickness
    Ht=zeros(ctr.imax,ctr.jmax);
    
    % Cold-temperate transition surface
    CTS=zeros(ctr.imax,ctr.jmax,ctr.kmax); CTSm=CTS; CTSp=CTS;
    
    % only keep the first value for the CTS where the condition apply
    % --> assumption that only one CTS exists for a given ice column
    first_one=false(size(CTS, 1), size(CTS, 2));
    second_one=false(size(CTS, 1), size(CTS, 2));
    for k=ctr.kmax-1 :-1:2
        % first grid point where E>=Epmp
        CTSm(:,:,k)=E(:,:,k)>=Epmp(:,:,k)&E(:,:,k-1)<Epmp(:,:,k-1) & ~first_one;    
        first_one=first_one | CTSm(:,:,k);
        % last grid point where E<Epmp
        CTSp(:,:,k)=E(:,:,k)<Epmp(:,:,k)&E(:,:,k+1)>=Epmp(:,:,k+1) & ~second_one;
        second_one = second_one | CTSp(:,:,k);
    end
    [idx_i,idx_j,idx_k]=ind2sub(size(CTSm), find(CTSm == 1));
    idx=zeros(ctr.imax,ctr.jmax);
    if ~isempty(idx_i)
        lin_idx = sub2ind(size(idx), idx_i, idx_j);
        idx(lin_idx) = idx_k;
    end
    for i=1:ctr.imax
        for j=1:ctr.jmax
            if MASK(i,j)>0 && H(i,j)>0
                % Find index of the CTS
                % idx1 = find(E(i,j,:)<Epmp(i,j,:), 1, 'last'); % CTSp
                % idx2 = find(E(i,j,:)>=Epmp(i,j,:), 1, 'first'); % CTSm
%                 idx = find(squeeze(CTSm(i,j,:))==1);
                if E(i,j,ctr.kmax)>=Epmp(i,j,ctr.kmax) && ...
                        E(i,j,ctr.kmax-1)>=Epmp(i,j,ctr.kmax-1)
                    % temperate layer thickness & CTS
                    if idx(i,j)~=0 
                        Ht(i,j) = H(i,j).*(1-zeta(idx(i,j))); 
                    end
                elseif E(i,j,ctr.kmax)>=Epmp(i,j,ctr.kmax) && ...
                        E(i,j,ctr.kmax-1)<Epmp(i,j,ctr.kmax-1)
                    CTSp(i,j,:)=0;
                    CTSm(i,j,:)=0;
                    Ht(i,j)=0;
                    % Correction when CTS is at the bed-ice interface
                    CTSm(i,j,ctr.kmax)=1;
                    E(i,j,ctr.kmax)=Epmp(i,j,ctr.kmax);
                else % no CTS if cold layer of ice above the bed
                    CTSp(i,j,:)=0;
                    CTSm(i,j,:)=0;
                    Ht(i,j)=0;
                end
            end
        end
    end
end


