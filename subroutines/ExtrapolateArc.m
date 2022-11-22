function Arc=ExtrapolateArc(MASK,oldMASK,Arc,Arc0,ctr)

% Extrapolate Arc between consecutive calls (to reduce
% calculation speed)

    dMASK=MASK-oldMASK; % difference between two consecutive
                        % masks (grounded - floated)
    Me=zeros(ctr.imax,ctr.jmax,8);
    OldArc=Arc;
    OldArc(MASK==1)=NaN;
    Me(:,:,1)=circshift(OldArc,[-1 1]);
    Me(:,:,2)=circshift(OldArc,[0 1]);
    Me(:,:,3)=circshift(OldArc,[1 1]);
    Me(:,:,4)=circshift(OldArc,[1 0]);
    Me(:,:,5)=circshift(OldArc,[-1 0]);
    Me(:,:,6)=circshift(OldArc,[-1 -1]);
    Me(:,:,7)=circshift(OldArc,[0 -1]);
    Me(:,:,8)=circshift(OldArc,[1 -1]);
    NewArc=nanmean(Me,3); % Take mean of neighbours for which Arc exists
    NewArc(isnan(NewArc))=0;
    Arc(dMASK==-1)=NewArc(dMASK==-1); % Apply Arc for points becoming floated
    % Apply Arc for points that were previously determined
    Arc(dMASK==-1 & Arc0~=0)=Arc0(dMASK==-1 & Arc0~=0);
    Arc(dMASK==1)=0; % remove melt for points becoming grounded
end


