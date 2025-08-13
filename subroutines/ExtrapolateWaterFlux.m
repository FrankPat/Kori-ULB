function Wd=ExtrapolateWaterFlux(MASK,oldMASK,Wd,Wd0,ctr,par)

% Kori-ULB
% Interpolate/extrapolate Wd between consecutive calls (to improve
% calculation speed)

    dMASK=MASK-oldMASK; % difference between two consecutive
                        % masks (grounded - floated)
    WD=zeros(ctr.imax,ctr.jmax,8);
    OldWD=Wd;
    OldWD(MASK==0)=NaN;
    WD(:,:,1)=circshift(OldWD,[-1 1]);
    WD(:,:,2)=circshift(OldWD,[0 1]);
    WD(:,:,3)=circshift(OldWD,[1 1]);
    WD(:,:,4)=circshift(OldWD,[1 0]);
    WD(:,:,5)=circshift(OldWD,[-1 0]);
    WD(:,:,6)=circshift(OldWD,[-1 -1]);
    WD(:,:,7)=circshift(OldWD,[0 -1]);
    WD(:,:,8)=circshift(OldWD,[1 -1]);
    NewWD=mean(WD,3,'omitnan'); % Take mean of neighbours for which Wd exists
    NewWD(isnan(NewWD))=par.Wdmin;
    Wd(dMASK==1)=NewWD(dMASK==1); % Apply depth for points becoming grounded
    % Apply flw depths for points that were previously determined
    %   with the subglacial water flow model
    Wd(dMASK==1 & Wd0~=par.Wdmin)=Wd0(dMASK==1 & Wd0~=par.Wdmin);
    Wd(dMASK==-1)=par.Wdmin; % set minimum Wd value for points becoming floated

end


