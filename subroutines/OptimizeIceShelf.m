function [Melt]=OptimizeIceShelf(ctr,MASKo,glMASK,H,Ho,Melt,bMASK)

% Kori-ULB
% Optimization of basal melt and accretion underneath ice shelves according
% to the method by Bernales

    Ftan=1.725;
    Melt=Melt+Ftan.*tan(max(-1.3,min(1.3,(H-Ho)./ctr.HinvMelt)));
    % smaller limit as Bernales (changes of Melt were too big for limit 
    % 1.5 --> created holes in Ross)
    Melt(MASKo==1)=0; % No Melt for grounded grid cells
    Melt=min(100,max(-100,Melt)); % limit melt
    Melt(glMASK==6)=0;
    if ctr.basin==1
        Melt(bMASK==1)=0;
    end
    
end


