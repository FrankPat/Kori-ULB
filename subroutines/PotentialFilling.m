function [pot]=PotentialFilling(pot,ctr)

% Kori-ULB
% Hollow filling for potential gradient

    for k=1:10
        pool=zeros(ctr.imax,ctr.jmax);
        p1=circshift(pot,[0 -1]);
        p2=circshift(pot,[0 1]);
        p3=circshift(pot,[-1 0]);
        p4=circshift(pot,[1 0]);
        pool(pot<p1 & pot<p2 & pot<p3 & pot<p4)=1;
        pot(pool==1)=(p1(pool==1)+p2(pool==1)+p3(pool==1)+p4(pool==1))/4;
    end

end


