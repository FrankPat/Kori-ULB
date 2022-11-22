function flux=DpareaWarGds(ewr,nsr,MASK,Bmelt,kegdsx,kegdsy, ...
    par,ctr,flux,gdsmag,funcnt)

% Kori-ULB
% Subglacial meltwater routing function from LeBrocq et al.

    if MASK(ewr,nsr)>=1 && flux(ewr,nsr)<0.
        flux(ewr,nsr)=Bmelt(ewr,nsr)*ctr.delta^2;
        count=1;
        for j=nsr-1:nsr+1
            for i=ewr-1:ewr+1
                pp=zeros(9,1);
                if count~=5
                    if i>=1 && i<ctr.imax && j>=1 && j<ctr.jmax
                        if kegdsy(i,j)>0.
                            pp(6)=kegdsy(i,j)/gdsmag(i,j);
                        elseif kegdsy(i,j)<0.
                            pp(4)=abs(kegdsy(i,j))/gdsmag(i,j);
                        end
                        if kegdsx(i,j)<0.
                            pp(2)=abs(kegdsx(i,j))/gdsmag(i,j);
                        elseif kegdsx(i,j)>0.
                            pp(8)=kegdsx(i,j)/gdsmag(i,j);
                        end
                        if pp(par.dirpp_war(count))>0 && funcnt<=5e4
                            funcnt=funcnt+1;
                            flux=DpareaWarGds(i,j,MASK,Bmelt,kegdsx,kegdsy, ...
                                par,ctr,flux,gdsmag,funcnt);
                            flux(ewr,nsr)=flux(ewr,nsr)+ ...
                                pp(par.dirpp_war(count))*flux(i,j);
                        end
                    end
                end
                count=count+1;
            end
        end
    end
    
end


