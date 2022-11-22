function [contshelfMASK]=ContinentalShelfMASK(MASK,H,B,par)

% Kori-ULB
% Definition of contshelfMASK

    [~,contshelf]=bwboundaries(MASK==0 & B>-2000,'noholes',4);
    MASK(MASK==0 & H>par.SeaIceThickness)=3; % consider ice shelves (MASK==3)
    contshelfMASK=B*0;% Initializing contshelfMASK
    MASK1=circshift(MASK,[0 -1]); % MASK(i,j+1)
    MASK2=circshift(MASK,[0 1]); % MASK(i,j-1)
    MASK3=circshift(MASK,[-1 0]); % MASK(i+1,j)
    MASK4=circshift(MASK,[1 0]); % MASK(i-1,j)

    % Iterate through contshelf items and check if they are somehow connected
    % to either shelves or grounded ice sheet. If so, set contshelfMASK=1 for
    % that item
    for i=1:max(max(contshelf))
        if sum(MASK1(contshelf==i)~=0)+sum(MASK2(contshelf==i)~=0)+ ...
                sum(MASK3(contshelf==i)~=0)+sum(MASK4(contshelf==i)~=0)~=0
            contshelfMASK(contshelf==i)=1;
        end
    end
    
end


