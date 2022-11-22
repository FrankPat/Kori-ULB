function u=vec2h(ux,uy) %VL: velocity on h-grid

    ux1=circshift(ux,[0 1]); % ux(i,j-1)
    uy1=circshift(uy,[1 0]); % uy(i-1,j)
    u=sqrt((0.5*(ux+ux1)).^2+(0.5*(uy+uy1)).^2);
    
end


