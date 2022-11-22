function [Ax,Ay,Ad]=StaggeredA(A)

% Kori-ULB
% Flow parameter A on different staggered grids. Originally, A is determined
% on the H-grid

    A1=circshift(A,[0 -1]); % A(i,j+1)
    A2=circshift(A,[-1 0]); % A(i+1,j)
    
    Ax=(A+A1)/2.; % A on u-grid
    Ay=(A+A2)/2.; % A on v-grid
    Ad=h2d(A); % A on d-grid

end


