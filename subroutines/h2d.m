function y=h2d(x)

% Kori-ULB
% Interpolate values initialiiy on d-grid to h-grid

    y=0.25*(x+circshift(x,[0 -1])+circshift(x,[-1 -1])+circshift(x,[-1 0]));
    
end


