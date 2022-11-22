function myfig(var,lim0,lim1)

    if nargin<1
        crameri;
    else
        FigHandle = figure('Position', [300, 100, 700, 600]);
        imagescn(var);
        if nargin>1
            caxis([lim0 lim1]);
        end
        axis xy;
        axis equal;
        axis tight;
        colormap(crameri('vik'));
        colorbar;
    end
end

