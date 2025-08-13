function misfig(var,lim0,lim1)

% plot MISMIP+ figure by doubling domain across symmetry axis

    if nargin<1
        crameri;
    else
        n=size(var,1);
        var1=zeros(2*n-3,size(var,2));
        var1(1:n,:)=flip(var);
        var1(n:2*n-3,:)=var(3:n,:);
        FigHandle = figure('Position', [300, 100, 1500, 300]);
        imagescn(var1);
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

