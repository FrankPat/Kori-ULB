function PlotMainFigure(ctr,par,x,y,sn,S0,H,u,B,MASK,glMASK,LSF)

% Kori-ULB
% Plot of the main figure during the model run with changes in ice
% thickness and total velocity field of the ice sheet

    if ctr.plotH==1
        sn1=sn;
    elseif ctr.plotH==2
        sn1=H;
        sn1(H==0)=NaN;
    else
        sn1=sn-S0;
        plim=max(max(abs(sn1)));
        plim(plim==0)=1;
    end
    u1=u;
    if ctr.glMASKexist==1
        sn1(glMASK==6)=NaN;
        u1(glMASK==6)=NaN;
    else
        sn1(MASK==0)=NaN;
        u1(MASK==0)=NaN;
    end
    ax1=subplot(1,2,1);
    if ctr.plotH>=1
        imagescn(x,y,sn1);
        if ctr.plotH==2
            title('Ice thickness (m)');
        else
            title('Surface elevation (m)');
        end
        colormap(ax1,crameri(par.color));
        colorbar;
    else
        imagescn(x,y,sn1,[-plim plim]);
        title('Surface elevation change (m)');
        colormap(ax1,crameri(par.dcolor));
        colorbar;
    end
    if ctr.glMASKexist==1 && ctr.plotGL==1
        hold on;
        contour(x,y,MASK,1,'LineColor','k','LineWidth',0.5);
        if ctr.calving>=1
            contour(x,y,LSF,[0 0],'LineColor','b','LineWidth',0.5);
        end
        hold off;
    end
    if ctr.mismip>=1
        hold on;
        contour(x,y,B,'LineColor','w','LineWidth',0.5);
        hold off;
    end
    axis xy;
    if ctr.mismip==0 || ctr.mismip==2
        axis equal;
    end
    axis tight;
    xlabel('x (km)');
    ylabel('y (km)');

    ax2=subplot(1,2,2);
    imagescn(x,y,log10(u1),[-1 3.5]);
    axis xy;
    if ctr.mismip==0 || ctr.mismip==2
        axis equal;
    end
    axis tight;
    colormap(ax2,crameri(par.color));
    colorbar;
    if ctr.glMASKexist==1 && ctr.plotGL==1
        hold on;
        contour(x,y,MASK,1,'LineColor','k','LineWidth',0.5);
        hold off;
    end
    if ctr.mismip>=1
        hold on;
        contour(x,y,B,'LineColor','w','LineWidth',0.5);
        hold off;
    end
    title('log_{10}(ice velocity) log_{10}(m a^{-1})');
    xlabel('x (km)');
    ylabel('y (km)');
    pause(0.0001);
end


