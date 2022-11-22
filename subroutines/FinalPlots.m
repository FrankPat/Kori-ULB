function FinalPlots(ctr,outfile,time,InvVol,SLC,u,v,MASK,bMASK,glMASK)

% Kori-ULB
% Produce basic plots at the end of the model run
% When inverse>0, plots of the quality of the inversion
% When velocity data is availbale, plot of observed versus modelled
% velocity.

    if ctr.inverse>=1
        figure('Name',['InvVol: ' outfile],'NumberTitle','off');
        plot(time,InvVol(:,1)*ctr.delta^2.*1e-9);
        grid on; hold on;
        plot(time,InvVol(:,2)*ctr.delta^2.*1e-9,'r');
        plot(time,InvVol(:,3)*ctr.delta^2.*1e-9,'k');
        xlabel('Time (year)');
        ylabel('|H-H_0| (km^3)');
        legend('Abs. Misfit','Misfit');
    end

    figure('Name',['SLC: ' outfile],'NumberTitle','off');
    plot(time,SLC); grid on; 
    xlabel('Time (year)');
    ylabel('SLR contribution (m)');
    
    if ctr.vexist==1
        figure('Name',['Velocity: ' outfile],'NumberTitle','off');
        if ctr.basin==1
            loglog(v(MASK==1 & bMASK==0),u(MASK==1 & bMASK==0),'.');
        else
            loglog(v(MASK==1),u(MASK==1),'.');
        end
        hold on;
        if ctr.glMASKexist==1
            loglog(v(glMASK==4),u(glMASK==4),'r.');
        end
        plot([1e-1 4e3],[1e-1 4e3],'k-','linewidth',2);
        axis([1e-1 4e3 1e-1 4e3]);
        grid on;
        xlabel('Observed velocity (m a^{-1})');
        ylabel('Modelled velocity (m a^{-1})');
    end
end


