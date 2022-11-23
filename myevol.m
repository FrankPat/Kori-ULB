function myevol(infile,var,vid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot time evolution of variable and create video (optional)
% Kori-ULB
%
% call:
%   myevol(infile, var, vid)
%
% infile:
%   name of output file (without toto extension)
%
% var:
%   h: ice thickness change
%   sn: surface elevation
%   sb: surface and ocean topography combined
%   b: bedrock elevation
%   u: ice velocity
%   t: basal temperature
%   n: effective pressure
%   sl: local sea level
%   m: sub-shelf melt rate
%   dh: dH/dt
%
% vid: (optional value)
%   v: make video of animation
%   s: variable color scale
%   w: both 'v' and 's'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all;

inputname=[infile,'_toto'];
if exist([inputname,'.mat'],'file')
    load(inputname);
end

if ctr.timeslice==0
    fprintf('The model run did not record\n');
    fprintf('time slices (ctr.timeslice=0)\n');
    return;
end

makevid=0;
fixedscale=1;

if nargin<3
    makevid=0;
    fixedscale=1;
else
    if vid=='v' || vid=='w'
        makevid=1;
        video=VideoWriter('KoriVideo','Uncompressed AVI');
        video.FrameRate=4;
        GF=[];
    end
    if vid=='s' || vid=='w'
        fixedscale=0;
    end
end

% link time to snapshot
plotst=floor(ctr.nsteps/(ctr.snapshot-1));

EVmask=readdata(ctr,outfile,'MASK');
EVmask1=EVmask;
EVH=readdata(ctr,outfile,'H');
EVmask1(EVmask==0 & EVH<=par.SeaIceThickness)=NaN;

maxvar=0;
minvar=0;

% determine variable to plot

switch var
    case 'n'
        plotvar=readdata(ctr,outfile,'Neff');
        var0=zeros(ctr.imax,ctr.jmax);
        for i=1:ctr.snapshot
            maxvar=max(maxvar,max(max(plotvar(:,:,i)-var0)));
            minvar=min(minvar,min(min(plotvar(:,:,i)-var0)));
        end
        plotvar=log10(plotvar);
        plim0=5;
        plim1=8;
        titlename='Effective pressure (Pa)';

    case 'h'
        plotvar=readdata(ctr,outfile,'H');
        var0=H0;
        for i=1:ctr.snapshot
            maxvar=max(maxvar,max(max(plotvar(:,:,i)-var0)));
            minvar=min(minvar,min(min(plotvar(:,:,i)-var0)));
        end
        plim0=-max(abs(maxvar),abs(minvar));
        plim1=max(abs(maxvar),abs(minvar));
        titlename='Ice thickness change (m)';
        
    case 'b'
        plotvar=readdata(ctr,outfile,'B');
        var0=B0;
        for i=1:ctr.snapshot
            maxvar=max(maxvar,max(max(plotvar(:,:,i)-var0)));
            minvar=min(minvar,min(min(plotvar(:,:,i)-var0)));
        end
        plim0=-max(abs(maxvar),abs(minvar));
        plim1=max(abs(maxvar),abs(minvar));
        titlename='Bedrock elevation change (m)';
        
    case 'm'
        plotvar=readdata(ctr,outfile,'Melt');
        var0=zeros(ctr.imax,ctr.jmax);
        if ctr.glMASKexist==1
            plotvar(EVmask==1)=NaN;
            plotvar(isnan(EVmask1))=NaN;
        end
        for i=1:ctr.snapshot
            maxvar=max(maxvar,max(max(plotvar(:,:,i)-var0)));
            minvar=min(minvar,min(min(plotvar(:,:,i)-var0)));
        end
        plim0=minvar;
        plim1=maxvar;
        titlename='Sub-shelf melt (m a^{-1})';
        
    case 'dh'
        plotvar=zeros(ctr.imax,ctr.jmax,ctr.snapshot);
        for i=2:ctr.snapshot
            plotvar(:,:,i)=(EVH(:,:,i)-EVH(:,:,i-1))/(ctr.dt*plotst);
        end
        var0=zeros(ctr.imax,ctr.jmax);
        if ctr.glMASKexist==1
            plotvar(isnan(EVmask1))=NaN;
        end
        for i=1:ctr.snapshot
            maxvar=max(maxvar,max(max(plotvar(:,:,i)-var0)));
            minvar=min(minvar,min(min(plotvar(:,:,i)-var0)));
        end
        plim0=-max(abs(maxvar),abs(minvar));
        plim1=max(abs(maxvar),abs(minvar));
        titlename='Ice thickness change (m a^{-1})';

    case 'u'
        EVUx=readdata(ctr,outfile,'ux');
        EVUy=readdata(ctr,outfile,'uy');
        plotvar=zeros(ctr.imax,ctr.jmax,ctr.snapshot);
        for i=1:ctr.snapshot
            plotvar(:,:,i)=velocity_on_grid(EVUx(:,:,i),EVUy(:,:,i));
        end
        var0=zeros(ctr.imax,ctr.jmax);
        if ctr.glMASKexist==1
            plotvar(isnan(EVmask1))=NaN;
        end
        plotvar=log10(plotvar);
        plim0=-1;
        plim1=3.5;
        titlename='log_{10}(ice velocity) log_{10}(m a^{-1})';
        
    case 'sn'
        EVB=readdata(ctr,outfile,'B');
        plotvar=zeros(ctr.imax,ctr.jmax,ctr.snapshot);
        for i=1:ctr.snapshot
            ti=(i-1)*plotst+1;
            HB=max(fc.DeltaSL(ti)-par.rho/par.rhow*EVH(:,:,i),EVB(:,:,i));
            plotvar(:,:,i)=EVH(:,:,i)+HB;
        end
        var0=zeros(ctr.imax,ctr.jmax);
        if ctr.glMASKexist==1
            plotvar(isnan(EVmask1))=NaN;
        end
        for i=1:ctr.snapshot
            maxvar=max(maxvar,max(max(plotvar(:,:,i)-var0)));
            minvar=min(minvar,min(min(plotvar(:,:,i)-var0)));
        end
        plim0=minvar;
        plim1=maxvar;
        titlename='Surface elevation (m)';

    case 'sb'
        EVB=readdata(ctr,outfile,'B');
        plotvar=zeros(ctr.imax,ctr.jmax,ctr.snapshot);
        for i=1:ctr.snapshot
            ti=(i-1)*plotst+1;
            B=EVB(:,:,i);
            HB=max(fc.DeltaSL(ti)-par.rho/par.rhow*EVH(:,:,i),B);
            HB(isnan(EVmask1(:,:,i)))=B(isnan(EVmask1(:,:,i)));
            plotvar(:,:,i)=EVH(:,:,i)+HB;
        end
        var0=zeros(ctr.imax,ctr.jmax);
        for i=1:ctr.snapshot
            maxvar=max(maxvar,max(max(plotvar(:,:,i)-var0)));
            minvar=min(minvar,min(min(plotvar(:,:,i)-var0)));
        end
        plim0=minvar;
        plim1=maxvar;
        titlename='Surface elevation (m)';

    case 't'
        plotvar=readdata(ctr,outfile,'Tbc');
        if ctr.glMASKexist==1
            plotvar(isnan(EVmask1))=NaN;
        end
        var0=zeros(ctr.imax,ctr.jmax);
        for i=1:ctr.snapshot
            maxvar=max(maxvar,max(max(plotvar(:,:,i)-var0)));
            minvar=min(minvar,min(min(plotvar(:,:,i)-var0)));
        end
        plim0=minvar;
        plim1=maxvar;
        titlename='Basal temperature (°C)';

    case 'sl'
        plotvar=readdata(ctr,outfile,'SLR');
        if ctr.glMASKexist==1
            plotvar(EVmask==1)=NaN;
        end
        var0=zeros(ctr.imax,ctr.jmax);
        for i=1:ctr.snapshot
            maxvar=max(maxvar,max(max(plotvar(:,:,i)-var0)));
            minvar=min(minvar,min(min(plotvar(:,:,i)-var0)));
        end
        plim0=-max(abs(maxvar),abs(minvar));
        plim1=max(abs(maxvar),abs(minvar));
        titlename='Sea level change (m)';
    otherwise
        fprintf('Variable does not exist\n');
        return;
end

scrsz=get(groot,'ScreenSize');
fig=figure('Position',[30 scrsz(4)/2.5 scrsz(3)/1.1 scrsz(4)/2]);
if makevid==1
    GF=[GF getframe(fig)];
end

yy=smooth(diff(SLC)/ctr.dt,20);
for i=1:ctr.snapshot
    subplot(1,2,1);
    yyaxis left;
    plot(time,SLC,'b','LineWidth',2); hold on;
    plot(time,VAF,'k','LineWidth',2);
    plot(time,(IVg(1)-IVg)*par.rho/(par.Aoc*par.rhow),'g','LineWidth',2);
    ti=(i-1)*plotst+1;
    plot(time(ti),SLC(ti),'bo','MarkerFaceColor','b');
    plot(time(ti),VAF(ti),'ko','MarkerFaceColor','k');
    plot(time(ti),(IVg(1)-IVg(ti))*par.rho/(par.Aoc*par.rhow),'go','MarkerFaceColor','g');
    title('Sea level contribution');
    xlabel('Time (years)');
    ylabel('Sea-level equivalent (m)');
    hold off;
    yyaxis right;
    plot(time(5:end-3),yy(4:end-3),'color',[1 0.6 0.6]);
%     plot(time,Ag/1e6,'r');
    legend('SLC','VAF','IVg','location','northwest');
    ylabel('Rate of SL change (m a^{-1})');
%     ylabel('Grounded Ice Area (km^2)');
    grid on;

    ax1=subplot(1,2,2);
    dvar=plotvar(:,:,i)-var0;
    if fixedscale==0
        plim1=nanmax(abs(dvar(:)));
        plim1(plim1==0)=1;
        plim0=-plim1;
        if var=='m'
            plim0=0;
        end
    end
    imagescn(x,y,dvar,[plim0 plim1]);
    title(titlename);
    if var=='sb'
        colormap(ax1,crameri('oleron','pivot',fc.DeltaSL(ti)));
    elseif var=='u'
        colormap(ax1,crameri('imola'));
    else
        colormap(ax1,crameri('broc'));
    end
    colorbar;
    if ctr.glMASKexist==1
        hold on;
        contour(x,y,EVmask(:,:,i),[1],'LineColor','k','LineWidth',0.5);
        hold off;
    end
    axis equal;
    axis tight;
    xlabel('x (km)');
    ylabel('y (km)');
    pause(0.0001);

    if makevid==1
        GF=[GF getframe(fig)];
    end
end

if makevid==1
    open(video);
    writeVideo(video,GF);
    close(video);
end

end


function EV=readdata(ctr,outfile,var)

    EV=zeros(ctr.imax,ctr.jmax,ctr.snapshot);
    for i=1:ctr.snapshot
        fname=strcat(outfile,'_',num2str(i-1,'%03i'));
        load(fname,var);
        EV(:,:,i)=eval(var);
    end

end


function u=velocity_on_grid(ux,uy)

% Determine magnitude of velocities based on velocities on the staggered u
% and v-grids. Only used for visualization purposes

    ux1=circshift(ux,[-1 0]); % ux(i+1,j)
    uy1=circshift(uy,[0 -1]); % uy(i,j+1)
    u=sqrt((0.5*(ux+ux1)).^2+(0.5*(uy+uy1)).^2);
    
end



