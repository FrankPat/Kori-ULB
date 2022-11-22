function myevol(infile,var,vid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot time evolution of variable and create video (optional)
% f.ETISh v1.8
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
        video=VideoWriter('fETIShVideo','Uncompressed AVI');
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
        if var=='m' || var=='cm'
            plim0=0;
        end
    end
    imagescn(x,y,dvar,[plim0 plim1]);
    title(titlename);
    if var=='sb'
        colormap(ax1,crameri('oleron','pivot',fc.DeltaSL(ti)));
    else
        colormap(ax1,crameri('vik'));
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


function cmap = crameri(ColormapName,varargin) 
% crameri returns perceptually-uniform scientific colormaps created
% by Fabio Crameri. 
% 
% Syntax 
% 
%  crameri 
%  cmap = crameri('ColormapName') 
%  cmap = crameri('-ColormapName') 
%  cmap = crameri(...,NLevels)
%  cmap = crameri(...,'pivot',PivotValue) 
%  crameri(...)
% 
% Description 
% 
% crameri without any inputs displays the options for colormaps. 
% 
% cmap = crameri('ColormapName') returns a 256x3 colormap.  For a visual
% depiction of valid colormap names, type |crameri|. 
%
% cmap = crameri('-ColormapName') a minus sign preceeding any ColormapName flips the
% order of the colormap. 
%
% cmap = crameri(...,NLevels) specifies a number of levels in the colormap.  Default
% value is 256. 
%
% cmap = crameri(...,'pivot',PivotValue) centers a diverging colormap such that white 
% corresponds to a given value and maximum extents are set using current caxis limits. 
% If no PivotValue is set, 0 is assumed. 
%
% crameri(...) without any outputs sets the current colormap to the current axes.  
% 
% Examples 
% For examples, type: 
% 
%  showdemo crameri_documentation
%
% Author Info 
% This function was written by Chad A. Greene of the University of Texas
% Institute for Geophysics (UTIG), August 2018, using Fabio Crameri's 
% scientific colormaps, version 4.0. http://www.fabiocrameri.ch/colourmaps.php
% 
% Citing this colormap: 
% Please acknowledge the free use of these colormaps by citing
% 
% Crameri, F. (2018). Scientific colour-maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
% 
% Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation and 
% StagLab 3.0, Geosci. Model Dev., 11, 2541-2562, doi:10.5194/gmd-11-2541-2018.
% 
% For more on choosing effective and accurate colormaps for science, be sure
% to enjoy this fine beach reading: 
% 
% Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True 
% colors of oceanography: Guidelines for effective and accurate colormap selection. 
% Oceanography 29(3):9-13, http://dx.doi.org/10.5670/oceanog.2016.66.
% 
% See also colormap and caxis.  

% Display colormap options: 

if nargin==0
   figure('menubar','none','numbertitle','off','Name','crameri options:')
   
   if license('test','image_toolbox')
      imshow(imread('crameri6.0.png')); 
   else
      axes('pos',[0 0 1 1])
      image(imread('crameri6.0.png')); 
      axis image off
   end
   
   return
end

% Error checks: 

assert(isnumeric(ColormapName)==0,'Input error: ColormapName must be a string.') 

% Set defaults: 

NLevels = 256; 
autopivot = false; 
PivotValue = 0; 
InvertedColormap = false; 

% Parse inputs: 

% Does user want to flip the colormap direction? 
dash = strncmp(ColormapName,'-',1); 
if any(dash) 
   InvertedColormap = true; 
   ColormapName(dash) = []; 
end

% Standardize all colormap names to lowercase: 
ColormapName = lower(ColormapName); 

% Oleron's too hard for me to remember, so I'm gonna use dem or topo. 
if ismember(ColormapName,{'dem','topo'})
   ColormapName = 'oleron'; 
end

% Does the user want to center a diverging colormap on a specific value? 
% This parsing support original 'zero' syntax and current 'pivot' syntax. 
tmp = strncmpi(varargin,'pivot',3); 
if any(tmp) 
   autopivot = true; 
   try
      if isscalar(varargin{find(tmp)+1})
         PivotValue = varargin{find(tmp)+1}; 
         tmp(find(tmp)+1) = 1; 
      end
   end
   varargin = varargin(~tmp); 
end

% Has user requested a specific number of levels? 
tmp = isscalar(varargin); 
if any(tmp) 
   NLevels = varargin{tmp}; 
end

% Load RGB values and interpolate to NLevels: 

try
   S = load('CrameriColourMaps6.0.mat',ColormapName); 
   cmap = S.(ColormapName); 
catch
   error(['Unknown colormap name ''',ColormapName,'''. Try typing crameri with no inputs to check the options and try again.'])
end

% Interpolate if necessary: 
if NLevels~=size(cmap,1) 
   cmap = interp1(1:size(cmap,1), cmap, linspace(1,size(cmap,1),NLevels),'linear');
end

% Invert the colormap if requested by user: 

if InvertedColormap
   cmap = flipud(cmap); 
end

% Adjust values to current caxis limits? 

if autopivot
   clim = caxis; 
   maxval = max(abs(clim-PivotValue)); 
   cmap = interp1(linspace(-maxval,maxval,size(cmap,1))+PivotValue, cmap, linspace(clim(1),clim(2),size(cmap,1)),'linear');
end

% Clean up 

if nargout==0
   colormap(gca,cmap) 
   clear cmap  
end

end


function h = imagescn(varargin) 
% imagescn behaves just like imagesc, but makes NaNs transparent, sets
% axis to xy (aka ydirection normal) if xdata and ydata are included, and has a little more 
% error checking than imagesc. 
% 
% Syntax 
% 
%  imagescn(C) 
%  imagescn(x,y,C) 
%  imagescn(x,y,C,clims) 
%  imagescn('PropertyName',PropertyValue,...) 
%  h = imagescn(...) 
% 
% Description 
% 
% imagescn(C) displays the data in array C as an image that uses the full range of colors in the colormap. 
% Each element of C specifies the color for 1 pixel of the image. The resulting image is an m-by-n grid of 
% pixels where m is the number of columns and n is the number of rows in C. The row and column indices of 
% the elements determine the centers of the corresponding pixels. NaN values in C appear transparent. 
% 
% imagescn(x,y,C) specifies x and y locations of the centers of the pixels in C. If x and y are two-element
% vectors, the outside rows and columns of C are centered on the values in x and y. Mimicking imagesc, if 
% x or y are vectors with more than two elements, only the first and last elements of of the vectors are 
% considered, and spacing is automatically scaled as if you entered two-element arrays. The imagescn function
% takes this one step further, and allows you to enter x and y as 2D grids the same size as C. If x and y
% are included, the imagescn function automatically sets axes to cartesian xy rather than the (reverse) ij axes. 
% 
% imagescn(x,y,C,clims) specifies the data values that map to the first and last elements of the colormap. 
% Specify clims as a two-element vector of the form [cmin cmax], where values less than or equal to cmin 
% map to the first color in the colormap and values greater than or equal to cmax map to the last color in 
% the colormap.
% 
% imagescn('PropertyName',PropertyValue,...) specifies image properties as name-value pairs. 
% 
% h = imagescn(...) retrns a handle of the object created. 
% 
% Differences between imagesc, imagescn, and pcolor
% The imagescn function plots data with imagesc, but after plotting, sets NaN pixels to an 
% alpha value of 0. The imagesc function allows input coordinates x and y to be grids, which 
% are assumed to be evenly-spaced and monotonic as if created by meshgrid. If x and y data 
% are included when calling imagescn, y axis direction is changed from reverse to normal. 
% 
% The imagescn function is faster than pcolor. Pcolor (nonsensically) deletes an outside row 
% and column of data, and pcolor also refuses to plot data points closest to any NaN holes. 
% The imagescn function does not delete any data.  However, you may still sometimes wish to 
% use pcolor if x,y coordinates are not evenly spaced or if you want interpolated shading. 
% 
% Examples 
% For examples, type 
% 
%  cdt imagescn
% 
% Author Info 
% 
% This function was written by Chad A. Greene of the University of 
% Texas Institute for Geophysics (UTIG), January 2017. 
% http://www.chadagreene.com 
% 
% See also imagesc, image, and pcolor.
% The imagesc function does not have error checking regarding number of elements
% in xdata, ydata versus number of elements in the input image, so I'm gonna add
% some error checking: 

% Check inputs: 

xydata = false; 
if nargin>2
   if all([isnumeric(varargin{1}) isnumeric(varargin{2})]) 
      % This is an assumption that should typically be safe to make: 
      xydata = true; 
      
      % Determine if input coordinates are meshgrid type and if so, convert to vector: 
      if isequal(size(varargin{1}),size(varargin{2}),size(varargin{3}))
         X = varargin{1}; 
         Y = varargin{2}; 
         
         varargin{1} = [X(1,1) X(end,end)]; 
         varargin{2} = [Y(1,1) Y(end,end)]; 
      end
   end
end

% Plot

% Plot imagesc: 
h = imagesc(varargin{:}); 

% Make NaNs transparent: 
cd = get(h,'CData'); 
set(h,'alphadata',isfinite(cd)); 

if xydata
   axis xy
end

if nargout==0
   clear h
end

end

