%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      %
%  Function to convert output of the Kori model into   %
%                    netcdf files                      %
%                 --- ANTARCTICA ---                   %
%                                                      %
%              Javier Blasco Navarro                   %
%                       ULB                            %
%                   January 2022                       %
%                                                      %
%              Input: Kori output .mat                 %
%                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat2nc(input,output,dx)

    % load mat file
    mat_file = load(input);
    [y,x] = size(mat_file.H);
    vars = fieldnames(mat_file);
    
    % --- Dimension variables ---
    % create .nc file
    % if variable exist overwrite, if not create
    if isfile(output)
        quest = 'Output file already exists. Do you wish to overwritte it?';
        answer = questdlg(quest,'Overwrite',...
                          'Yes','No','No');
        switch answer
            case 'No'
                return; % Or break or continue
        end
        delete(output);
    end
    nccreate(output,'y','Dimension',{'y' y});
    nccreate(output,'x','Dimension',{'x' x});
    
    % dx is an optional parameter (resolution in km)
    if exist('dx','var')
        % y dimensions
        ncwrite(output,'y',linspace(-0.5*(y-1)*dx,0.5*(y-1)*dx,y));
        % x dimensions
        ncwrite(output,'x',linspace(-0.5*(x-1)*dx,0.5*(x-1)*dx,x));
    else
        % y dimensions
        ncwrite(output,'y',linspace(0,y,y));
        % x dimensions
        ncwrite(output,'x',linspace(0,x,x));
    end
    
    % t dimension
    % nccreate(output,'time','Dimension',{'time' t});
    % ncwrite(output,'time',??:??:??);
    % create attributes
    ncwriteatt(output, 'x', 'units', 'km');
    ncwriteatt(output, 'x', '_CoordinateAxisType', 'x');
    ncwriteatt(output, 'y', 'units', 'km');
    ncwriteatt(output, 'y', '_CoordinateAxisType', 'y');
    %ncwriteatt(output, 'time', 'units', 'years');
    %ncwriteatt(output, 'time', '_CoordinateAxisType', 'time');
    
    % --- General variables ----
    variables = ["lat" "lon" "v" "vx" "vy" "uxssa" "uyssa"...
	         "u" "ux" "uy" ...
                 "H" "Ho" "B" "Bo" "Bor"...
                 "stdB" "MASK" "MASKo" "Db" "G"...
                 "Btau" "Mb" "Ts" "Pr" "Evp"...
                 "To" "So" "ZB" "Bmelt" "Melt"...
		 "A" "eta" "dudx" "dudy" "dvdx" "dvdy"  "scale_eta"];
    
    stnd_name = ["Latitude coordinates" "Longitude coordinates" "Obs. Surface velocity" "Obs. Horizontal surface velocity" "Obs. Vertical surface velocity" "Horizontal integrated velocity (SSA)" "Vertical integrated velocity (SSA)"...
                 "Surface velocity" "Horizontal surface velocity" "Vertical surface velocity" ... 
	         "Ice thickness" "Observed ice thickness" "Bedrock elevation (saturated)" "Observed bedrock elevation" "Bedrock elevation"...
                  "Bedrock variability" "Ice mask" "Observed ice mask" "Lithospheric thickness" "Geothermal heatflow"...
                  "Asthenosphere relaxation time" "Surface mass balance" "Mean annual surface temperature" "Annual precipitation" "Evaporation"...
                  "Ocean temperature" "Ocean salinity" "Zwally basins" "Basal melt" "Basal melt (grounded)"...
		  "" "" "" "" "" "" "Damage percent"];
    
    long_name = ["Latitude coordinates" "Longitude coordinates" "Obs. Surface velocity" "Obs. Horizontal surface velocity" "Obs. Vertical surface velocity" "Horizontal integrated velocity (SSA)" "Vertical integrated velocity (SSA)"...
                 "Surface velocity" "Horizontal surface velocity" "Vertical surface velocity" ... 
	         "Ice thickness" "Observed ice thickness" "Bedrock elevation (saturated)" "Observed bedrock elevation" "Bedrock elevation"...
                 "Bedrock variability (std dev representing bed roughness)" "Ice mask (grounded = 1; floating = 0)" "Observed ice mask (invariable in time)" "Flexural rigidity of the lithosphere" "Geothermal heatflow"...
                  "Asthenosphere relaxation time" "Surface mass balance" "Mean annual surface temperature" "Annual precipitation" "Evaporation"...
                  "Ocean temperature" "Ocean salinity" "Zwally basins" "Basal melt" "Basal melt (grounded)"...
		  "" "" "" "" "" "" "Damage percent (eta scaling)"];
    
    units      = ["Degrees north" "Degrees east" "m/yr" "m/yr" "m/yr" "m/yr" "m/yr"...
                  "m/yr" "m/yr" "m/yr" ...
	          "m" "m" "m" "m" "m"...
                  "m" "#" "#" "?" "W m −2"...
                  "?" "m a-1" "degC" "?" "?"...
                  "degC" "PSU" "#" "m/yr" "m/yr"...
		  "Pa-n yr-1" "Pa yr m" "m/yr" "m/yr" "m/yr" "m/yr" "%"];
    
    for i=1:length(vars)
        index = 0;
        % create variable if available in the list
        index = find(ismember(variables,vars(i)));
        type = isa(mat_file.(vars{i}),'logical');
        if (index>0 & type==0)
            nccreate(output,variables(index),'datatype' ,'double','Dimension',{'x' x 'y' y}); 
            ncwrite(output,variables(index),transpose(mat_file.(vars{i})));
            % create attributes
            ncwriteatt(output, variables(index), 'standard_name', stnd_name(index));
            ncwriteatt(output, variables(index), 'long_name', long_name(index));
            ncwriteatt(output, variables(index), 'units',units(index));
        end
    end
