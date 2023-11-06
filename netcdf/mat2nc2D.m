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

function mat2nc2D(input,output,dx)

    % load mat file
    mat_file = load(input);
    
    % load dimension variables
    [y,x] = size(mat_file.H);
    %[y,x] = size(mat_file.H);
    
    % load file names 
    vars = fieldnames(mat_file);
    
    % --- Dimension variables ---
    % create .nc file
    % if variable exist overwrite, if not create
    if isfile(output)
        quest = 'Output file already exists. Do you wish to overwritte it?';
        answer = questdlg(quest,'Overwrite',...
                          'Yes','No','Yes');
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
    
    for i=1:length(vars)
        % create variable if size coincides with dimensions
        [sy,sx] = size(mat_file.(vars{i}));
        if (sy==y && sx==x)
            nccreate(output,vars{i},'datatype' ,'double','Dimension',{'x' x 'y' y}); 
            ncwrite(output,vars{i},transpose(mat_file.(vars{i})));
            % create attributes
            %ncwriteatt(output, variables(index), 'standard_name', stnd_name(index));
            %ncwriteatt(output, variables(index), 'long_name', long_name(index));
            %ncwriteatt(output, variables(index), 'units',units(index));
        end
    end
