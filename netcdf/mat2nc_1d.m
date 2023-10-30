%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      %
% Function to convert output of the f.ETISh model into %
%                   netcdf files                       %
%                 --- ANTARCTICA ---                   %
%                                                      %
%              Javier Blasco Navarro                   %
%                       ULB                            %
%                   January 2022                       %
%                                                      %
%              Input: Kori output .mat                 %
%                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat2nc_1d(input,output)

    % load mat file
    mat_file = load(input);
    [t,num1] = size(mat_file.time);
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
    nccreate(output,'t','Dimension',{'t' t});
    
    ncwrite(output,'t',linspace(0,mat_file.time(t,1),t));
    
    % t dimension
    % nccreate(output,'time','Dimension',{'time' t});
    % ncwrite(output,'time',??:??:??);
    % create attributes
    ncwriteatt(output, 't', 'units', 'years');
    ncwriteatt(output, 't', '_CoordinateAxisType', 'time');
    %ncwriteatt(output, 'time', 'units', 'years');
    %ncwriteatt(output, 'time', '_CoordinateAxisType', 'time');
    
    % --- General variables ----
    variables = ["VAF","Ag","SLC"];
    
    stnd_name = ["Volume above flotation","Grounded ice area","Sea-level contribution"];
    
    long_name = ["Volume above flotation","Grounded ice area","Sea-level contribution"];
    
    units      = ["?","km2","m"];
    
    for i=1:length(vars)
        index = 0;
        % create variable if available in the list
        index = find(ismember(variables,vars(i)));
        type = isa(mat_file.(vars{i}),'logical');
        if (index>0 & type==0)
            nccreate(output,variables(index),'datatype' ,'double','Dimension',{'t' t}); 
            ncwrite(output,variables(index),mat_file.(vars{i}));
            % create attributes
            ncwriteatt(output, variables(index), 'standard_name', stnd_name(index));
            ncwriteatt(output, variables(index), 'long_name', long_name(index));
            ncwriteatt(output, variables(index), 'units',units(index));
        end
    end
