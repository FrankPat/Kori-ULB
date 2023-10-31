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

function mat2nc1D(input,output)

    % load mat file
    mat_file = load(input);
    t    = length(mat_file.time);
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
    nccreate(output,'t','Dimension',{'t' t});
    
    ncwrite(output,'t',linspace(0,mat_file.time(t,1),t));
    
    % create attributes
    ncwriteatt(output, 't', 'units', 'years');
    ncwriteatt(output, 't', '_CoordinateAxisType', 'time');
    
    for i=1:length(vars)
        index = 0;
        % create variable if available in the list
        [sy,sx] = size(mat_file.(vars{i}));
        if (sy==t && sx==1)
            nccreate(output,vars{i},'datatype' ,'double','Dimension',{'t' t}); 
            ncwrite(output,vars{i},mat_file.(vars{i}));
            % create attributes
            %ncwriteatt(output, variables(index), 'standard_name', stnd_name(index));
            %ncwriteatt(output, variables(index), 'long_name', long_name(index));
            %ncwriteatt(output, variables(index), 'units',units(index));
        end
    end
