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
%              Input: fetish output .mat               %
%                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat2nc_3d(input,output)

    % load mat file
    mat_file = load(input);
    [x,y] = size(mat_file.lon)

    % create .nc file
    % first the basic dimensions
    nccreate(output,'y','Dimension',{'y' y});
    nccreate(output,'x','Dimension',{'x' x});
    % nccreate(output,'time','Dimension',{'time' t});
    % create 2D/3D variables
    nccreate(output,'lat2D','datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'lon2D','datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'B'    ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'stdB' ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'H'    ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    %nccreate(output,'MASK' ,'datatype' ,'int'   ,'Dimension',{'x' x 'y' y});
    nccreate(output,'MASK' ,'datatype' ,'double'   ,'Dimension',{'x' x 'y' y});
    nccreate(output,'Bor'  ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'Db'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'Btau' ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'To'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'So'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'ZB'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'Mb'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'Ts'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'G'    ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'vx'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'vy'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'v'    ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'Ho'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'Bo'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    %nccreate(output,'MASKo','datatype' ,'int'   ,'Dimension',{'x' x 'y' y});
    nccreate(output,'MASKo','datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'Pr'   ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    nccreate(output,'Evp'  ,'datatype' ,'double','Dimension',{'x' x 'y' y});
    % save
    ncdisp(output);

    % x dimensions
    ncwrite(output,'y',-3040:32:3040);
    % y dimensions
    ncwrite(output,'x',-3040:32:3040);
    % t dimension
    %ncwrite(output,'time',??:??:??);

    % write 1D data

    % write 2D data
    ncwrite(output,'lat2D',transpose(mat_file.lat));
    ncwrite(output,'lon2D',transpose(mat_file.lon));
    ncwrite(output,'B'    ,transpose(mat_file.B));
    ncwrite(output,'stdB' ,transpose(mat_file.stdB));
    ncwrite(output,'H'    ,transpose(mat_file.H));
    ncwrite(output,'MASK' ,transpose(mat_file.MASK));
    ncwrite(output,'Bor'  ,transpose(mat_file.Bor));
    ncwrite(output,'Db'   ,transpose(mat_file.Db));
    ncwrite(output,'Btau' ,transpose(mat_file.Btau));
    ncwrite(output,'To'   ,transpose(mat_file.To));
    ncwrite(output,'So'   ,transpose(mat_file.So));
    ncwrite(output,'ZB'   ,transpose(mat_file.ZB));
    ncwrite(output,'Mb'   ,transpose(mat_file.Mb));
    ncwrite(output,'Ts'   ,transpose(mat_file.Ts));
    ncwrite(output,'G'    ,transpose(mat_file.G));
    ncwrite(output,'vx'   ,transpose(mat_file.vx));
    ncwrite(output,'vy'   ,transpose(mat_file.vy));
    ncwrite(output,'v'    ,transpose(mat_file.v));
    ncwrite(output,'Ho'   ,transpose(mat_file.Ho));
    ncwrite(output,'Bo'   ,transpose(mat_file.Bo));
    ncwrite(output,'MASKo',transpose(mat_file.MASKo));
    ncwrite(output,'Pr'   ,transpose(mat_file.Pr));
    ncwrite(output,'Evp'  ,transpose(mat_file.Evp));
    % create attributes
    %ncwriteatt(output, 'y', 'Vertical resolution','km');
    %ncwriteatt(output, 'x', 'Horizontal resolution','km');
    ncwriteatt(output, 'x', 'units', 'km');
    ncwriteatt(output, 'x', '_CoordinateAxisType', 'x');
    ncwriteatt(output, 'y', 'units', 'km');
    ncwriteatt(output, 'y', '_CoordinateAxisType', 'y');
    ncwriteatt(output, 'H', 'standard_name', 'Ice thickness');
    ncwriteatt(output, 'H', 'long_name', 'Ice thickness');
    ncwriteatt(output, 'H', 'units','m');    
    
    
end
