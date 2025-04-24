function DS = prep_schism_init(Mobj, dst)
% Prepare the initial data for SCHISM.
%
%% Syntax
% DS = prep_schism_init(Mobj, dst)
%
%% Description
% DS = prep_schism_init(Mobj, dst) prepares the initial data for SCHISM
%
%% Example
% DS = prep_schism_init(Mobj, 'hycom')
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% dst - dateset; char
%       the dateset name.
%
%% Output Arguments
% DS - data object; datastruct
%       the datastruct used to store the initial data.
%
%% Notes
% All variables are stored in the datastruct 'DS' as individual sub-structures,
% allowing each variable to have its own coordinate system. This design
% supports input from multiple data sources or products.
%
% This is a simple wrapper function, allowing users to integrate additional
% datasets easily, as long as they follow the same output format.
%
% Important format requirements:
% (1) Longitude and latitude vectors must be in strictly ascending order.
% (2) Variable dimensions must be either (lon × lat × depth) or (lon × lat).
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 22 Apr 2025.
% Email: wwu@vims.edu
%
% See also: prep_schism_bdry

%% Select the database
switch dst
    case 'hycom_online' % directly download hycom data from the internet
        varList = {'ssh','temp','salt'};
        C = get_hycom_online(Mobj.aimpath, Mobj.region, Mobj.time(1), varList);
        nVars = numel(varList);
        DS(nVars, 1) = struct('Variable', [], 'Data', [], 'Lon', [], 'Lat', [], 'Depth', [], 'Time', []);
        for iVar = 1:nVars
            clear D; 
            varName = varList{iVar};
            D.Variable = varName; D.Data = C.(varName);
            D.Lon = C.lon(:); D.Lat = C.lat(:); D.Depth = C.depth(:); D.Time = C.time;
            DS(iVar) = D;
        end

    case 'hycom_bys_clim'
        filepath = 'example_hycom_clim_1995_2020.mat'; % regional dataset for the BYS
        C = load(filepath); iMon = month(Mobj.time(1));
        D.Lon = C.lon; D.Lat = C.lat; D.Depth = abs(C.depth); D.Time = Mobj.time(1);
        varList = {'ssh', 'temp', 'salt'}; nVars = numel(varList);
        DS(nVars, 1) = struct('Variable', [], 'Data', [], 'Lon', [], 'Lat', [], 'Depth', [], 'Time', []);
        for ii = 1:nVars
            varName = varList{ii};
            DS(ii).Lon = D.Lon;  DS(ii).Lat = D.Lat; DS(ii).Time = D.Time;
            if strcmpi(varName, 'ssh'); DS(ii).Depth = 0;  else; DS(ii).Depth = D.Depth; end
            DS(ii).Variable = varName; DS(ii).Data = squeeze(C.(varName)(:,:,:,iMon));
        end

    case 'hycom_bys'
        meta_data.src_file = 'C:\Users\15641\OneDrive - vims.edu\GitHub_Projects\schism-toolbox\data\hycom\W117E127S33N41_****.mat';
        meta_data.raw_list = {'ssh','salt','temp'};      % original variable name in the file
        meta_data.new_list = {'ssh','salt','temp'};      % standard variable name for output
        meta_data.dim_vars = {'lon', 'lat', 'depth'};    % dimensional variables in the file
        meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
        DS = get_hycom_init(Mobj, meta_data);

    case 'hycom_chesbay'
        meta_data.src_file = 'E:\HYCOM\HYCOM_ChesBay\Wn79En72S33N41_****.mat';
        meta_data.raw_list = {'ssh','salt','temp'};      % original variable name in the file
        meta_data.new_list = {'ssh','salt','temp'};      % standard variable name for output
        meta_data.dim_vars = {'lon', 'lat', 'depth'};    % dimensional variables in the file
        meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
        DS = get_hycom_init(Mobj, meta_data);

    case 'glorys_chesbay'
        DS = get_glorys_init(Mobj);
end
end
