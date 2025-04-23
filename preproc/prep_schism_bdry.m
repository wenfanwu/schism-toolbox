function DS = prep_schism_bdry(Mobj, dst, obc_bnds)
% Prepare the boundary iniputs for SCHISM.
%
%% Syntax
% DS = prep_schism_bdry(Mobj, dst)
% DS = prep_schism_bdry(Mobj, dst, obc_bnds)
%
%% Description
% DS = prep_schism_bdry(Mobj, dst) prepares the boundry inputs for SCHISM.
% DS = prep_schism_bdry(Mobj, dst, obc_bnds) specifies the open boundaries.
%
%% Example
% obc_bnds = 1:Mobj.obc_counts;
% DS = prep_schism_bdry(Mobj, 'hycom_bys', obc_bnds);
%
% obc_bnds = {1, [1 2], [1 2], [1 2], [1 2]};
% DS = prep_schism_bdry(Mobj, 'hycom_chesbay', obc_bnds);
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% dst - dateset; char
%       the dateset name.
% obc_bnds - the open boundary index; numeric/cell
%       the index of open boundary segments, e.g., obc_bnds =
%       1:Mobj.obc_counts, indicating  that all open boundaries are used.
%       obc_bnds can also be specified as "cell", which means that the
%       variables use different open boundaries.
%
%% Output Arguments
% DS - data object; datastruct
%       the datastruct used to store the boundary inputs.
%
%% Notes
% This is a simple wrapper function, allowing users to integrate additional
% datasets easily, as long as they follow the same output format.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 22 Apr 2025.
% Email: wwu@vims.edu
%
% See also: prep_schism_init

%% Parse inputs
if nargin<3; obc_bnds = 1:Mobj.obc_counts; end

%% Select the database
switch dst
    case 'hycom_bys'
        meta_data.src_file = 'C:\Users\15641\OneDrive - vims.edu\GitHub_Projects\schism-toolbox\data\hycom\W117E127S33N41_****.mat';   % the source files
        meta_data.src_time = unique(dateshift(Mobj.time, 'start', 'days')); % timestamps of the boundary data
        meta_data.raw_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};  % original variable name in the files
        meta_data.new_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; % standard variable name for output
        meta_data.dim_vars = {'lon', 'lat', 'depth'};    % dimensional variables in the file
        meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
        meta_data.par_flag = 0;    % 1: parallel mode; 0: serial mode.
        DS = get_hycom_bdry(Mobj, obc_bnds, meta_data);

    case 'hycom_chesbay'
        meta_data.src_file = 'E:\HYCOM\HYCOM_ChesBay\Wn79En72S33N41_****.mat';
        meta_data.src_time = unique(dateshift(Mobj.time, 'start', 'days')); % timestamps of the boundary data
        meta_data.raw_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};  % original variable name in the files
        meta_data.new_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; % standard variable name for output
        meta_data.dim_vars = {'lon', 'lat', 'depth'};    % dimensional variables in the file
        meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
        meta_data.par_flag = 0;    % 1: parallel mode; 0: serial mode.
        DS = get_hycom_bdry(Mobj, obc_bnds, meta_data);

    case 'hycom_bys_nc'
        meta_data.src_file = 'E:\HYCOM\HYCOM_BYS\W117E129S31N42_****.nc';
        meta_data.src_time = unique(dateshift(Mobj.time, 'start', 'days')); % timestamps of the boundary data
        meta_data.raw_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};  % original variable name in the files
        meta_data.new_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; % standard variable name for output
        meta_data.dim_vars = {'lon', 'lat', 'depth'};    % dimensional variables in the file
        meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
        meta_data.par_flag = 1;    % 1: parallel mode; 0: serial mode.
        DS = get_hycom_bdry_nc(Mobj, obc_bnds, meta_data);

    case 'glorys_chesbay'
        meta_data.src_file = 'E:\CMEMS\GLORYS_Reanalysis\ChesBay\L0_data\CMEMS_GLORYS_ChesBay_****.nc';
        meta_data.src_time = unique(dateshift(Mobj.time, 'start', 'days')); % timestamps of the boundary data
        meta_data.raw_list = {'zos','thetao', 'so', 'uo', 'vo'};      % original variable name in the file
        meta_data.new_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};      % standard variable name for output
        meta_data.dim_vars = {'longitude', 'latitude', 'depth'};    % dimensional variables in the file
        meta_data.date_fmt = 'yyyymmdd';  % the date format of filename
        meta_data.par_flag = 1;    % 1: parallel mode; 0: serial mode.
        DS = get_glorys_bdry(Mobj, obc_bnds, meta_data);

    case 'hycom_bys_clim'
        C = load('example_hycom_clim_1995_2020.mat');
        varList = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; nVars = numel(varList);
        if isnumeric(obc_bnds); obc_bnds = repmat({obc_bnds(:)}, [nVars, 1]); end

        DS(nVars, 1) = struct('Variable', [], 'Data', [],'Depth', [], 'Nodes', [], 'Time', []);
        for iVar = 1:numel(varList)
            obc_nodes = Mobj.obc_nodes(:, obc_bnds{iVar}); obc_nodes(obc_nodes == 0) = [];
            lon_obc = Mobj.lon(obc_nodes); lat_obc = Mobj.lat(obc_nodes);

            varName = varList{iVar}; varData = C.(varName);
            vz = max(1, ~strcmp(varName, 'ssh')*numel(C.depth));
            varTmp = nan(vz, numel(obc_nodes), 12);
            for iMon = 1:12
                varTmp(:, :, iMon) = interp_tri(lon_obc, lat_obc, C.lon, C.lat, squeeze(varData(:,:,:,iMon)), [1 1]);  % nps*depth*time
            end

            DS(iVar).Variable = varName;
            DS(iVar).Depth = abs(C.depth(1:vz));
            DS(iVar).Nodes = obc_nodes(:);
            DS(iVar).Time = Mobj.time(:);
            DS(iVar).Data = varTmp(:, :, month(Mobj.time)); % depth*nodes*time
        end
end

end