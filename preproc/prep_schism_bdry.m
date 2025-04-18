function DS = prep_schism_bdry(Mobj, dst, obc_bnds)
% Prepare the boundary iniputs for SCHISM.
%
%% Syntax
% DS = prep_schism_bdry(Mobj, dst)
%
%% Description
% DS = prep_schism_bdry(Mobj, dst) prepares the boundary inputs for SCHISM
%
%% Example
% DS = prep_schism_bdry(Mobj, 'hycom_bys')
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% dst - dateset; char
%       the dateset name.
% obc_bnds - the open boundaries; numeric
%       the open boundary segments with boundary inputs applied, e.g.,
%       obc_bnds = 1:Mobj.obc_counts or "all" means that all open
%       boundaries are used.
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
% Last Updated on 18 Apr 2025.
% Email: wwu@vims.edu
%
% See also: prep_schism_init

%% Parse inputs
if nargin==1; obc_bnds = 'all'; end
if strcmpi(obc_bnds, 'all'); obc_bnds = 1:Mobj.obc_counts; end

%% Select the database
switch dst
    case 'hycom_bys'
        meta_data.src_file = 'C:\Users\wwu\OneDrive - vims.edu\GitHub_Projects\schism-toolbox\data\hycom\W117E127S33N41_****.mat';   % the source files
        meta_data.src_time = unique(dateshift(Mobj.time, 'start', 'days')); % timestamps of the boundary data
        meta_data.raw_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};  % original variable name in the files
        meta_data.new_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; % standard variable name for output
        meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
        meta_data.par_flag = 1;    % 1: parallel mode; 0: serial mode.
        DS = get_hycom_bdry(Mobj, obc_bnds, meta_data);

    case 'hycom_chesbay'
        meta_data.src_file = 'E:\HYCOM\HYCOM_ChesBay\Wn79En72S33N41_****.mat';
        meta_data.src_time = unique(dateshift(Mobj.time, 'start', 'days')); % timestamps of the boundary data
        meta_data.raw_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};  % original variable name in the files
        meta_data.new_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; % standard variable name for output
        meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
        meta_data.par_flag = 1;    % 1: parallel mode; 0: serial mode.
        DS = get_hycom_bdry(Mobj, obc_bnds, meta_data);

    case 'hycom_chesbay_nc'
        meta_data.src_file = 'E:\HYCOM\HYCOM_ChesBay\Wn79En72S33N41_****.nc';
        meta_data.src_time = unique(dateshift(Mobj.time, 'start', 'days')); % timestamps of the boundary data
        meta_data.raw_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};  % original variable name in the files
        meta_data.new_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; % standard variable name for output
        meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
        meta_data.par_flag = 1;    % 1: parallel mode; 0: serial mode.
        DS = get_hycom_bdry_nc(Mobj, obc_bnds, meta_data);

    case 'hycom_bys_clim'
        C = load('example_hycom_clim_1995_2020.mat');
        varList = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};

        obc_nodes = Mobj.obc_nodes(:, obc_bnds); obc_nodes(obc_nodes == 0) = [];
        lon_obc = Mobj.lon(obc_nodes); lat_obc = Mobj.lat(obc_nodes);

        D.ind = obc_nodes; D.lon = Mobj.lon(D.ind); D.lat = Mobj.lat(D.ind);
        D.depth = abs(C.depth); D.time = Mobj.time;
        for iVar = 1:numel(varList)
            varName = varList{iVar}; varData = C.(varName);
            nz = max(1, ~strcmp(varName, 'ssh')*numel(D.depth));
            varTmp = nan(numel(obc_nodes), nz, 12);
            for iMon = 1:12
                varTmp(:, :, iMon) = interp_tri(lon_obc, lat_obc, C.lon, C.lat, squeeze(varData(:,:,:,iMon)), [1 0]);  % nps*depth*time
            end
            D.var = varTmp(:, :, month(Mobj.time));
            DS.(varName) = D;
        end
end

end