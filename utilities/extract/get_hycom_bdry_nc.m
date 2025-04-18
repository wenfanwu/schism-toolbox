function DS = get_hycom_bdry_nc(Mobj, obc_bnds, meta_data)
% Parallel extraction of boundary inputs from NetCDF database.
%
%% Requirements
% 1) All NetCDF files should be named and saved by timestamp.  
% 2) Each file must include lon and lat (can be different), and a consistent depth.
% 3) Variables in "raw_list" should have size nLons*nLats or nLons*nLats*nDeps.  
% 4) Parallel Computing Toolbox is required (used by parfor).
% 
% "get_hycom_online" can download HYCOM data that meets the above requirements.
%
%% This function accelerates data extraction by:
% 1) Utilizing parallel for-loops (parfor) to speed up file processing.
% 2) Partial loading data within the sub-region that merely cover open boundary nodes.
% 3) Pre-loading grid info and subregion indices before the loop.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025.
% Last Updated on 18 Apr 2025.
% Email: wwu@vims.edu
%
% See also: get_hycom_bdry

%% Data specification   
if nargin < 3
    meta_data.src_file = 'E:\HYCOM\HYCOM_BYS\W117E129S31N42_****.nc';   % the source files
    meta_data.src_time = unique(dateshift(Mobj.time, 'start', 'days')); % timestamps of the boundary data
    meta_data.raw_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};  % original variable name in the files
    meta_data.new_list = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; % standard variable name for output
    meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
    meta_data.par_flag = 1;    % 1: parallel mode; 0: serial mode.
end
%% Essential variables
src_file = meta_data.src_file; src_time = meta_data.src_time;
raw_list = meta_data.raw_list; new_list = meta_data.new_list;
date_fmt = meta_data.date_fmt; par_flag = meta_data.par_flag;

if nargin==1; obc_bnds = 'all'; end
if strcmpi(obc_bnds, 'all'); obc_bnds = 1:Mobj.obc_counts; end

obc_nodes = Mobj.obc_nodes(:, obc_bnds); obc_nodes = obc_nodes(:);
obc_nodes(obc_nodes == 0) = [];
nps = numel(obc_nodes);

lon_obc = Mobj.lon(obc_nodes);
lat_obc = Mobj.lat(obc_nodes);
nt = numel(src_time);

bbox = [min(lon_obc) max(lon_obc) min(lat_obc) max(lat_obc)];
%% Pre-load grid info and subregion index
lon_all = cell(nt,1); lat_all = cell(nt,1); depth_all = cell(nt,1);
ind_lons = cell(nt,1); ind_lats = cell(nt,1); 

for iTime = 1:nt
    filepath = strrep(src_file, '****', datestr(src_time(iTime), date_fmt)); %#ok<DATST>
    info = dir(filepath);
    if iTime==1
        sz = info.bytes;
        C.lon = ncread(filepath, 'lon');
        C.lat = ncread(filepath, 'lat');
        C.depth = ncread(filepath, 'depth');
    end
    if abs(info.bytes-sz)/sz > 0.10  % may use different coordinates
        sz = info.bytes;
        C.lon = ncread(filepath, 'lon');
        C.lat = ncread(filepath, 'lat');
        C.depth = ncread(filepath, 'depth');
    end
    [idx_x, idx_y, lon_reg, lat_reg] = sub_region(C.lon, C.lat, bbox, 2);

    ind_lons{iTime} = idx_x;
    ind_lats{iTime} = idx_y;
    lon_all{iTime} = lon_reg;
    lat_all{iTime} = lat_reg;
    depth_all{iTime} = abs(C.depth);
end
nz = numel(depth_all{1});  % depth must be the same for all files

%% Check the parallel pool
if par_flag == 1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        num_workers = parcluster('local').NumWorkers;
        disp('parallel pool is being activated ...');
        parpool('local', num_workers);
    else
        disp(['parallel pool is running with ', num2str(poolobj.NumWorkers), ' logical threads']);
    end
end
%% Loop over each variable
for iVar = 1:numel(raw_list)
    raw_name = raw_list{iVar};
    new_name = new_list{iVar};
    disp([new_name, ' is being extracted ...'])

    varAll = zeros(nps, max(1, ~strcmp(new_name, 'ssh')*nz), nt);
    tic
    if par_flag == 1
        parfor iTime = 1:nt
            filepath = strrep(src_file, '****', datestr(src_time(iTime), date_fmt)); %#ok<DATST>
            idx_x = ind_lons{iTime}; idx_y = ind_lats{iTime};
            switch new_name
                case 'ssh'
                    varTmp = ncread(filepath, raw_name, [min(idx_x) min(idx_y)], [abs(diff(idx_x))+1 abs(diff(idx_y))+1]);
                otherwise
                    varTmp = ncread(filepath, raw_name, [min(idx_x) min(idx_y) 1], [abs(diff(idx_x))+1 abs(diff(idx_y))+1 nz]);
            end
            varAll(:,:, iTime) = interp_tri(lon_obc, lat_obc, lon_all{iTime}, lat_all{iTime}, varTmp);
        end
    else
        for iTime = 1:nt
            progressbar(iTime/nt)
            filepath = strrep(src_file, '****', datestr(src_time(iTime), date_fmt)); %#ok<DATST>
            idx_x = ind_lons{iTime}; idx_y = ind_lats{iTime};
            switch new_name
                case 'ssh'
                    varTmp = ncread(filepath, raw_name, [min(idx_x) min(idx_y)], [abs(diff(idx_x))+1 abs(diff(idx_y))+1]);
                otherwise
                    varTmp = ncread(filepath, raw_name, [min(idx_x) min(idx_y) 1], [abs(diff(idx_x))+1 abs(diff(idx_y))+1 nz]);
            end
            varAll(:,:, iTime) = interp_tri(lon_obc, lat_obc, lon_all{iTime}, lat_all{iTime}, varTmp);
        end
    end
    cst = toc;

    % store in output struct
    D.ind = obc_nodes;
    D.depth = depth_all{1};  % depth is same across time steps
    D.time = src_time;
    D.var = varAll;
    DS.(new_name) = D;
    
    disp(['It takes ', num2str(cst,'%.2f'),' secs to extract ', new_name])
end

end
