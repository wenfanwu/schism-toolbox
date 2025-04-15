function TideForc = get_fes2014_tide(Mobj, tideList, obc_bnds)
% Extract the FES2014 tidal data
%
%% Syntax
% TideForc = get_fes2014_tide(Mobj, tideList)
% TideForc = get_fes2014_tide(Mobj, tideList, obc_bnds)
%
%% Description
% TideForc = get_fes2014_tide(Mobj, tideList) extracts FES2014 tidal forcing.
% TideForc = get_fes2014_tide(Mobj, tideList, obc_bnds) specifies the open boundary segments.
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct containing mesh info.
% tideList - tidal constituent list; cell
%       the tidal constituents to be extracted. by default: tideList = {'S2','M2','N2','K2','K1','P1','O1','Q1'}; 
% obc_bnds - the open boundary index; numeric
%       the index for the open boundary segments to be extracted. 
%       obc_bnds = 1 means only the first open boundary will be used.
%       by default: obc_bnds = 'all', or 1:Mobj.obc_counts, which means all
%       open boundaries are included.
%
%% Output Arguments
% TideForc -  tidal forcing; datastruct
%       the returned datastruct containing all tidal forcing.
%
%% Notes
% The FES2014 data was obtained from the AVISO ftp site
% (ftp://ftp-access.aviso.altimetry.fr/), it is user-authorized. 
%
% Detailed information about FES2014 can be found at 
% https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 31 Mar 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_bctides

%% Parse inputs
if nargin == 1
    tideList = {'S2','M2','N2','K2','K1','P1','O1','Q1'};
end
if nargin < 3
    obc_bnds = 'all';
end
%% Specify the open boundary segments
if strcmpi(obc_bnds, 'all')
    obc_bnds = 1:Mobj.obc_counts;
end

obc_bnds = sort(obc_bnds);
obc_nodes = Mobj.obc_nodes(:, obc_bnds);
obc_nodes(obc_nodes==0) = [];
%% Paths to FES2014 NetCDF files
% Make sure you have downloaded the fes2014 data set from AVISO website!!!
% https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html
elev_path = 'E:\Tidal-Models\FES2014\aviso\L1_data\fes2014b_ocean_tide\';
ubar_path = 'E:\Tidal-Models\FES2014\aviso\L1_data\fes2014a_eastward_velocity\';
vbar_path = 'E:\Tidal-Models\FES2014\aviso\L1_data\fes2014a_northward_velocity\';

% Load lon/lat grid
lonTide = ncread(fullfile(elev_path, 'm2.nc'), 'lon');   % [0, 360]
latTide = ncread(fullfile(elev_path, 'm2.nc'), 'lat');   % [-90, 90]

% Adjust lon if needed
[lon_adjust, lon_flag] = check_lons(Mobj.lon, lonTide);
if lon_flag ~= 0
    disp('the longitude coordinate system is inconsistent with FES2014 (0-360)')
end

% Get boundary node coordinates
lon_nodes = lon_adjust(obc_nodes);
lat_nodes = Mobj.lat(obc_nodes);

%% Determine bounding box in lon/lat indices
min_lon = min(lon_nodes); max_lon = max(lon_nodes);
min_lat = min(lat_nodes); max_lat = max(lat_nodes);

% Expand the bounding box slightly
margin = 2;
i1 = max(1, find(lonTide >= min_lon, 1, 'first') - margin);
i2 = min(length(lonTide), find(lonTide <= max_lon, 1, 'last') + margin);
j1 = max(1, find(latTide >= min_lat, 1, 'first') - margin);
j2 = min(length(latTide), find(latTide <= max_lat, 1, 'last') + margin);

% Subregion grid
lon_sub = lonTide(i1:i2);
lat_sub = latTide(j1:j2);

% Compute index of each boundary node within subgrid
ind_x = minfind(lon_sub, lon_nodes);
ind_y = minfind(lat_sub, lat_nodes);

% Linear index into subregion
lin_idx = sub2ind([length(lon_sub), length(lat_sub)], ind_x, ind_y);

%% Extract data
% Preallocate outputs
nps = numel(lon_nodes); nTides = length(tideList);

amp  = nan(nps, nTides); pha  = nan(nps, nTides);
uamp = nan(nps, nTides); upha = nan(nps, nTides);
vamp = nan(nps, nTides); vpha = nan(nps, nTides);

tic
for iTide = 1:nTides
    tideName = lower(tideList{iTide});

    % Read subregion only
    elev_amp_data = ncread([elev_path, tideName, '.nc'], 'amplitude', [i1, j1], [i2-i1+1, j2-j1+1]);
    elev_pha_data = ncread([elev_path, tideName, '.nc'], 'phase',     [i1, j1], [i2-i1+1, j2-j1+1]);

    u_amp_data = ncread([ubar_path, tideName, '.nc'], 'Ua', [i1, j1], [i2-i1+1, j2-j1+1]);
    u_pha_data = ncread([ubar_path, tideName, '.nc'], 'Ug', [i1, j1], [i2-i1+1, j2-j1+1]);

    v_amp_data = ncread([vbar_path, tideName, '.nc'], 'Va', [i1, j1], [i2-i1+1, j2-j1+1]);
    v_pha_data = ncread([vbar_path, tideName, '.nc'], 'Vg', [i1, j1], [i2-i1+1, j2-j1+1]);

    % Extract using subregion linear index
    amp(:,  iTide) = elev_amp_data(lin_idx);
    pha(:,  iTide) = elev_pha_data(lin_idx);
    uamp(:, iTide) = u_amp_data(lin_idx);
    upha(:, iTide) = u_pha_data(lin_idx);
    vamp(:, iTide) = v_amp_data(lin_idx);
    vpha(:, iTide) = v_pha_data(lin_idx);
end
cst = toc;

% Output structure
TideForc.tide_list = tideList;
TideForc.elev_amp = amp/100;               % convert cm to meters
TideForc.elev_pha = mod(pha, 360);
TideForc.u_amp = uamp/100;
TideForc.u_pha = mod(upha, 360);
TideForc.v_amp = vamp/100;
TideForc.v_pha = mod(vpha, 360);

disp(['It takes ', num2str(cst, '%.2f'), ' secs to extract tidal forcing from the FES2014 dataset.'])

%% Kill potential NaNs adjacent to the coast
% Warn if any NaNs present
nan_locs = find(isnan(sum(amp,2)));
if ~isempty(nan_locs)
    warning on
    for ii = 1:numel(nan_locs)
        warning(['NaN values found at node #', num2str(nan_locs(ii)),' have been killed!'])
    end

    field_list = fieldnames(TideForc);
    for ii = 2:numel(field_list)
        tide_var = field_list{ii};
        TideForc.(tide_var) = fillmissing(TideForc.(tide_var), 'previous', 1);
    end
end
end