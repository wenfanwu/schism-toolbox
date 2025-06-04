function SAL = get_fes2014_SAL(Mobj, tideList)
% Extract self-attracting and loading tide (SAL) from FES2014
%
%% Syntax
% SAL = get_fes2014_SAL(Mobj, tideList)
%
%% Description
% SAL = get_fes2014_SAL(Mobj, tideList) extracts SAL data from FES2014.
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct containing mesh info.
% tideList - tidal constituent list; cell
%       the tidal constituents to be extracted. Default: tideList = {'S2','M2','N2','K2','K1','P1','O1','Q1'}; 
%
%% Output Arguments
% SAL -  SAL data; datastruct
%       the datastruct containing SAL data (sal_amp, sal_pha).
%
%% Notes
% The FES2014 data was obtained from the AVISO ftp site
% (ftp://ftp-access.aviso.altimetry.fr/), it is user-authorized. 
%
% Detailed information about FES2014 can be found at 
% https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 15 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_SAL

%% Parse inputs
if nargin == 1; tideList = {'S2','M2','N2','K2','K1','P1','O1','Q1'}; end

%% Paths to FES2014 NetCDF files
% Make sure you have downloaded the fes2014 data set from AVISO website!!!
% https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html
load_path = 'E:\Tidal-Models\FES2014\aviso\L1_data\fes2014a_load_tide\';

% Load lon/lat grid
lonTide = ncread(fullfile(load_path, 'm2.nc'), 'lon');   % [0, 360]
latTide = ncread(fullfile(load_path, 'm2.nc'), 'lat');   % [-90, 90]

% Adjust lon if needed
[lon_adjust, lon_flag] = check_lons(Mobj.lon, lonTide);
if lon_flag ~= 0
    disp('the longitude coordinate system is inconsistent with FES2014 (0-360)')
end

% Get boundary node coordinates
lon_nodes = lon_adjust;
lat_nodes = Mobj.lat;

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

tic
for iTide = 1:nTides
    tideName = lower(tideList{iTide});

    % Read subregion only
    amp_data = ncread([load_path, tideName, '.nc'], 'amplitude', [i1, j1], [i2-i1+1, j2-j1+1]);
    pha_data = ncread([load_path, tideName, '.nc'], 'phase',     [i1, j1], [i2-i1+1, j2-j1+1]);

    % Extract using subregion linear index
    amp(:,  iTide) = amp_data(lin_idx);
    pha(:,  iTide) = pha_data(lin_idx);
end
cst = toc;

% Remove junk values
amp(abs(amp)>1.e10) = 0;
pha(abs(pha)>1.e10) = 0;

% Output structure
SAL.tide_list = tideList;
SAL.sal_amp = amp/100;     % convert cm to meters
SAL.sal_pha = mod(pha, 360);  % [0, 360]

disp(['It takes ', num2str(cst, '%.2f'), ' secs to extract SAL data from the FES2014 dataset.'])
%% Kill potential NaNs
% Warn if any NaNs present
if any(isnan(amp(:)))
    warning('NaN values have been killed!')

    SAL.sal_amp(isnan(SAL.sal_amp)) = 0;
    SAL.sal_pha(isnan(SAL.sal_pha)) = 0;
end

% Double-check NaN values
if any(isnan(SAL.sal_amp(:))) || any(isnan(SAL.sal_pha(:)))
    warning on
    warning('NaN values still exist in the SAL data!!!')
end
end