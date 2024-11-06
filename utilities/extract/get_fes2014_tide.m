function TideForc = get_fes2014_tide(Mobj, tideList)
% Extract the FES2014 tidal data
%
%% Syntax
% TideForc = get_fes2014_tide(Mobj, tideList)
%
%% Description
% TideForc = get_fes2014_tide(Mobj, tideList) returns the datastruct
% containing tidal forcing from FES2014.
%
%% Input Arguments
% Mobj --- the mesh object
% tideList --- tide constituent list.
%
%% Output Arguments
% TideForc --- the tide forcing datastruct.
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
% Last Updated on 23 Sep 2024. 
% Email: wwu@vims.edu
% 
% See also: add_nodal_factors and write_schism_bctides

%% Parse inputs
if nargin == 1
    tideList = {'S2','M2','N2','K2', 'K1','P1','O1','Q1'};
end
nTides = length(tideList);

% Make sure you have downloaded the fes2014 data set from AVISO website
% https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html
elev_path = 'E:\Tidal-Models\FES2014\aviso\L1_data\fes2014b_ocean_tide\';
ubar_path = 'E:\Tidal-Models\FES2014\aviso\L1_data\fes2014a_eastward_velocity\';
vbar_path = 'E:\Tidal-Models\FES2014\aviso\L1_data\fes2014a_northward_velocity\';

%% Extract
lonTide = ncread([elev_path, 'm2.nc'], 'lon');   % [0, 360]
latTide = ncread([elev_path, 'm2.nc'], 'lat');    % [-90, 90]

% Adjust the coordinate system if necessary
[lon_adjust, lon_flag] = check_lons(Mobj.lon, lonTide);
if lon_flag~=0
    disp('longitude coordinate system is inconsistent with FES2014 (0-360)')
end

indLon = minfind(lonTide, lon_adjust(Mobj.obc_nodes_tot));
indLat = minfind(latTide, Mobj.lat(Mobj.obc_nodes_tot));

nNodes_obc = numel(Mobj.obc_nodes_tot);

amp = nan(nNodes_obc, nTides);
pha = nan(nNodes_obc, nTides);
uamp = nan(nNodes_obc, nTides);
vamp = nan(nNodes_obc, nTides);
upha = nan(nNodes_obc, nTides);
vpha = nan(nNodes_obc, nTides);
tic
for iTide = 1:nTides
    tideName = tideList{iTide};
    amp(:,iTide) = arrayfun(@(x, y) ncread([elev_path, lower(tideName), '.nc'], 'amplitude', [x y], [1 1]), indLon, indLat);
    pha(:,iTide) = arrayfun(@(x, y) ncread([elev_path, lower(tideName), '.nc'], 'phase', [x y], [1 1]), indLon, indLat);
    uamp(:,iTide) = arrayfun(@(x, y) ncread([ubar_path, lower(tideName), '.nc'], 'Ua', [x y], [1 1]), indLon, indLat);
    upha(:,iTide) = arrayfun(@(x, y) ncread([ubar_path, lower(tideName), '.nc'], 'Ug', [x y], [1 1]), indLon, indLat);
    vamp(:,iTide) = arrayfun(@(x, y) ncread([vbar_path, lower(tideName), '.nc'], 'Va', [x y], [1 1]), indLon, indLat);
    vpha(:,iTide) = arrayfun(@(x, y) ncread([vbar_path, lower(tideName), '.nc'], 'Vg', [x y], [1 1]), indLon, indLat);
end
cst = toc;

if numel(find(isnan(amp(:)))) ~= 0 || numel(find(isnan(pha(:)))) ~= 0
    warning('there are NaN values in the extracted tide data!')
end

TideForc.tide_list = tideList;
TideForc.elev_amp = amp/100;
TideForc.elev_pha = pha;
TideForc.u_amp = uamp/100;
TideForc.u_pha = upha;
TideForc.v_amp = vamp/100;
TideForc.v_pha = vpha;

%% Put phases into [0,360)
TideForc.elev_pha(TideForc.elev_pha<0) = TideForc.elev_pha(TideForc.elev_pha<0)+360;
TideForc.u_pha(TideForc.u_pha<0) = TideForc.u_pha(TideForc.u_pha<0)+360;
TideForc.v_pha(TideForc.v_pha<0) = TideForc.v_pha(TideForc.v_pha<0)+360;

disp(['It takes ', num2str(cst, '%.2f'), ' secs to extract tidal forcing from the FES2014 data set.'])
end









