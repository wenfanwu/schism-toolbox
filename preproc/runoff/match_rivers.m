function river_info = match_rivers(lon_elems, lat_elems, ind_elems)
% Match the source elements with rivers
%
%% Syntax:
%   river_info = match_rivers(lon_elems, lat_elems, ind_elems)
%
%% Description:
%   Matches source elements (specified by longitude and latitude) with
%   pre-defined estuarine zones representing major rivers.
%
%% Input Arguments:
%   lon_elems --- vector of longitudes
%   lat_elems --- vector of latitudes
%   ind_elems --- indices of the source elements
%
%% Output Arguments:
%   river_info.River  --- river names of matched elements
%   river_info.Lon    --- longitudes of matched elements
%   river_info.Lat    --- latitudes of matched elements
%   river_info.Elements --- indices of matched elements
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2023.
% Last Updated on 6 May 2025.
% Email: wwu@vims.edu
%
% See also: add_river_runoff

%% Define estuarine zones for rivers (CAN BE CHANGED)
river_definitions = {'Huanghe',       [118.8 119.6 37.5 38.1];
    'Huaihe',        [120.2 120.45 34 34.2];
    'Haihe',         [117.6 117.88 38.9 39.1];
    'Luanhe',        [119.2 119.4 39.2 39.55];
    'Liaohe',        [121.6 122 40.73 41.2];
    'Changjiang',    [119.5 122 31.5 32.5];
    'Qiantangjiang', [120 122 30 31];
    'Yalujiang',     [124.1693 124.4999 39.7008 40.0907];
};

%% Begin to match rivers
% Initialize river info
river_info(numel(lon_elems), 1) = struct('River', [], 'Lon', [], 'Lat', [], 'Elements', []);

% Match each river
for ii = 1:size(river_definitions, 1)
    river_name = river_definitions{ii,1};
    river_zone = river_definitions{ii,2};
    ind = is_river(lon_elems, lat_elems, river_zone);
    for jj = 1:numel(ind)
        river_info(ind(jj)).River = river_name;
        river_info(ind(jj)).Lon = lon_elems(ind(jj));
        river_info(ind(jj)).Lat = lat_elems(ind(jj));
        river_info(ind(jj)).Elements = ind_elems(ind(jj));       
    end
end

% Display summary
ind_empty = arrayfun(@(x) isempty(x.River), river_info);
if any(ind_empty)
    warning('%d unvalid source elements have been removed', sum(ind_unvalid));
else
    disp('All selected elements are matched!');
end
river_info = river_info(~ind_empty);  % Remove empty rows.

end

function ind_river = is_river(lon, lat, zone)
% Check if coordinates fall within the estuarine zone
ind_river = find(lon > zone(1) & lon < zone(2) & lat > zone(3) & lat < zone(4));
end
