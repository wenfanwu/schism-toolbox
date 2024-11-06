function river_info = match_rivers(lon_elems, lat_elems, ind_elems)
% Match the source elements with rivers
%
%% Syntax
% river_list = match_rivers(lon_elems, lat_elems)
%
%% Description
% river_list = match_rivers(lon_elems, lat_elems) matches the source
% elements with rivers.
%
%% Examples
% river_list = match_rivers(lon_elems, lat_elems)
%
%% Input Arguments
% lon_elems --- a vector containing the longitude of source elements.
% lat_elems --- a vector containing the latitude of source elements.
%
%% Output Arguments
% river_list --- corresponding rivers for all given elements. There can be
% dulplicate rivers, meaning that more than one outlets for the same river.
%
%% Notes
% User can easily add rivers in the code, and the only variable you need to
% provide is the 'estuary_zone', namely the region covering the river
% estuary, e.g. estuary_zone = [119.5 122 31.5 32.5] denotes the Changjiang
% River. 
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023.
% Last Updated on 2023-11-24.
% Email: wenfanwu@stu.ouc.edu.cn
%
% See also: def_schism_source and add_river_runoff

%% Parse inputs
river_list = cell(size(lon_elems));
river_list(:) = {'unknown'};

estuary_zone = [118.8 119.6 37.5 38.1];
ind_river  = is_river(lon_elems, lat_elems, estuary_zone);
river_list(ind_river) = {'Huanghe'};

estuary_zone = [117.6 117.88 38.9 39.1];
ind_river  = is_river(lon_elems, lat_elems, estuary_zone);
river_list(ind_river) = {'Haihe'};

estuary_zone = [119.2 119.4 39.2 39.55];
ind_river  = is_river(lon_elems, lat_elems, estuary_zone);
river_list(ind_river) = {'Luanhe'};

estuary_zone = [121.6 122 40.73 41.2];
ind_river  = is_river(lon_elems, lat_elems, estuary_zone);
river_list(ind_river) = {'Liaohe'};

estuary_zone = [119.5 122 31.5 32.5];
ind_river  = is_river(lon_elems, lat_elems, estuary_zone);
river_list(ind_river) = {'Changjiang'};

estuary_zone = [120 122 30 31];
ind_river  = is_river(lon_elems, lat_elems, estuary_zone);
river_list(ind_river) = {'Qiantangjiang'};

% remove unvalid source elements
ind_unvalid = contains(river_list, 'unknown');

if nargin>=3
    river_info.river_elems = ind_elems;
end
river_info.lon_list = lon_elems(~ind_unvalid);
river_info.lat_list = lat_elems(~ind_unvalid);
river_info.river_list = river_list(~ind_unvalid);

if numel(find(ind_unvalid))>0
    warning([num2str(numel(find(ind_unvalid))), ' unvalid source elements have been removed'])
else
    disp('all selected elements are matched!')
end

end

function ind_river  = is_river(lon_elems, lat_elems, estuary_zone)
% Determine whether the given point is inside a certain region

ind_river = lon_elems > estuary_zone(1) & lon_elems < estuary_zone(2) & lat_elems > estuary_zone(3) & lat_elems < estuary_zone(4);
end

















