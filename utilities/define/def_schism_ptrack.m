function [lon_list, lat_list] = def_schism_ptrack(Mobj, N, method_type)
% Define the locations for particle tracking.
% 
%% Syntax
% [lon_list, lat_list] = def_schism_ptrack(Mobj)
% [lon_list, lat_list] = def_schism_ptrack(Mobj, N)
% 
%% Description 
% [lon_list, lat_list] = def_schism_ptrack(Mobj) uses polygon to define
% a released region for the particles on the map.
% [lon_list, lat_list] = def_schism_ptrack(Mobj, N) allows you to define
% multiple regions.
% [lon_list, lat_list] = def_schism_ptrack(Mobj, N, method_type) selects a
% defining method.
%
%% Example
% [lon_list, lat_list] = def_schism_ptrack(Mobj)
% [lon_list, lat_list] = def_schism_ptrack(Mobj, 3, 'polygon)
% 
%% Input Arguments
% Mobj --- the mesh object.
% N --- the # of released regions. Default: N = 1;
% method_type --- the method to select the released regions.
%
%% Output Arguments
% lon_list --- longitude vector of the particles
% lat_list --- latitude vector of the particles
% 
%% Notes
% All the 'define' functions (def_schism_*.m) follows the same idea, that is,
% a basemap should be activated first before you use these functions. The
% so-called basemap is often a bathymetry map generated from
% 'disp_schism_hgrid.m', or it can be a velocity map, and any other maps.
% It depends on your needs.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2023-11-27. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: write_schism_ptrack

%% Parse inputs
if nargin < 2
    N = 1;
end
if nargin < 3
    method_type = 'polygon';
end
%% Define on the map
lon_list = cell(N,1);
lat_list = cell(N,1);

hold on
for ii = 1:N
    disp('Draw a polygon on the map to select the released region, press ENTER to end...')
    switch lower(method_type)
        case 'polygon'
            geo_handle = drawpolygon;
            lonRoi = geo_handle.Position(:,1)';
            latRoi = geo_handle.Position(:,2)';
            in_flag = inpolygon(Mobj.lon, Mobj.lat, lonRoi, latRoi);
            lon_tmp = Mobj.lon(in_flag);
            lat_tmp = Mobj.lat(in_flag);
        case 'points'
            pause;
            dataTips = findobj(gcf,'Type','datatip');
            if ~isempty(dataTips)
                select_coords = cell2array(arrayfun(@(x) [dataTips(x).X; dataTips(x).Y], 1:length(dataTips), 'UniformOutput', false));
                lon_tmp = select_coords(:,1);
                lat_tmp =  select_coords(:,2);
            else
                lon_tmp = [];
                lat_tmp = [];
            end
    end
    lon_list{ii} = lon_tmp(:);
    lat_list{ii} = lat_tmp(:);
    scatter(lon_tmp, lat_tmp, 10, 'filled', 'm')
end
lon_list = cell2mat(lon_list);
lat_list = cell2mat(lat_list);

end





















