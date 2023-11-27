function sect_info = def_schism_transect(Mobj, method_flag, ref_val)
% Define a transect on the mesh grid
%
%% Syntax
% sect_info = def_schism_transect(Mobj)
% sect_info = def_schism_transect(Mobj, method_flag)
% sect_info = def_schism_transect(Mobj, method_flag, ref_val)
%
%% Description
% sect_info = def_schism_transect(Mobj) defines a transect on the map.
% sect_info = def_schism_transect(Mobj, method_flag) specifies the
% defineing method.
% sect_info = def_schism_transect(Mobj, method_flag, ref_val) specifies the
% refining value on the transect, with lower values indicating denser
% points
%
%% Input Arguments
% Mobj --- mesh object (vertical grid info must be added)
% method_flag --- the flags used to decide the method to define the
% transect. Avaliable flags are given below. Default: method_flag = 1; 
%       -1: draw a line on the map and get the interpolated vertical layers
%            on these points.
%        1: draw a line on the map and return the vertical lalyers on the
%            nearset nodes to these points. 
%        2: tap some nodes on the map using datatip and return the
%            vertical layers on these nodes.
%        3: draw multiple lines on the map using polyline and return the
%            vertical layers on these points (Not Work Now)
% 
%% Output Arguments
% sect_Info --- the datastruct containing transect data.
%
%% Examples
% figure('Color', 'w')
% disp_schism_hgrid(Mobj, [1 0], 'EdgeAlpha', 0.05, 'LineWidth', 0.5);
% sect_info = def_schism_transect(Mobj, -1, 0.01);
% 
%% Notes
% All the 'define' functions (def_schism_*.m) follows the same usage. A
% basemap should be activated first before you use these functions. The
% so-called basemap is often a bathymetry map generated from
% 'disp_schism_hgrid.m', or it can be a velocity map, and any other maps.
% It depends on your needs. 
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2022-05-19.
% Email: wenfanwu@stu.ouc.edu.cn
%
% See also: drawline and findobj

%% Parse inputs
if nargin < 2
    method_flag = 1;
end
if nargin<3
    if method_flag == 2
        ref_val = 1;
    else
        ref_val = 0.02;
    end
end
switch method_flag
    case -1 % drawline (will not find nearest nodes)
        disp('Please draw a line on the map and press Enter to end to continue')
        sect_handle = drawline;
        pause;
        lonRoi = sect_handle.Position(:,1)';
        latRoi = sect_handle.Position(:,2)';
        close gcf
        lon_list = interp1(1:2, lonRoi, 1:ref_val:2);
        lat_list = interp1(1:2, latRoi, 1:ref_val:2);

        F = scatteredInterpolant(Mobj.lon, Mobj.lat, Mobj.depth);
        depList = F(lon_list, lat_list);
        indNodes = minfind(Mobj.depth, depList);
    case 1 % drawline (find the nearest nodes on the mesh)
        disp('Please draw a line on the map and press Enter to end to continue')
        sect_handle = drawline;
        pause;
        lonRoi = sect_handle.Position(:,1)';
        latRoi = sect_handle.Position(:,2)';
        close gcf
        lonRoi2 = interp1(1:2, lonRoi, 1:ref_val:2);
        latRoi2 = interp1(1:2, latRoi, 1:ref_val:2);
        indNodes = geomin(Mobj.lon, Mobj.lat, lonRoi2, latRoi2);
        indNodes = unique(indNodes, 'stable');
        lon_list = Mobj.lon(indNodes);
        lat_list = Mobj.lat(indNodes);
    case 2  % datatip
        disp('Please point out the section on the nodes and press Enter to end to continue')
        pause
        dataTips = findobj(gcf,'Type','datatip');
        if ~isempty(dataTips)
            select_coords = cell2array(arrayfun(@(x) [dataTips(x).X; dataTips(x).Y], 1:length(dataTips), 'UniformOutput', false));
            close gcf
            siteNums = size(select_coords,1);
            lon_list = select_coords(:,1);
            lat_list =  select_coords(:,2);
        end
        if ref_val == 1
            [lat_list, ind] = sort(lat_list, 'descend');
            lon_list = lon_list(ind);
        else
            [lat_list, ind] = sort(lat_list, 'ascend');
            lon_list = lon_list(ind);
        end
        indNodes = arrayfun(@(x) geomin(Mobj.lon, Mobj.lat, lon_list(x), lat_list(x)), 1:siteNums);
    case 3  % polyline (may not work)
        geo_handle = drawpolyline;
        lonRoi = geo_handle.Position(:,1)';
        latRoi = geo_handle.Position(:,2)';
        close gcf
        clear geo_handle
        F = scatteredInterpolant(Mobj.lon, Mobj.lat, Mobj.depth);

        lon_list = []; lat_list = [];
        for iSeg = 1:numel(lonRoi)-1
            lonSegs = interp1(1:2, lonRoi(iSeg:iSeg+1), 1:ref_val:2);
            latSegs = interp1(1:2, latRoi(iSeg:iSeg+1), 1:ref_val:2);
            lon_list = [lon_list; lonSegs(:)]; %#ok<*AGROW>
            lat_list = [lat_list; latSegs(:)];
        end
        [~, ia, ~] = unique(hypot(lon_list, lat_list), 'stable');
        lon_list = lon_list(ia);
        lat_list = lat_list(ia);

        depList = F(lon_list, lat_list);
        indNodes = minfind(Mobj.depth, depList);
end

switch lower(Mobj.coord)
    case 'geographic'
        distRoi = distance(lat_list(1), lon_list(1), lat_list, lon_list, [6378.137 0.0818191910428158]);
    case 'cartesian'
        distRoi = hypot(lon_list-lon_list(1), lat_list-lat_list(1));
end

sect_info.lon = lon_list(:);
sect_info.lat = lat_list(:);
sect_info.dist = distRoi(:);

if isfield(Mobj, 'depLayers')
    sect_info.depth = Mobj.depLayers(:, indNodes(:));
end

% Returns the ANTI-clockwise angle from A to B
if method_flag == -1 || method_flag == 1
    A = [diff(lonRoi) diff(latRoi)];
    B = [1, 0];
    sect_info.theta = vector_angle(A, B);
    sect_info.vector = A;
end
end
