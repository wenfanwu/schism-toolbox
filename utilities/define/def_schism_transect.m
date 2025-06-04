function sect_info = def_schism_transect(Mobj, mtype, dx)
% Define a transect on the mesh grid
%
%% Syntax
% sect_info = def_schism_transect(Mobj)
% sect_info = def_schism_transect(Mobj, mtype)
% sect_info = def_schism_transect(Mobj, mtype, dx)
%
%% Description
% sect_info = def_schism_transect(Mobj) defines a transect on the map.
% sect_info = def_schism_transect(Mobj, mtype) specifies the defining
%       method of transect. 
% sect_info = def_schism_transect(Mobj, mtype, dx) specifies the density of
%       transect nodes.
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the mesh object with vertical grid info added.
% mtype - method type; numeric
%       the flag used to specify the defning method of transect. 
%       Four different options are available. Default: mtype = 1;   
%       -1: draw a line on the map and return the vertical layers strictly
%            along the transect (transect points are not on the node center).
%        1: draw a line on the map and return the vertical lalyers on the
%            nearset nodes to the transect points.
%        2: click some node centers on the map using datatip and return the
%            vertical layers on these nodes.
%        3: draw broken lines on the map using polyline and return the
%            vertical layers strickly along the transect (transect points
%            are not on the node center). 
% dx - the density of transect points; numeric
%       if mtype =-1, 1, or 3: dx can be used to determine the density of
%       transect points, ranging from 0 to 1, and lower values mean high
%       density. Default: dx = 0.02; 
%       if mtype = 2: dx =1 means the transect points will be arranged in
%       a descending order based on the latitudes, while -1 means an
%       ascending order. Default: dx = 1;
% 
%% Output Arguments
% sect_Info - transect info; datastruct
%       the datastruct containing transect data.
%
%% Examples
% figure('Color', 'w')
% disp_schism_hgrid(Mobj, [1 0], 'EdgeAlpha', 0.05, 'LineWidth', 0.5);
% sect_info = def_schism_transect(Mobj, -1, 0.01);
% 
%% Notes
% All the 'define' functions (def_schism_*.m) follows the same style. A
% basemap should be activated first before you use these functions. The
% so-called basemap is often a bathymetry map generated from
% 'disp_schism_hgrid.m', or it can be a velocity map, and any other maps.
% It depends on your needs. 
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 1 Jun 2025.
% Email: wwu@vims.edu
%
% See also: disp_schism_vgrid

%% Parse inputs
if nargin < 2; mtype = 1; end
if nargin<3
    if mtype == 2; dx = 1; else; dx = 0.02; end
end
%% Define transect
switch mtype
    case -1 % drawline (will not find nearest nodes)
        disp('Please draw a line on the map and press ENTER to end')
        sect_handle = drawline; pause;
        lon_roi = sect_handle.Position(:,1)'; lat_roi = sect_handle.Position(:,2)'; close gcf
        lon_list = interp1(1:2, lon_roi, 1:dx:2); lat_list = interp1(1:2, lat_roi, 1:dx:2);
        F = scatteredInterpolant(Mobj.lon, Mobj.lat, Mobj.depth);
        dep_list = F(lon_list, lat_list);
        ind_nodes = minfind(Mobj.depth, dep_list);  % use depth layers from the nearest depth.

    case 1 % drawline (find the nearest nodes on the mesh)
        disp('Please draw a line on the map and press ENTER to end')
        sect_handle = drawline; 
        pause;
        lon_roi = sect_handle.Position(:,1)'; lat_roi = sect_handle.Position(:,2)';
        close gcf
        lon_roi_ref = interp1(1:2, lon_roi, 1:dx:2); lat_roi_ref = interp1(1:2, lat_roi, 1:dx:2);  % refine the transect points
        ind_nodes = unique(geomin(Mobj.lon, Mobj.lat, lon_roi_ref, lat_roi_ref), 'stable');
        lon_list = Mobj.lon(ind_nodes); lat_list = Mobj.lat(ind_nodes);  % find the nearest nodes

    case 2  % datatip
        disp('Please click all transect nodes on the map using datatips and press ENTER to end')
        pause
        dataTips = findobj(gcf,'Type','datatip');
        if ~isempty(dataTips)
            xy_data = cell2array(arrayfun(@(x) [dataTips(x).X; dataTips(x).Y], 1:length(dataTips), 'UniformOutput', false));
            close gcf
            lon_list = xy_data(:,1); lat_list =  xy_data(:,2);
        end
        if dx == 1  % latitude-descending order 
            [lat_list, ind] = sort(lat_list, 'descend');
            lon_list = lon_list(ind);
        elseif dx==-1 % latitude-ascending order 
            [lat_list, ind] = sort(lat_list, 'ascend');
            lon_list = lon_list(ind);
        end
        ind_nodes = arrayfun(@(x) geomin(Mobj.lon, Mobj.lat, lon_list(x), lat_list(x)), 1:nps);

    case 3  % polyline (may not work)
        geo_handle = drawpolyline;
        lon_roi = geo_handle.Position(:,1)'; lat_roi = geo_handle.Position(:,2)';
        close gcf; clear geo_handle
        lon_list = []; lat_list = [];
        for s = 1:numel(lon_roi)-1
            lon_seg = interp1(1:2, lon_roi(s:s+1), 1:dx:2);
            lat_seg = interp1(1:2, lat_roi(s:s+1), 1:dx:2);
            lon_list = [lon_list; lon_seg(:)]; lat_list = [lat_list; lat_seg(:)]; %#ok<AGROW>
        end
        [~, ia, ~] = unique(hypot(lon_list, lat_list), 'stable');
        lon_list = lon_list(ia); lat_list = lat_list(ia);

        F = scatteredInterpolant(Mobj.lon, Mobj.lat, Mobj.depth);
        dep_list = F(lon_list, lat_list);
        ind_nodes = minfind(Mobj.depth, dep_list);
end
%% Distance
switch lower(Mobj.coord)
    case 'geographic'
        dist_roi = distance(lat_list(1), lon_list(1), lat_list, lon_list, [6378137 0.0818191910428158]); % meters
    case 'cartesian'
        dist_roi = hypot(lon_list-lon_list(1), lat_list-lat_list(1));
end
%% Save transect info
sect_info.lon = lon_list(:);
sect_info.lat = lat_list(:);
sect_info.dist = dist_roi(:);

% Find the depth layers.
if isfield(Mobj, 'depLayers')
    sect_info.depth = Mobj.depLayers(:, ind_nodes(:));
end

% Returns the ANTI-clockwise angle from A to B
if mtype == -1 || mtype == 1  % only for straight line.
    A = [diff(lon_roi) diff(lat_roi)]; 
    B = [1, 0];
    sect_info.theta = vector_angle(A, B);
    sect_info.vector = A;
end
end
