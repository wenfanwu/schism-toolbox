function [min_dist, nodes] = coast_distance(Mobj, with_island)
% Calculate the minimum distance of all nodes to the coastline.
% 
%% Syntax
% [min_dist, nodes] = coast_distance(Mobj)
% [min_dist, nodes] = coast_distance(Mobj, with_island)
%
%% Description
% [min_dist, nodes] = coast_distance(Mobj) calculates the minimum distance to the coastline 
% [min_dist, nodes] = coast_distance(Mobj, with_island) includes the island nodes or not.
%
%% Examples
% [min_dist, nodes] = coast_distance(Mobj);
% figure; disp_schism_var(Mobj, min_dist/1e3)
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store the mesh info.
% with_island - island flags; char
%       the flag used to include the island nodes or not (yes/no). 
%       Default: with_island = 'yes'.
%
%% Output Arguments
% min_dist - minimum distance; numeric
%       the minimum distance (in meters) of all nodes to the coastline.
% nodes - nearest coast nodes; numeric
%       the nearest coastline nodes for each node.
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025.
% Last Updated on 30 Sep 2025.
% Email: wwu@vims.edu
%
% See also: 

%% Parse inputs
if nargin<2; with_island = 'yes'; end

%% Distance to the coastline
% Index of shore nodes (w/ or w/o islands)
if strcmpi(with_island, 'yes')
    coast_nodes = [Mobj.land_nodes(:); Mobj.island_nodes(:)];
else
    coast_nodes = Mobj.land_nodes(:);
end
coast_nodes(coast_nodes==0) = []; coast_nodes = coast_nodes(:);

% Coordinate systems
if strncmpi(Mobj.coord, 'geographic', 3)
    method = @haversine_pdist2; 
else
    method = 'euclidean';
end
[min_dist, index] = min(pdist2([Mobj.lat(:) Mobj.lon(:)], [Mobj.lat(coast_nodes) Mobj.lon(coast_nodes)], method), [], 2);
nodes = coast_nodes(index);

end

function d = haversine_pdist2(Xi, Y)
% Great-circle distance between one point and a set of points.
%
%   INPUTS:
%       Xi - 1×2 vector [lat lon] in degrees
%       Y  - Q×2 array of [lat lon] in degrees
%
%   OUTPUT:
%       d  - 1×Q row vector of distances (in meters)

R = 6378.137*1e3;  % earth radius (m)
dlat = deg2rad(Y(:,1) - Xi(1));
dlon = deg2rad(Y(:,2) - Xi(2));
a = sin(dlat/2).^2 + cos(deg2rad(Xi(1))).*cos(deg2rad(Y(:,1))).*sin(dlon/2).^2;
d = (R * 2 * atan2(sqrt(a), sqrt(1-a)))';  % m

end
