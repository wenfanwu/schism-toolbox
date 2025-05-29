function Mobj = add_grid_metrics(Mobj, lon, lat, tri, depth)
% Add grid geometry metrics (work for mixed triangular/quadrangular grid).
%
%% Syntax
% Mobj = add_grid_metrics(Mobj, lon, lat, tri, depth)
%
%% Description
% Mobj = add_grid_metrics(Mobj, lon, lat, tri) add grid geomery metrics
% Mobj = add_grid_metrics(Mobj, lon, lat, tri, depth) provides depth info.
%
%% Examples 
% Mobj = add_grid_metrics(Mobj, lon, lat, tri) 
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
% lon - longitude vector@nodes; double
%       longitude vector.
% lat - latitude vector@nodes; double
%       latitude vector.
% tri - connectivity matrix; double
%       connectivity matrix (N*4).
% depth - depth vector@nodes (optional); double;
%       depth vector.
%
%% Output Arguments
% Mobj - mesh object; datastruct
%       the updated datastruct containing mesh info.
%
%% Notes
% This function was generated with the help of ChatGPT
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 24 Oct 2024. 
% Email: wwu@vims.edu
% 
% See also: mesh2schism and read_schism_hgrid

%% Parse inputs
% Ensure all cells/elements are aligned anti-clockwise
i34 = ~isnan(tri(:,4))+3;
tri3 = tri(i34==3,1:3);
tri4 = tri(i34==4,1:4);

% calculate the signed areas of triangles
x3 = lon(tri3);  y3 = lat(tri3); 
signed_areas = (x3(:,2) - x3(:,1)).*(y3(:,3) - y3(:,1)) - (y3(:,2) - y3(:,1)).*(x3(:,3) - x3(:,1));
ind3_cw = signed_areas<0;
tri3(ind3_cw,:) = flip(tri3(ind3_cw,:),2);

% split quad into two triangles to calculate the directed areas
x4_p1 = lon(tri4(:, [1 2 3])); y4_p1 = lat(tri4(:, [1 2 3]));
x4_p2 =lon(tri4(:, [3 4 1])); y4_p2 = lat(tri4(:, [3 4 1]));
direct_areas1 = (x4_p1(:,2) - x4_p1(:,1)).*(y4_p1(:,3) - y4_p1(:,1)) - (y4_p1(:,2) - y4_p1(:,1)).*(x4_p1(:,3) - x4_p1(:,1));
direct_areas2 = (x4_p2(:,2) - x4_p2(:,1)).*(y4_p2(:,3) - y4_p2(:,1)) - (y4_p2(:,2) - y4_p2(:,1)).*(x4_p2(:,3) - x4_p2(:,1));
ind4_cw = (direct_areas1+direct_areas2)<0;
tri4(ind4_cw,:) = flip(tri4(ind4_cw,:),2);

tri(i34==3,1:3) = tri3;
tri(i34==4,1:4) = tri4;

if sum(ind3_cw)==0
    disp('• all triangular cells are aligned anti-clockwise (√)')
else
    disp(['• ', num2str(sum(ind3_cw)), ' triangular cells was adjusted from clockwise to anti-clockwise'])
end
if sum(ind4_cw)==0
    disp('• all quadrangular cells are aligned anti-clockwise (√)')
else
    disp(['•', sum(ind4_cw), ' quadrangular cells was adjusted from clockwise to anti-clockwise'])
end

%% Node metrics
Mobj.region = round([min(lon)-0.5 max(lon)+0.5 min(lat)-0.5 max(lat)+0.5], 2);
Mobj.nNodes = length(lon);
Mobj.lon = lon;
Mobj.lat = lat;

if nargin==5
    Mobj.depth = depth;
else
    warning on
    warning('• depth info is not provided')
end
%% Element/Cell metrics
Mobj.nElems = size(tri,1);
Mobj.lonc = nan(Mobj.nElems, 1);
Mobj.lonc(i34==3) = mean(Mobj.lon(tri(i34==3, 1:3)), 2);
Mobj.lonc(i34==4) = mean(Mobj.lon(tri(i34==4, 1:4)), 2);

Mobj.latc = nan(Mobj.nElems, 1);
Mobj.latc(i34==3) = mean(Mobj.lat(tri(i34==3, 1:3)), 2);
Mobj.latc(i34==4) = mean(Mobj.lat(tri(i34==4, 1:4)), 2);

if isfield(Mobj, 'depth')
    Mobj.depthc = nan(Mobj.nElems, 1);
    Mobj.depthc(i34==3) = mean(Mobj.depth(tri(i34==3, 1:3)), 2);
    Mobj.depthc(i34==4) = mean(Mobj.depth(tri(i34==4, 1:4)), 2);
end

Mobj.tri = tri;
Mobj.i34 = i34;
%% Edge metrics 
% index for the 3rd edge of triangular cells
ind_e3 = ones(size(tri'));  
ind_e3(3,i34==3) = nan; 
ind_e3 = ind_e3(:);

tri(i34==3,4) = tri(i34==3,1);

edges = nan(4*Mobj.nElems, 2);
edges(1:4:end, :) = tri(:, [2, 3]); % start from the 2nd edge 
edges(2:4:end, :) = tri(:, [3, 4]);
edges(3:4:end, :) = tri(:, [4, 1]);  % invalid edge for triangular cells
edges(4:4:end, :) = tri(:, [1, 2]);

edges(isnan(ind_e3),:) = [];

% find all the independent edges
sorted_edges = sort(edges, 2);
edges = unique(sorted_edges, 'rows','stable');

Mobj.nEdges = size(edges, 1);

% calculate the coordinates@side centers
Mobj.lons = (lon(edges(:, 1)) + lon(edges(:, 2))) / 2;
Mobj.lats = (lat(edges(:, 1)) + lat(edges(:, 2))) / 2;

if isfield(Mobj, 'depth')
    Mobj.depths = mean(Mobj.depth(edges), 2);
end

Mobj.edg = edges;

disp('grid geometry metrics have been added')
end



