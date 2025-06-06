function Mobj = def_schism_obc(Mobj, obc_counts)
% Define open boundary nodes on the SCHISM grid.
%
%% Syntax
% Mobj = def_schism_obc(Mobj)
% Mobj = def_schism_obc(Mobj, obc_counts)
%
%% Description
% Mobj = def_schism_obc(Mobj) define the open boundary segments on the map.
% Mobj = def_schism_obc(Mobj, obc_counts) specifies the # of open boundary segments.
%
%% Examples
% Mobj2 = def_schism_obc(Mobj, 2);
% figure; disp_schism_hgrid(Mobj2, [0 1])
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
%       Required fields: lon, lat, edg, tri.
% obc_counts - the # of open boundaries; numeric
%       the # of open boundary segments. Default: obc_counts = 1.
%
%% Output Arguments
% Mobj - mesh object; datastruct
%       the updated mesh object with boundary info added.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025.
% Last Updated on 5 Jun 2025.
% Email: wwu@vims.edu
%
% See also: def_schism_mask

%% Parse inputs
if nargin < 2; obc_counts = 1; end
ux = Mobj.lon; uy = Mobj.lat;

%% Find the land-sea loop
loop_nodes = find_loop_nodes(Mobj.tri, Mobj.edg);
land_sea_nodes = loop_nodes(:,1);
head_node = land_sea_nodes(1);
tail_node = land_sea_nodes(end);

%% Check the basemap
if isempty(findall(0, 'Type', 'figure'))
    figure
    disp_schism_hgrid(Mobj, [0 0])
    axis image
    auto_center
end
%% Define open boundary segments
nps_max = numel(land_sea_nodes);
obc_nodes = zeros(nps_max, obc_counts);

hold on
for ii = 1:obc_counts
    disp('draw a polygon on the map and press ENTER')
    geo_handle = drawpolygon; 
    sx = geo_handle.Position(:,1)'; sy = geo_handle.Position(:,2)';

    idx_nodes = find(inpolygon(ux, uy, sx, sy));
    tmp_nodes = intersect(land_sea_nodes(:), idx_nodes, 'stable');

    % Exceptional case: both head and tail nodes of the land-sea loop are included.
    ind_head = find(tmp_nodes==head_node, 1);
    ind_tail = find(tmp_nodes==tail_node, 1);
    if ~isempty(ind_head) && ~isempty(ind_tail)
        [~, idx] = ismember(tmp_nodes, land_sea_nodes);
        ind_cut = find(diff(idx)~=1);
        tmp_nodes = [tmp_nodes(ind_cut+1:end); tmp_nodes(1:ind_cut)];
    end
    nps = numel(tmp_nodes);
    obc_nodes(1:nps, ii) = tmp_nodes;
    
    delete(geo_handle)
    scatter(ux(tmp_nodes), uy(tmp_nodes), 20, 's','filled','b')  % highlight the selected boundary nodes
end

%% Find land/island nodes
[land_nodes, island_nodes] = find_land_island(Mobj.tri, Mobj.edg, obc_nodes);

%% Update the mesh object
Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes);

end