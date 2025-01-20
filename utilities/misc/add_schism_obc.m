function Mobj = add_schism_obc(Mobj, obc_nodes)
% Add new open boundary segments for SCHISM
%
%% Syntax
% Mobj = add_schism_obc(Mobj, obc_nodes)
%
%% Description
% Mobj = add_schism_obc(Mobj, obc_nodes) adds new open boundary nodes
%
%% Examples 
% mesh_file = 'BYS_20814.2dm';
% Mobj = mesh2schism(mesh_file);
% 
% obc_nodes = [4576 4845 0; 4431 4155 4154]';  % added open boundary nodes (Luanhe River and Yellow River)
% Mobj = add_schism_obc(Mobj, obc_nodes);
% 
% figure('Position', [666, 42, 672, 867], 'Color', 'w')
% disp_schism_hgrid(Mobj, [0 1 1])
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
% obc_nodes - added open boundary nodes (n_nodes*n_counts); double
%       added open boundary nodes; note that each column represents a new
%       open boundary segment, and the nodes in each new segment must be
%       adjacent on the map.
%
%% Output Arguments
% Mobj - mesh object; datastruct
%       the updated mesh object.
%
%% Notes
% This function can be used to add river boundaries to your grid.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institue of Marine Science in 2025. 
% Last Updated on 19 Jan 2025. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
[n_nodes, n_counts] = size(obc_nodes);
n_max = size(Mobj.obc_nodes ,1);

obc_nodes2 = zeros(n_max, n_counts);
obc_nodes2(1:n_nodes, :) = obc_nodes;

obc_nodes_plus = [Mobj.obc_nodes obc_nodes2];

[land_nodes, island_nodes] = find_land_island(Mobj, obc_nodes_plus);
Mobj = add_bnd_metrics(Mobj, obc_nodes_plus, land_nodes, island_nodes);

end




