function [land_nodes, island_nodes] = find_land_island(tri, edg, obc_nodes)
% Find land and island nodes (work for mixed triangular/quadrangular grid).
%
%% Syntax
% [land_nodes, island_nodes] = find_land_island(tri, edg, obc_nodes)
%
%% Description
% [land_nodes, island_nodes] = find_land_island(tri, edg, obc_nodes) finds the
%       land and island nodes for a set of unstructured grid.
%
%% Input Arguments
% tri - connectivity table of elements; numeric
%       the connectivity table between elements and nodes.
% edg - connectivity table of edges; numeric
%       the connectivity table between edges and nodes.
% obc_nodes - open boundary nodes; double 
%       obc_nodes is a matrix (M*N) providing the open or ocean boundary
%       nodes, with M indicating the node number of the longest open
%       boundary, N indicating the # of open boundaries. N = 1 means there
%       is only one open boundary for your mesh grid. 
%
%% Output Arguments
% land_nodes - land boundary nodes; numeric
%       land_nodes is a matrix (M*N) providing the land boundary nodes,
%       with M indicating the node number of the longest land boundary, N
%       indicating the # of land boundaries.
% island_nodes - island boundary nodes; numeric
%       island_nodes is a matrix (M*N) providing the island boundary nodes,
%       with M indicating the node number of the longest island boundary, N
%       indicating the # of islands.
%
%% Notes
% If you have already obtained the open boundary nodes for a set of
% unstructured gird, this function can help you get the land/island nodes.  
% 
% This function was created with the help of ChatGPT.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024.
% Last Updated on 26 Mar 2025.
% Email: wwu@vims.edu
%
% See also: add_bnd_metrics and add_grid_metrics

%% Parse inputs
% remove zero lines in the obc_nodes
obc_counts = size(obc_nodes, 2);
tmp_nodes = obc_nodes(:); tmp_nodes(tmp_nodes==0) = [];
obc_nodes_tot = tmp_nodes;

% the # of land and obc segments should be the same by default, this means
% that they can form a big loop with island loops located inside.
land_counts = obc_counts; 

% find all loop nodes (island or land-sea loop)
loop_nodes = find_loop_nodes(tri, edg);

%% Identify land and island nodes from all loop nodes
loop_lens = sum(loop_nodes~=0,1);
[~, ind_sorted] = sort(loop_lens, 'descend');
loop_nodes = loop_nodes(:, ind_sorted); % sort by the # of loop lens

outer_nodes = loop_nodes(:,1);  % the longest loop should be the land-sea boundary
island_nodes = loop_nodes(:,2:end);

max_len = max(sum(island_nodes~=0,1));
island_nodes = island_nodes(1:max_len,:);

ind_obc = ismember(outer_nodes, obc_nodes_tot);
ind_cut = find(diff(ind_obc)==1, 1);
if isempty(ind_cut)
    outer_nodes2 = outer_nodes;  % it already starts from obc nodes
else
    outer_nodes2 = [outer_nodes(ind_cut+1:end); outer_nodes(1:ind_cut)];  % ensure it starts from obc nodes
end

ind_land = ~ismember(outer_nodes2, obc_nodes_tot);
start_locs = find(diff(ind_land)==1);
end_locs = [find(diff(ind_land)==-1); numel(outer_nodes2)]+1;

land_nodes = zeros(numel(outer_nodes), land_counts); % land-sea adjacent nodes must be included in both land and open boundaries
for ii = 1:land_counts
    if end_locs(ii)>numel(outer_nodes2)
        land_tmp = [outer_nodes2(start_locs(ii):end); outer_nodes2(1)];
    else
        land_tmp = outer_nodes2(start_locs(ii):end_locs(ii));
    end
    land_nodes(1:numel(land_tmp), ii) = land_tmp(:);
end

max_len = max(sum(land_nodes~=0,1));
land_nodes = land_nodes(1:max_len,:);
%% All loops are arranged in a descending order
bnd_lens = sum(land_nodes~=0,1);
[~, ind_sorted] = sort(bnd_lens, 'descend');
land_nodes = double(land_nodes(:, ind_sorted)); % sort by the # of loop lens

bnd_lens = sum(island_nodes~=0,1);
[~, ind_sorted] = sort(bnd_lens, 'descend');
island_nodes = double(island_nodes(:, ind_sorted)); % sort by the # of loop lens

disp('land/island nodes are identified automatically')
end

