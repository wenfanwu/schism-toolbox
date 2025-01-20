function [land_nodes, island_nodes] = find_land_island(Mobj, obc_nodes)
% Find land and island nodes (work for mixed triangular/quadrangular grid).
%
%% Syntax
% [land_nodes, island_nodes] = find_land_island(Mobj, obc_nodes)
%
%% Description
% [land_nodes, island_nodes] = find_land_island(Mobj, obc_nodes) finds the
% land and island nodes for a set of unstructured grid.
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info. Reuiqred fields: tri and edg.
% obc_nodes - open boundary nodes; double 
%       obc_nodes is a matrix (M*N) providing the open or ocean boundary
%       nodes, with M indicating the node number of the longest open
%       boundary, N indicating the # of open boundaries. N = 1 means there
%       is only one open boundary for your mesh grid. 
%
%% Output Arguments
% land_nodes - land boundary nodes; double 
%       land_nodes is a matrix (M*N) providing the land boundary nodes,
%       with M indicating the node number of the longest land boundary, N
%       indicating the # of land boundaries.
% island_nodes - island boundary nodes; double 
%       island_nodes is a matrix (M*N) providing the island boundary nodes,
%       with M indicating the node number of the longest island boundary, N
%       indicating the # of islands.
%
%% Notes
% If you have already obtained the open boundary nodes for a set of
% unstructured gird (lon/lat/tri), this function can help you get the land
% and island nodes. 
% 
% This function was created with the help of ChatGPT.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024.
% Last Updated on 03 Nov 2024.
% Email: wwu@vims.edu
%
% See also: add_bnd_metrics and add_grid_metrics

%% Parse inputs
i34 = ~isnan(Mobj.tri(:,4))+3;
tri3 = Mobj.tri(i34==3,1:3);
tri4 = Mobj.tri(i34==4,1:4);

obc_counts = size(obc_nodes, 2);
tmp_nodes = obc_nodes(:);
tmp_nodes(tmp_nodes==0) = [];
obc_nodes_tot = tmp_nodes;

land_counts = obc_counts;
%% Find edges belonging to only one element (viz. boundary edges)
edges_tri = [tri3(:,[1,2]); tri3(:,[2,3]); tri3(:,[3,1])];   % all edges used by triangular elements (including shared edges)
edges_quad = [tri4(:,[1,2]); tri4(:,[2,3]); tri4(:,[3,4]); tri4(:,[4,1])];  % all edges used by quad elements (including shared edges)
edges = [edges_tri; edges_quad];   % all edges used by elements (including shared edges)

e0 = Mobj.edg(:,1)+Mobj.edg(:,2)*1i;  % complex edge list (unique)
e1 = edges(:,1)+edges(:,2)*1i;   % compex used edges list (position direction)
e2 = edges(:,2)+edges(:,1)*1i; % compex used edges list (negative direction)

[e1_unique, ~, ind_obc_nodes] = unique(e1);
e1_counts = accumarray(ind_obc_nodes, 1);

[is_in, locs] = ismember(e0, e1_unique);
counts_e1 = zeros(size(e0));
counts_e1(is_in) = e1_counts(locs(is_in));

[e2_unique, ~, ind_obc_nodes] = unique(e2);
e2_counts = accumarray(ind_obc_nodes, 1);

[is_in, locs] = ismember(e0, e2_unique);
counts_e2 = zeros(size(e0));
counts_e2(is_in) = e2_counts(locs(is_in));

edge_counts = counts_e1 + counts_e2;  %  the # of times each edge is used
ind_bdry_edges = find(edge_counts == 1); % the edges used by only one times are island or land-sea boundary nodes

bdry_edges = Mobj.edg(ind_bdry_edges,:);
bdry_nodes = unique(bdry_edges);

[~, bdry_edges(:,1)] = ismember(bdry_edges(:,1), bdry_nodes);
[~, bdry_edges(:,2)] = ismember(bdry_edges(:,2), bdry_nodes);

G = graph(bdry_edges(:,1), bdry_edges(:,2));
[ind_loops, loop_lens] = conncomp(G);  % identify all loop components (islands or land-sea boundaries)

%% Make sure all loop edges meet end to end
num_loops = max(ind_loops);
max_loop_len =  max(loop_lens);

ind_bdry_nodes = zeros(max_loop_len, num_loops);
for ii = 1:num_loops
    nodes_in_loop = bdry_nodes(ind_loops == ii);
    ind_loop_edges = ismember(Mobj.edg(:,1), nodes_in_loop) & ismember(Mobj.edg(:,2), nodes_in_loop) & edge_counts==1;
    edges_in_loop = Mobj.edg(ind_loop_edges,:);

    nps = size(edges_in_loop,1);
    ind_sorted = nan(nps,1); ind_sorted(1) = 1; % start from the first node
    edges_in_loop_sorted = edges_in_loop';
    for jj = 1:nps-1
        if jj==1
            ind_next = 1;
        end
        d1 = edges_in_loop(ind_next,1);
        d2 = edges_in_loop(ind_next,2);

        ind_next = find(edges_in_loop(:,1) == d2);
        if isempty(ind_next)
            ind_next = find(edges_in_loop(:,2) == d2 & edges_in_loop(:,1) ~= d1);
            edges_in_loop(ind_next,:) = flip(edges_in_loop(ind_next,:),2);
        end
        ind_sorted(jj+1) = ind_next;
        edges_in_loop_sorted(:,jj+1) = edges_in_loop(ind_next,:);
    end

    ind_bdry_nodes(1:nps, ii) = unique(edges_in_loop_sorted(:), 'stable');
end
%% Identify land and island nodes from all loop nodes
loop_lens = sum(ind_bdry_nodes~=0,1);
[~, ind_sorted] = sort(loop_lens, 'descend');
ind_bdry_nodes = ind_bdry_nodes(:, ind_sorted); % sort by the # of loop lens

outer_nodes = ind_bdry_nodes(:,1);  % the longest loop should be the land-sea boudary
island_nodes = ind_bdry_nodes(:,2:end);

max_len = max(sum(island_nodes~=0,1));
island_nodes = island_nodes(1:max_len,:);

ind_obc = ismember(outer_nodes, obc_nodes_tot);
ind_cut = find(diff(ind_obc)==1, 1);
outer_nodes2 = [outer_nodes(ind_cut+1:end); outer_nodes(1:ind_cut)];  % ensure it starts from obc nodes

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

disp('land/island nodes have been identified automatically')
end

