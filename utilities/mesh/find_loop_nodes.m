function loop_nodes = find_loop_nodes(tri, edg)
% Find the loop nodes from an unstructured grid
%
%% Syntax
% loop_nodes = find_loop_nodes(tri, edg)
%
%% Description
% loop_nodes = find_loop_nodes(tri, edg)
%
%% Input Arguments
% tri - connectivity table of elements; numeric
%       the connectivity table between elements and nodes.
% edg - connectivity table of edges; numeric
%       the connectivity table between edges and nodes.
%
%% Output Arguments
% loop_nodes - loop nodes; numeric
%       loop_nodes includes all island loops and the land-sea loop. The
%       land-sea loop should be the biggest loop with all island loops
%       located inside it. This is the most common case for 2-D unstructured
%       grids.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 26 Mar 2025. 
% Email: wwu@vims.edu
% 
% See also: find_land_island

%% Parse inputs
i34 = ~isnan(tri(:,4))+3;
tri3 = tri(i34==3,1:3);
tri4 = tri(i34==4,1:4);

%% Find edges that belong to only one element (boundary edges)
edges_tri = [tri3(:,[1,2]); tri3(:,[2,3]); tri3(:,[3,1])];   % all edges used by triangular elements (including shared edges)
edges_quad = [tri4(:,[1,2]); tri4(:,[2,3]); tri4(:,[3,4]); tri4(:,[4,1])];  % all edges used by quad elements (including shared edges)
edges = [edges_tri; edges_quad];   % all edges used by elements (including shared edges)

e0 = edg(:,1) + edg(:,2)*1i;  % complex edge list (unique)
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
ind_bdry_edges = edge_counts == 1; % the edges used by only one times are island or land-sea boundary nodes

bdry_edges = edg(ind_bdry_edges,:);
bdry_nodes = unique(bdry_edges);

[~, bdry_edges(:,1)] = ismember(bdry_edges(:,1), bdry_nodes);
[~, bdry_edges(:,2)] = ismember(bdry_edges(:,2), bdry_nodes);

G = graph(bdry_edges(:,1), bdry_edges(:,2));
[ind_loops, loop_lens] = conncomp(G);  % identify all loop components (islands or land-sea boundaries)

%% Ensure that the edges within the same loop are connected end-to-end
num_loops = max(ind_loops);
max_loop_len =  max(loop_lens);

loop_nodes = zeros(max_loop_len, num_loops);
for ii = 1:num_loops
    nodes_in_loop = bdry_nodes(ind_loops == ii);
    ind_loop_edges = ismember(edg(:,1), nodes_in_loop) & ismember(edg(:,2), nodes_in_loop) & edge_counts==1;
    edges_in_loop = edg(ind_loop_edges,:);

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
        edges_in_loop_sorted(:, jj+1) = edges_in_loop(ind_next,:);
    end

    loop_nodes(1:nps, ii) = unique(edges_in_loop_sorted(:), 'stable');
end
%% All loops are arranged in a descending order
[~, ind_sorted] = sort(loop_lens, 'descend');
loop_nodes = double(loop_nodes(:, ind_sorted)); % sort by the # of loop lens

end