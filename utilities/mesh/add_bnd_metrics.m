function Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes, sort_flag)
% Add open/land/island boundary info.
%
%% Syntax
% Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes)
% Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes, sort_flag)
%
%% Description
% Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes) add
%       boundary info into the mesh object
% Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes, sort_flag) 
%       sorts the nodes in descending order or not. 
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       a datastruct containing the mesh info.
% obc_nodes - open boundary nodes; double
%       a matrix (M*N) containing the open boundary info, with M indicating
%       the node number of the longest open boundary, while N indicating
%       the # of open boundaries.
% land_nodes - open boundary nodes; double
%       a matrix (M*N) containing the land boundary info, with M indicating
%       the node number of the longest land boundary, while N indicating
%       the # of land boundaries.
% island_nodes - open boundary nodes; double
%       a matrix (M*N) containing the island boundary info, with M indicating
%       the node number of the longest island boundary, while N indicating
%       the # of island boundaries.
% sort_flag - sort flags; double (0/1)
%       the flag used to determine whether sort the nodes or not. default: sort_flag = 1;
%
%% Output Arguments
% Mobj - mesh object; datastruct
%       Mobj with updated boundary info.
%
%% Notes
% The outer boundaies (land/sea) should be aligned anti-clockwise while the
% inner boundaries (island) should be clockwise.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 11 Nov 2024. 
% Email: wwu@vims.edu
% 
% See also: add_grid_metrics

%% Parse inputs
if nargin<5
    sort_flag = 1;
end
% 1) arrange all boundaries in descending order based on the # of nodes (if sort_flag=1)
% 2) trim redundant lines that are all zeros
obc_nodes = trim_bnd_nodes(obc_nodes, sort_flag);
land_nodes = trim_bnd_nodes(land_nodes, sort_flag);
island_nodes = trim_bnd_nodes(island_nodes, sort_flag);

%% Check the island loops (inner loops)
% the islands don't need to form a loop in hgrid.gr3 or hgrid.ll
ind_btm = sub2ind(size(island_nodes), sum(island_nodes~=0, 1), 1:size(island_nodes,2));
if sum(island_nodes(1,:) ~= island_nodes(ind_btm))==0
    disp('removed ending points of islands, since all islands start/end at the same node')
    island_nodes(ind_btm) = 0;
    island_nodes(end, :) = [];
end

% make sure each island loop is aligned clockwise
note_flag = 0;
island_counts = size(island_nodes, 2);
for ii = 1:island_counts
    tmp_nodes = island_nodes(:, ii);
    ind_valid = tmp_nodes~=0;

    if ispolycw(Mobj.lon(tmp_nodes(ind_valid)), Mobj.lat(tmp_nodes(ind_valid)))==0
        disp(['the island boundary (#', num2str(ii),') is aligned anti-clockwise and has been adjusted to clockwise'])
        island_nodes(ind_valid, ii) = flip(tmp_nodes(ind_valid),1);
        note_flag = 1;
    end
end

if note_flag==0
    disp('all island boundaries are aligned clockwise (great)')
end
%% Check the land-sea loop (outer loop)
% make sure land/open boundary are aligned anti-clockwise
obc_counts = size(obc_nodes,2);
land_counts = size(land_nodes,2);

% only work for the case that land and ocean boudaries forms a big loop,
% and all islands are inside it so far. 
outer_bnd_cells = cell(obc_counts+land_counts,1);  
tmp_nodes = obc_nodes(:,1);
tmp_nodes(tmp_nodes==0) = [];
outer_bnd_cells{1} = tmp_nodes(:);
ind_flags = ones(obc_counts+land_counts,1);
for ii = 2:obc_counts+land_counts
    head_nodes = outer_bnd_cells{ii-1};
    switch mod(ii,2)
        case 0
            next_nodes = land_nodes;
        case 1
            next_nodes = obc_nodes;
    end
    [idx1, idx2] = find(next_nodes==head_nodes(end));
    tmp_nodes = next_nodes(:,idx2); 
    tmp_nodes(tmp_nodes==0) = []; tmp_nodes = tmp_nodes(:);
    if idx1~=1
        tmp_nodes = flip(tmp_nodes,1);
        ind_flags(ii) = -1;
    end
    outer_bnd_cells{ii} = tmp_nodes(:);
end
outer_loop = unique(cell2mat(outer_bnd_cells), 'stable');

% check the land-sea loop direction
if ispolycw(Mobj.lon(outer_loop), Mobj.lat(outer_loop)) == 1  % if the land-sea loop is aligned clockwise
    outer_bnd_cells = cellfun(@(x) flip(x,1), outer_bnd_cells, 'UniformOutput',false);
    ind_flags = ind_flags*-1;
end

% arrange the land/open boundary nodes in descending order
obc_flags = ind_flags(1:2:end);
obc_cells = outer_bnd_cells(1:2:end);
obc_lens = cellfun(@(x) length(x), obc_cells);
[~, ind_sorted] = sort(obc_lens, 'descend');
obc_cells = obc_cells(ind_sorted); % sort by the obc  lens
obc_flags = obc_flags(ind_sorted);

land_flags = ind_flags(2:2:end);
land_cells = outer_bnd_cells(2:2:end);
land_lens = cellfun(@(x) length(x), land_cells);
[~, ind_sorted] = sort(land_lens, 'descend');
land_cells = land_cells(ind_sorted); % sort by the land lens
land_flags = land_flags(ind_sorted);

note_flag = 0;
obc_nodes = zeros(max(obc_lens), obc_counts);
for ii = 1:obc_counts
    tmp_nodes = obc_cells{ii};
    obc_nodes(1:numel(tmp_nodes), ii) = tmp_nodes;
    if obc_flags(ii)==-1
        disp(['the open boundary (#', num2str(ii),') is aligned clockwise and has been adjusted to anti-clockwise'])
        note_flag = 1;
    end
end
if note_flag==0
    disp('all open boundaries are aligned anti-clockwise (great)')
end

note_flag = 0;
land_nodes = zeros(max(land_lens), land_counts);
for ii = 1:land_counts
    tmp_nodes = land_cells{ii};
    land_nodes(1:numel(tmp_nodes), ii) = tmp_nodes;
    if land_flags(ii)==-1
        disp(['the land boundary (#', num2str(ii),') is aligned clockwise and has been adjusted to anti-clockwise'])
        note_flag = 1;
    end
end
if note_flag==0
    disp('all land boundaries are aligned anti-clockwise (great)')
end
%% Add boundary info
%========= open boundary=========
Mobj.obc_nodes = obc_nodes;
Mobj.obc_counts = size(Mobj.obc_nodes, 2);
Mobj.obc_lens = sum(Mobj.obc_nodes~=0)';
Mobj.obc_nodes_tot = Mobj.obc_nodes(:);
Mobj.obc_nodes_tot(Mobj.obc_nodes_tot==0) = [];

%========= land boundary=========
Mobj.land_nodes = land_nodes;
Mobj.land_counts = size(Mobj.land_nodes,2);
Mobj.land_lens = sum(Mobj.land_nodes~=0)';
Mobj.land_nodes_tot = Mobj.land_nodes(:);
Mobj.land_nodes_tot(Mobj.land_nodes_tot==0) = [];

%========= island boundary=========
Mobj.island_nodes = island_nodes;
Mobj.island_counts = size(Mobj.island_nodes,2);
Mobj.island_lens = sum(Mobj.island_nodes~=0)';
Mobj.island_nodes_tot = Mobj.island_nodes(:);
Mobj.island_nodes_tot(Mobj.island_nodes_tot==0) = [];

disp('open/land/island boundary metrics have been added')
end

function bnd_nodes = trim_bnd_nodes(bnd_nodes, sort_flag)
% 1) arrange all boundaries in descending order based on the # of nodes
% 2) trim redundant lines that are all zeros

bnd_nodes = abs(bnd_nodes);
max_len = max(sum(bnd_nodes~=0,1));
bnd_nodes = bnd_nodes(1:max_len,:);

if sort_flag == 1
    bnd_lens = sum(bnd_nodes~=0,1);
    [~, ind_sorted] = sort(bnd_lens, 'descend');
    bnd_nodes = double(bnd_nodes(:, ind_sorted)); % sort by the # of loop lens
end
end

%% Debug (check the direction of boundaries)
% tmp_nodes = Mobj.land_nodes(:,1);
% tmp_nodes(tmp_nodes==0) = [];
% figure
% sz = rescale(1:numel(tmp_nodes))*60+1;
% scatter(Mobj.lon(tmp_nodes), Mobj.lat(tmp_nodes), sz)
% hold on
% scatter(Mobj.lon(tmp_nodes(1)), Mobj.lat(tmp_nodes(1)), 20, 'red')
% plot(Mobj.lon(tmp_nodes), Mobj.lat(tmp_nodes))
% plot_schism_bnds(Mobj)
