function Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes)
% Add ocean/land/island boundary info.
%
%% Syntax
%
%
%% Description
%
%
%% Examples 
% 
%
%% Input Arguments
%
%
%% Output Arguments
%
%% Notes
% The outer boundaies (land/ocean) should be aligned anti-clockwise while the
% inner boundaries (island) should be clockwise.
% 
% This function was generated with the help of ChatGPT
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 1 Nov 2024. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
% 1) the islands don't need to form a loop in hgrid.gr3 or hgrid.ll
ind_btm = sub2ind(size(island_nodes), sum(island_nodes~=0, 1), 1:size(island_nodes,2));
if sum(island_nodes(1,:) ~= island_nodes(ind_btm))==0
    disp('removed ending points of islands, since all islands start/end at the same node')
    island_nodes(ind_btm) = 0;
end

% 2) arrange all boundaries in descending order based on the # of nodes
% 3) trim redundant lines that are all zeros
% 4) ensure outer/inner boundary nodes are aligned anti-clockwise/clockwise
obc_nodes = trim_bnd_nodes(Mobj, obc_nodes, 'open');
land_nodes = trim_bnd_nodes(Mobj, land_nodes, 'land');
island_nodes = trim_bnd_nodes(Mobj, island_nodes, 'island');

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
% Mobj.nNodes_land = sum(Mobj.land_lens);
Mobj.land_nodes_tot = Mobj.land_nodes(:);
Mobj.land_nodes_tot(Mobj.land_nodes_tot==0) = [];

%========= island boundary=========
Mobj.island_nodes = island_nodes;
Mobj.island_counts = size(Mobj.island_nodes,2);
Mobj.island_lens = sum(Mobj.island_nodes~=0)';
% Mobj.nNodes_island = sum(Mobj.island_lens);
Mobj.island_nodes_tot = Mobj.island_nodes(:);
Mobj.island_nodes_tot(Mobj.island_nodes_tot==0) = [];

disp('land/island/open boundary metrics have been added')
end

function bnd_nodes = trim_bnd_nodes(Mobj, bnd_nodes, bnd_str)
% This function has three goals:
% 1) arrange all boundaries in descending order based on the # of nodes
% 2) trim redundant lines that are all zeros
% 3) ensure all boundary nodes are placed anti-clockwise

max_len = max(sum(bnd_nodes~=0,1));
bnd_nodes = bnd_nodes(1:max_len,:);

bnd_lens = sum(bnd_nodes~=0,1);
[~, ind_sorted] = sort(bnd_lens, 'descend');
bnd_nodes = bnd_nodes(:, ind_sorted); % sort by the # of loop lens

bnd_counts = size(bnd_nodes, 2);
for ii = 1:bnd_counts
    bnd_tmp = bnd_nodes(:, ii);
    ind_valid = bnd_tmp~=0;
    if strcmp(bnd_str, 'island')
        cw_flag = 0; old_str = 'anti-clockwise'; new_str = 'clockwise';
    else
        cw_flag = 1; old_str = 'clockwise'; new_str = 'anti-clockwise';
    end
    if ispolycw(Mobj.lon(bnd_tmp(ind_valid)), Mobj.lat(bnd_tmp(ind_valid)))==cw_flag
        disp(['the ',bnd_str,' boundary (#', num2str(ii),') is aligned ', old_str,' and has been adjusted to ', new_str])
        bnd_nodes(ind_valid, ii) = flip(bnd_tmp(ind_valid),1);
    end
end

% debug
% sz = rescale(1:numel(bnd_tmp(ind_valid)))*20+1;
% scatter(Mobj.lon(bnd_tmp(ind_valid)), Mobj.lat(bnd_tmp(ind_valid)), sz)

end

% function Mobj = get_obc_elems(Mobj)
% % calculate the index of open boundary nodes.
% 
% obc_nodes = Mobj.obc_nodes;
% Mobj.obc_elems = zeros(size(obc_nodes));
% nBnds = size(obc_nodes,2);
% for ibnd = 1:nBnds
%     indVert1 = Mobj.tri(:,1)==obc_nodes(:,ibnd)';
%     indVert2 = Mobj.tri(:,2)==obc_nodes(:,ibnd)';
%     indVert3 = Mobj.tri(:,3)==obc_nodes(:,ibnd)';
%     indVert = indVert1+indVert2+indVert3;
%     tmp = sum(indVert,2);
%     tmp = find(tmp>1);
%     Mobj.obc_elems(1:length(tmp),ibnd) = tmp;
% end
% 
% Mobj.obc_elems(sum(Mobj.obc_elems, 2)==0, :) = [];
% end
