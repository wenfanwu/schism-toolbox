function Mobj = read_schism_hgrid(Mobj, hgrid_file)
% Read the horizontal grids from hgrid.gr3/hgrid.ll (works for triangular/quad)
% 
%% Syntax
% Mobj = read_schism_hgrid(Mobj, hgrid_file)
% 
%% Description 
% Mobj = read_schism_hgrid(Mobj, hgrid_file) reads horizontal grid 
% information from the hgrid file 
% 
%% Example
% Mobj.expName = 'test';
% hgrid_file = 'E:/Exp1/inputs/hgrid.gr3';
% Mobj = read_schism_hgrid(Mobj, hgrid_file)
%
%% Input Arguments
% Mobj --- the mesh object
% hgrid_file --- the absolute filepath of the hgrid.gr3 or hgrid.ll file.
% 
%% Output Arguments
% Mobj --- the mesh object with the information of horizontal grids loaded
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 17 Sep 2024. 
% Email: wwu@vims.edu
% 
% See also: read_schism_vgrid and importdata

%% Parse inputs
Mobj.aimpath = [fileparts(hgrid_file), '\'];
D = importdata(hgrid_file, '%/s', inf);
D = cellfun(@(x) strtrim(x), D, 'UniformOutput',false);  % Adapt to earlier MATLAB versions

%% Basic mesh info.
head_info = strsplit(D{2});
nElems = str2double(head_info{1});
nNodes = str2double(head_info(2));

node_part = double(split(string(D(3:3+nNodes-1))));
node_part(:, isnan(sum(node_part, 1))) = [ ];

lon = node_part(:,2);
lat = node_part(:,3);
depth = node_part(:, 4);

elem_info = D(3+nNodes:3+nNodes+nElems-1);
i34 = cellfun(@(x) numel(strsplit(x)), elem_info)-2;

elem3_info = double(split(string(elem_info(i34==3))));
elem4_info = double(split(string(elem_info(i34==4))));

% works for triangular and quad mesh now.
elem_part = nan(nElems, 6);
elem_part(i34==3, 1:5) = elem3_info;
elem_part(i34==4, :) = elem4_info;

tri = elem_part(:,3:6);

Mobj = add_grid_metrics(Mobj, lon, lat, tri, depth);
%% Read the open boundary nodes
ind_cut = 3+Mobj.nNodes+Mobj.nElems; % open boundary part begins from this line
hgrid_str = double(split(string(D{ind_cut})));

obc_counts = hgrid_str(1);
obc_lens = nan(1, obc_counts);
obc_nodes = [];
for ii = 1:obc_counts
    switch ii
        case 1
            begind = ind_cut+2;
        otherwise
            begind = begind+N+1;
    end

    bdry_str = double(split(string(D{begind})));
    N = bdry_str(1);
    ind_nodes = begind+1:begind+N;
    tmp_nodes = double(string(D(ind_nodes)));
    obc_nodes = wisecat(obc_nodes, tmp_nodes(:));

    obc_lens(ii) = N;
end

nNodes_obc= sum(obc_lens);
%% Read the land and island nodes
ind_cut = 3+nNodes+nElems+nNodes_obc+obc_counts+2; % land boundary part begins from this line
hgrid_str = double(split(string(D{ind_cut})));

land_island_counts = hgrid_str(1);
land_island_lens = nan(1, land_island_counts);
bdry_flags = nan(1, land_island_counts);
land_island_nodes = [];
for ii = 1:land_island_counts
    switch ii
        case 1
            begind = ind_cut+2;
        otherwise
            begind = begind+N+1;
    end

    bdry_str = double(split(string(D{begind})));
    N = bdry_str(1);
    bdry_flags(ii) = bdry_str(2);  % used to distinguish land and island nodes

    ind_nodes = begind+1:begind+N;
    tmp_nodes = double(string(D(ind_nodes)));
    land_island_nodes = wisecat(land_island_nodes, tmp_nodes(:));

    land_island_lens(ii) = N;
end

ind_island = logical(bdry_flags);
land_nodes = land_island_nodes(:, ~ind_island);
island_nodes = land_island_nodes(:, ind_island);

max_len = max(sum(land_nodes~=0,1));
land_nodes = land_nodes(1:max_len,:);

max_len = max(sum(island_nodes~=0,1));
island_nodes = island_nodes(1:max_len,:);

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

function C = wisecat(A, B, maxLen, padval)

B = B(:); 
if dimnum(A) == 1
    A = A(:);
end
if nargin < 3
    maxLen = max([size(A,1), size(B,1)]);
end
if nargin < 4
    padval = 0;
end
C = [padarray(A, [maxLen-size(A,1) 0], padval, 'post') padarray(B, [maxLen-size(B,1) 0], padval, 'post')];
end

% function Mobj = get_obc_elems(Mobj)
% % Calculate the index of open boundary nodes.
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




















