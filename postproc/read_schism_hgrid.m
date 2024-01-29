function Mobj = read_schism_hgrid(Mobj, hgrid_file)
% Read the horizontal grids from hgrid.gr3                            
% 
%% Syntax
% Mobj = read_schism_hgrid(Mobj, hgrid_file)
% 
%% Description 
% Mobj = read_schism_hgrid(Mobj, hgrid_file) reads horizontal grid 
% information from the hgrid.gr3 file 
% 
%% Example
% Mobj.expName = 'test';
% hgrid_file = 'E:/Exp1/inputs/hgrid.gr3';
% Mobj = read_schism_hgrid(Mobj, hgrid_file)
%
%% Input Arguments
% Mobj --- the mesh object
% hgrid_file --- the absolute filepath of the hgrid.gr3 file.
% 
%% Output Arguments
% Mobj --- the mesh object with the information of horizontal grids loaded
% 
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2023-12-08. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: read_schism_vgrid and importdata

%% Parse inputs
Mobj.aimpath = fileparts(hgrid_file);
D = importdata(hgrid_file, '%/s', inf);
D = cellfun(@(x) strtrim(x), D, 'UniformOutput',false);  % Adapt to earlier versions of MATLAB

%% Basic mesh info.
head_info = strsplit(D{2});
Mobj.nElems = str2double(head_info{1});
Mobj.nNodes = str2double(head_info(2));

NodePart = double(split(string(D(3:3+Mobj.nNodes-1))));
NodePart(:, isnan(sum(NodePart, 1))) = [ ];

Mobj.lon = NodePart(:,2);
Mobj.lat = NodePart(:,3);
Mobj.depth = NodePart(:, 4);

ElemPart = double(split(string(D(3+Mobj.nNodes:3+Mobj.nNodes+Mobj.nElems-1))));
ElemPart(:, isnan(sum(ElemPart, 1))) = [ ];

Mobj.tri = ElemPart(:,3:5);
Mobj.region = round([min(Mobj.lon)-0.5 max(Mobj.lon)+0.5 min(Mobj.lat)-0.5 max(Mobj.lat)+0.5], 2);

Mobj.lonc = mean(Mobj.lon(Mobj.tri), 2);
Mobj.latc = mean(Mobj.lat(Mobj.tri), 2);
Mobj.depthc = mean(Mobj.depth(Mobj.tri), 2);
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

Mobj.obc_nodes = obc_nodes;
Mobj.obc_counts = obc_counts;
Mobj.obc_lens = obc_lens;

tmp = Mobj.obc_nodes(:);
tmp(tmp==0) = [];
Mobj.obc_nodes_tot = tmp;

ind = Mobj.obc_nodes;
ind_nan = ind==0;
ind(ind==0) = 1;
lon_tmp = Mobj.lon(ind);
lat_tmp = Mobj.lat(ind);

lon_tmp(ind_nan) = nan;
lat_tmp(ind_nan) = nan;
Mobj.lon_obc = lon_tmp;
Mobj.lat_obc = lat_tmp;

Mobj = get_obc_elems(Mobj);
Mobj.nElems_obc = numel(find(Mobj.obc_elems(:)~=0));
Mobj.nNodes_obc= sum(Mobj.obc_lens);

disp('open boundary nodes have been added')
%% Read the land nodes (island part is included now)
ind_cut = 3+Mobj.nNodes+Mobj.nElems+Mobj.nNodes_obc+Mobj.obc_counts+2; % land boundary part begins from this line
hgrid_str = double(split(string(D{ind_cut})));

land_counts = hgrid_str(1);
land_lens = nan(1, land_counts);
bdry_flags = nan(1, land_counts);
land_nodes = [];
for ii = 1:land_counts
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
    land_nodes = wisecat(land_nodes, tmp_nodes(:));

    land_lens(ii) = N;
end

Mobj.land_nodes = land_nodes;
Mobj.land_counts = land_counts;
Mobj.land_lens = land_lens;

tmp = Mobj.land_nodes(:);
tmp(tmp==0) = [];
Mobj.land_nodes_tot = tmp;
Mobj.nNodes_land = sum(Mobj.land_lens);

disp('land nodes have been added')
%% Speprate the island nodes from the land part
ind_island = logical(bdry_flags);
Mobj.island_nodes = Mobj.land_nodes(:, ind_island);
Mobj.island_counts = numel(find((ind_island)));
Mobj.island_lens = Mobj.land_lens(ind_island);

max_island_len = max(Mobj.island_lens);
Mobj.island_nodes = Mobj.island_nodes(1:max_island_len, :);

tmp = Mobj.island_nodes(:);
tmp(tmp==0) = [];
Mobj.island_nodes_tot = tmp;
Mobj.nNodes_island = sum(Mobj.island_lens);

disp('island nodes have been added')
%% Modify the land outputs
Mobj.land_nodes(:, ind_island) = [ ];
Mobj.land_counts = numel(find((~ind_island)));
Mobj.land_lens(ind_island) = [];

tmp = Mobj.land_nodes(:);
tmp(tmp==0) = [];
Mobj.land_nodes_tot = tmp;
Mobj.nNodes_land = sum(Mobj.land_lens);

end

function C = wisecat(A, B, maxLen, padval)
B = B(:); 
if dimnum(A) == 1
    A = A(:);
end
if nargin < 3
    maxLen = max([length(A), length(B)]);
end
if nargin < 4
    padval = 0;
end
C = [padarray(A, [maxLen-length(A) 0], padval, 'post') padarray(B,[maxLen-length(B) 0], padval, 'post')];
end

function Mobj = get_obc_elems(Mobj)
% Calculate the index of open boundary nodes.

obc_nodes = Mobj.obc_nodes;
Mobj.obc_elems = zeros(size(obc_nodes));
nBnds = size(obc_nodes,2);
for ibnd = 1:nBnds
    indVert1 = Mobj.tri(:,1)==obc_nodes(:,ibnd)';
    indVert2 = Mobj.tri(:,2)==obc_nodes(:,ibnd)';
    indVert3 = Mobj.tri(:,3)==obc_nodes(:,ibnd)';
    indVert = indVert1+indVert2+indVert3;
    tmp = sum(indVert,2);
    tmp = find(tmp>1);
    Mobj.obc_elems(1:length(tmp),ibnd) = tmp;
end
end




















