function Mobj = read_schism_hgrid(Mobj, hgrid_file)
% Read the horizontal grids from hgrid.gr3/hgrid.ll file.
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
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
% hgrid_file - hgrid filepath; char
%       the absolute filepath of the hgrid.gr3 or hgrid.ll file.
% 
%% Output Arguments
% Mobj - mesh object; datastruct
%       the mesh object with the updated grid info.
% 
%% Notes
% This function will check the rotation directions of land, island, or
% ocean boundaries in the provided hgrid file, however, it will not change
% the index of boundaries. If the rotation direction is problematic, please
% re-generate a correct hgrid file, otherwise the resulting input files
% cannot match the old grid!
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 28 May 2025.
% Email: wwu@vims.edu
% 
% See also: read_schism_vgrid

%% Parse inputs
Mobj.aimpath = [fileparts(hgrid_file), '\'];
fid = fopen(hgrid_file); D = textscan(fid, '%s', 'Delimiter', '\n'); fclose(fid);
D = strtrim(D{1});
disp('read horizontal grid from the hgrid file')

%% Basic mesh info.
head_info = strsplit(D{2});
nElems = str2double(head_info{1});
nNodes = str2double(head_info(2));

node_part = reshape(sscanf(strjoin(D(3:3+nNodes-1)'), '%f'), [], nNodes)';
lon = node_part(:,2); lat = node_part(:,3); depth = node_part(:, 4);

elem_info = D(3+nNodes:3+nNodes+nElems-1);
i34 = cellfun(@(s) numel(regexp(s, '\s+')), elem_info) - 1;  % consider consecutive spaces

elem3_info = double(split(string(elem_info(i34==3))));
elem4_info = double(split(string(elem_info(i34==4))));

% works for triangular and quad mesh now.
elem_part = nan(nElems, 6);
elem_part(i34==3, 1:5) = elem3_info;
elem_part(i34==4, :) = elem4_info;

tri = elem_part(:,3:6);
Mobj = add_grid_metrics(Mobj, lon, lat, tri, depth);  % add grid info.

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

ind_island = mod(bdry_flags,10)==1;   % "1" or "11" stands for island.
land_nodes = land_island_nodes(:, ~ind_island);
island_nodes = land_island_nodes(:, ind_island);

% Trim redundant zeros
land_nodes = land_nodes(1:max(sum(land_nodes~=0,1)), :);
island_nodes = island_nodes(1: max(sum(island_nodes~=0,1)), :);

%% Add boundary info
% the boundary index will not be changed, but the rotation directions will
% be checked and corrected if problems were detected.
Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes, 0);
if ~isequal(Mobj.obc_nodes, obc_nodes)
    warning on;
    warning('the provided hgrid file has rotation issues, please re-generate new hgrid file!')
end
end

function C = wisecat(A, B, maxLen, padval)

B = B(:); 
if dimnum(A) == 1; A = A(:); end
if nargin < 3; maxLen = max([size(A,1), size(B,1)]); end
if nargin < 4; padval = 0; end

C = [padarray(A, [maxLen-size(A,1) 0], padval, 'post') padarray(B, [maxLen-size(B,1) 0], padval, 'post')];

end




















