function Mobj = mesh2schism(mesh_file, sname)
% Load the unstructured grid from OceanMesh2D or SMS.
%
%% Syntax
% Mobj = mesh2schism(mesh_file)
% Mobj = mesh2schism(mesh_file, sname)
%
%% Description
% Mobj = mesh2schism(mesh_file) read the mesh info from 2dm/mat file into a
%       datastruct (Mobj)
% Mobj = mesh2schism(mesh_file, sname) specifies the mesh generation
%       software (SMS or OceanMesh2D)
%
%% Example
% mesh_file  = 'E:\Exp2_BYES\inputs\BYES_64676.mat';
% Mobj = mesh2schism(mesh_file);
%
% mesh_file = 'D:\WorkDisk\test.2dm';
% Mobj = mesh2schism(mesh_file);
%
%% Input Arguments
% mesh_file - mesh file; char
%       the absolute filepath of mesh file from OceanMesh2D or SMS.
% sname - software name; char
%       the mesh generation software; By default: the function will
%       automatically decide based on the file suffix
%
%% Output Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
%
%% Notes
% This function aims to load the mesh grid generated from OceanMesh2D/SMS
% and save it as a datastruct named "Mobj", meaning "Mesh object".
% The obtained datastruct contains essential information of mesh grid,
% including the # of nodes and elements, depth and so on.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 04 Nov 2024.
% Email: wwu@vims.edu
%
% See also: add_grid_metrics and add_bnd_metrics

%% Parse inputs
% if contains(mesh_file, ' ') == 1
%     error('there must be no spaces in the filepath!')
% end

if nargin<2
    switch mesh_file(end-3:end)
        case '.2dm'
            sname = 'SMS';
        case '.mat'
            sname = 'OceanMesh2D';
        otherwise
            error('unrecognized mesh file!')
    end
end
%% Load mesh
switch lower(sname)
    case 'sms'
        disp('read horizontal grid from the 2dm file (SMS)')
        Mobj = read_2dm_info(mesh_file);
    case 'oceanmesh2d'
        disp('read horizontal grid from the mat file (OceanMesh2D)')
        Mobj = read_mat_info(mesh_file);
    otherwise
        error('unrecognized software name!')
end

end

function Mobj = read_2dm_info(mesh_file)
%% Load the mesh info from 2dm file created by SMS.
fid = fopen(mesh_file, 'r'); D = textscan(fid, '%s', 'Delimiter', '\n'); fclose(fid);
D = strtrim(D{1}); D = D(3:end);

ind_node = startsWith(D, 'ND');  % startsWith was introduced in R2016b
ind_elem3 = startsWith(D, 'E3T');
ind_elem4 = startsWith(D, 'E4Q');
ind_elem = ind_elem3 | ind_elem4;

node_part = D(ind_node); elem_part = D(ind_elem);
i34 = startsWith(elem_part, 'E4Q')+3;
elem3_info = double(split(string(elem_part(i34==3))));
elem4_info = double(split(string(elem_part(i34==4))));
node_info = double(split(string(node_part)));

nElems = size(elem_part,1);
tri = nan(nElems, 4);
tri(i34==3, 1:3) = elem3_info(:,3:5);
if ~isempty(elem4_info)
    tri(i34==4, 1:4) = elem4_info(:,3:6);
end

lon = node_info(:,3); lat = node_info(:,4); depth = node_info(:,5);
if sum(depth~=0)==0
    warning on
    warning('No depth info in the 2dm file')
    depth = depth+nan;
end

Mobj.aimpath = [fileparts(mesh_file), '\'];
Mobj = add_grid_metrics(Mobj, lon, lat, tri, depth);

% boundary info
ind_obc = startsWith(D, 'NS');
if any(ind_obc)
    obc_part = cellfun(@(x) strtrim(x), D(ind_obc), 'UniformOutput',false);  % adapt to earlier versions of MATLAB
    end_flags = cellfun(@(x) count(x, '-'), obc_part);
    end_locs = find(end_flags); end_locs = [1; end_locs(:)]; end_locs(1) = 0;

    nNodes_est = length(obc_part)*10;
    obc_counts = sum(end_flags);

    obc_nodes = zeros(nNodes_est, obc_counts);
    for ii = 1:obc_counts
        loc = end_locs(ii)+1:end_locs(ii+1)-1;

        obc_line = double(split(string(obc_part(loc,:))));
        obc_line = obc_line(:, 2:end)'; obc_line = obc_line(:);

        obc_end = double(split(string(obc_part(end_locs(ii+1),:))));
        obc_end = obc_end(2:end)'; 
        obc_end = abs(obc_end(1:find(obc_end<=0)));

        obc_tot = [obc_line(:); abs(obc_end(:))];
        obc_nodes(1:numel(obc_tot), ii) = obc_tot(:);
    end
    if isempty(obc_nodes)
        warning on
        warning('open boundary nodes are not found in the 2dm file!')
    end

    [land_nodes, island_nodes] = find_land_island(Mobj.tri, Mobj.edg, obc_nodes);
    Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes);
else
    warning on
    warning('No open boundary info in the 2dm file')
end
end

function Mobj = read_mat_info(mesh_file)
%% Load the mesh info from MAT file created by OceanMesh2D
% essential information
load(mesh_file, 'm')

Mobj.aimpath = [fileparts(mesh_file), '\'];
lon = m.p(:,1);  lat = m.p(:,2);  tri = m.t; tri(:,4) = nan; depth = m.b(:);

Mobj = add_grid_metrics(Mobj, lon, lat, tri, depth);
% Mobj.slope =  hypot(m.bx, m.by);  % may be removed in the future for consistency

% load the land/island/open nodes
obc_nodes = m.op.nbdv;

land_island_nodes = m.bd.nbvv;
ind_land = m.bd.ibtype==20;  % the land ID is 20
land_nodes = land_island_nodes(:,ind_land);

ind_island = m.bd.ibtype==21;  % the island ID is 21
island_nodes = land_island_nodes(:,ind_island);

Mobj = add_bnd_metrics(Mobj, obc_nodes, land_nodes, island_nodes);
end


