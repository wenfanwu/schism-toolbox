function Mobj = mesh2schism(meshfile)
% Load the horizontal grids
% 
%% Syntax
% Mobj = mesh2schism(meshfile)
%
%% Description 
% Mobj = mesh2schism(meshfile) returns a datastruct containing the basic
% information of the horizontal grids.
% 
%% Example
% mesh_file  = 'E:\Exp2_BYES\inputs\BYES_64676.mat';
% Mobj = mesh2schism(mesh_file);
%
%% Input Arguments
% meshfile --- the absolute filepath of the mesh file generated from the 
% OceanMesh2D toolbox.
% 
%% Output Arguments
% Mobj --- a datastruct containing the basic information of horizontal grids.
% 
%% Notes
% this function aims to load the mesh grid generated from OceanMesh2D
% toolbox, and save it as a datastruct named "Mobj", meaning "Mesh object".
% The obtained datastruct contains the basic information of the mesh grid,
% including the # of nodes and elements, depth and so on. In addition, this
% function will check whether the open/land boundary nodes are
% anti-clockwise. If not, it will modify them automatically. 
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 9 Nov. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: ispolycw

%% Parse inputs
if contains(meshfile, ' ') == 1
    error('There must be no Spaces in the filepath!')
end

Mobj.aimpath = [fileparts(meshfile), '\'];
load(meshfile, 'm')

% node metrics
Mobj.nNodes = length(m.p(:,1));
Mobj.lon = m.p(:,1);
Mobj.lat = m.p(:,2);
Mobj.depth = m.b;
Mobj.slope =  hypot(m.bx, m.by);

% element metrics
Mobj.nElems = size(m.t,1);
Mobj.tri = m.t;  % The vertex order of the triangle mesh must be counterclockwise in SCHISM
Mobj.lonc = mean(Mobj.lon(Mobj.tri), 2);
Mobj.latc = mean(Mobj.lat(Mobj.tri), 2);
Mobj.depthc = mean(Mobj.depth(Mobj.tri), 2);

% side metrics
T = Mobj.tri;
P = [Mobj.lon, Mobj.lat];
TR = triangulation(T,P);
E = edges(TR);
Mobj.nSides = size(E,1);
Mobj.edges = E;
% Mobj.lons = mean(Mobj.lon(Mobj.edges), 2);  % lons and lonc should be obtained from the model ouputs
% Mobj.lats = mean(Mobj.lat(Mobj.edges), 2);

Mobj.region = round([min(Mobj.lon)-0.5 max(Mobj.lon)+0.5 min(Mobj.lat)-0.5 max(Mobj.lat)+0.5], 2);

%% Check the rotation direction
% open boundary nodes
op_bnds = m.op.nbdv;
bnd_lens = m.op.nvdll;
nBnds = size(op_bnds, 2);
disp([num2str(nBnds), ' open boundaries were considered'])
bnd_tmp = op_bnds(:,1);
bnd_tmp(bnd_tmp==0) = [];
if ispolycw(Mobj.lon(bnd_tmp), Mobj.lat(bnd_tmp)) == 1
    disp('The open boundary nodes are aligned clockwise and has been adjusted to anti-clockwise')
    tmp = cell2mat(arrayfun(@(x) circshift(op_bnds(:,x), -bnd_lens(x), 1), 1:nBnds,'UniformOutput',false));
    obcNodes_cali = flipud(tmp);
else
    disp('The open boundary nodes are aligned anti-clockwise')
    obcNodes_cali = op_bnds;
end

% land boundary nodes
landLens = m.bd.nvell;
landBnds = m.bd.nbvv;  % land ID is 20; island ID is 21
nLands = size(landBnds, 2);
disp([num2str(nLands), ' land and island boundaries were considered'])
land_tmp = landBnds(:,1);
land_tmp(land_tmp==0) = [];
if ispolycw(Mobj.lon(land_tmp), Mobj.lat(land_tmp)) == 1
    disp('The land boundary nodes are aligned clockwise and has been adjusted to anti-clockwise')
    tmp = cell2mat(arrayfun(@(x) circshift(landBnds(:,x), -landLens(x), 1), 1:nLands,'UniformOutput',false));
    landNodes_cali = flipud(tmp);
else
    disp('The land boundary nodes are aligned anti-clockwise')
    landNodes_cali = landBnds;
end

%% Load the ocean/land/island info.
Mobj.obc_nodes = obcNodes_cali;  
Mobj.obc_counts = size(Mobj.obc_nodes, 2);
Mobj.obc_lens = sum(Mobj.obc_nodes~=0);
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

% the land ID is 20
indLand = m.bd.ibtype==20;
Mobj.land_nodes = landNodes_cali(:,indLand);
Mobj.land_counts = size(Mobj.land_nodes,2);
Mobj.land_lens = sum(Mobj.land_nodes~=0);
Mobj.nNodes_land = sum(Mobj.land_lens);
tmp = Mobj.land_nodes(:);
tmp(tmp==0) = [];
Mobj.land_nodes_tot = tmp;

% the island ID is 21
indLand2 = m.bd.ibtype==21;
Mobj.island_nodes = landNodes_cali(:,indLand2);
Mobj.island_counts = size(Mobj.island_nodes,2);
Mobj.island_lens = sum(Mobj.island_nodes~=0);
Mobj.nNodes_island = sum(Mobj.island_lens);
tmp = Mobj.island_nodes(:);
tmp(tmp==0) = [];
Mobj.island_nodes_tot = tmp;

% outer boundary nodes
bnd_nodes = [Mobj.land_nodes_tot; Mobj.obc_nodes_tot];
Mobj.lon_bnd = Mobj.lon(bnd_nodes);
Mobj.lat_bnd = Mobj.lat(bnd_nodes);

end

function Mobj = get_obc_elems(Mobj)
% calculate the index of open boundary nodes.

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

% !...  Count # of tracer models and tracers
% !...  Each trace module has a pre-defined ID as follows:
% !     1: T
% !     2: S
% !     3: GEN
% !     4: AGE
% !     5: SED3D
% !     6: EcoSim
% !     7: ICM
% !     8: CoSINE
% !     9: Feco
% !    10: TIMOR
% !    11: FABM
% !    12: DVD numerical mixing analysis of Klingbeit
