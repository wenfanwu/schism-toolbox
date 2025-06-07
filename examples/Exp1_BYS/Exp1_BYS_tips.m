%% This program shows some TIPS when using the schism-toolbox
%% Tested platform: Matlab 2024a (Windows)
%% Matlab add-on: Image Processing Toolbox; Mapping Toolbox
%% Public package: OceanMesh2D
%% Model verision: SCHISM v5.10
%% Author: Wenfan Wu, CCRM, Virginia Institute of Marine Science. 2024
%% Load the mesh grid
clc;clearvars
mesh_file = 'Exp1_BYS\inputs\BYS_20814.mat';  % NEED TO BE CHANGED!!!

Mobj = mesh2schism(mesh_file); 
Mobj.expname = 'Exp1_BYS';      
Mobj.time = (datetime(2020,6,1):hours(1):datetime(2020,6,10))'; 
Mobj.rundays = days(Mobj.time(end)-Mobj.time(1)+1); 
Mobj.dt = 150; % dt (secs), the same as in param.nml
Mobj.coord = 'geographic'; % geographic or Cartesian coordinate

%% Load the vertical grids
dep_edges = [10, 20, 30, 45, 55, 65, 75, 90];
dep_nums =  [20 21 22 23 24 25 27 28];
Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_nums, [4 5 3 5], 0.25);

%% Tip-1: display the variables@Nodes
figure('Color', 'w')
subplot(221)
disp_schism_var(Mobj, Mobj.depth)
axis image
subplot(222)
disp_schism_var(Mobj, Mobj.depth,'EdgeColor', 'k')
axis image
subplot(223)
disp_schism_var(Mobj, Mobj.depth,'EdgeColor', 'r', 'EdgeAlpha', 0.5, 'LineWidth', 0.15)
axis image
subplot(224)
disp_schism_var(Mobj, Mobj.depth,'EdgeColor', 'b', 'FaceColor', 'none')
axis image

%% Tip-2: display the variables@Elems or @Sides
figure('Color', 'w')
disp_schism_var(Mobj, Mobj.depthc,'EdgeColor', 'k')
axis image
auto_center

figure('Color', 'w')
disp_schism_var(Mobj, Mobj.depths)
axis image
auto_center

%% Tip-3: find the values at the seabed
varTest = Mobj.vgrids;

varBtm = get_schout_btm(Mobj, varTest, 0);
figure('Color', 'w')
disp_schism_var(Mobj, varBtm)
axis image
title('sigma depth (0-1)')

varBtm = get_schout_btm(Mobj, varTest, 1);
figure('Color', 'w')
disp_schism_var(Mobj, varBtm)
axis image
title('sigma depth (0-1)')

%% Tip-4: define a mask on the triangular mesh
varTest = Mobj.depth;

figure('Color', 'w')
disp_schism_var(Mobj, varTest)
axis image
hold on
msk = def_schism_mask(Mobj, 2, 'test', 'rebuild');  % "2" means two separate regions

% draw a polygon on the map and press ENTER, and the mask info. will be
% saved into the aimpath as a MAT file automatically 

varTest(msk) = nan;
figure('Color', 'w')
disp_schism_var(Mobj, varTest)
axis image

%% Tip-5: interplolate variables from nodes onto elemes
varNode = Mobj.depth;
varElem = convert_schism_var(Mobj, varNode, 'node2elem');

figure('Color', 'w')
subplot(211)
disp_schism_var(Mobj, varNode)
axis image

subplot(212)
disp_schism_var(Mobj, varElem)
axis image

%% Tip-6: interplolate variables from elems onto nodes
varElem = Mobj.depthc;
varNode = convert_schism_var(Mobj, varNode, 'elem2node');

figure('Color', 'w')
subplot(211)
disp_schism_var(Mobj, varNode)
axis image

subplot(212)
disp_schism_var(Mobj, varElem)
axis image

%% Tip-7: extract data along the transect
figure
disp_schism_hgrid(Mobj, [1 0])
axis image
auto_center
hold on
sect_info = def_schism_transect(Mobj, 1);

var_tri = Mobj.depLayers; % test data

[var2d, dist2d, dep2d] = read_schism_transect(Mobj, sect_info, var_tri);

figure
pcolor(dist2d/1e3, dep2d, var2d)
shading interp
colorbar
colormap(jet)
xlabel('Along transect distance (km)')
ylabel('Depth (m)')

%% Tip-8: prepare for particle tracking
figure('Color', 'w')
disp_schism_var(Mobj, Mobj.depth)
axis image
hold on
[lon_list, lat_list] = def_schism_ptrack(Mobj, 1, 'polygon'); % draw a polygon on the map and the selected nodes will emerge

dep_list = repmat(-3, [nps 1]);  % all the particles are released at the 3-m depth
xyz_data = [lon_list(:), lat_list(:), dep_list(:)];
nps = numel(lon_list);
drop_time_list = repmat(datetime(2020,6,3), [nps 1]);  % all the particles are released on 3 Jun 2020.
life_day = 3;
ptrack_vars = [1 0 1 0];

% particle.bp file can be found in 'Exp1_BYS/inputs'
write_schism_ptrack(Mobj, xyz_data, drop_time_list, life_day, ptrack_vars)
%% Tip-9: prepare station.in file
figure('Color', 'w')
disp_schism_var(Mobj, Mobj.depth)
axis image
hold on
[lon_list, lat_list] = def_schism_ptrack(Mobj, 1, 'polygon');  % draw a polygon on the map and the selected nodes will emerge

nps = numel(lon_list);
dep_list = zeros(nps, 1);  % all the particles are released at the 0-m depth
xyz_data = [lon_list(:), lat_list(:), dep_list(:)];
switch_flags = [1 1 1 1 1 1 1 1 1];

% station.in file can be found in 'Exp1_BYS/inputs'
write_schism_station_in(Mobj, xyz_data, switch_flags)

%% END



















