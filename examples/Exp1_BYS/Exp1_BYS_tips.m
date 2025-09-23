%% This program shows some TIPS when using the schism-toolbox
%% Tested platform: Matlab 2024a (Windows)
%% Matlab add-on: Image Processing Toolbox; Mapping Toolbox
%% Public package: OceanMesh2D
%% Model verision: SCHISM v5.10
%% Author: Wenfan Wu, CCRM, Virginia Institute of Marine Science. 2024
%% Horizontal grid
clc;clearvars
mesh_file = 'Exp1_BYS\inputs\BYS_20814.mat';  % NEED TO BE CHANGED

Mobj = mesh2schism(mesh_file); 
Mobj.expname = 'Exp1_BYS';      
Mobj.time = (datetime(2020,6,1):hours(1):datetime(2020,6,10))'; 
Mobj.rundays = days(Mobj.time(end)-Mobj.time(1)+1); 
Mobj.dt = 150; % dt (secs), the same as in param.nml
Mobj.coord = 'geographic'; % geographic or Cartesian coordinate

%% Vertical grid
dep_edges = [10, 20, 30, 45, 55, 65, 75, 90];
dep_nums =  [20 21 22 23 24 25 27 28];
Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_nums, [4 5 3 5], 0.25);

%% Visualize data defined on different grid centers
% Variable defined at node centers
figure('Color', 'w')
subplot(221)
disp_schism_var(Mobj, Mobj.depth)
axis image

% Customize visualization styles
subplot(222)
disp_schism_var(Mobj, Mobj.depth,'EdgeColor', 'k')
axis image
subplot(223)
disp_schism_var(Mobj, Mobj.depth,'EdgeColor', 'r', 'EdgeAlpha', 0.5, 'LineWidth', 0.15)
axis image
subplot(224)
disp_schism_var(Mobj, Mobj.depth,'EdgeColor', 'b', 'FaceColor', 'none')
axis image

% Variable defined at element centers
figure('Color', 'w')
disp_schism_var(Mobj, Mobj.depthc)
axis image

% Variable defined at side centers
figure('Color', 'w')
disp_schism_var(Mobj, Mobj.depths)
axis image

%% Find values at the bottom layer
var_test = Mobj.vgrids;

% the lowest level
var_btm = get_schout_btm(Mobj, var_test, 0);
figure('Color', 'w')
disp_schism_var(Mobj, var_btm)
axis image
title('sigma depth (0-1)')

% the second-to-last level
var_btm = get_schout_btm(Mobj, var_test, 1);
figure('Color', 'w')
disp_schism_var(Mobj, var_btm)
axis image
title('sigma depth (0-1)')

%% Define masks on the model grid
var_test = Mobj.depth;

figure('Color', 'w')
disp_schism_var(Mobj, var_test)
axis image
hold on
msk = def_schism_mask(Mobj, 2, 'test', 'rebuild');  % "2" means two separate regions
% draw two polygons on the map and press ENTER, and the mask info. will be
% saved into the aimpath as a MAT file automatically 

var_test(msk) = nan;
figure('Color', 'w')
disp_schism_var(Mobj, var_test)
axis image

%% Convert variables between different grid centers
% from nodes to elements
v1 = Mobj.depth;
v2 = convert_schism_var(Mobj, v1, 'node2elem');

figure('Color', 'w')
subplot(211)
disp_schism_var(Mobj, v1)
axis image

subplot(212)
disp_schism_var(Mobj, v2)
axis image

% from nodes to sides
v1 = Mobj.depth;
v2 = convert_schism_var(Mobj, v1, 'node2side');

figure('Color', 'w')
subplot(211)
disp_schism_var(Mobj, v1)
axis image
subplot(212)
disp_schism_var(Mobj, v2)
axis image

% from elements to nodes
v1 = Mobj.depthc;
v2 = convert_schism_var(Mobj, v1, 'elem2node');

figure('Color', 'w')
subplot(211)
disp_schism_var(Mobj, v1)
axis image
subplot(212)
disp_schism_var(Mobj, v2)
axis image

% from sides to nodes
v1 = Mobj.depths;
v2 = convert_schism_var(Mobj, v1, 'side2node');

figure('Color', 'w')
subplot(211)
disp_schism_var(Mobj, v1)
axis image
subplot(212)
disp_schism_var(Mobj, v2)
axis image

%% Define transect and extract data
% extract data along a straight transect
figure
disp_schism_hgrid(Mobj, [1 0])
axis image; auto_center
hold on
sect_info = def_schism_transect(Mobj, -1);

disp_schism_vgrid(Mobj, sect_info)

var_tri = Mobj.depLayers; % test data
[var2d, dist2d, dep2d] = read_schism_transect(Mobj, sect_info, var_tri);

figure
pcolor(dist2d/1e3, dep2d, var2d)
shading interp
colorbar
colormap(jet)
xlabel('Along transect distance (km)')
ylabel('Depth (m)')

% define a curved transect along the isobaths
figure
disp_schism_hgrid(Mobj)
clim([-60 0])
sect_info = def_schism_transect(Mobj, -3, 30);  % curved transect

x = sect_info.lon; y = sect_info.lat;
tvec = sect_info.tvec; % tangential unit vector
nvec = sect_info.nvec; % normal unit vector

figure
disp_schism_hgrid(Mobj)
hold on
colormap(turbo(25))
plot_schism_bnds(Mobj)
h1 = plot(x, y, 'LineWidth', 3, 'Marker', '.', 'Color', 'g');
h2 = quiver(x, y, tvec(:,1), tvec(:,2), 0.5, 'r');
h3 = quiver(x, y, nvec(:,1), nvec(:,2), 0.5, 'b');
legend([h1,h2,h3], {'Transect', 'Tangent','Normal'});

%% Extract contour lines and export as shapefiles
levels = [5, 10, 30, 60]; % depth levels
S = calc_schism_contour(Mobj, Mobj.depth, levels);

% bathymetry contour lines
figure
hold on
for ii = 1:length(S)
    plot(S(ii).X, S(ii).Y)
end
axis image
legend(string(levels))

shapewrite(S, 'D:\test')  % save as shapefiles if necessary

%% Calculate gradient on unstructured grids
[Fx, Fy, Fxy] = calc_schism_grad(Mobj, Mobj.depth);

% bathymetry gradient
figure
disp_schism_var(Mobj, Fxy)
caxis([0 0.0025]) %#ok<*CAXIS>
axis image; auto_center

%% Prepare for particle tracking
figure('Color', 'w')
disp_schism_var(Mobj, Mobj.depth)
axis image
hold on
[lon_list, lat_list] = def_schism_ptrack(Mobj, 1, 'polygon'); % draw a polygon on the map and the selected nodes will emerge

nps = numel(lon_list);
dep_list = repmat(-3, [nps 1]);  % all the particles are released at the 3-m depth
xyz_data = [lon_list(:), lat_list(:), dep_list(:)];
drop_time_list = repmat(datetime(2020,6,3), [nps 1]);  % all the particles are released on 3 Jun 2020.
life_day = 3;
ptrack_vars = [1 0 1 0];

% particle.bp file can be found in 'Exp1_BYS/inputs'
write_schism_ptrack(Mobj, xyz_data, drop_time_list, life_day, ptrack_vars)

%% Generate station.in file
figure('Color', 'w')
disp_schism_var(Mobj, Mobj.depth)
axis image
hold on
[lon_list, lat_list] = def_schism_ptrack(Mobj, 1, 'polygon'); 

nps = numel(lon_list);
dep_list = zeros(nps, 1); 
xyz_data = [lon_list(:), lat_list(:), dep_list(:)];
switch_flags = [1 1 1 1 1 1 1 1 1];

% station.in file can be found in 'Exp1_BYS/inputs'
write_schism_station_in(Mobj, xyz_data, switch_flags)

%% Re-define the open boundaries
obc_counts = 2;  % two open boundaries will be defined.
Mobj2 = def_schism_obc(Mobj, obc_counts);

% draw two polygons on the map, ensuring that each polygon encloses a
% contiguous boundary segment.

figure
disp_schism_hgrid(Mobj2, [0 1])
%% END



















