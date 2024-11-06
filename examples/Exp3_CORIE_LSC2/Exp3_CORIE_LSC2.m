%% This program is an EXAMPLE (Exp3_CORIE_LSC2) in the schism-toolbox
%% Tested platform: Matlab 2022a (Windows)
%% Matlab add-on: Image Processing Toolbox; Mapping Toolbox
%% Public package: M_Map; OceanMesh2D
%% Author: Wenfan Wu, COAS, Ocean Univ. of China. 2023
% =============================================================== 
% ========== A purely triangular mesh for the Columbia River Estaury ==========
% ==== hgrid.gr3/vgrid.in is obtained from the schism_verification_test website =====
% =============================================================== 
%% Load the mesh grid (Change current folder first)
% This example will shows you to create the mesh object (Mobj) directly
% from an existed hgrid.gr3 file.
clc;clearvars
Mobj.time = (datetime(2020,6,1):hours(1):datetime(2020,6,10))'; 
Mobj.rundays = days(Mobj.time(end)-Mobj.time(1)); 
Mobj.dt = 150;           % dt (secs), the same as in param.nml
Mobj.coord = 'geographic';             % geographic or Cartesian coordinate

hgrid_file = 'Exp3_CORIE_LSC2\inputs\hgrid.gr3';  % NEED TO BE CHANGED
vgrid_file = 'Exp3_CORIE_LSC2\inputs\vgrid.in';     % NEED TO BE CHANGED

Mobj = read_schism_hgrid(Mobj, hgrid_file);
Mobj = read_schism_vgrid(Mobj, vgrid_file, 'v5.10'); 

%% Horizontal grids
figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 1])
axis image
colormap(gray(25))
title('Locally refined mesh grid for the Columbia River Estuary')

% The resulting figure (CORIE_model_domain.png) can be found in the working
% path. In this figure, blue rectangle: open boundaries; red points: land
% nodes; green points: island nodes; 

% It is seen that there are four segments of open boundaries in this
% hgrid.gr3, including three river boundaries and one ocean boundary. 

%% Vertical grids
% draw a line on the map and you'll get the vertical layers on the transect
figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 0], 'EdgeAlpha', 0.05, 'LineWidth', 0.5);
auto_center
sect_info = def_schism_transect(Mobj, -1, 0.01);

disp_schism_vgrid(Mobj, sect_info)

% the resulting figures (CORIE_transect_location.png and
% CORIE_vertical_layers.png) can be found in the working path. 

%% END