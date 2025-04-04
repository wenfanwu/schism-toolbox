%% This program is an EXAMPLE (Exp4_ChesBay) in the schism-toolbox
%% Tested platform: Matlab 2024a (Windows)
%% Matlab add-on: Image Processing Toolbox; Mapping Toolbox
%% Author: Wenfan Wu, CCRM, Virginia Institute of Marine Science. 2024
% ==========================================================
% ====== mixed triangular/quadrangular grid for the Chesapeake Bay =====
% ==== hgrid.ll/vgrid.in is obtained at schism_verification_tests==========
% ==========================================================
%% Load the mesh grid
clc;clearvars
Mobj.time = (datetime(2020,6,1):hours(1):datetime(2020,6,10))'; 
Mobj.rundays = days(Mobj.time(end)-Mobj.time(1)); 
Mobj.dt = 150;           
Mobj.coord = 'geographic';          

hgrid_file = 'Exp4_ChesBay\inputs\hgrid.ll';  % NEED TO BE CHANGED
vgrid_file = 'Exp4_ChesBay\inputs\vgrid.in';     % NEED TO BE CHANGED

Mobj = read_schism_hgrid(Mobj, hgrid_file);
Mobj = read_schism_vgrid(Mobj, vgrid_file, 'v5.10'); 

%% Check the grid quality
check_schism_CFL(Mobj);
check_schism_hydrostatic(Mobj);

%% Horizontal grids
figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 1])
axis image
colormap(gray(25))
caxis([-60 0]) %#ok<*CAXIS>
title('Mixed triangular/quad grid for the Chesapeake Bay')

%% Vertical grids
% draw a line on the map and you'll get the vertical layers on the transect
figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 0], 'EdgeAlpha', 0.05, 'LineWidth', 0.5);
auto_center
sect_info = def_schism_transect(Mobj, -1, 0.01);

disp_schism_vgrid(Mobj, sect_info)

%% END