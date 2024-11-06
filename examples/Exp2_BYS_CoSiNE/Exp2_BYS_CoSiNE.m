%% This program is an EXAMPLE (Exp2_BYS_CoSiNE) in the schism-toolbox
%% Tested platform: Matlab 2022a (Windows)
%% Matlab add-on: Image Processing Toolbox; Mapping Toolbox
%% Public package: OceanMesh2D
%% Model verision: SCHISM v5.10
%% Author: Wenfan Wu, CCRM, Virginia Institute of Marine Science. 2024
% ================================================================
% ======== This is a purely hydrological simulation on the lon/lat coordinate =======
% ======== in the Bohai, Yellow Seas (BYS) ===============================
% ================================================================
%% Step-1: Load the mesh grid
% This is a supplement to the first example (Exp1_BYS). It will show you
% how to prepare input files when another tracer module such as CoSiNE is
% activated. You can complete and extend this example to other tracer
% modules in a smiliar vein.
clc;clearvars
workpath = pwd; 
mesh_file = 'Exp2_BYS_CoSiNE\inputs\BYS_20814.mat';  % NEED TO BE CHANGED

Mobj = mesh2schism(mesh_file); 
Mobj.expname = 'Exp2_BYS_CoSiNE';      
Mobj.time = (datetime(2020,6,1):hours(1):datetime(2020,6,10))'; 
Mobj.rundays = days(Mobj.time(end)-Mobj.time(1)); 
Mobj.dt = 150;         
Mobj.coord = 'geographic';         

%% Step-2: Activated modules 
% A total of 15 tracers are considered since the CoSiNE module is activated. 
Mobj.use_cosine = 'yes';
Mobj = call_schism_tracers(Mobj);

%% Step-3: Horizontal grids
figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 0])
axis image
hold on
plot_schism_bnds(Mobj, [1 1 1], 'Color', 'k')

write_schism_hgrid(Mobj)
%% Step-4: Check the inverse CFL constraints and Hydrostatic
check_schism_metrics(Mobj);
calc_schism_CFL(Mobj)
check_schism_hydrostatic(Mobj);

%% Step-5: Vertical grids
% LSC2 coordinates
dep_edges = [10, 20, 30, 45, 55, 65, 75, 90];
dep_nums =  [20 21 22 23 24 25 27 28];
Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_nums, [4 5 3 5], 0.25);

figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 0], 'EdgeAlpha', 0.05, 'LineWidth', 0.5);
auto_center
sect_info = def_schism_transect(Mobj, -1, 0.01);
disp_schism_vgrid(Mobj, sect_info)

write_schism_vgrid(Mobj);
%% Step-6: River Inputs@Element Sources
SS = def_schism_source(Mobj, [1 0], 'load', 'on');
river_info = match_rivers(SS.source.lonc, SS.source.latc, SS.source.elems);

river_info = add_river_runoff(river_info, Mobj.time, 'real_time');

% There are 13 bgc tracers in the CoSiNE module
tracer_list = {'temp', 'salt', 'no3', 'sio4', 'nh4', 's1', 's2', 'z1', 'z2', 'dn', 'dsi', 'po4', 'dox', 'co2', 'alk'};  
river_info = add_river_tracer(river_info, tracer_list, 'real_time');

D = prep_river_source(river_info, tracer_list); 
write_schism_source_nc(Mobj, D,  tracer_list) % wrong tracer order will report errors

% Make sure all the tracer variables are accessible in the
% 'example_river_data.mat' for your selected rivers.

%% Step-7: Initial Conditions (elev/temp/salinity) (time-depedent)
DS1 = prep_schism_init(Mobj, 'hycom_bys'); 
varList = {'ssh', 'temp', 'salt'};
InitCnd1= interp_schism_init(Mobj, DS1, varList);

% Take an example with WOA18 data here, which provides partial bgc data for
% the BYS with coarse resolution.
DS2 = prep_cosine_init(Mobj, 'test_data');   
varList = {'no3', 'sio4', 'po4', 'dox'};  
InitCnd2 = interp_schism_init(Mobj, DS2, varList);

DS = add_structs(DS1, DS2); 
InitCnd = add_structs(InitCnd1, InitCnd2); 

% check the no3 interpolation
check_schism_init(Mobj, DS, InitCnd, 'no3')

% option-1: space-varying but vertically uniform initial field (temp.ic&salt.ic)
write_schism_ic(Mobj, 'elev', InitCnd.ssh)
write_schism_ic(Mobj, 'temp', InitCnd.temp(:,1))
write_schism_ic(Mobj, 'salt', InitCnd.salt(:,1))

write_cosine_ic(Mobj, InitCnd)  % COS_hvar_*.ic

% option-2: 3D initial fields (hotstart.nc)
% start_time = Mobj.time(1);
% hotstart_data = write_schism_hotstart(Mobj, InitCnd, start_time);
%% Step-8: Boundary Conditions (elev/temp/salinity/velocity/module-tracers)
DS1 = prep_schism_bdry(Mobj, 'hycom_bys');
varList = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};
BdryCnd1 = interp_schism_bdry(Mobj, DS1, varList);

% Take an example with a constructed data set here
DS2 = prep_cosine_bdry(Mobj, 'test_data');
varList = {'no3','sio4','nh4','s1'	's2','z1','z2','dn','dsi','po4','dox','co2','alk'}; % DO NOT change the order
BdryCnd2 = interp_schism_bdry(Mobj, DS2, varList);

DS = add_structs(DS1, DS2); 
BdryCnd = add_structs(BdryCnd1, BdryCnd2);

% check the no3 interpolation on the open boundary
check_schism_bdry(Mobj, DS, BdryCnd, 'no3', 1)

write_schism_th_nc(Mobj, 'elev2D', BdryCnd)
write_schism_th_nc(Mobj, 'TEM_3D', BdryCnd)
write_schism_th_nc(Mobj, 'SAL_3D', BdryCnd)
write_schism_th_nc(Mobj, 'uv3D', BdryCnd)

write_schism_th_nc(Mobj, 'COS_3D', BdryCnd)

% check the consistency between initial fields and boundary inputs
% check_schism_icbc(Mobj, 'temp', Mobj.maxLev)

%% END












