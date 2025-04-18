%% This program is an EXAMPLE (Exp1_BYS) in the schism-toolbox
%% Tested platform: Matlab 2024a (Windows)
%% Matlab add-on: Image Processing Toolbox; Mapping Toolbox
%% Public package: OceanMesh2D
%% Model verision: SCHISM v5.10
%% Author: Wenfan Wu, CCRM, Virginia Institute of Marine Science. 2024
% ================================================================
% ====== This is a purely hydrological simulation on geographic coordinate =======
% ================ in the Bohai and Yellow Seas (BYS) ====================
% ================================================================
%% Step-1: Load the mesh grid
% the mesh grid can be generated from OceanMesh2D or SMS. If your grid is
% from other softwares, please use 'read_schism_hgrid.m' and refer to Exp3_CORIE_LSC2.  
clc;clearvars
% option-1: Load mesh grid created by OceanMesh2D
mesh_file = 'Exp1_BYS\inputs\BYS_20814.mat';  % NEED TO BE CHANGED

% option-2: Load mesh grid created by SMS
% mesh_file = 'Exp1_BYS\inputs\BYS_20814.2dm';  % NEED TO BE CHANGED

Mobj = mesh2schism(mesh_file); 
Mobj.expname = 'Exp1_BYS';      
Mobj.time = (datetime(2020,6,1):hours(1):datetime(2020,6,10))';
Mobj.rundays = days(Mobj.time(end)-Mobj.time(1)); 
Mobj.dt = 150; % dt (secs), the same as in param.nml
Mobj.coord = 'geographic'; % geographic or Cartesian coordinates

% All the input files generated afterwards wiil be placed in the directory
% where the meshfile is located (Exp1_BYS\inputs).
%% Step-2: Activated modules 
% only hydrological module is activated in this run and thus there are only two activated tracers (TEM&SAL).
Mobj = call_schism_tracers(Mobj);

%% Step-3: Horizontal grids
% visualize the horizontal grids
figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 0])
axis image
hold on
plot_schism_bnds(Mobj, [1 1 1], 'Color', 'k')

% write the hgrid.gr3 and hgrid.ll files 
write_schism_hgrid(Mobj)
%% Step-4: Check the grid quality
% check the invese CFL constraints
check_schism_CFL(Mobj);

% check the hydrostatic assumption
check_schism_hydrostatic(Mobj);

%% Step-5: Vertical grids
% option-1: LSC2 coordinates
dep_edges = [10, 20, 30, 45, 55, 65, 75, 90];
dep_nums =  [20 21 22 23 24 25 27 28];
Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_nums, [4 5 3 5], 0.25);

% option-2: SZ coordinates
% s_consts = [10, 0.7, 5, 20];
% zcors = 20:2:(fix(max(Mobj.depth))+10);
% Mobj = gen_schism_SZ(Mobj, s_consts, zcors);

% check the quality of vertical grids
% draw a line on the map and press ENTER
figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 0], 'EdgeAlpha', 0.05, 'LineWidth', 0.5);
auto_center
sect_info = def_schism_transect(Mobj, -1, 0.01);

disp_schism_vgrid(Mobj, sect_info) % display the vertical layers on your selected transect

% write the vgrid.in file. Note that the format of vgrid.in has changed
% since v5.10, and thus you need to specify the version number here (v5.10
% or v5.9). v5.10 is default. This option is invalid for SZ coordinate.
write_schism_vgrid(Mobj, 'v5.10');
%% Step-6: River Inputs@Element Sources
% Left-click the points at the center of elements to select river sources
% (avtivate the datatips mode first), and press SHIFT to select multiple
% points simultaneously. The selected river source points will be saved as
% a MAT file named 'source_sink.mat'. 

% load the source_sink.mat existed in the aimpath; change the 'load' to
% 'rebuild' to select your own source/sink elements if needed.
SS = def_schism_source(Mobj, [1 0], 'rebuild', 'on');   % select Yellow River Mouth here
river_info = match_rivers(SS.source.lonc, SS.source.latc, SS.source.elems);

river_info = add_river_runoff(river_info, Mobj.time, 'real_time');

tracer_list = {'temp', 'salt'};
river_info = add_river_tracer(river_info, tracer_list, 'real_time');

D = prep_river_source(river_info, tracer_list);
write_schism_source_nc(Mobj, D,  tracer_list)

% Two things should be done before preparing your own application.
% 1) prepare your own 'example_river_data.mat' file in 'add_river_runoff'
%     and 'add_river_tracer' according to your needs. 
% 2) add corresponding rivers in the 'match_rivers.m' function.

% In this case, the provided 'example_river_data.mat' is just a sample,
% and the tracer data inside it are set to be constant for simplicity.

% If you want to add rivers in the form of open boundaries, please refer to
% add_schism_obc.m function. 

%% Step-7: Initial Conditions (elev/temp/salinity) (time-depedent)
% DS contains raw data in a standardized format:
% 1) 'lon', 'lat', and 'depth' vectors must be in ascending order;
% 2) 'depth' must be positive; ensure the 'lon' and 'lat' ranges fully cover your model domain;
% 3) variable matrix ('var') must have dimensions of either [lon x lat] or [lon x lat x depth].

% prep_schism_init is a simple wrapper function, so you can add more data
% sources in it as needed, just make sure the format of DS complies with
% the requirements above.

% option-1: real-time hycom data
DS = prep_schism_init(Mobj, 'hycom_bys'); 

% option-2: monthly climatology hycom data
% DS = prep_schism_init(Mobj, 'hycom_bys_clim');

% option-3: directly download real-time hycom data from the internet
% DS = prep_schism_init(Mobj, 'hycom_online'); 

varList = {'ssh', 'temp', 'salt'};  % it can be changed if you only want to interpolate for partial variables.
InitCnd = interp_schism_init(Mobj, DS, varList);

% check the interpolation
check_schism_init(Mobj, DS, InitCnd, 'temp')

% option-1: space-varying but vertically uniform initial field (temp.ic&salt.ic)
write_schism_ic(Mobj, 'elev', InitCnd.ssh)
write_schism_ic(Mobj, 'temp', InitCnd.temp(:,1))
write_schism_ic(Mobj, 'salt', InitCnd.salt(:,1))

% option-2: 3D initial fields (hotstart.nc)
start_time = Mobj.time(1);
hst_data = write_schism_hotstart(Mobj, InitCnd, start_time);
%% Step-8: Boundary Conditions (elev/temp/salinity/velocity/module-tracers)
% This step can be time-consuming when using high-resolution, real-time boundary
% inputs, particularly over long time periods.

% prep_schism_bdry is also a wrapper function. It integrates two general-purpose 
% functions for handling HYCOM data: get_hycom_bdry and get_hycom_bdry_nc.
% These functions support both serial and parallel data extraction.
obc_bnds = 1:Mobj.obc_counts;  % extract data for all open boundaries.

% option-1: real-time boundary inputs from hycom.
DS = prep_schism_bdry(Mobj, 'hycom_bys', obc_bnds);  % supoort parallel extraction

% option-2: monthly climatology boundary inputs from hycom.
% DS = prep_schism_bdry(Mobj, 'hycom_bys_clim', obc_bnds);

varList = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; 
BdryCnd = interp_schism_bdry(Mobj, DS, varList);

write_schism_th_nc(Mobj, 'elev2D', BdryCnd)
write_schism_th_nc(Mobj, 'TEM_3D', BdryCnd)
write_schism_th_nc(Mobj, 'SAL_3D', BdryCnd)
write_schism_th_nc(Mobj, 'uv3D', BdryCnd)

% check the interpolation
check_schism_bdry(Mobj, DS, BdryCnd, 'temp', Mobj.time(1))

% check the consistency between initial fields and boundary inputs
check_schism_icbc(Mobj, 'temp', 1)

%% Step-9: Tide Forcing (bctides.in)
% download the fes2014 tidal products first, and change the directory in
% 'get_fes2014_tide.m' (Lines 58–60).
tideList = {'S2','M2','N2','K2', 'K1','P1','O1','Q1'};
obc_bnds = 1:Mobj.obc_counts;
TideForc = get_fes2014_tide(Mobj, tideList, obc_bnds);   

TideForc.cutoff_depth = 10;
TideForc.nf_temp = 0.8;
TideForc.nf_salt = 0.8;

bc_flags = [5 5 4 4];
write_schism_bctides(Mobj, TideForc, bc_flags)

% Add self-attracting and loading tide (optional)
% download fes2014 data and change the directory in 'get_fes2014_SAL,m' (Line 40) 
SAL = get_fes2014_SAL(Mobj, tideList);
write_schism_SAL(Mobj, SAL)
%% Step-10: Bottom friction
% Type-1: roughness
z0 = 0.001;  % set constant roughness in the model domain
write_schism_gr3(Mobj, 'rough', z0)

% Type-2: drag
Cd = calc_schism_bfric(Mobj, 1, [0.07 3], 'on');
write_schism_gr3(Mobj, 'drag', Cd)

% Type-3: manning
fmc = 0.025;
write_schism_gr3(Mobj, 'manning', fmc)

%% Step-11: Misc. files ending in gr3
% shapiro.gr3
shapiro_val = calc_schism_shapiro(Mobj, [0.001, 0.05], 0.5, 'on');
write_schism_gr3(Mobj, 'shapiro', shapiro_val)

% windrot_geo2proj.gr3
write_schism_gr3(Mobj, 'windrot_geo2proj', 0)

% albedo.gr3
albedo_val = calc_schism_albedo(Mobj, 1, 'on');
write_schism_gr3(Mobj, 'albedo', albedo_val)

% watertype.gr3
wtypes = 5;
write_schism_gr3(Mobj, 'watertype', wtypes)

% diffmax.gr3 & diffmin.gr3
write_schism_gr3(Mobj, 'diffmax', 1)
write_schism_gr3(Mobj, 'diffmin', 1.0e-6)

% bdef.gr3
bdef_factor = zeros(size(Mobj.depth));
bdef_factor(Mobj.depth<5) = -2;
write_schism_gr3(Mobj, 'bdef', bdef_factor)

% hdif.gr3
hdif = calc_schism_hdif(Mobj, 0.25, 5e-3, 'on');
write_schism_gr3(Mobj, 'hdif', hdif)

%% Step-12: Misc. files ending in prop
tvd_flags = ones(Mobj.nElems, 1);
tvd_flags(Mobj.depthc<5) = 0;
write_schism_prop(Mobj, 'tvd', tvd_flags)

% define flux_flags if needed
figure('Color', 'w')
disp_schism_hgrid(Mobj, [0 1])
hold on
flux_flags = def_schism_fluxflag(Mobj, 2);
write_schism_prop(Mobj, 'fluxflag', flux_flags)

%% Step-13: Atmospheric forcing
% AtmForc inculdes the extracted atmospheric forcing data with a fixed format:
% 1) lon/lat matrix (nLons*nLats) are generated from the meshgrid function, and in ascending order
% 2) variable matrix should be of nLons*nLats*nTimes
% 3) the 'region' and 'time' fields must cover your simulation period, model domain, respectively.

time_steps = 30; % time steps in each netcdf file

load('example_AtmForc_era5.mat')
AtmForc.aimpath = Mobj.aimpath;

write_schism_sflux(AtmForc, 'prc', time_steps)
write_schism_sflux(AtmForc, 'rad', time_steps)
write_schism_sflux(AtmForc, 'air', time_steps)

% Note the "time_steps" of each sflux nc file can not exceed 1000 in the
% model. In addition, the 'hour' component of 'base_date' attribute is
% unused in each nc file.

% AtmForc was created by the function get_era5_forcing.m. However, you need to
% modify/replace this function carefully based on your own data sources,
% just make sure the obtained AtmForc meets the required format above.    

%% Step-14: Boundary nudging (optional)
% define a boundary nuding zone (90-km width)
% 20 km is the width of max-nudging zone adjacent to the boundary
[nudge_factor, nudge_nodes] = calc_schism_nudge(Mobj, [20, 90, 4e-5], 'all', 'on');

nudge_time = Mobj.time(1):Mobj.time(end); % daily inputs
D.time = seconds(nudge_time-nudge_time(1));
D.map_to_global_node = nudge_nodes;
D.tracer_concentration = 25*ones(1,Mobj.maxLev, numel(nudge_nodes), numel(nudge_time));  % constant temperature

write_schism_gr3(Mobj, 'TEM_nudge', nudge_factor)
write_schism_nu_nc(Mobj, 'TEM', D)
%% END

























































