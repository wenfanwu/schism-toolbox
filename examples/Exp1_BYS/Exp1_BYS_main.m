%% This program is an EXAMPLE (Exp1_BYS) in the schism-toolbox
%% Tested platform: Matlab 2022a (Windows)
%% Matlab add-on: Image Processing Toolbox; Mapping Toolbox
%% Public package: M_Map; OceanMesh2D
%% Model verision: SCHISM v5.10
%% Author: Wenfan Wu, COAS, Ocean Univ. of China. 2023
% ================================================================
% ======== This is a purely hydrological simulation on the lon/lat coordinate =======
% ======== in the Bohai, Yellow Seas (BYS) ===============================
% ================================================================
%% Step-1: Load the mesh grid
% Load mesh grid from the MAT file, in this example, the mesh grid is
% generated from the OceanMesh2D. If your grid is from SMS or other
% softwares, please use 'read_schism_hgrid.m' and refer to Exp3_CORIE.
clc;clearvars
mesh_file = 'E:\Code-repository\Matlab-codes\functions-test\schism-toolbox-v1.0-beta\examples\Exp1_BYS\BYS_20814.mat';  % NEED TO BE CHANGED!!!

Mobj = mesh2schism(mesh_file); 
Mobj.expname = 'Exp1_BYS';      
Mobj.time = (datetime(2020,6,1):hours(1):datetime(2020,6,10))';
Mobj.rundays = days(Mobj.time(end)-Mobj.time(1)); 
Mobj.dt = 150; % dt (secs), the same as in param.nml
Mobj.coord = 'geographic'; % geographic or Cartesian coordinate

% All the input files generated afterwards wiil be placed in the directory
% where the meshfile is located (Exp1_BYS\inputs).
%% Step-2: Activated modules 
% in this case, only hydrological module is activated and thus there are
% only two activated tracers (TEM&SAL).
% This step is reserved for future extension to other modules
Mobj = call_schism_tracers(Mobj);

%% Step-3: Horizontal grids
% visualize the horizontal grids
figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 0])
axis image
hold on
plot_schism_bnds(Mobj, [1 1], 'Color', 'k')

% write the hgrid.gr3 and hgrid.ll files 
write_schism_hgrid(Mobj)
%% Step-4: Check the grid quality
% check the invese CFL constraints
check_schism_metrics(Mobj);

% display the Max. acceptable resolutions as a function of water depth 
calc_schism_CFL(Mobj)

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

% Write the vgrid.in file. Note that the format of vgrid.in has changed
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
SS = def_schism_source(Mobj, [1 0], 'load', 'on');  
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

% River can also be added in the form of open boundaries but it is not
% supported in this toolbox so far.

%% Step-7: Initial Conditions (elev/temp/salinity) (time-depedent)
% DS contains the original initial fields with a fixed format:
% 1) 'lon', 'lat', 'depth' vectors must be in ascending order;
% 2) 'depth' vector should be positive; and ensure the range of lon/latcovers you model domain
% 3) 'var' must in the dimension of lon*lat or lon&lat*depth.

% prep_schism_init is just a packing function, and thus you can easily add
% more options in it according to your needs. Just make sure the format of
% DS meets the requirements above.

% option-1: real-time initial field from hycom.
DS = prep_schism_init(Mobj, 'hycom_bys'); 

% option-2: monthly clim. initial field from hycom.
% DS = prep_schism_init(Mobj, 'hycom_clim'); 

% option-3: directly download the real-time hycom data from the internet
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
% This step can be quite time-consuming if you choose real-time boundary
% inputs with a high time resolution, especially when the model timespan is long.
% option-1: real-time boundary inputs from hycom.
DS = prep_schism_bdry(Mobj, 'hycom_bys'); 

% option-2: monthly clim. boundary inputs from hycom.
% DS = prep_schism_bdry(Mobj, 'hycom_clim');

varList = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; 
BdryCnd = interp_schism_bdry(Mobj, DS, varList);

write_schism_th_nc(Mobj, BdryCnd, 'elev2D')
write_schism_th_nc(Mobj, BdryCnd, 'TEM_3D')
write_schism_th_nc(Mobj, BdryCnd, 'SAL_3D')
write_schism_th_nc(Mobj, BdryCnd, 'uv3D')

% check the interpolation
check_schism_bdry(Mobj, DS, BdryCnd, 'temp', 1)

% check the consistency between initial fields and boundary inputs
check_schism_icbc(Mobj, 'temp', Mobj.maxLev)

%% Step-9: Tide Forcing (bctides.in)
% download the fes2014 tidal products first, and change the directory in
% 'get_fes2014_tide.m' (Lines 40–42).
tideList = {'S2','M2','N2','K2', 'K1','P1','O1','Q1'};
TideForc = get_fes2014_tide(Mobj, tideList);   

% kill the potential NaN values adjacent to the coast
field_list = fieldnames(TideForc);
for ii = 2:numel(field_list)
    tide_var = field_list{ii};
    TideForc.(tide_var) = fillmissing(TideForc.(tide_var), 'previous', 1);
end

TideForc = add_nodal_factors(Mobj, TideForc);
TideForc.cutoff_depth = 10;

bc_flags = [5 5 4 4];
TideForc.nf_temp = 0.8;
TideForc.nf_salt = 0.8;

write_schism_bctides(Mobj, TideForc, bc_flags)
%% Step-10: Bottom friction
% Type-1: roughness
z0 = 0.001;  % set constant roughness in the model domain
write_schism_gr3(Mobj, 'rough', z0)

% Type-2: drag
Cd = calc_schism_bfric(Mobj, 1, [0.07 3], 'on');

write_schism_gr3(Mobj, 'drag', Cd)

% Type-3: manning
fmc = 0.025;
write_schism_gr3(Mobj, 'rough', fmc)

%% Step-11: Misc. files ending in gr3
% shapiro.gr3
shapiro_val = gen_slope_filter2(Mobj, [0.001, 0.05], 0.5, 'on');
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
hdif = calc_schism_hdif(Mobj, 0.25, 5e-5, 'on');
write_schism_gr3(Mobj, 'hdif', hdif)

%% Step-12: Misc. files ending in prop
tvd_flags = ones(Mobj.nElems, 1);
tvd_flags(Mobj.depthc<5) = 0;
write_schism_prop(Mobj, 'tvd', tvd_flags)

% define flux_flags if needed
flux_flags = def_schism_fluxflag(Mobj, 2);
write_schism_prop(Mobj, 'fluxflag', flux_flags)

%% Step-13: Atmospheric forcing
% AtmForc inculdes the extracted atmospheric forcing data with a fixed format:
% 1) lon/lat matrix (nLons*nLats) are generated from the meshgrid function, and in ascending order
% 2) variable matrix should be of nLons*nLats*nTimes
% 3) the 'region' and 'time' fields must cover your simulation period, model domain, respectively.

nFiles = 5; 

load('test_era5_AtmForc.mat')
AtmForc.aimpath = Mobj.aimpath;

write_schism_sflux(AtmForc, 'prc', nFiles)
write_schism_sflux(AtmForc, 'rad', nFiles)
write_schism_sflux(AtmForc, 'air', nFiles)

% nFiles is the estimated # of sflux_prc/air/rad_*.nc files. Note that the
% nFiles can not be two small, since the time_steps of each sflux nc file
% can not exceed 1000 by default. In addition, write_schism_sflux.m will
% slightly adjust this value to ensure that each file starts at 0 o'clock
% in certain day, since the 'hour' component of 'base_date' property is
% unused in each nc file.     

% Considering that the raw data of atmospheric forcing is too large, this
% toolbox did not offer the function to create 'AtmForc', but directly
% uploads the result. 

% You need to prepare the AtmForc variable by yourself according to the
% required format, and then use 'write_schism_sflux.m' to create the nc files. 

%% END

























































