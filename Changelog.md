# Changelog

All notable changes to the project will be documented in this file. 

## [v1.2-beta] - 2025-04-23

### New features!

- Use structure arrays in IC/BC functions to allow variable-specific open boundaries and time series.
- Add functions for the serial/parallel extraction of open boundary or nudging data.
- Improve the efficiency of many file writing functions.
- Add SAL-related functions (self-attracting and loading tide).

### Added

- **.\utilities\general\merge_structs**: replace *add_structs*
- **.\utilities\extract\get_hycom_nudge**: extract nudging data from the database (Mat)
- **.\utilities\extract\get_hycom_bdry**: extract boundary inputs from the database (Mat)
- **.\utilities\extract\get_hycom_bdry_nc**: extract boundary inputs from the database (NetCDF)
- **.\utilities\extract\get_fes2014_SAL**: extract SAL data from FES2014
- **.\utilities\mesh\interactive_mesh_edit**: edit the unstructured grid manually.
- **.\utilities\calculate\calc_schism_contour**: calculate contour lines.
- **.\utilities\calculate\calc_schism_grad**: calculate the gradient of variable
- **.\utilities\general\sub_region**: extract the index of sub-region.
- **.\utilities\mesh\def_schism_obc**: define open boundaries for the grid.
- **.\preproc\write_schism_SAL**: write load tide files.
- **.\preproc\write_schism_HA**: write the harm.in file.
- **.\preproc\write_schism_nml**: write the namelist file.
- **.\preproc\check\check_schism_CFL**: check the inverse CFL constraints.
- **.\preproc\runoff\add_river_inputs**: add river inputs for source elements.
- **.\postproc\read_schism_bp**: read the bp file.
- **.\postproc\read_schism_nml**: read the namelist file.
- **.\postproc\read_schism_transect**: read variable on a given transect.
- **.\CONTENT**: show the list of all functions.

### Changed

- **.\preproc\prep_schism_init**: updated.
- **.\preproc\prep_schism_bdry**: updated.
- **.\preproc\write_schism_hotstart**: updated.
- **.\preproc\write_schism_gr3**: adopt faster writing methods.
- **.\preproc\write_schism_th**: adopt faster writing methods.
- **.\preproc\write_schism_th_nc**: adopt faster writing methods.
- **.\preproc\write_schism_nu_nc**: keep consistent with *write_schism_th_nc*
- **.\preproc\write_schism_prop**: adopt faster writing methods.
- **.\preproc\write_schism_ic**: support both hvar and vvar-type ic files.
- **.\preproc\write_schism_sflux**: use time_steps rather than the # of files.
- **.\preproc\check\check_schism_icbc**: add a time series plot.
- **.\postproc\read_schism_ic**: support both hvar- and vvar-type ic files.
- **.\utilities\general\auto_center**: rewritten.
- **.\utilities\general\any2time**: updated.
- **.\utilities\general\geomin**: updated.
- **.\utilities\extract\get_hycom_online**: support the latest hycom product (ESPC-D-V02).
- **.\utilities\interp\interp_tri**: improved.
- **.\utilities\interp\interp_deps**: improved.
- **.\utilities\interp\interp_schism_init**: improved.
- **.\utilities\interp\interp_schism_bdry**: improved.
- **.\utilities\define\def_schism_mask**: can specify grid centers.

### Removed

- **.\utilities\general\add_structs**: replaced by *merge_structs*
- **.\utilities\extract\get_hycom_bdry_clim**: redundant functions
- **.\utilities\extract\get_hycom_init_clim**: redundant functions
- **.\utilities\extract\get_hycom_bdry_bys**: redundant functions
- **.\utilities\calculate\calc_schism_CFL**: replaced by *check_schism_CFL*
- **.\preproc\check\check_schism_metrics**: replaced by *check_schism_CFL*
- **.\preproc\runoff\add_river_runoff**: replaced by *add_river_inputs*
- **.\preproc\runoff\add_river_tracer**: replaced by *add_river_inputs*

## [v1.1-beta] - 2024-11-06

### New features!

- work for mixed triangular/quadrangular grid
- work for cartesian or geographic coordinates
- map projection is allowed for variable visualization
- nudging functions are provided.

### Added

- **.\utilities\misc\check_lons**: check the consistency of longitude coordinates
- **.\utilities\misc\add_grid_metrics**: add grid geometry metrics
- **.\utilities\misc\add_bnd_metrics**: add boundary geometry metrics
- **.\utilities\misc\add_schism_obc**: add new open boundaries
- **.\utilities\misc\find_land_island**: find land and island nodes for unstructured grid.
- **.\utilities\misc\schism_datatips**: enhance datatips when visualizing data.
- **.\utilities\extract\get_era5_forcing**: extract atmospheric forcing from era5
- **.\utilities\extract\get_hycom_init**: extract initial data from hycom
- **.\utilities\calculate\calc_schism_nudge**: calculate boundary nudging factors (TEM/SAL/ICM_nudge.gr3)
- **.\utilities\calculate\calc_schism_reso**: calculate grid resolution with dfferent methods
- **.\utilities\calculate\calc_schism_edge**: calculate edge length info
- **.\utilities\calculate\calc_schism_shapiro**: renamed from *gen_slope_filter2*
- **.\utilities\calculate\calc_schism_angles**: calculate interior angles of each cell
- **.\utilities\calculate\calc_schism_skew**: calculate the skewness of each cell
- **.\utilities\interp\convert_schism_var**: convert variables between different mesh centers
- **.\preproc\write_schism_nu_nc**: write nudging files (TEM/SAL/ICM_nu.nc)
- **.\postproc\read_schism_gr3**: read "depth" info from the *gr3 file
- **.\postproc\read_schism_prop**: read the *prop file
- **.\examples\Exp4_ChesBay**: new example with mixed triangular/quad grid
- **.\mesh_object.png**: shows details of the mesh object (*Mobj*)

### Fixed

- **.\preproc\write_schism_sflux**: minor fix on potential time lag

### Changed

- **.\preproc\mesh2schism**: work for 2dm files created by SMS
- **.\preproc\write_schism_hotstart**: updated 
- **.\preproc\write_schism_ic**: updated 
- **.\preproc\write_schism_gr3**: support mixed triangular/quad grid
- **.\preproc\write_schism_hgrid**: supoort cartesian or geographic coordinates
- **.\preproc\write_schism_th_nc**: minor improvements
- **.\preproc\call_schism_tracers**: minor changes
- **.\preproc\check\check_schism_hydrostatic**: updated 
- **.\preproc\check\check_schism_icbc**: updated 
- **.\preproc\check\check_schism_metrics**: updated 
- **.\preproc\runoff\match_rivers**: location change
- **.\postproc\read_schism_vgrid**: updated 
- **.\postproc\read_schism_hgrid**: support mixed triangular/quad grid
- **.\utilities\define\def_schism_source**: updated 
- **.\utilities\define\def_schism_transect**: updated 
- **.\utilities\general\minfind**: updated 
- **.\utilities\general\tripcolor**: updated 
- **.\utilities\general\geomin**: more rigorous searching algorithms
- **.\utilities\general\auto_center**: minor improvements
- **.\utilities\generate\gen_schism_LSC2**: updated 
- **.\utilities\calculate\calc_schism_area**: support mixed triangular/quad grid
- **.\utilities\calculate\calc_schism_bfric**: minor changes
- **.\utilities\calculate\calc_schism_hdif**: minor changes
- **.\utilities\visualize\disp_schism_var**: support mixed triangular/quad grid
- **.\utilities\visualize\plot_schism_bnds**: allow map projection
- **.\utilities\extract\get_fes2014_tide**: add coordinate system checking
- **.\utilities\extract\get_hycom_online**: minor improvements
- **.\examples\Exp1_BYS**: minor changes
- **.\examples\Exp2_BYS_CoSiNE**: minor changes
- **.\examples\Exp3_CORIE_LSC2**: minor changes

### Removed

- **.\utilities\generate\gen_slope_filter2**: renamed to *calc_schism_shapiro*
- **.\utilities\calculate\calc_schism_cradius**: replaced by *calc_schism_reso*
- **.\utilities\interp\interp_schism_elem2node**: replaced by *convert_schism_var*
- **.\utilities\interp\interp_schism_node2elem**: replaced by *convert_schism_var*
- **.\utilities\calculate\calc_schism_sidelen**: replaced by *calc_schism_edge*
