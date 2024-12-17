# Changelog

All notable changes to the project will be documented in this file. 

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
- **.\utilities\misc\find_land_island**: find land and island nodes for unstructured grid.
- **.\utilities\misc\schism_datatips**: enhance datatips when visualizing data.
- **.\utilities\extract\get_era5_forcing**: extract atmospheric forcing from era5
- **.\utilities\calculate\calc_schism_nudge**: calculate boundary nudging factors (TEM/SAL/ICM_nudge.gr3)
- **.\utilities\calculate\calc_schism_reso**: calculate grid resolution with dfferent methods
- **.\utilities\calculate\calc_schism_edge**: calculate edge length info
- **.\utilities\calculate\calc_schism_shapiro**: renamed from *gen_slope_filter2*
- **.\utilities\calculate\calc_schism_angles**: calculate interior angles of each cell
- **.\utilities\calculate\calc_schism_skew**: calculate the skewness of each cell
- **.\utilities\calculate\calc_schism_grad**: calculate the gradient of variable
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
