# CONTENT

This file lists all functions in the schism-toolbox.

## postproc

- [`read_schism_bp.m`](../postproc/read_schism_bp.m) - Read the *.bp file.
- [`read_schism_gr3.m`](../postproc/read_schism_gr3.m) - Read the *.gr3 file.
- [`read_schism_hgrid.m`](../postproc/read_schism_hgrid.m) - Read the horizontal grids from hgrid.gr3/hgrid.ll file.
- [`read_schism_ic.m`](../postproc/read_schism_ic.m) - Read the *.ic file (hvar or vvar).
- [`read_schism_nml.m`](../postproc/read_schism_nml.m) - Read the *.nml file for SCHISM.
- [`read_schism_prop.m`](../postproc/read_schism_prop.m) - Read the *.prop file.
- [`read_schism_transect.m`](../postproc/read_schism_transect.m) - Read SCHISM variable along a given transect.
- [`read_schism_vgrid.m`](../postproc/read_schism_vgrid.m) - Read the vertical grids from vgrid.in file.

## preproc

- [`call_schism_tracers.m`](../preproc/call_schism_tracers.m) - Load the info of activated tracers.
- [`mesh2schism.m`](../preproc/mesh2schism.m) - Load the unstructured grid from OceanMesh2D or SMS.
- [`prep_schism_bdry.m`](../preproc/prep_schism_bdry.m) - Prepare the boundary iniputs for SCHISM.
- [`prep_schism_init.m`](../preproc/prep_schism_init.m) - Prepare the initial data for SCHISM.
- [`write_schism_HA.m`](../preproc/write_schism_HA.m) - Write the harm.in file for SCHISM
- [`write_schism_SAL.m`](../preproc/write_schism_SAL.m) - Write loadtide_[FREQ].gr3 files for SCHISM
- [`write_schism_bctides.m`](../preproc/write_schism_bctides.m) - Write the bctides.in file for SCHISM.
- [`write_schism_gr3.m`](../preproc/write_schism_gr3.m) - Write the *.gr3 file for SCHISM
- [`write_schism_hgrid.m`](../preproc/write_schism_hgrid.m) - Write the hgrid.gr3 & hgrid.ll files for SCHISM.
- [`write_schism_hotstart.m`](../preproc/write_schism_hotstart.m) - Write hotstart.nc for SCHISM model.
- [`write_schism_ic.m`](../preproc/write_schism_ic.m) - Write the *ic file for SCHISM (hvar or vvar)
- [`write_schism_nml.m`](../preproc/write_schism_nml.m) - Write*.nml file for SCHISM (no comments).
- [`write_schism_nu_nc.m`](../preproc/write_schism_nu_nc.m) - Write the *_nu.nc files for SCHISM.
- [`write_schism_prop.m`](../preproc/write_schism_prop.m) - Write the *.prop file for SCHISM.
- [`write_schism_ptrack.m`](../preproc/write_schism_ptrack.m) - Write the particle.bp file for SCHISM (Not Completed Now)
- [`write_schism_sflux.m`](../preproc/write_schism_sflux.m) - Write sflux nc files for SCHISM
- [`write_schism_source.m`](../preproc/write_schism_source.m) - Write source_sink.in and relevant files for SCHISM
- [`write_schism_source_nc.m`](../preproc/write_schism_source_nc.m) - Write the source.nc file for SHCISM
- [`write_schism_station_in.m`](../preproc/write_schism_station_in.m) - Write station.in file for SCHISM model
- [`write_schism_th.m`](../preproc/write_schism_th.m) - Write the *.th file for SCHISM
- [`write_schism_th_nc.m`](../preproc/write_schism_th_nc.m) - Write *th.nc files for SCHISM.
- [`write_schism_vgrid.m`](../preproc/write_schism_vgrid.m) - Write the vgrid.in file for SCHISM.

## preproc/check

- [`check_schism_CFL.m`](../preproc/check/check_schism_CFL.m) - Check the inverse CFL constraints for SCHISM grid.
- [`check_schism_bdry.m`](../preproc/check/check_schism_bdry.m) - Check the interpolation at the boundary nodes
- [`check_schism_hydrostatic.m`](../preproc/check/check_schism_hydrostatic.m) - Check the hydrostatic assumption
- [`check_schism_icbc.m`](../preproc/check/check_schism_icbc.m) - Check the consistency between hotstart.nc and *.th.nc files
- [`check_schism_init.m`](../preproc/check/check_schism_init.m) - Check the interpolation of initial fields (surface/bottom)

## preproc/runoff

- [`add_river_inputs.m`](../preproc/runoff/add_river_inputs.m) - Add river runoff at the source elements
- [`match_rivers.m`](../preproc/runoff/match_rivers.m) - Match the source elements with rivers
- [`prep_river_source.m`](../preproc/runoff/prep_river_source.m) - Prepare river source data based on "river_info"

## utilities/calculate

- [`calc_schism_albedo.m`](../utilities/calculate/calc_schism_albedo.m) - Calculate the albedo based on different empirical formulae
- [`calc_schism_angles.m`](../utilities/calculate/calc_schism_angles.m) - Calculate the interior angles (deg) of each cell
- [`calc_schism_area.m`](../utilities/calculate/calc_schism_area.m) - Calculate the element area (m^2)
- [`calc_schism_bfric.m`](../utilities/calculate/calc_schism_bfric.m) - Calculate the bottom friction using different empirical formulae
- [`calc_schism_contour.m`](../utilities/calculate/calc_schism_contour.m) - Calculate contour lines based on SCHISM grid
- [`calc_schism_edge.m`](../utilities/calculate/calc_schism_edge.m) - Calculate side lengths (m) of triangular cells
- [`calc_schism_grad.m`](../utilities/calculate/calc_schism_grad.m) - Calculate the gradient (units m-1) on an unstructured grid.
- [`calc_schism_hdif.m`](../utilities/calculate/calc_schism_hdif.m) - Calculate the horizontal diffusivity
- [`calc_schism_nudge.m`](../utilities/calculate/calc_schism_nudge.m) - Calculate open boundary nudging factors
- [`calc_schism_reso.m`](../utilities/calculate/calc_schism_reso.m) - Calculate horizontal resolution of grid cells
- [`calc_schism_shapiro.m`](../utilities/calculate/calc_schism_shapiro.m) - Calculate the shapiro filter based on bathymetric slope.
- [`calc_schism_skew.m`](../utilities/calculate/calc_schism_skew.m) - Calculate the skewness (0-1) of each cell

## utilities/define

- [`def_schism_fluxflag.m`](../utilities/define/def_schism_fluxflag.m) - Define fluxflags on the SCHISM grid
- [`def_schism_mask.m`](../utilities/define/def_schism_mask.m) - Define a mask on the SCHISM grid
- [`def_schism_ptrack.m`](../utilities/define/def_schism_ptrack.m) - Define the locations for particle tracking.
- [`def_schism_source.m`](../utilities/define/def_schism_source.m) - Define the source/sink elements on the map.
- [`def_schism_transect.m`](../utilities/define/def_schism_transect.m) - Define a transect on the mesh grid

## utilities/extract

- [`get_era5_forcing.m`](../utilities/extract/get_era5_forcing.m) - Extract ERA5 atmopheric forcing
- [`get_fes2014_SAL.m`](../utilities/extract/get_fes2014_SAL.m) - Extract self-attracting and loading tide (SAL) from FES2014
- [`get_fes2014_tide.m`](../utilities/extract/get_fes2014_tide.m) - Extract the FES2014 tidal data
- [`get_hycom_bdry.m`](../utilities/extract/get_hycom_bdry.m) - Parallel extraction of boundary inputs from MAT (v7.3) database.
- [`get_hycom_bdry_nc.m`](../utilities/extract/get_hycom_bdry_nc.m) - Parallel extraction of boundary inputs from NetCDF database.
- [`get_hycom_init.m`](../utilities/extract/get_hycom_init.m) - Extract real-time hycom data as the initial field.
- [`get_hycom_nudge.m`](../utilities/extract/get_hycom_nudge.m) - Parallel extraction of boundary inputs from MAT (v7.3) database.
- [`get_hycom_online.m`](../utilities/extract/get_hycom_online.m) - Download the hycom data in a flexible way

## utilities/general

- [`any2time.m`](../utilities/general/any2time.m) - Convert various time formats to datetime
- [`auto_center.m`](../utilities/general/auto_center.m) - Automatically center the figure on the specified or current mouse screen
- [`cell2array.m`](../utilities/general/cell2array.m) - Convert the numeric vectors inside a cell to an array.
- [`dimnum.m`](../utilities/general/dimnum.m) - Calculate the # of valid dimensions
- [`geomin.m`](../utilities/general/geomin.m) - Find the indices of closet points
- [`merge_structs.m`](../utilities/general/merge_structs.m) - Merges multiple datastructs.
- [`minfind.m`](../utilities/general/minfind.m) - Find the position of the closest value
- [`multi_interp1.m`](../utilities/general/multi_interp1.m) - Interpolate along a specific dimension of an N-D array.
- [`nannum.m`](../utilities/general/nannum.m) - Calculate the # of NaNs
- [`sub_region.m`](../utilities/general/sub_region.m) - Extract the index and coordinates of sub-region.
- [`tripcolor.m`](../utilities/general/tripcolor.m) - Visualize the variable on an unstructured grid.
- [`ut_FUV.m`](../utilities/general/ut_FUV.m) - Compute nodal/satellite correction factors and astronomical argument

## utilities/generate

- [`gen_schism_LSC2.m`](../utilities/generate/gen_schism_LSC2.m) - Generate LSC2-type vertical grids
- [`gen_schism_SZ.m`](../utilities/generate/gen_schism_SZ.m) - Generate SZ-type (Sigma-Z) vertical grids

## utilities/interp

- [`convert_schism_var.m`](../utilities/interp/convert_schism_var.m) - Convert variables between different mesh centeres
- [`interp_deps.m`](../utilities/interp/interp_deps.m) - Interpolate data from z-layers onto sigma-layers
- [`interp_schism_bdry.m`](../utilities/interp/interp_schism_bdry.m) - Interpolate boundary data onto SCHISM vertical layers.
- [`interp_schism_init.m`](../utilities/interp/interp_schism_init.m) - Interpolate initial data onto SCHISM grids
- [`interp_tri.m`](../utilities/interp/interp_tri.m) - Interpolate gridded data onto scattered points

## utilities/mesh

- [`add_bnd_metrics.m`](../utilities/mesh/add_bnd_metrics.m) - Add open/land/island boundary info.
- [`add_grid_metrics.m`](../utilities/mesh/add_grid_metrics.m) - Add grid geometry metrics (work for mixed triangular/quadrangular grid).
- [`def_schism_obc.m`](../utilities/mesh/def_schism_obc.m) - Define open boundary nodes on the SCHISM grid.
- [`find_land_island.m`](../utilities/mesh/find_land_island.m) - Find land and island nodes (work for mixed triangular/quadrangular grid).
- [`find_loop_nodes.m`](../utilities/mesh/find_loop_nodes.m) - Find the loop nodes from an unstructured grid
- [`interactive_mesh_edit.m`](../utilities/mesh/interactive_mesh_edit.m) - Interactive mesh editing tool for unstructured grids

## utilities/misc

- [`add_nodal_factors.m`](../utilities/misc/add_nodal_factors.m) - Add nodal factors for the tidal simulation
- [`check_lons.m`](../utilities/misc/check_lons.m) - Check the consistency of longitude vectors
- [`schism_datatips.m`](../utilities/misc/schism_datatips.m) - Display data cursor position in a data tip
- [`vector_angle.m`](../utilities/misc/vector_angle.m) - Calculate the anti-clockwise angle from vector A to B.

## utilities/visualize

- [`disp_schism_hgrid.m`](../utilities/visualize/disp_schism_hgrid.m) - Visualize horizontal grids.
- [`disp_schism_var.m`](../utilities/visualize/disp_schism_var.m) - Visualize variables on an unstructured grid.
- [`disp_schism_vgrid.m`](../utilities/visualize/disp_schism_vgrid.m) - Visualize the vertical layers of transect
- [`plot_schism_bnds.m`](../utilities/visualize/plot_schism_bnds.m) - Plot the land/island/ocean boundaries for SCHISM

## modules/CoSiNE

- [`prep_cosine_bdry.m`](../modules/CoSiNE/prep_cosine_bdry.m) - Prepare boundary inputs for the CoSiNE module
- [`prep_cosine_init.m`](../modules/CoSiNE/prep_cosine_init.m) - Prepare initial data for the CoSiNE module
- [`write_cosine_ic.m`](../modules/CoSiNE/write_cosine_ic.m) - Write the COS_hvar_*.ic files for the CoSiNE module

## experimental

- [`drag2rough.m`](../experimental/drag2rough.m) - Convert from Cd to roughness [m] in SCHISM model
- [`find_schism_inside.m`](../experimental/find_schism_inside.m) - Return the index of points located within the model domain.
- [`get_bnd_vars.m`](../experimental/get_bnd_vars.m) - Find the variable matrix along the boundary nodes
- [`get_schout_btm.m`](../experimental/get_schout_btm.m) - Get the variable at the near-bottom layer
- [`interp_zcors.m`](../experimental/interp_zcors.m) - Interp SCHISM outputs onto standard z levels
