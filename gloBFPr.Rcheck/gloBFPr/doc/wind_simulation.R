## ----setup, message=FALSE, eval=FALSE-----------------------------------------
# library(gloBFPr)
# library(sf)
# library(terra)
# library(ggplot2)

## ----qs-data, eval=FALSE------------------------------------------------------
# data(globfp_example)
# buildings_list <- list(poly = globfp_example, binary = NULL, graduated = NULL)
# foam_case <- file.path(path.expand("~"), "openfoam_quickstart")
# 
# foam_inputs <- prepare_openfoam_inputs(
#   case_dir            = foam_case,
#   buildings_list      = buildings_list,
#   height_col          = "Height",
#   include_buildings   = TRUE,
#   include_fused_dsm   = FALSE,
#   include_tree_canopy = FALSE,
#   include_morphology  = FALSE,
#   include_neighbors   = FALSE,
#   include_greenspace  = FALSE,
#   overwrite = TRUE, quiet = TRUE
# )
# 
# case_files <- prepare_foam_case(
#   case_dir = foam_inputs$case_dir,
#   stl_file = foam_inputs$files$building_stl,
#   domain   = foam_inputs$domain,
#   inlet_velocity = c(5, 0, 0), z_ref = 10, z0 = 0.1,
#   base_cell_size = 5, sim_hours = 0.5, overwrite = TRUE
# )
# 
# run_openfoam_docker(case_dir = foam_inputs$case_dir,
#                     image = "opencfd/openfoam-run:2506", wait = TRUE)
# 
# maps <- read_foam_pedestrian_slice(case_dir = foam_inputs$case_dir,
#                                    base_cell_size = 5, resolution = 2)
# 
# plot_foam_map(maps, layer = "U_mag",
#               palette = "YlOrRd", legend_title = "Speed (m/s)")

## ----inputs-spatial, eval=FALSE-----------------------------------------------
# foam_inputs <- prepare_openfoam_inputs(
#   case_dir            = "~/openfoam_cases/detroit",
#   bbox                = c(-83.065644, 42.333792, -83.045217, 42.346988),
#   # target_crs is optional — the local UTM zone is chosen automatically
#   include_buildings   = TRUE,
#   include_fused_dsm   = TRUE,           # OpenTopography DEM + buildings + canopy
#   include_tree_canopy = TRUE,           # Meta CHM → canopy drag
#   canopy_source       = "metachm",
#   min_tree_height     = 2,
#   include_morphology  = TRUE,
#   include_neighbors   = TRUE,
#   include_greenspace  = TRUE,           # land cover → z0 raster
#   mask_tree_cover     = TRUE,
#   landcover_source    = "esri",         # "esa" (2020-2021) or "esri" (2017-2025)
#   landcover_year      = 2022,
#   opentopo_key        = Sys.getenv('OPENTOPO_API'),
#   overwrite           = TRUE
# )

## ----geometry, eval=FALSE-----------------------------------------------------
# geo <- prepare_foam_geometry(
#   case_dir        = foam_inputs$case_dir,
#   fused_dsm       = foam_inputs$files$fused_dsm,
#   building_height = foam_inputs$files$building_height_raster,
#   canopy_height   = foam_inputs$files$canopy_height_raster,
#   buildings       = readRDS(foam_inputs$files$buildings_rds),
#   height_col      = "Height",
#   domain          = foam_inputs$domain,
#   base_cell_size  = 10
# )
# # geo$terrain_stl   → constant/triSurface/terrain.stl
# # geo$canopy_stl    → constant/triSurface/canopy.stl
# # geo$building_stl  → buildings.stl, re-based onto the terrain
# # geo$dem           → recovered bare-earth SpatRaster

## ----inputs-z0, eval=FALSE----------------------------------------------------
# z0_rast <- rast(foam_inputs$files$roughness_raster)
# z0_mean <- global(z0_rast, "mean", na.rm = TRUE)[[1]]
# cat(sprintf("Mean z0 = %.4f m\n", z0_mean))

## ----inputs-era5-wind, eval=FALSE---------------------------------------------
# # Daytime — fetches 10-m wind components and 2-m temperature
# met <- get_era5_met(
#   lon      = -83.05,
#   lat      =  42.34,
#   datetime = "2023-07-15 14:00",
#   cds_key  = Sys.getenv("CDS_API_KEY")
# )
# # met$inlet_velocity  →  c(u10, v10, 0) m/s, ready for prepare_foam_case()
# # met$T_ref           →  2-m temperature in K

## ----wind-case, eval=FALSE----------------------------------------------------
# case_files <- prepare_foam_case(
#   case_dir            = foam_inputs$case_dir,
#   # geo$building_stl when you ran Section 2.2, otherwise the flat-ground STL
#   stl_file            = if (exists("geo") && !is.null(geo$building_stl)) {
#                           geo$building_stl
#                         } else {
#                           foam_inputs$files$building_stl
#                         },
#   domain              = foam_inputs$domain,
#   inlet_velocity      = met$inlet_velocity,   # ERA5 u10, v10 — any direction
#   z_ref               = met$z_ref,            # 10 m
#   T_ref               = met$T_ref,            # uniform: buoyancy inactive
#   z0                  = z0_mean,
#   terrain_stl         = if (exists("geo")) geo$terrain_stl else NULL,
#   terrain_dem         = if (exists("geo")) geo$dem else NULL,
#   canopy_stl          = if (exists("geo")) geo$canopy_stl else NULL,
#   base_cell_size      = 15,
#   building_refinement = 2L,
#   sim_hours           = 1,
#   overwrite           = TRUE
# )
# case_files$params$patch_roles
# #> xMin      xMax      yMin      yMax
# #> "inlet"   "outlet"  "lateral" "lateral"

## ----wind-oblique, eval=FALSE-------------------------------------------------
# prepare_foam_case(..., inlet_velocity = c(3, 3, 0))
# #> Wind: 4.24 m/s at 10 m, dir (0.707 0.707)
# #> patches: xMin=inlet xMax=outlet yMin=inlet yMax=outlet

## ----wind-run, eval=FALSE-----------------------------------------------------
# run_openfoam_docker(
#   case_dir = foam_inputs$case_dir,
#   image    = "opencfd/openfoam-run:2506",
#   wait     = TRUE
# )

## ----wind-viz, eval=FALSE-----------------------------------------------------
# maps <- read_foam_pedestrian_slice(
#   case_dir       = foam_inputs$case_dir,
#   T_ref          = case_files$params$T_ref,
#   base_cell_size = case_files$params$base_cell_size,
#   resolution     = 2
# )
# 
# # Wind speed — buildings overlaid automatically
# plot_foam_map(maps, layer = "U_mag",
#               palette = "YlOrRd", legend_title = "Speed (m/s)",
#               title = "Pedestrian wind speed at 1.5 m AGL")
# 
# # Speed with flow vectors
# p <- plot_foam_map(maps, layer = "U_mag",
#                    palette = "Blues", reverse = TRUE,
#                    legend_title = "Speed (m/s)")
# add_flow_vectors(p, maps, spacing = 20, colour = "grey20",
#                  alpha = 0.7, linewidth = 0.25, arrow_size = 0.07)
# 
# # Wind speed ratio (> 1.3 = potentially uncomfortable)
# plot_foam_map(maps, layer = "U_mag", palette = "RdYlGn", reverse = TRUE,
#               legend_title = "U / U_ref",
#               max_u_ref = sqrt(sum(met$inlet_velocity^2)))

