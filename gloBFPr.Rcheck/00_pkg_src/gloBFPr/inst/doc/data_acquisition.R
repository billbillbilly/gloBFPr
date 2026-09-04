## ----eval=FALSE---------------------------------------------------------------
# buildings_list <- gloBFPr::search_3dglobdf(bbox = c(-83.065644,42.333792,-83.045217,42.346988),
#                                            out_type = "all",
#                                            cell_size = 1)

## ----eval=FALSE---------------------------------------------------------------
# buildings_rast <- gloBFPr::search_3dglobdf(bbox = c(-83.065644,42.333792,-83.045217,42.346988),
#                                            out_type = "graduated_rast",
#                                            mask = TRUE,
#                                            cell_size = 1)

## ----eval=FALSE---------------------------------------------------------------
# buildings_gba <- gloBFPr::search_3dglobdf(bbox = c(-83.065644,42.333792,-83.045217,42.346988),
#                                            out_type = "poly",
#                                            data_source = "GBA")

## ----eval=FALSE---------------------------------------------------------------
# dsm <- gloBFPr::get_fused_dsm(
#   x = buildings_list$poly,
#   datasource_canopy_height = "metachm",
#   min_tree_height = 2,
#   opentopo_key = "YOUR_OPENTOPOGRAPHY_API_KEY",
#   quiet = FALSE
# )
# 
# terra::plot(dsm)

## ----eval=FALSE---------------------------------------------------------------
# dsm_no_canopy <- gloBFPr::get_fused_dsm(
#   x = buildings_list$poly,
#   datasource_canopy_height = NULL,
#   opentopo_key = "YOUR_OPENTOPOGRAPHY_API_KEY"
# )

## ----eval=FALSE---------------------------------------------------------------
# dsm_1m <- gloBFPr::get_fused_dsm(
#   x = buildings_list$poly,
#   datasource_canopy_height = "metachm",
#   resolution = 1,
#   opentopo_key = "YOUR_OPENTOPOGRAPHY_API_KEY"
# )

## ----eval=FALSE---------------------------------------------------------------
# dsm_with_canopy <- gloBFPr::get_fused_dsm(
#   x = buildings_list$poly,
#   datasource_canopy_height = "metachm",
#   min_tree_height = 1,
#   opentopo_key = "YOUR_OPENTOPOGRAPHY_API_KEY",
#   resolution = 1
# )

## ----setup, message=FALSE-----------------------------------------------------
library(gloBFPr)
library(sf)
library(terra)

data(globfp_example)
buildings <- globfp_example

## ----quickstart---------------------------------------------------------------
out_dir <- file.path(tempdir(), "world_flat")
world <- get_3d_world(
  x       = buildings,
  terrain = FALSE,
  canopy  = NULL,
  out_dir = out_dir
)

list.files(out_dir)
world$n_buildings

## ----full, eval=FALSE---------------------------------------------------------
# data(globfp_example_dem)
# data(globfp_example_canopy_height)
# 
# world <- get_3d_world(
#   x              = buildings,
#   terrain        = TRUE,
#   dem            = rast(globfp_example_dem),
#   canopy_height  = rast(globfp_example_canopy_height),
#   canopy         = NULL,
#   roads          = "overture",
#   water          = "overture",
#   greenspace     = TRUE,
#   basemap        = TRUE,
#   facade_palette = TRUE,
#   out_dir        = file.path(tempdir(), "world_textured"),
#   quiet          = FALSE
# )

## ----voxel, eval=FALSE--------------------------------------------------------
# world_vox <- get_3d_world(
#   x              = buildings,
#   terrain        = TRUE,
#   dem            = rast(globfp_example_dem),
#   canopy_height  = rast(globfp_example_canopy_height),
#   canopy         = NULL,
#   roads          = "overture",
#   water          = "overture",
#   greenspace     = TRUE,
#   facade_palette = TRUE,
#   all_vox        = TRUE,
#   vox_size       = 1,
#   out_dir        = file.path(tempdir(), "world_vox")
# )

