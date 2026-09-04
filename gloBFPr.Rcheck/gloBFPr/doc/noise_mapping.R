## ----setup, message=FALSE, eval=FALSE-----------------------------------------
# library(gloBFPr)
# library(sf)
# library(terra)

## ----eval=FALSE---------------------------------------------------------------
# data(globfp_example)
# data(globfp_example_dem)
# data(globfp_example_canopy_height)
# 
# buildings <- globfp_example
# dem <- rast(globfp_example_dem)
# canopy_height <- rast(globfp_example_canopy_height)
# 
# names(buildings)

## ----eval=FALSE---------------------------------------------------------------
# noise_inputs <- prepare_noisemodelling_inputs(
#   x = buildings,
#   height_field = "Height",
#   datasource_greenspace = "esri",
#   greenspace_zoom = 14,
#   canopy_height = canopy_height,
#   dem = dem,
#   receiver = "grid",
#   resolution = 25,
#   quiet = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# noise_inputs <- prepare_noisemodelling_inputs(
#   x = buildings,
#   height_field = "Height",
#   canopy_height = canopy_height,
#   dem = dem,
#   receiver = "grid",
#   resolution = 25,
#   quiet = TRUE
# )
# 
# names(noise_inputs)
# nrow(noise_inputs$receivers)

## ----eval=FALSE---------------------------------------------------------------
# noise_inputs <- prepare_noisemodelling_inputs(
#   x = buildings,
#   roads = roads,
#   canopy_height = canopy_height,
#   dem = dem,
#   out_dir = tempdir(),
#   write = TRUE,
#   quiet = TRUE
# )
# 
# noise_inputs$gpkg

## ----eval=FALSE---------------------------------------------------------------
# noise_inputs <- prepare_noisemodelling_inputs(
#   x = buildings,
#   height_field = "Height",
#   min_tree_height = 2,
#   datasource_canopy_height = "metachm",
#   datasource_greenspace = "esri",
#   opentopo_key = Sys.getenv("OPENTOPOGRAPHY_KEY"),
#   receiver = "grid",
#   resolution = 25,
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# install_noisemodelling(version = "5.0.1")

## ----eval=FALSE---------------------------------------------------------------
# noise_result <- get_noise_map(
#   x = buildings,
#   height_field = "Height",
#   datasource_canopy_height = "metachm",
#   datasource_greenspace = "esri",
#   dem = dem,
#   receiver = "grid",
#   resolution = 25,
#   run = TRUE,
#   keep_files = TRUE,
#   quiet = FALSE,
#   java = 17
# )
# plot_noise_map(noise_result, period = "DEN", scalebar = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# plot_noise_map(noise_result, period = "DEN")

## ----eval=FALSE---------------------------------------------------------------
# noise_result <- get_noise_map(
#   x = buildings,
#   height_field = "Height",
#   canopy_height = canopy_height,
#   dem = dem,
#   receiver = "grid",
#   resolution = 25,
#   run = TRUE,
#   java = 17,
#   reflection_order = 1,
#   max_src_distance = 500,
#   max_reflection_distance = 350,
#   diffraction_horizontal = TRUE,
#   diffraction_vertical = FALSE,
#   wall_alpha = 0.1,
#   humidity = 75,
#   temperature = 31,
#   favourable_occurrences = rep(0.5, 16),
#   max_error = 0.1,
#   export_source_id = FALSE,
#   frequency_field_prepend = "HZ"
# )

## ----eval=FALSE---------------------------------------------------------------
# plot_noise_map(noise_result, period = "DEN", scalebar = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# install.packages("osmextract")
# 
# osm_file <- osmextract::oe_get(
#   place = "Detroit, Michigan",
#   provider = "geofabrik",
#   download_directory = tempdir(),
#   force_download = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# osm_file <- file.path(tempdir(), "michigan-latest.osm.pbf")
# 
# utils::download.file(
#   "https://download.geofabrik.de/north-america/us/michigan-latest.osm.pbf",
#   osm_file,
#   mode = "wb"
# )

## ----eval=FALSE---------------------------------------------------------------
# noise_result <- get_noise_map(
#   x = buildings,
#   height_field = "Height",
#   canopy_height = canopy_height,
#   dem = dem,
#   osm_file = osm_file,
#   receiver = "grid",
#   resolution = 25,
#   run = TRUE,
#   java = 17
# )

