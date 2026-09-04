## ----setup, message=FALSE-----------------------------------------------------
library(gloBFPr)
library(sf)
library(dplyr)

## -----------------------------------------------------------------------------
data(globfp_example)

buildings <- globfp_example
buildings <- buildings[seq_len(min(10, nrow(buildings))), ]

names(buildings)

## ----eval=FALSE---------------------------------------------------------------
# buildings <- search_3dglobdf(
#   bbox = c(-83.065644, 42.333792, -83.045217, 42.346988),
#   out_type = "poly",
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# morphology <- get_morphology(buildings, quiet = TRUE)
# 
# morphology |>
#   st_drop_geometry() |>
#   select(id, Height, g_area, pmeter, vol, rec, cnv, elo_z) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# basic_shape <- get_morphology(
#   buildings,
#   metrics = c("g_area", "pmeter", "vol", "rec"),
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# neighbors <- get_neighbors(
#   morphology,
#   radius = 500,
#   quiet = TRUE
# )
# 
# neighbors |>
#   st_drop_geometry() |>
#   select(id, n_count, m_ndist, min_ndist, max_ndist, sd_ndist) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# neighbors_100m <- get_neighbors(buildings, radius = 100, quiet = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# dng <- get_dng(
#   neighbors,
#   datasource = "metachm",
#   radius = 800,
#   min_tree_height = 2,
#   min_area = 500,
#   unit = "m2",
#   quiet = TRUE
# )
# 
# dng |>
#   st_drop_geometry() |>
#   select(id, dng, dng_method) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# dng_esri <- get_dng(
#   buildings,
#   datasource = "esri",
#   zoom = 17,
#   radius = 800,
#   min_area = 500,
#   unit = "m2",
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# dng_net <- get_dng(
#   neighbors,
#   datasource = "metachm",
#   radius = 800,
#   min_area = 500,
#   unit = "m2",
#   network = "osm",
#   quiet = TRUE
# )
# 
# dng_net |>
#   st_drop_geometry() |>
#   select(id, dng, dng_method) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# # roads is an sf LINESTRING layer you already have
# dng_custom <- get_dng(
#   neighbors,
#   datasource = "metachm",
#   radius = 800,
#   min_area = 500,
#   unit = "m2",
#   network = roads,
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# table(dng_net$dng_method)
# 
# comparison <- data.frame(
#   id        = dng$id,
#   euclidean = dng$dng,
#   network   = dng_net$dng
# )
# comparison$detour_ratio <- comparison$network / comparison$euclidean
# 
# # Ratios near 1 mean the green space is straight ahead; large ratios flag
# # buildings cut off by barriers, dead ends, or missing crossings.
# head(comparison[order(-comparison$detour_ratio), ])

## ----eval=FALSE---------------------------------------------------------------
# bgvi <- get_bgvi(
#   buildings[1:10, ],
#   datasource_canopy_height = "metachm",
#   datasource_greenspace = "esri",
#   min_tree_height = 2,
#   radius = 800,
#   floor = FALSE,
#   workers = 1,
#   key = Sys.getenv('OPENTOPO_API'),
#   quiet = TRUE
# )
# 
# bgvi |>
#   st_drop_geometry() |>
#   select(id, mean_gvi, bottom_gvi, top_gvi) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# bgvi_2d_green <- get_bgvi(
#   buildings[1:3, ],
#   datasource_canopy_height = NULL,
#   datasource_greenspace = "esri",
#   radius = 800,
#   short_building_threshold = 6,
#   workers = 1,
#   key = Sys.getenv('OPENTOPO_API'),
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# bgvi_directional <- get_bgvi(
#   buildings[1:3, ],
#   datasource_canopy_height = "metachm",
#   datasource_greenspace = "esri",
#   radius = 800,
#   directions = c("southwest", "south", "southeast"),
#   field_of_view = 45,
#   workers = 1,
#   key = Sys.getenv('OPENTOPO_API'),
#   quiet = TRUE
# )
# 
# bgvi_directional |>
#   st_drop_geometry() |>
#   select(id, mean_gvi_south, gvi_bottom_south, gvi_top_south) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# top_direction_view <- plot_bgvi_viewshed(
#   buildings,
#   building = 195,
#   level = "top",
#   orientation = "north",
#   field_of_view = 90,
#   datasource_canopy_height = "metachm",
#   datasource_greenspace = "esri",
#   radius = 500,
#   key = Sys.getenv('OPENTOPO_API'),
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# fifth_floor_view <- plot_bgvi_viewshed(
#   buildings,
#   building = 1,
#   floor = 5,
#   orientation = 135,
#   field_of_view = 90,
#   datasource_canopy_height = "metachm",
#   datasource_greenspace = "esri",
#   radius = 500,
#   key = Sys.getenv('OPENTOPO_API'),
#   quiet = TRUE
# )
# 
# fifth_floor_view$gvi
# fifth_floor_view$visible_green

## ----eval=FALSE---------------------------------------------------------------
# bgvi_by_floor <- get_bgvi(
#   buildings[1:3, ],
#   datasource_canopy_height = "metachm",
#   radius = 800,
#   floor = TRUE,
#   floor_step = 3,
#   workers = 2,
#   key = Sys.getenv('OPENTOPO_API'),
#   quiet = TRUE
# )
# 
# bgvi_by_floor |>
#   st_drop_geometry() |>
#   select(id, estimated_floors, mean_gvi, bottom_gvi, top_gvi, min_gvi, max_gvi, sd_gvi) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# result <- buildings |>
#   get_morphology(quiet = TRUE) |>
#   get_neighbors(radius = 500, quiet = TRUE) |>
#   get_dng(
#     datasource = "metachm",
#     radius = 800,
#     min_area = 500,
#     unit = "m2",
#     network = "osm",
#     quiet = TRUE
#   )
# 
# metric_table <- result |>
#   st_drop_geometry() |>
#   select(
#     id, Height, g_area, vol, rec,
#     n_count, m_ndist, min_ndist, max_ndist, sd_ndist,
#     dng, dng_method
#   )
# 
# head(metric_table)

