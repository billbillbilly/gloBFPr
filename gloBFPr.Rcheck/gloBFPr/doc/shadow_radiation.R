## ----setup, message=FALSE-----------------------------------------------------
library(gloBFPr)
library(sf)
library(terra)

## ----eval=FALSE---------------------------------------------------------------
# data(globfp_example)
# data(globfp_example_dem)
# data(globfp_example_canopy_height)
# 
# buildings <- globfp_example
# # buildings <- buildings[seq_len(min(10, nrow(buildings))), ]
# dem <- rast(globfp_example_dem)
# canopy_height <- rast(globfp_example_canopy_height)

## ----eval=FALSE---------------------------------------------------------------
# solar_time <- "2026-06-21 15:00:00"
# time_zone <- "America/Detroit"
# 
# azimuth <- 135
# elevation <- 35

## ----eval=FALSE---------------------------------------------------------------
# shadow_footprints <- get_shadow_footprint(
#   buildings,
#   solar_time = solar_time,
#   time_zone = time_zone,
#   plot = TRUE,
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# solar_times <- data.frame(
#   sun_id = c("morning", "midday", "afternoon"),
#   solar_time = c(
#     "2026-06-21 09:00:00",
#     "2026-06-21 12:00:00",
#     "2026-06-21 16:00:00"
#   )
# )
# 
# multi_shadow_footprints <- get_shadow_footprint(
#   buildings,
#   solar_time = solar_times$solar_time,
#   time_zone = time_zone,
#   plot = TRUE,
#   plot_overlap_gradient = TRUE,
#   quiet = TRUE
# )
# 
# unique(multi_shadow_footprints$sun_id)

## ----eval=FALSE---------------------------------------------------------------
# combined_shadow_footprints <- get_shadow_footprint(
#   buildings,
#   solar_time = solar_times$solar_time,
#   time_zone = time_zone,
#   overlap_shadow = TRUE,
#   plot = TRUE,
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# template <- rast(
#   xmin = st_bbox(shadow_footprints)[["xmin"]],
#   xmax = st_bbox(shadow_footprints)[["xmax"]],
#   ymin = st_bbox(shadow_footprints)[["ymin"]],
#   ymax = st_bbox(shadow_footprints)[["ymax"]],
#   resolution = 4,
#   crs = st_crs(buildings)$wkt
# )
# 
# shadow_height_surface <- get_shadow_height(
#   buildings,
#   shadow_locations = template,
#   solar_time = solar_time,
#   time_zone = time_zone,
#   quiet = TRUE
# )
# 
# plot(shadow_height_surface, main = "Shadow height")
# plot(st_geometry(buildings), col = NA, border = "black", add = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# shadow_height_with_downloaded_trees <- get_shadow_height(
#   buildings,
#   shadow_locations = template,
#   solar_time = solar_time,
#   time_zone = time_zone,
#   datasource_canopy_height = "metachm",
#   key = "YOUR_OPENTOPOGRAPHY_API_KEY",
#   min_tree_height = 2,
#   quiet = TRUE
# )
# 
# radiation_with_downloaded_trees <- get_radiation(
#   buildings,
#   solar_time = solar_time,
#   time_zone = time_zone,
#   solar_normal = 850,
#   solar_diffuse = 120,
#   datasource_canopy_height = "metachm",
#   key = "YOUR_OPENTOPOGRAPHY_API_KEY",
#   min_tree_height = 2,
#   canopy_transmissivity = 0.15,
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# plot(dem, main = "Sample DEM")
# plot(st_geometry(buildings), col = NA, border = "black", add = TRUE)
# plot(canopy_height, main = "Canopy height")
# plot(st_geometry(buildings), col = NA, border = "black", add = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# shadow_footprints_with_trees <- get_shadow_footprint(
#   buildings,
#   solar_time = solar_time,
#   time_zone = time_zone,
#   canopy_height = canopy_height,
#   dem = dem,
#   min_tree_height = 1.5,
#   quiet = TRUE,
#   plot = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# shadow_height_with_trees <- get_shadow_height(
#   buildings,
#   shadow_locations = template,
#   solar_time = solar_time,
#   time_zone = time_zone,
#   canopy_height = canopy_height,
#   dem = dem,
#   min_tree_height = 2,
#   quiet = TRUE
# )
# 
# plot(shadow_height_with_trees, main = "Building and canopy shadow height")
# plot(st_geometry(buildings), col = NA, border = "black", add = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# radiation <- get_radiation(
#   buildings,
#   solar_time = solar_time,
#   time_zone = time_zone,
#   solar_normal = 850,
#   solar_diffuse = 120,
#   grid_res = 4,
#   plot_3d = TRUE
# )
# 
# head(st_drop_geometry(radiation))

## ----eval=FALSE---------------------------------------------------------------
# radiation_with_trees <- get_radiation(
#   buildings,
#   solar_time = solar_time,
#   time_zone = time_zone,
#   solar_normal = 850,
#   solar_diffuse = 120,
#   canopy_height = canopy_height,
#   dem = dem,
#   min_tree_height = 2,
#   canopy_transmissivity = 0.15,
#   grid_res = 4,
#   plot_3d = TRUE,
#   quiet = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# roof_building_only <- aggregate(
#   total ~ building_id,
#   data = st_drop_geometry(radiation[radiation$surface == "roof", ]),
#   FUN = mean
# )
# roof_with_trees <- aggregate(
#   total ~ building_id,
#   data = st_drop_geometry(radiation_with_trees[radiation_with_trees$surface == "roof", ]),
#   FUN = mean
# )
# 
# roof_compare <- merge(
#   roof_building_only,
#   roof_with_trees,
#   by = "building_id",
#   suffixes = c("_building_only", "_with_trees")
# )
# 
# head(roof_compare)

## ----eval=FALSE---------------------------------------------------------------
# radiation_diff <- radiation
# radiation_diff$direct_diff <- radiation_with_trees$direct - radiation$direct
# radiation_diff$total_diff  <- radiation_with_trees$total  - radiation$total
# 
# # Summary by surface type
# aggregate(
#   cbind(direct_diff, total_diff) ~ surface,
#   data = st_drop_geometry(radiation_diff),
#   FUN = function(x) round(mean(x), 1)
# )

## ----eval=FALSE---------------------------------------------------------------
# diff_pal <- hcl.colors(100, "Blue-Red 3")
# max_abs  <- max(abs(radiation_diff$total_diff), na.rm = TRUE)
# diff_breaks <- cut(radiation_diff$total_diff,
#                    breaks = seq(-max_abs, max_abs, length.out = 101),
#                    include.lowest = TRUE, labels = FALSE)
# 
# plot(st_geometry(buildings), col = "grey95", border = "grey45",
#      main = "Total radiation change due to tree canopy (W/m²)")
# plot(st_geometry(radiation_diff), pch = 16, cex = 0.45,
#      col = diff_pal[diff_breaks], add = TRUE)
# legend("topright",
#        legend = c(paste0("≤ ", -round(max_abs, 0)), "0",
#                   paste0("≥ +", round(max_abs, 0))),
#        pch = 16,
#        col = c(diff_pal[1], diff_pal[50], diff_pal[100]),
#        bty = "n")

## ----eval=FALSE---------------------------------------------------------------
# roof_radiation   <- radiation[radiation$surface == "roof",   ]
# facade_radiation <- radiation[radiation$surface == "facade", ]
# 
# plot(roof_radiation["total"],   pch = 16, cex = 0.8,  key.pos = 4)
# plot(st_geometry(buildings), col = NA, border = "grey30", add = TRUE)
# 
# plot(facade_radiation["total"], pch = 16, cex = 0.45, key.pos = 4)
# plot(st_geometry(buildings), col = NA, border = "grey30", add = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# solar_day <- c(
#   "2026-06-21 08:00:00",
#   "2026-06-21 11:00:00",
#   "2026-06-21 14:00:00",
#   "2026-06-21 17:00:00"
# )
# 
# radiation_day <- get_radiation(
#   buildings,
#   solar_time    = solar_day,
#   time_zone     = time_zone,
#   solar_normal  = c(500, 850, 900, 650),
#   solar_diffuse = c(180, 120, 110, 150),
#   grid_res = 4,
#   plot_3d  = TRUE,
#   quiet    = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# roof_day   <- radiation_day[radiation_day$surface == "roof",   ]
# facade_day <- radiation_day[radiation_day$surface == "facade", ]
# 
# roof_day_breaks <- cut(roof_day$total, breaks = 100,
#                        include.lowest = TRUE, labels = FALSE)
# 
# plot(st_geometry(buildings), col = "grey95", border = "grey45")
# plot(st_geometry(facade_day), pch = 16, cex = 0.25,
#      col = "#2563EB55", add = TRUE)
# plot(st_geometry(roof_day), pch = 16, cex = 0.75,
#      col = roof_cols[roof_day_breaks], add = TRUE)
# legend("topright",
#        legend = c("lower roof total", "higher roof total", "facade samples"),
#        pch = 16,
#        col = c(roof_cols[1], roof_cols[100], "#2563EB55"),
#        bty = "n")

## ----eval=FALSE---------------------------------------------------------------
# radiation_ground <- get_radiation(
#   buildings,
#   solar_time    = c("2026-06-21 08:00:00", "2026-06-21 11:00:00",
#                     "2026-06-21 14:00:00", "2026-06-21 17:00:00"),
#   time_zone     = time_zone,
#   solar_normal  = c(500, 850, 900, 650),
#   solar_diffuse = c(180, 120, 110, 150),
#   grid_res      = 5,
#   ground        = TRUE,
#   ground_res    = 5,
#   plot          = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# ground_rad <- radiation_ground[radiation_ground$surface == "ground", ]
# 
# summary(st_drop_geometry(ground_rad[, c("svf", "direct", "diffuse", "total")]))

## ----eval=FALSE---------------------------------------------------------------
# radiation_ground_trees <- get_radiation(
#   buildings,
#   solar_time    = c("2026-06-21 08:00:00", "2026-06-21 11:00:00",
#                     "2026-06-21 14:00:00", "2026-06-21 17:00:00"),
#   time_zone             = time_zone,
#   solar_normal  = c(500, 850, 900, 650),
#   solar_diffuse = c(180, 120, 110, 150),
#   canopy_height         = canopy_height,
#   dem                   = dem,
#   min_tree_height       = 2,
#   canopy_transmissivity = 0.15,
#   grid_res              = 5,
#   ground                = TRUE,
#   ground_res            = 5,
#   quiet                 = TRUE,
#   plot = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# ground_trees <- radiation_ground_trees[
#   radiation_ground_trees$surface == "ground", ]
# 
# # Mean radiation reduction at ground level due to tree canopy
# mean(ground_trees$total - ground_rad$total, na.rm = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# roof_ground <- radiation_ground[radiation_ground$surface %in% c("roof", "ground"), ]
# 
# aggregate(
#   cbind(direct, diffuse, total) ~ surface,
#   data = st_drop_geometry(roof_ground),
#   FUN  = function(x) round(mean(x, na.rm = TRUE), 1)
# )

