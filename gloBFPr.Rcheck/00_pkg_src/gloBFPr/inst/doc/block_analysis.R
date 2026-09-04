## ----setup, message=FALSE-----------------------------------------------------
library(gloBFPr)
library(sf)
library(dplyr)

## -----------------------------------------------------------------------------
data(globfp_example)
buildings <- globfp_example

## ----eval=FALSE---------------------------------------------------------------
# block_result <- generate_block(buildings, quiet = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# blocks    <- block_result$blocks
# buildings <- block_result$buildings
# 
# # How many blocks were found?
# nrow(blocks)
# 
# # How many buildings were assigned?
# sum(!is.na(buildings$block_id))

## ----eval=FALSE---------------------------------------------------------------
# library(osmdata)
# net <- opq(bbox = sf::st_bbox(buildings)) |>
#   add_osm_feature("highway",
#     value = c("motorway", "trunk", "primary", "secondary",
#               "tertiary", "residential", "unclassified")) |>
#   osmdata_sf()
# net_lines <- net$osm_lines
# 
# block_result <- generate_block(buildings, network = net_lines, quiet = FALSE)
# 
# plot(block_result$blocks)

## ----eval=FALSE---------------------------------------------------------------
# # Include primary roads in dual-carriageway simplification
# block_result <- generate_block(
#   buildings,
#   dc_highway_types    = c("motorway", "trunk", "primary"),
#   dc_overlap_threshold = 0.7,
#   quiet = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# block_result <- generate_block(
#   buildings,
#   dc_highway_types = character(0),
#   quiet = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# block_result <- generate_block(
#   buildings,
#   min_block_area = 1000,
#   quiet = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Compute building metrics, then generate blocks
# buildings_with_metrics <- buildings |>
#   get_morphology(quiet = TRUE)
# 
# block_result <- generate_block(buildings_with_metrics, quiet = FALSE)
# 
# # Aggregate to block level
# block_metrics <- aggregate_block(block_result, quiet = FALSE)
# 
# block_metrics |>
#   st_drop_geometry() |>
#   select(block_id, n_buildings, coverage_ratio, g_area, vol) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# # Use maximum volume instead of sum; use median height instead of mean
# block_metrics_custom <- aggregate_block(
#   block_result,
#   .fns = list(vol = max, Height = median),
#   quiet = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# block_metrics <- aggregate_block(
#   block_result,
#   population      = TRUE,
#   population_year = 2020,
#   quiet = FALSE
# )
# 
# block_metrics |>
#   st_drop_geometry() |>
#   select(block_id, n_buildings, pop_total) |>
#   arrange(desc(pop_total)) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# block_metrics <- aggregate_block(
#   block_result,
#   residential      = TRUE,
#   residential_year = 2020,
#   quiet = FALSE
# )
# 
# block_metrics |>
#   st_drop_geometry() |>
#   select(block_id, n_buildings, res_prop) |>
#   arrange(desc(res_prop)) |>
#   head()

## ----eval=FALSE---------------------------------------------------------------
# library(ggplot2)
# 
# ggplot(block_metrics) +
#   geom_sf(aes(fill = coverage_ratio), color = "white", linewidth = 0.15) +
#   scale_fill_viridis_c(
#     option = "magma",
#     labels = function(x) paste0(round(100 * x), "%"),
#     na.value = "grey90"
#   ) +
#   coord_sf(datum = NA) +
#   labs(
#     fill = "Coverage",
#     title = "Building coverage ratio by block"
#   ) +
#   theme_minimal()

## ----eval=FALSE---------------------------------------------------------------
# ggplot(block_metrics) +
#   geom_sf(aes(fill = n_buildings), color = "grey85", linewidth = 0.1) +
#   scale_fill_viridis_c(
#     option = "plasma",
#     trans = "sqrt",
#     na.value = "grey90"
#   ) +
#   coord_sf(datum = NA) +
#   labs(
#     fill = "Buildings",
#     title = "Number of buildings per block"
#   ) +
#   theme_minimal()
# 
# ggplot(block_metrics) +
#   geom_sf(aes(fill = vol), color = NA) +
#   scale_fill_viridis_c(
#     option = "cividis",
#     trans = "log10",
#     labels = function(x) format(x, scientific = TRUE),
#     na.value = "grey90"
#   ) +
#   coord_sf(datum = NA) +
#   labs(
#     fill = "Volume",
#     title = "Total building volume by block"
#   ) +
#   theme_minimal()

## ----eval=FALSE---------------------------------------------------------------
# ggplot(block_metrics) +
#   geom_sf(aes(fill = pop_total), color = "white", linewidth = 0.1) +
#   scale_fill_viridis_c(
#     option = "inferno",
#     trans = "sqrt",
#     na.value = "grey90"
#   ) +
#   coord_sf(datum = NA) +
#   labs(
#     fill = "Population",
#     title = "Estimated population by block"
#   ) +
#   theme_minimal()
# 
# ggplot(block_metrics) +
#   geom_sf(aes(fill = res_prop), color = "white", linewidth = 0.1) +
#   scale_fill_viridis_c(
#     option = "viridis",
#     limits = c(0, 1),
#     labels = function(x) paste0(round(100 * x), "%"),
#     na.value = "grey90"
#   ) +
#   coord_sf(datum = NA) +
#   labs(
#     fill = "Residential",
#     title = "Residential built-up surface proportion by block"
#   ) +
#   theme_minimal()

## ----eval=FALSE---------------------------------------------------------------
# buildings_for_map <- block_result$buildings |>
#   left_join(
#     block_metrics |>
#       st_drop_geometry() |>
#       select(block_id, block_coverage_ratio = coverage_ratio),
#     by = "block_id"
#   )
# 
# ggplot() +
#   geom_sf(data = block_metrics, fill = "grey95", color = "white",
#           linewidth = 0.15) +
#   geom_sf(data = buildings_for_map,
#           aes(fill = block_coverage_ratio),
#           color = "grey25", linewidth = 0.05) +
#   scale_fill_viridis_c(
#     option = "magma",
#     labels = function(x) paste0(round(100 * x), "%"),
#     na.value = "grey90"
#   ) +
#   coord_sf(datum = NA) +
#   labs(
#     fill = "Block coverage",
#     title = "Buildings colored by their block-level coverage ratio"
#   ) +
#   theme_minimal()

## ----eval=FALSE---------------------------------------------------------------
# metric_labels <- c(
#   coverage_ratio = "Coverage ratio",
#   n_buildings = "Building count",
#   vol = "Building volume"
# )
# 
# metrics_long <- bind_rows(lapply(names(metric_labels), function(metric_col) {
#   out <- block_metrics[, c("block_id", metric_col, "geometry")]
#   names(out)[names(out) == metric_col] <- "value"
#   out$metric <- metric_labels[[metric_col]]
#   out
# }))
# 
# ggplot(metrics_long) +
#   geom_sf(aes(fill = value), color = "white", linewidth = 0.1) +
#   scale_fill_viridis_c(
#     option = "mako",
#     trans = "sqrt",
#     na.value = "grey90"
#   ) +
#   facet_wrap(vars(metric), nrow = 1) +
#   coord_sf(datum = NA) +
#   labs(fill = "Value", title = "Block-level metric comparison") +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     strip.text = element_text(face = "bold")
#   )

## ----eval=FALSE---------------------------------------------------------------
# bbox <- c(-83.065644, 42.333792, -83.045217, 42.346988)
# 
# # 1. Retrieve building footprints
# buildings <- search_3dglobdf(bbox = bbox, out_type = "poly", quiet = TRUE)
# 
# # 2. Compute building-level metrics
# buildings <- buildings |>
#   get_morphology(quiet = TRUE)
# 
# # 3. Generate blocks
# block_result <- generate_block(buildings, quiet = FALSE)
# 
# # 4. Aggregate metrics to block level, including population and residential
# block_metrics <- aggregate_block(
#   block_result,
#   population       = TRUE,
#   population_year  = 2025,
#   residential      = TRUE,
#   residential_year = 2025,
#   quiet = FALSE
# )
# 
# # 5. Inspect
# block_metrics |>
#   st_drop_geometry() |>
#   select(block_id, n_buildings, coverage_ratio, g_area, vol,
#          pop_total, res_prop) |>
#   arrange(desc(n_buildings)) |>
#   head(10)

