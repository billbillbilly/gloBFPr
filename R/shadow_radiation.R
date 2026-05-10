#' Building shadow and radiation calculations
#'
#' @name get_shadows
#' @param x An `sf` polygon object with building footprints and a height field.
#' @param height_field Character. Name of the building height column. Defaults to
#'   `"Height"`, matching [search_3dglobdf()] output.
#' @param azimuth Numeric vector or list. Solar azimuth in decimal degrees,
#'   measured clockwise from north. Must have the same length as `elevation`.
#' @param elevation Numeric vector or list. Solar elevation in decimal degrees
#'   above the horizon. Must have the same length as `azimuth`.
#' @param solar_time Character vector or list of character strings. Local solar
#'   times such as `"2026-06-21 15:00:00"`. If `solar_time` and `time_zone` are
#'   supplied, `azimuth` and `elevation` are ignored and solar position is
#'   estimated from time and the building-layer centroid.
#' @param time_zone Character. A single time zone used to interpret
#'   `solar_time`, for example `"America/Denver"` or `"UTC"`.
#' @param b Numeric buffer tolerance used when cleaning footprint unions.
#' @param overlap_shadow Logical. For `get_shadow_footprint()`, if `TRUE`,
#'   dissolve overlapping shadows across all supplied solar positions by shadow
#'   source.
#' @param plot Logical. For `get_shadow_footprint()`, draw a base R map of the
#'   building footprints and shadow polygons before returning the `sf` result.
#'   For `get_radiation()`, draw a base R map of radiation sample points colored
#'   by `total`.
#' @param plot_overlap_gradient Logical. For `get_shadow_footprint()` plots with
#'   multiple `solar_time` values, if `TRUE`, draw all shadows in transparent
#'   gray so overlapping areas appear darker.
#' @param shadow_locations Optional query locations for shadow height, as an
#'   `sf` point layer or a `terra::SpatRaster`.
#' @param ... Reserved for compatibility. `location` is accepted as a deprecated
#'   alias of `shadow_locations`; `solar_pos`, `time`, and `combine` are accepted
#'   as deprecated aliases.
#' @param cell_size Numeric cell resolution in CRS units when
#'   `shadow_locations` is omitted.
#' @param extent_buffer Optional numeric buffer around `x` used when creating an
#'   automatic `terra::SpatRaster` template. If omitted, a buffer is estimated
#'   from building heights and solar elevation.
#' @param parallel Ignored. Kept for API compatibility.
#' @param filter_footprint Ignored. Shadow footprints are always used to limit
#'   height calculations.
#' @param min_tree_height Numeric. Minimum canopy height, in map units, used as
#'   a tree obstacle.
#' @param datasource_canopy_height Character or `NULL`. Canopy height source to
#'   retrieve internally when `canopy_height` is not supplied. Currently supports
#'   `"metachm"`, `"ethCHM"`, or `NULL`.
#' @param key Character or `NULL`. OpenTopography API key used to retrieve DEM
#'   internally when `dem` is not supplied.
#' @param raster_buffer Numeric or `NULL`. Buffer distance in CRS units around
#'   buildings used when retrieving CHM/DEM internally. If `NULL`, a buffer is
#'   estimated from building height and solar elevation.
#' @param canopy_transmissivity Numeric from 0 to 1. Fraction of direct
#'   irradiance transmitted through canopy shadows in `get_radiation()`.
#' @param canopy_height Optional `terra::SpatRaster` canopy height map. Values
#'   are interpreted as height above ground.
#' @param dem Optional `terra::SpatRaster` digital elevation model. When
#'   supplied, canopy and building shadows are compared in absolute elevation
#'   and shadow-height outputs are returned above local ground.
#' @param grid Optional 3D `sf` point surface grid. If omitted, it is created
#'   from building roofs and facades.
#' @param grid_res Numeric surface-grid resolution in CRS units.
#' @param offset Numeric vertical offset added to generated surface-grid points.
#' @param solar_normal Direct Normal Irradiance vector, one value per solar
#'   position.
#' @param solar_diffuse Diffuse Horizontal Irradiance vector, one value per solar
#'   position.
#' @param radius Ignored. Kept for API compatibility.
#' @param return_list Logical. If `TRUE`, return per-timestep radiation matrices
#'   instead of a summed `sf` surface grid.
#' @param quiet Logical. If `FALSE`, emit progress messages.
#' @return
#' `get_shadow_footprint()` returns an `sf` polygon layer.
#'
#' `get_shadow_height()` returns a `terra::SpatRaster` for `terra` locations or
#' a numeric matrix for point locations.
#'
#' `get_radiation()` returns an `sf` point layer with `svf`, `direct`,
#' `diffuse`, and `total` columns, unless `return_list = TRUE`.
#'
#' @details
#' These functions are implemented directly with `sf` and `terra` using a
#' projected 2.5D building model.
#'
#' @references
#' Dorman, M. et al. `shadow`: Geometric Shadow Calculations.
#' <https://github.com/michaeldorman/shadow>
NULL

#' @description
#' `get_shadow_footprint()` computes ground shadow footprints for extruded
#' building polygons.
#'
#' @export
#' @rdname get_shadows
get_shadow_footprint <- function(x = NULL,
                                 solar_time = NULL,
                                 time_zone = NULL,
                                 azimuth = NULL,
                                 elevation = NULL,
                                 height_field = "Height",
                                 min_tree_height = 2,
                                 datasource_canopy_height = NULL,
                                 key = NULL,
                                 canopy_height = NULL,
                                 dem = NULL,
                                 raster_buffer = NULL,
                                 b = 0.01,
                                 overlap_shadow = FALSE,
                                 plot = FALSE,
                                 plot_overlap_gradient = FALSE,
                                 quiet = TRUE,
                                 ...) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygons.")
    return(NULL)
  }
  dots <- list(...)
  solar_args <- extract_deprecated_solar_args(dots)
  dots <- solar_args$dots
  if ("combine" %in% names(dots)) {
    if (!identical(overlap_shadow, FALSE)) {
      stop("Use only one of `overlap_shadow` or deprecated `combine`.", call. = FALSE)
    }
    warning(
      "`combine` is deprecated; use `overlap_shadow` instead.",
      call. = FALSE
    )
    overlap_shadow <- dots$combine
  }
  unknown_args <- setdiff(names(dots), "combine")
  if (length(unknown_args) > 0) {
    stop("Unused argument(s): ", paste(unknown_args, collapse = ", "), call. = FALSE)
  }
  buildings <- prepare_shadow_buildings(x, height_field)
  solar_pos <- resolve_solar_inputs(
    buildings,
    azimuth = azimuth,
    elevation = elevation,
    solar_time = solar_time,
    time_zone = time_zone,
    solar_pos = solar_args$solar_pos,
    time = solar_args$time
  )
  sun_ids <- make_shadow_sun_ids(solar_time, nrow(solar_pos))
  raster_inputs <- resolve_shadow_raster_inputs(
    buildings = buildings,
    solar_pos = solar_pos,
    height_field = height_field,
    canopy_height = canopy_height,
    dem = dem,
    datasource_canopy_height = datasource_canopy_height,
    key = key,
    min_tree_height = min_tree_height,
    raster_buffer = raster_buffer,
    quiet = quiet
  )
  canopy <- prepare_canopy_obstacles(
    raster_inputs$canopy_height,
    raster_inputs$dem,
    buildings,
    min_tree_height
  )
  if (!quiet) cli::cli_alert_info("Computing building shadow footprints ...")
  out <- lapply(seq_len(nrow(solar_pos)), function(i) {
    building_shadows <- compute_shadow_footprints(buildings, solar_pos[i, ], height_field, b)
    building_shadows$shadow_source <- "building"
    building_shadows$solar_index <- i
    building_shadows$sun_id <- sun_ids[i]
    building_shadows$azimuth <- solar_pos[i, 1]
    building_shadows$elevation <- solar_pos[i, 2]
    if (is.null(canopy)) {
      return(building_shadows)
    }
    canopy_shadows <- compute_canopy_shadow_footprints(canopy, buildings, solar_pos[i, ], b)
    if (nrow(canopy_shadows) == 0) {
      return(building_shadows)
    }
    canopy_shadows$solar_index <- i
    canopy_shadows$sun_id <- sun_ids[i]
    canopy_shadows$azimuth <- solar_pos[i, 1]
    canopy_shadows$elevation <- solar_pos[i, 2]
    combine_shadow_sources(building_shadows, canopy_shadows)
  })
  shadows <- do.call(rbind, out)
  if (isTRUE(overlap_shadow)) {
    shadows <- combine_shadow_footprint_overlaps(shadows)
    if (isTRUE(plot)) {
      plot_shadow_footprints(
        buildings,
        shadows,
        canopy = canopy,
        plot_overlap_gradient = plot_overlap_gradient
      )
    }
    return(shadows)
  }
  if (isTRUE(plot)) {
    plot_shadow_footprints(
      buildings,
      shadows,
      canopy = canopy,
      plot_overlap_gradient = plot_overlap_gradient
    )
  }
  shadows
}

#' @description
#' `get_shadow_height()` computes shadow height at points or across a `terra`
#' surface. If `shadow_locations` is omitted, a `terra` template is generated
#' around the buildings.
#'
#' @export
#' @rdname get_shadows
get_shadow_height <- function(x = NULL,
                              shadow_locations = NULL,
                              solar_time = NULL,
                              time_zone = NULL,
                              azimuth = NULL,
                              elevation = NULL,
                              height_field = "Height",
                              min_tree_height = 2,
                              datasource_canopy_height = NULL,
                              key = NULL,
                              raster_buffer = NULL,
                              canopy_height = NULL,
                              dem = NULL,
                              cell_size = 2,
                              extent_buffer = NULL,
                              b = 0.01,
                              parallel = 1,
                              filter_footprint = FALSE,
                              quiet = TRUE,
                              ...) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygons.")
    return(NULL)
  }
  dots <- list(...)
  solar_args <- extract_deprecated_solar_args(dots)
  dots <- solar_args$dots
  if ("location" %in% names(dots)) {
    if (!is.null(shadow_locations)) {
      stop("Use only one of `shadow_locations` or deprecated `location`.", call. = FALSE)
    }
    warning(
      "`location` is deprecated; use `shadow_locations` instead.",
      call. = FALSE
    )
    shadow_locations <- dots$location
  }
  unknown_args <- setdiff(names(dots), "location")
  if (length(unknown_args) > 0) {
    stop("Unused argument(s): ", paste(unknown_args, collapse = ", "), call. = FALSE)
  }
  buildings <- prepare_shadow_buildings(x, height_field)
  solar_pos <- resolve_solar_inputs(
    buildings,
    azimuth = azimuth,
    elevation = elevation,
    solar_time = solar_time,
    time_zone = time_zone,
    solar_pos = solar_args$solar_pos,
    time = solar_args$time
  )
  raster_inputs <- resolve_shadow_raster_inputs(
    buildings = buildings,
    solar_pos = solar_pos,
    height_field = height_field,
    canopy_height = canopy_height,
    dem = dem,
    datasource_canopy_height = datasource_canopy_height,
    key = key,
    min_tree_height = min_tree_height,
    raster_buffer = if (is.null(raster_buffer)) extent_buffer else raster_buffer,
    quiet = quiet
  )
  canopy <- prepare_canopy_obstacles(
    raster_inputs$canopy_height,
    raster_inputs$dem,
    buildings,
    min_tree_height
  )
  loc <- prepare_shadow_location(
    shadow_locations = shadow_locations,
    buildings = buildings,
    height_field = height_field,
    solar_pos = solar_pos,
    cell_size = cell_size,
    extent_buffer = extent_buffer
  )
  if (!quiet) cli::cli_alert_info("Computing shadow height ...")
  if (inherits(loc, "SpatRaster")) {
    return(compute_shadow_height_spatraster(loc, buildings, solar_pos, height_field, b, canopy))
  }
  compute_shadow_height_points(loc, buildings, solar_pos, height_field, b, canopy)
}

#' @description
#' `get_radiation()` estimates direct, diffuse, and total radiation load on
#' roofs and facades represented by a 3D `sf` surface grid.
#'
#' @export
#' @rdname get_shadows
get_radiation <- function(x = NULL,
                          grid = NULL,
                          solar_time = NULL,
                          time_zone = NULL,
                          azimuth = NULL,
                          elevation = NULL,
                          solar_normal,
                          solar_diffuse,
                          height_field = "Height",
                          min_tree_height = 2,
                          datasource_canopy_height = NULL,
                          key = NULL,
                          raster_buffer = NULL,
                          canopy_transmissivity = 0.15,
                          canopy_height = NULL,
                          dem = NULL,
                          grid_res = 2,
                          offset = 0.01,
                          radius = Inf,
                          return_list = FALSE,
                          parallel = 1,
                          plot = FALSE,
                          quiet = TRUE,
                          ...) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygons.")
    return(NULL)
  }
  if (missing(solar_normal) || missing(solar_diffuse)) {
    stop("`solar_normal` and `solar_diffuse` are required.", call. = FALSE)
  }
  dots <- list(...)
  solar_args <- extract_deprecated_solar_args(dots)
  dots <- solar_args$dots
  if (length(dots) > 0) {
    stop("Unused argument(s): ", paste(names(dots), collapse = ", "), call. = FALSE)
  }
  buildings <- prepare_shadow_buildings(x, height_field)
  solar_pos <- resolve_solar_inputs(
    buildings,
    azimuth = azimuth,
    elevation = elevation,
    solar_time = solar_time,
    time_zone = time_zone,
    solar_pos = solar_args$solar_pos,
    time = solar_args$time
  )
  check_radiation_vectors(solar_pos, solar_normal, solar_diffuse)
  check_canopy_transmissivity(canopy_transmissivity)
  raster_inputs <- resolve_shadow_raster_inputs(
    buildings = buildings,
    solar_pos = solar_pos,
    height_field = height_field,
    canopy_height = canopy_height,
    dem = dem,
    datasource_canopy_height = datasource_canopy_height,
    key = key,
    min_tree_height = min_tree_height,
    raster_buffer = raster_buffer,
    quiet = quiet
  )
  canopy <- prepare_canopy_obstacles(
    raster_inputs$canopy_height,
    raster_inputs$dem,
    buildings,
    min_tree_height
  )
  surface <- prepare_radiation_grid(grid, buildings, height_field, grid_res, offset)
  if (!quiet) cli::cli_alert_info("Computing building surface radiation ...")
  rad <- compute_surface_radiation(surface, buildings, solar_pos, solar_normal,
                                   solar_diffuse, height_field, canopy,
                                   canopy_transmissivity)
  if (isTRUE(return_list)) {
    return(rad$by_time)
  }
  out <- surface
  out$svf <- rad$svf
  out$direct <- rad$direct
  out$diffuse <- rad$diffuse
  out$total <- rad$total
  if (isTRUE(plot)) {
    plot_radiation_surface(buildings, out)
  }
  out
}
