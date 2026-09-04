#' Building shadow and radiation calculations

# The functions of shadow and radiation were developed based on
# the archived package "shadow" by Dorman et al. (2019).
# The credit goes to Dorman et al.

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
#'   For `get_radiation()`, draw the default 2D base R radiation map colored by
#'   `total`. When ground samples are included, the 2D layout shows separate
#'   ground, facade, and roof maps with one shared legend. When canopy data are
#'   supplied, a second 2D map shows canopy impact as `canopy - no_canopy`
#'   total-radiation difference.
#' @param plot_3d Logical. For `get_radiation()`, draw a base R 3D-style view
#'   with separate panels for direct, diffuse, and total radiation. This is
#'   opt-in; `plot = TRUE` uses the 2D map layout by default.
#' @param plot_overlap_gradient Logical. For `get_shadow_footprint()` plots with
#'   multiple `solar_time` values, if `TRUE`, draw all shadows in transparent
#'   gray so overlapping areas appear darker.
#' @param scalebar Logical. Draw a distance scale bar using the shared map
#'   layout. Defaults to `TRUE` for plotted maps.
#' @param scalebar_unit Scale bar unit: `"auto"` (default), `"km"`, or `"m"`.
#' @param scalebar_cex Scale bar label size. Defaults to `0.7`.
#' @param north_arrow Logical. Draw a north arrow in the map panel. Defaults to
#'   `TRUE`.
#' @param shadow_locations Optional query locations for shadow height, as an
#'   `sf` point layer or a `terra::SpatRaster`.
#' @param cell_size Numeric cell resolution in CRS units when
#'   `shadow_locations` is omitted.
#' @param extent_buffer Optional numeric buffer around `x` used when creating an
#'   automatic `terra::SpatRaster` template. If omitted, a buffer is estimated
#'   from building heights and solar elevation.
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
#'   from building roofs and facades (and optionally the ground).
#' @param grid_res Numeric surface-grid resolution in CRS units.
#' @param ground Logical. If `TRUE`, add a regular grid of ground-level sample
#'   points over the study-area bounding box (excluding building footprints).
#'   Ground points have an upward normal and receive direct radiation whenever
#'   they are not in a building or canopy shadow, and diffuse radiation scaled
#'   by their Sky View Factor. Returned rows have `surface = "ground"` and
#'   `building_id = NA`.
#' @param ground_res Numeric resolution for the ground sample grid in CRS
#'   units. If `NULL`, defaults to `grid_res`.
#' @param offset Numeric vertical offset added to generated surface-grid points.
#' @param solar_normal Direct Normal Irradiance vector, one value per solar
#'   position.
#' @param solar_diffuse Diffuse Horizontal Irradiance vector, one value per solar
#'   position.
#' @param radius Maximum obstacle search distance in CRS units for radiation
#'   Sky View Factor calculations. Defaults to `500`. Obstacles beyond this
#'   distance contribute negligibly to Sky View Factor but dominate runtime,
#'   so a finite radius enables spatial culling and is typically many times
#'   faster. Use `Inf` to consider all obstacles regardless of distance.
#' @param svf_res_angle Numeric. Azimuth sampling interval in decimal degrees
#'   used when estimating Sky View Factor inside `get_radiation()`.
#' @param return_list Logical. If `TRUE`, return per-timestep radiation matrices
#'   instead of a summed `sf` surface grid.
#' @param res_angle Numeric. Azimuth sampling interval in decimal degrees for
#'   `svf()`. Smaller values are slower and more detailed.
#' @param observer_height Numeric. Height above local ground for SVF query
#'   locations. Defaults to `1.7`, representing pedestrian eye level.
#' @param max_distance Numeric. Maximum obstacle search distance in CRS units
#'   for `svf()`.
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
#' `svf()` returns a `terra::SpatRaster` with Sky View Factor values from 0 to 1.
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
#' `svf()` computes a Sky View Factor raster from building and optional canopy
#' obstacles.
#'
#' @export
#' @rdname get_shadows
svf <- function(x = NULL,
                height_field = "Height",
                min_tree_height = 2,
                datasource_canopy_height = NULL,
                key = NULL,
                canopy_height = NULL,
                dem = NULL,
                raster_buffer = NULL,
                grid_res = 2,
                extent_buffer = NULL,
                res_angle = 5,
                observer_height = 1.7,
                max_distance = NULL,
                plot = FALSE,
                scalebar = TRUE,
                scalebar_unit = c("auto", "km", "m"),
                scalebar_cex = 0.7,
                north_arrow = TRUE,
                quiet = TRUE) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygons.")
    return(NULL)
  }
  scalebar_unit <- match.arg(scalebar_unit)
  buildings <- prepare_shadow_buildings(x, height_field)
  if (!is.numeric(res_angle) || length(res_angle) != 1 || is.na(res_angle) ||
      res_angle <= 0 || res_angle > 360) {
    stop("`res_angle` must be a single positive value up to 360.", call. = FALSE)
  }
  if (!is.numeric(observer_height) || length(observer_height) != 1 || is.na(observer_height)) {
    stop("`observer_height` must be a single numeric value.", call. = FALSE)
  }
  if (is.null(extent_buffer)) {
    extent_buffer <- estimate_svf_buffer(buildings, height_field, canopy_height, max_distance)
  }
  raster_inputs <- resolve_shadow_raster_inputs(
    buildings = buildings,
    solar_pos = matrix(c(0, 45), ncol = 2, dimnames = list(NULL, c("azimuth", "elevation"))),
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
  template <- make_svf_template(buildings, grid_res, extent_buffer)
  if (!quiet) cli::cli_alert_info("Computing Sky View Factor raster ...")
  out <- compute_svf_raster(
    template = template,
    buildings = buildings,
    height_field = height_field,
    canopy = canopy,
    dem = raster_inputs$dem,
    res_angle = res_angle,
    observer_height = observer_height,
    max_distance = max_distance
  )
  if (isTRUE(plot)) {
    plot_svf_raster(
      buildings,
      out,
      scalebar = scalebar,
      scalebar_unit = scalebar_unit,
      scalebar_cex = scalebar_cex,
      north_arrow = north_arrow
    )
  }
  out
}

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
                                 scalebar = TRUE,
                                 scalebar_unit = c("auto", "km", "m"),
                                 scalebar_cex = 0.7,
                                 north_arrow = TRUE,
                                 quiet = TRUE) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygons.")
    return(NULL)
  }
  scalebar_unit <- match.arg(scalebar_unit)
  buildings <- prepare_shadow_buildings(x, height_field)
  solar_pos <- resolve_solar_inputs(
    buildings,
    azimuth = azimuth,
    elevation = elevation,
    solar_time = solar_time,
    time_zone = time_zone
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
        plot_overlap_gradient = plot_overlap_gradient,
        scalebar = scalebar,
        scalebar_unit = scalebar_unit,
        scalebar_cex = scalebar_cex,
        north_arrow = north_arrow
      )
    }
    return(shadows)
  }
  if (isTRUE(plot)) {
    plot_shadow_footprints(
      buildings,
      shadows,
      canopy = canopy,
      plot_overlap_gradient = plot_overlap_gradient,
      scalebar = scalebar,
      scalebar_unit = scalebar_unit,
      scalebar_cex = scalebar_cex,
      north_arrow = north_arrow
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
                              filter_footprint = FALSE,
                              quiet = TRUE) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygons.")
    return(NULL)
  }
  buildings <- prepare_shadow_buildings(x, height_field)
  solar_pos <- resolve_solar_inputs(
    buildings,
    azimuth = azimuth,
    elevation = elevation,
    solar_time = solar_time,
    time_zone = time_zone
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
                          ground = FALSE,
                          ground_res = NULL,
                          offset = 0.01,
                          radius = 500,
                          svf_res_angle = 15,
                          return_list = FALSE,
                          plot = FALSE,
                          plot_3d = FALSE,
                          scalebar = TRUE,
                          scalebar_unit = c("auto", "km", "m"),
                          scalebar_cex = 0.7,
                          north_arrow = TRUE,
                          quiet = TRUE) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygons.")
    return(NULL)
  }
  scalebar_unit <- match.arg(scalebar_unit)
  if (missing(solar_normal) || missing(solar_diffuse)) {
    stop("`solar_normal` and `solar_diffuse` are required.", call. = FALSE)
  }
  buildings <- prepare_shadow_buildings(x, height_field)
  solar_pos <- resolve_solar_inputs(
    buildings,
    azimuth = azimuth,
    elevation = elevation,
    solar_time = solar_time,
    time_zone = time_zone
  )
  check_radiation_vectors(solar_pos, solar_normal, solar_diffuse)
  check_canopy_transmissivity(canopy_transmissivity)
  if (!is.numeric(svf_res_angle) || length(svf_res_angle) != 1 ||
      is.na(svf_res_angle) || svf_res_angle <= 0 || svf_res_angle > 360) {
    stop("`svf_res_angle` must be a single positive value up to 360.", call. = FALSE)
  }
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
  surface <- prepare_radiation_grid(
    grid, buildings, height_field, grid_res, offset,
    ground = ground, ground_res = ground_res,
    dem = raster_inputs$dem
  )
  if (!quiet) cli::cli_alert_info("Computing building surface radiation ...")
  rad <- compute_surface_radiation(surface, buildings, solar_pos, solar_normal,
                                   solar_diffuse, height_field, canopy,
                                   canopy_transmissivity,
                                   dem = raster_inputs$dem,
                                   svf_res_angle = svf_res_angle,
                                   radius = radius)
  if (isTRUE(return_list)) {
    return(rad$by_time)
  }
  out <- surface
  out$svf <- rad$svf
  out$direct <- rad$direct
  out$diffuse <- rad$diffuse
  out$total <- rad$total
  if (isTRUE(plot)) {
    no_canopy <- NULL
    if (!is.null(canopy)) {
      if (!quiet) cli::cli_alert_info("Computing no-canopy radiation baseline for impact plot ...")
      no_canopy <- compute_surface_radiation(
        surface, buildings, solar_pos, solar_normal, solar_diffuse,
        height_field, canopy = NULL, canopy_transmissivity = 1,
        dem = raster_inputs$dem,
        svf_res_angle = svf_res_angle,
        radius = radius
      )
    }
    plot_radiation_surface(
      buildings,
      out,
      scalebar = scalebar,
      scalebar_unit = scalebar_unit,
      scalebar_cex = scalebar_cex,
      north_arrow = north_arrow
    )
    if (!is.null(no_canopy)) {
      plot_radiation_canopy_impact(
        buildings,
        out,
        no_canopy,
        scalebar = scalebar,
        scalebar_unit = scalebar_unit,
        scalebar_cex = scalebar_cex,
        north_arrow = north_arrow
      )
    }
  }
  if (isTRUE(plot_3d)) {
    plot_radiation_surface_3d(buildings, out, height_field)
  }
  out
}
