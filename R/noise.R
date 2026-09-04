#' Default OSM traffic assumptions for screening-level road noise
#'
#' @description
#' Returns the built-in lookup table used by `infer_osm_traffic()` when observed
#' traffic counts are unavailable. The defaults mirror the road-category
#' assumptions embedded in NoiseModelling's `Import_OSM.groovy`, which cites the
#' Good Practice Guide for Strategic Noise Mapping and the Production of
#' Associated Data on Noise Exposure, version 2. They should still be replaced
#' with local traffic counts whenever available.
#'
#' @return A data frame with OSM highway class, speed, and hourly light/heavy
#' vehicle assumptions.
#' @export
osm_noise_traffic_defaults <- function() {
  highways <- list(
    c("motorway", "motorway_link"),
    c("trunk", "trunk_link"),
    c("primary", "primary_link"),
    c("secondary", "secondary_link"),
    c("tertiary", "tertiary_link", "unclassified"),
    "residential",
    c("service", "living_street")
  )
  category <- rep(seq_along(highways) - 1L, lengths(highways))
  highway <- unlist(highways, use.names = FALSE)
  aadf_d <- c(26103, 17936, 7124, 1400, 700, 350, 175)
  aadf_e <- c(7458, 3826, 1069, 400, 200, 100, 50)
  aadf_n <- c(3729, 2152, 712, 200, 100, 50, 25)
  hv_d <- c(0.25, 0.2, 0.2, 0.15, 0.10, 0.05, 0.02)
  hv_e <- c(0.35, 0.2, 0.15, 0.10, 0.06, 0.02, 0.01)
  hv_n <- c(0.45, 0.2, 0.1, 0.05, 0.03, 0.01, 0.0)
  speed <- c(130, 110, 80, 80, 50, 30, 30)
  index <- category + 1L
  data.frame(
    highway = highway,
    category = category,
    speed_kmh = speed[index],
    light_veh_h = (1 - hv_d[index]) * aadf_d[index] / 12,
    heavy_veh_h = hv_d[index] * aadf_d[index] / 12,
    LV_D = (1 - hv_d[index]) * aadf_d[index] / 12,
    LV_E = (1 - hv_e[index]) * aadf_e[index] / 4,
    LV_N = (1 - hv_n[index]) * aadf_n[index] / 8,
    HGV_D = hv_d[index] * aadf_d[index] / 12,
    HGV_E = hv_e[index] * aadf_e[index] / 4,
    HGV_N = hv_n[index] * aadf_n[index] / 8,
    LV_SPD_D = speed[index],
    LV_SPD_E = speed[index],
    LV_SPD_N = speed[index],
    HGV_SPD_D = pmin(90, speed[index]),
    HGV_SPD_E = pmin(90, speed[index]),
    HGV_SPD_N = pmin(90, speed[index]),
    PVMT = "NL08",
    stringsAsFactors = FALSE
  )
}

#' Infer screening-level road traffic inputs from OSM road attributes
#'
#' @description
#' Adds traffic volume and speed columns to an OSM road `sf` object based on the
#' `highway` class and optional `maxspeed`, `lanes`, or `oneway` values. These
#' values are intended for relative/scenario noise mapping when measured traffic
#' counts are unavailable, not calibrated regulatory noise maps.
#'
#' @param roads `sf` line object containing at least a `highway` column.
#' @param defaults Data frame of road-class assumptions. Defaults to
#' `osm_noise_traffic_defaults()`.
#' @param use_maxspeed Logical. If `TRUE`, parse numeric speeds from `maxspeed`
#' and use them when available. Defaults to `FALSE` to match NoiseModelling's
#' OSM import defaults.
#' @param use_lanes Logical. If `TRUE`, scale vehicle counts by lane count
#' relative to a two-lane road. Defaults to `FALSE` to match NoiseModelling's
#' OSM import defaults.
#' @param use_oneway Logical. If `TRUE`, halve traffic counts on OSM one-way
#' roads. Defaults to `FALSE`.
#' @param quiet Logical. If `TRUE`, suppress informational messages.
#'
#' @return An `sf` object with inferred traffic columns, including
#' `speed_kmh`, `light_veh_h`, `heavy_veh_h`, and NoiseModelling-style
#' `LV_*`, `MV_*`, `HGV_*`, `WAV_*`, `WBV_*`, and speed columns for
#' day/evening/night.
#' @export
infer_osm_traffic <- function(roads,
                              defaults = osm_noise_traffic_defaults(),
                              use_maxspeed = FALSE,
                              use_lanes = FALSE,
                              use_oneway = FALSE,
                              quiet = TRUE) {
  if (!inherits(roads, "sf")) {
    stop("`roads` must be an sf object.", call. = FALSE)
  }
  if (!"highway" %in% names(roads)) {
    stop("`roads` must include an OSM `highway` column.", call. = FALSE)
  }
  required <- c(
    "highway", "speed_kmh", "light_veh_h", "heavy_veh_h",
    "LV_D", "LV_E", "LV_N", "HGV_D", "HGV_E", "HGV_N", "PVMT"
  )
  missing_defaults <- setdiff(required, names(defaults))
  if (length(missing_defaults) > 0) {
    stop("`defaults` is missing required columns: ", paste(missing_defaults, collapse = ", "), call. = FALSE)
  }

  out <- roads
  out$highway_class <- normalize_osm_highway(out$highway)
  default_index <- match(out$highway_class, defaults$highway)

  fallback <- match("residential", defaults$highway)
  if (is.na(fallback)) fallback <- 1L
  default_index[is.na(default_index)] <- fallback

  out$speed_kmh <- as.numeric(defaults$speed_kmh[default_index])
  out$light_veh_h <- as.numeric(defaults$light_veh_h[default_index])
  out$heavy_veh_h <- as.numeric(defaults$heavy_veh_h[default_index])

  if (isTRUE(use_maxspeed) && "maxspeed" %in% names(out)) {
    parsed_speed <- parse_osm_speed_kmh(out$maxspeed)
    out$speed_kmh[!is.na(parsed_speed)] <- parsed_speed[!is.na(parsed_speed)]
  }

  if (isTRUE(use_lanes) && "lanes" %in% names(out)) {
    lanes <- parse_osm_number(out$lanes)
    lane_factor <- pmax(lanes / 2, 1, na.rm = TRUE)
    lane_factor[is.na(lanes)] <- 1
    out$light_veh_h <- out$light_veh_h * lane_factor
    out$heavy_veh_h <- out$heavy_veh_h * lane_factor
  }

  out$traffic_source <- "osm_inferred"
  out$LV_D <- as.numeric(defaults$LV_D[default_index])
  out$LV_E <- as.numeric(defaults$LV_E[default_index])
  out$LV_N <- as.numeric(defaults$LV_N[default_index])
  out$HGV_D <- as.numeric(defaults$HGV_D[default_index])
  out$HGV_E <- as.numeric(defaults$HGV_E[default_index])
  out$HGV_N <- as.numeric(defaults$HGV_N[default_index])

  if (isTRUE(use_lanes) && "lanes" %in% names(out)) {
    out$LV_D <- out$LV_D * lane_factor
    out$LV_E <- out$LV_E * lane_factor
    out$LV_N <- out$LV_N * lane_factor
    out$HGV_D <- out$HGV_D * lane_factor
    out$HGV_E <- out$HGV_E * lane_factor
    out$HGV_N <- out$HGV_N * lane_factor
  }
  if (isTRUE(use_oneway) && "oneway" %in% names(out)) {
    one_way <- parse_osm_oneway(out$oneway)
    out$LV_D[one_way] <- out$LV_D[one_way] / 2
    out$LV_E[one_way] <- out$LV_E[one_way] / 2
    out$LV_N[one_way] <- out$LV_N[one_way] / 2
    out$HGV_D[one_way] <- out$HGV_D[one_way] / 2
    out$HGV_E[one_way] <- out$HGV_E[one_way] / 2
    out$HGV_N[one_way] <- out$HGV_N[one_way] / 2
  }
  out$light_veh_h <- out$LV_D
  out$heavy_veh_h <- out$HGV_D
  out$MV_D <- 0
  out$MV_E <- 0
  out$MV_N <- 0
  out$WAV_D <- 0
  out$WAV_E <- 0
  out$WAV_N <- 0
  out$WBV_D <- 0
  out$WBV_E <- 0
  out$WBV_N <- 0
  out$LV_SPD_D <- out$speed_kmh
  out$LV_SPD_E <- out$speed_kmh
  out$LV_SPD_N <- out$speed_kmh
  out$MV_SPD_D <- out$speed_kmh
  out$MV_SPD_E <- out$speed_kmh
  out$MV_SPD_N <- out$speed_kmh
  out$HGV_SPD_D <- pmin(90, out$speed_kmh)
  out$HGV_SPD_E <- pmin(90, out$speed_kmh)
  out$HGV_SPD_N <- pmin(90, out$speed_kmh)
  out$WAV_SPD_D <- out$speed_kmh
  out$WAV_SPD_E <- out$speed_kmh
  out$WAV_SPD_N <- out$speed_kmh
  out$WBV_SPD_D <- out$speed_kmh
  out$WBV_SPD_E <- out$speed_kmh
  out$WBV_SPD_N <- out$speed_kmh
  out$PVMT <- as.character(defaults$PVMT[default_index])
  out$TEMP_D <- 20
  out$TEMP_E <- 20
  out$TEMP_N <- 20
  out$WAY <- 3

  if (!quiet) {
    cli::cli_alert_info("Traffic volumes and speeds were inferred from OSM road class; replace with observed counts when available.")
  }
  out
}

#' Prepare NoiseModelling-style input layers
#'
#' @description
#' Builds a compact set of layers for an integrated or external NoiseModelling
#' workflow: `BUILDINGS`, `ROADS`, `GROUND`, `RECEIVERS`, and optional `DEM`.
#' The returned layers can be inspected directly in R and optionally written to
#' a GeoPackage.
#'
#' @param x `sf` polygon object with building footprints and a height field.
#' @param roads Optional `sf` line object, usually from OSM. If `NULL`, OSM
#' roads are downloaded from the bounding box of `x`. If inferred traffic
#' columns are absent, `infer_osm_traffic()` is applied.
#' @param download_roads Logical. If `TRUE` and `roads = NULL`, download roads
#' from OSM using the building extent. Set to `FALSE` only when a later workflow
#' step supplies roads, for example `get_noise_map(osm_file = ...)`.
#' @param population Logical. If `TRUE`, assign GHSL population to buildings
#' with the package population helper (or `get_pop()` when available) before writing
#' NoiseModelling `BUILDINGS`. Existing `POP`
#' or `population_field` values are used when available.
#' @param population_field Optional column containing building population. The
#' value is copied to NoiseModelling's `POP` field.
#' @param population_year GHSL population year passed to the package population
#' function when
#' `population = TRUE` and no usable population field is already present.
#' @param greenspace Optional `sf` polygons or `terra::SpatRaster` marking green
#' areas. Green areas are translated to ground absorption `G = 1`.
#' @param height_field Building height column name. Defaults to `"Height"`.
#' @param min_tree_height Minimum canopy height treated as green cover.
#' @param datasource_canopy_height Character or `NULL`. Canopy height source to
#' retrieve internally when `canopy_height` is not supplied. Currently supports
#' `"metachm"` and `"ethCHM"`.
#' @param datasource_greenspace character or `NULL`. Optional 2D greenspace map
#' tile source for the visible-green feature layer. Supports \code{"esri"} and
#' \code{"sentinel2"}. If both canopy height and greenspace sources are supplied,
#' the visible-green layer is the union of height-filtered canopy and 2D
#' greenspace.
#' @param greenspace_year numeric. The desired year for Sentinel-2 cloudless mosaic
#' tiles. (This has to be specified when `datasource_greenspace = "sentinel2"`)
#' @param greenspace_zoom numeric. Zoom level of map tile when
#' `datasource_greenspace = "esri"` or `"sentinel2"`.
#' @param opentopo_key OpenTopography API key used to retrieve DEM data internally when
#' `dem` is not supplied.
#' @param canopy_height Optional `terra::SpatRaster` canopy height map. Cells
#' greater than or equal to `min_tree_height` are treated as green ground.
#' @param dem Optional `terra::SpatRaster` terrain layer.
#' @param receiver One of `"grid"` or `"none"`.
#' @param resolution Receiver grid resolution in map units.
#' @param ground_default Default ground absorption for the analysis extent.
#' @param green_ground Ground absorption assigned to greenspace/canopy polygons.
#' @param out_dir Optional output directory for a GeoPackage.
#' @param write Logical. If `TRUE`, write layers to `noise_inputs.gpkg`.
#' @param quiet Logical. If `TRUE`, suppress informational messages.
#' @return A list with prepared `sf`/`terra` layers and optional GeoPackage path.
#' @export
prepare_noisemodelling_inputs <- function(x = NULL,
                                          height_field = "Height",
                                          datasource_canopy_height = NULL,
                                          datasource_greenspace = NULL,
                                          greenspace_year = NULL,
                                          greenspace_zoom = 17,
                                          opentopo_key = NULL,
                                          canopy_height = NULL,
                                          roads = NULL,
                                          download_roads = TRUE,
                                          greenspace = NULL,
                                          dem = NULL,
                                          population = FALSE,
                                          population_field = NULL,
                                          population_year = 2025,
                                          min_tree_height = 2,
                                          receiver = c("grid", "none"),
                                          resolution = 10,
                                          ground_default = 0,
                                          green_ground = 1,
                                          out_dir = NULL,
                                          write = FALSE,
                                          quiet = TRUE) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygons.")
    return(NULL)
  }
  receiver <- match.arg(receiver)
  if (!inherits(x, "sf")) stop("`x` must be an sf object.", call. = FALSE)
  if (!height_field %in% names(x)) {
    stop("`x` must include the height column `", height_field, "`.", call. = FALSE)
  }
  if (!is.numeric(resolution) || length(resolution) != 1 || is.na(resolution) || resolution <= 0) {
    stop("`resolution` must be a positive numeric value.", call. = FALSE)
  }

  target_crs <- sf::st_crs(x)
  if (is.na(target_crs)) {
    stop("`x` must have a coordinate reference system.", call. = FALSE)
  }

  raster_inputs <- resolve_noise_raster_inputs(
    buildings = x,
    roads = roads,
    greenspace = greenspace,
    canopy_height = canopy_height,
    dem = dem,
    datasource_canopy_height = datasource_canopy_height,
    datasource_greenspace = datasource_greenspace,
    greenspace_year = greenspace_year,
    greenspace_zoom = greenspace_zoom,
    opentopo_key = opentopo_key,
    min_tree_height = min_tree_height,
    quiet = quiet
  )

  buildings_nm <- prepare_noise_buildings(
    x,
    height_field,
    population = population,
    population_field = population_field,
    population_year = population_year,
    quiet = quiet
  )
  if (is.null(roads) && isTRUE(download_roads)) {
    if (!quiet) cli::cli_alert_info("Downloading OSM roads for noise analysis ...")
    roads <- get_osm_roads_for_buildings(x)
  }
  if (!is.null(roads) && !inherits(roads, "sf")) stop("`roads` must be an sf object.", call. = FALSE)
  roads_nm <- NULL
  if (!is.null(roads)) {
    roads <- sf::st_transform(roads, target_crs)
    roads_nm <- prepare_noise_roads(roads, quiet = quiet)
  }
  analysis_extent <- if (is.null(roads_nm)) {
    sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(buildings_nm)))
  } else {
    make_noise_analysis_extent(buildings_nm, roads_nm)
  }
  ground_nm <- prepare_noise_ground(
    analysis_extent,
    greenspace = raster_inputs$greenspace,
    canopy_height = raster_inputs$canopy_height,
    ground_default = ground_default,
    green_ground = green_ground,
    min_tree_height = min_tree_height
  )
  receivers_nm <- if (receiver == "grid" && !is.null(roads_nm)) {
    make_noise_receiver_grid(buildings_nm, roads_nm, resolution)
  } else {
    NULL
  }

  dem_out <- raster_inputs$dem
  if (!is.null(dem_out)) {
    if (!inherits(dem_out, "SpatRaster")) {
      stop("`dem` must be a terra SpatRaster when supplied.", call. = FALSE)
    }
    dem_out <- terra::project(dem_out, target_crs$wkt, method = "bilinear")
  }

  gpkg <- NULL
  if (isTRUE(write)) {
    if (is.null(out_dir)) out_dir <- tempdir()
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    gpkg <- file.path(out_dir, "noise_inputs.gpkg")
    write_noise_gpkg(
      gpkg,
      buildings = buildings_nm,
      roads = roads_nm,
      ground = ground_nm,
      receivers = receivers_nm,
      quiet = quiet
    )
  }

  list(
    buildings = buildings_nm,
    roads = roads_nm,
    ground = ground_nm,
    receivers = receivers_nm,
    dem = dem_out,
    gpkg = gpkg,
    traffic_note = "Traffic volume and speed are screening-level OSM inferences unless replaced with observed counts."
  )
}

#' Install the headless NoiseModelling runner
#'
#' @description
#' Downloads and unzips the official headless NoiseModelling release into an R
#' user cache directory. This is used by [get_noise_map()] when `run = TRUE` and
#' no `nm_path` is supplied.
#'
#' @param version NoiseModelling release version. Defaults to `"5.0.1"`.
#' @param destdir Destination directory. Defaults to the `gloBFPr` user cache.
#' @param quiet Logical. If `TRUE`, suppress download messages.
#'
#' @return Path to the installed NoiseModelling directory.
#' @export
install_noisemodelling <- function(version = "5.0.1",
                                   destdir = noisemodelling_cache_dir(),
                                   quiet = TRUE) {
  version <- as.character(version[1])
  target <- file.path(destdir, paste0("NoiseModelling_without_gui-", version))
  wps <- noisemodelling_wps_script(target)
  if (file.exists(wps)) {
    return(target)
  }
  if (!dir.exists(destdir)) dir.create(destdir, recursive = TRUE)
  url <- sprintf(
    "https://github.com/Universite-Gustave-Eiffel/NoiseModelling/releases/download/v%s/NoiseModelling_without_gui-%s.zip",
    version,
    version
  )
  zipfile <- tempfile(fileext = ".zip")
  on.exit(unlink(zipfile), add = TRUE)
  if (!quiet) cli::cli_alert_info("Downloading NoiseModelling {version} headless runner ...")
  utils::download.file(url, zipfile, mode = "wb", quiet = quiet)
  utils::unzip(zipfile, exdir = destdir)
  wps <- noisemodelling_wps_script(target)
  if (!file.exists(wps)) {
    stop("NoiseModelling was downloaded, but the `bin/wps_scripts` runner was not found.", call. = FALSE)
  }
  if (.Platform$OS.type != "windows") {
    Sys.chmod(wps, mode = "0755")
  }
  target
}

#' Create a screening-level urban road-noise workflow object
#'
#' @description
#' Convenience wrapper around `prepare_noisemodelling_inputs()`. With
#' `run = FALSE`, it returns reproducible input layers. With `run = TRUE`, it
#' runs the official headless NoiseModelling WPS scripts and reads the
#' `RECEIVERS_LEVEL` output table back into R. The returned `noise_map` element
#' joins the calculated levels to receiver geometries for point mapping, and
#' `isophones` contains NoiseModelling's `CONTOURING_NOISE_MAP` polygons.
#'
#' @inheritParams prepare_noisemodelling_inputs
#' @param run Logical. If `TRUE`, run NoiseModelling after preparing inputs.
#' @param nm_path Optional path to a `NoiseModelling_without_gui-*` directory.
#' If `NULL`, the package uses `options(gloBFPr.noisemodelling.path)`,
#' `NOISEMODELLING_HOME`, or a cached download.
#' @param nm_version NoiseModelling version to download when needed.
#' @param download_nm Logical. If `TRUE`, download the headless NoiseModelling
#' runner when `nm_path` is not available.
#' @param osm_file Optional path to a local `.osm`, `.osm.gz`, or `.osm.pbf`
#' extract. When supplied with `roads = NULL` and `run = TRUE`, NoiseModelling's
#' `Import_OSM.groovy` creates the `ROADS` table internally using its OSM road
#' defaults. This is different from the default R workflow, which downloads
#' roads from the building bounding box.
#' @param osm_remove_tunnels Logical passed to NoiseModelling `Import_OSM` when
#' `osm_file` is used. If `TRUE`, OSM roads tagged `tunnel=yes` are removed.
#' @param osm_eliminate_no_traffic_roads Logical passed to NoiseModelling
#' `Import_OSM` when `osm_file` is used. If `TRUE`, keeps only road classes that
#' NoiseModelling treats as traffic-bearing roads.
#' @param java Optional Java executable, Java home directory, or installed Java
#' major version such as `17`. If `NULL`, the function uses `JAVA_HOME` or
#' `java` on `PATH`. NoiseModelling 5.0.1 works with Java 11-21; Java 17 is
#' recommended.
#' @param keep_files Logical. If `TRUE`, keep the temporary NoiseModelling work
#' directory and include it in the result.
#' @param wall_alpha Wall absorption coefficient passed to NoiseModelling.
#' @param reflection_order Reflection order passed to NoiseModelling.
#' @param max_src_distance Maximum source-receiver distance in meters.
#' @param max_reflection_distance Maximum reflection distance in meters.
#' @param thread_count Number of NoiseModelling worker threads. `0` lets
#' NoiseModelling choose.
#' @param diffraction_vertical Logical. Passes `confDiffVertical`; enables
#' diffraction around vertical edges. NoiseModelling notes that CNOSSOS-EU uses
#' this mainly for rail and industrial sources.
#' @param diffraction_horizontal Logical. Passes `confDiffHorizontal`; enables
#' diffraction over horizontal building/terrain edges.
#' @param export_source_id Logical. Passes `confExportSourceId`; if `TRUE`,
#' receiver levels are kept by source identifier instead of merged across all
#' sources. This is useful for source-contribution diagnostics and can greatly
#' enlarge the output.
#' @param humidity Relative humidity percentage for atmospheric absorption.
#' Passes `confHumidity`. Defaults to NoiseModelling's script default of `70`.
#' @param temperature Air temperature in degrees Celsius for atmospheric
#' absorption. Passes `confTemperature`. Defaults to NoiseModelling's script
#' default of `15`.
#' @param favourable_occurrences Probability of favourable propagation
#' conditions by 16 wind-direction sectors, clockwise, where the north sector is
#' the last value. Supply one value to recycle to all sectors or 16 values. The
#' NoiseModelling default is sixteen `0.5` values.
#' @param rays_name Optional table name or file URL passed as `confRaysName`.
#' When supplied, NoiseModelling exports propagation rays/attenuation details for
#' advanced diagnostics. This can create very large outputs.
#' @param max_error Maximum allowed error in dB for pruning negligible source
#' contributions. Passes `confMaxError`; NoiseModelling's default is `0.1`.
#' @param frequency_field_prepend Prefix for source spectral columns, passed as
#' `frequencyFieldPrepend`. Defaults to `"HZ"` for columns such as `HZ1000`.
#' @param noise_wps_args Optional named list of additional raw arguments passed
#' to `Noise_level_from_source.groovy`. Use this only for advanced
#' NoiseModelling options not yet exposed directly.
#' @param delaunay_max_area Maximum Delaunay triangle area in square map units.
#' Defaults to `resolution^2`. Smaller values create denser receiver/contour
#' meshes and slower runs.
#' @param road_width Receiver exclusion distance around roads in meters for the
#' NoiseModelling Delaunay grid.
#' @param building_buffer Receiver exclusion distance around buildings in meters
#' for the NoiseModelling Delaunay grid.
#' @param iso_levels Numeric vector of isosurface breakpoints in dB passed to
#' `Create_Isosurface.groovy`.
#' @param iso_field Result field used to build isosurfaces. Defaults to `"LAEQ"`.
#' @param iso_smooth Smoothing coefficient passed to `Create_Isosurface.groovy`.
#'
#' @return Prepared noise input layers when `run = FALSE`; otherwise a list with
#' prepared inputs, receiver-level results, output file paths, and logs.
#' @export
get_noise_map <- function(x = NULL,
                          height_field = "Height",
                          datasource_canopy_height = NULL,
                          datasource_greenspace = NULL,
                          greenspace_year = NULL,
                          greenspace_zoom = 17,
                          opentopo_key = NULL,
                          roads = NULL,
                          greenspace = NULL,
                          canopy_height = NULL,
                          dem = NULL,
                          population = FALSE,
                          population_field = NULL,
                          population_year = 2025,
                          min_tree_height = 2,
                          receiver = c("grid", "none"),
                          resolution = 10,
                          out_dir = NULL,
                          write = FALSE,
                          run = FALSE,
                          nm_path = NULL,
                          nm_version = "5.0.1",
                          download_nm = TRUE,
                          osm_file = NULL,
                          osm_remove_tunnels = TRUE,
                          osm_eliminate_no_traffic_roads = TRUE,
                          java = NULL,
                          keep_files = FALSE,
                          wall_alpha = 0.1,
                          reflection_order = 0,
                          max_src_distance = 150,
                          max_reflection_distance = 50,
                          thread_count = 0,
                          diffraction_vertical = FALSE,
                          diffraction_horizontal = FALSE,
                          export_source_id = FALSE,
                          humidity = NULL,
                          temperature = NULL,
                          favourable_occurrences = NULL,
                          rays_name = NULL,
                          max_error = NULL,
                          frequency_field_prepend = "HZ",
                          noise_wps_args = NULL,
                          delaunay_max_area = NULL,
                          road_width = 2,
                          building_buffer = 2,
                          iso_levels = c(35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 200),
                          iso_field = "LAEQ",
                          iso_smooth = 0.5,
                          quiet = TRUE) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygons.")
    return(NULL)
  }
  if (!is.null(osm_file) && !is.null(roads)) {
    stop("Use either `roads` or `osm_file`, not both.", call. = FALSE)
  }
  inputs <- prepare_noisemodelling_inputs(
    x = x,
    roads = roads,
    greenspace = greenspace,
    height_field = height_field,
    min_tree_height = min_tree_height,
    datasource_canopy_height = datasource_canopy_height,
    datasource_greenspace = datasource_greenspace,
    greenspace_year = greenspace_year,
    greenspace_zoom = greenspace_zoom,
    opentopo_key = opentopo_key,
    canopy_height = canopy_height,
    dem = dem,
    population = population,
    population_field = population_field,
    population_year = population_year,
    download_roads = is.null(osm_file),
    receiver = receiver,
    resolution = resolution,
    out_dir = out_dir,
    write = write,
    quiet = quiet
  )
  if (!isTRUE(run)) {
    return(inputs)
  }
  run_noisemodelling(
    inputs = inputs,
    nm_path = nm_path,
    nm_version = nm_version,
    download_nm = download_nm,
    osm_file = osm_file,
    osm_remove_tunnels = osm_remove_tunnels,
    osm_eliminate_no_traffic_roads = osm_eliminate_no_traffic_roads,
    java = java,
    out_dir = out_dir,
    keep_files = keep_files,
    wall_alpha = wall_alpha,
    reflection_order = reflection_order,
    max_src_distance = max_src_distance,
    max_reflection_distance = max_reflection_distance,
    thread_count = thread_count,
    diffraction_vertical = diffraction_vertical,
    diffraction_horizontal = diffraction_horizontal,
    export_source_id = export_source_id,
    humidity = humidity,
    temperature = temperature,
    favourable_occurrences = favourable_occurrences,
    rays_name = rays_name,
    max_error = max_error,
    frequency_field_prepend = frequency_field_prepend,
    noise_wps_args = noise_wps_args,
    delaunay_max_area = if (is.null(delaunay_max_area)) resolution^2 else delaunay_max_area,
    road_width = road_width,
    building_buffer = building_buffer,
    iso_levels = iso_levels,
    iso_field = iso_field,
    iso_smooth = iso_smooth,
    surface_resolution = resolution,
    quiet = quiet
  )
}

#' @noRd
run_noisemodelling <- function(inputs,
                               nm_path = NULL,
                               nm_version = "5.0.1",
                               download_nm = TRUE,
                               osm_file = NULL,
                               osm_remove_tunnels = TRUE,
                               osm_eliminate_no_traffic_roads = TRUE,
                               java = NULL,
                               out_dir = NULL,
                               keep_files = FALSE,
                               wall_alpha = 0.1,
                               reflection_order = 0,
                               max_src_distance = 150,
                               max_reflection_distance = 50,
                               thread_count = 0,
                               diffraction_vertical = FALSE,
                               diffraction_horizontal = FALSE,
                               export_source_id = FALSE,
                               humidity = NULL,
                               temperature = NULL,
                               favourable_occurrences = NULL,
                               rays_name = NULL,
                               max_error = NULL,
                               frequency_field_prepend = "HZ",
                               noise_wps_args = NULL,
                               delaunay_max_area = 100,
                               road_width = 2,
                               building_buffer = 2,
                               iso_levels = c(35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 200),
                               iso_field = "LAEQ",
                               iso_smooth = 0.5,
                               surface_resolution = NULL,
                               quiet = TRUE) {
  nm_path <- resolve_noisemodelling_path(nm_path, nm_version, download_nm, quiet)
  wps <- noisemodelling_wps_script(nm_path)
  check_noisemodelling_java(wps, java)

  if (is.null(out_dir)) out_dir <- tempdir()
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  work_dir <- tempfile("gloBFPr_noisemodelling_", tmpdir = out_dir)
  dir.create(work_dir, recursive = TRUE)
  if (!isTRUE(keep_files)) {
    on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  }
  input_dir <- file.path(work_dir, "input")
  output_dir <- file.path(work_dir, "output")
  dir.create(input_dir)
  dir.create(output_dir)

  srid <- noisemodelling_srid(inputs$buildings)
  paths <- write_noisemodelling_shapefiles(inputs, input_dir, quiet = quiet)
  db_name <- "h2gisdb"

  import_script <- nm_wps_file(nm_path, "Import_and_Export", "Import_File.groovy")
  import_osm_script <- nm_wps_file(nm_path, "Import_and_Export", "Import_OSM.groovy")
  delaunay_script <- nm_wps_file(nm_path, "Receivers", "Delaunay_Grid.groovy")
  emission_script <- nm_wps_file(nm_path, "NoiseModelling", "Road_Emission_from_Traffic.groovy")
  emission_script <- patch_road_emission_wps_script(emission_script, work_dir)
  noise_script <- nm_wps_file(nm_path, "NoiseModelling", "Noise_level_from_source.groovy")
  isosurface_script <- nm_wps_file(nm_path, "Acoustic_Tools", "Create_Isosurface.groovy")
  export_script <- nm_wps_file(nm_path, "Import_and_Export", "Export_Table.groovy")

  logs <- list()
  import_paths <- paths[setdiff(names(paths), "RECEIVERS")]
  for (table in names(import_paths)) {
    logs[[paste0("import_", table)]] <- run_wps_script(
      wps,
      work_dir,
      db_name,
      import_script,
      c("-pathFile", import_paths[[table]], "-inputSRID", as.character(srid), "-tableName", table),
      java = java,
      quiet = quiet
    )
  }
  if (!is.null(osm_file)) {
    if (!file.exists(osm_file)) {
      stop("`osm_file` was not found: ", osm_file, call. = FALSE)
    }
    logs$import_osm_roads <- run_wps_script(
      wps,
      work_dir,
      db_name,
      import_osm_script,
      c(
        "-pathFile", normalizePath(osm_file, mustWork = TRUE),
        "-targetSRID", as.character(srid),
        "-ignoreBuilding", "true",
        "-ignoreGround", "true",
        "-ignoreRoads", "false",
        "-removeTunnels", noise_wps_value(osm_remove_tunnels),
        "-eliminateNoTrafficRoads", noise_wps_value(osm_eliminate_no_traffic_roads)
      ),
      java = java,
      quiet = quiet
    )
  }

  logs$delaunay <- run_wps_script(
    wps,
    work_dir,
    db_name,
    delaunay_script,
    c(
      "-tableBuilding", "BUILDINGS",
      "-sourcesTableName", "ROADS",
      "-outputTableName", "RECEIVERS",
      "-maxArea", as.character(delaunay_max_area),
      "-maxCellDist", as.character(max_src_distance),
      "-roadWidth", as.character(road_width),
      "-buildingBuffer", as.character(building_buffer),
      "-height", "4"
    ),
    java = java,
    quiet = quiet
  )

  logs$emission <- run_wps_script(
    wps,
    work_dir,
    db_name,
    emission_script,
    c("-tableRoads", "ROADS"),
    java = java,
    quiet = quiet
  )

  noise_args <- noise_wps_args_vector(
    tableBuilding = "BUILDINGS",
    tableSources = "LW_ROADS",
    tableReceivers = "RECEIVERS",
    tableGroundAbs = "GROUND",
    paramWallAlpha = wall_alpha,
    confReflOrder = reflection_order,
    confMaxSrcDist = max_src_distance,
    confMaxReflDist = max_reflection_distance,
    confThreadNumber = thread_count,
    confDiffVertical = diffraction_vertical,
    confDiffHorizontal = diffraction_horizontal,
    confExportSourceId = export_source_id,
    confHumidity = humidity,
    confTemperature = temperature,
    confFavourableOccurrencesDefault = format_favourable_occurrences(favourable_occurrences),
    confRaysName = rays_name,
    confMaxError = max_error,
    frequencyFieldPrepend = frequency_field_prepend
  )
  if ("DEM" %in% names(paths)) {
    noise_args <- c(noise_args, noise_wps_args_vector(tableDEM = "DEM"))
  }
  if (!is.null(noise_wps_args)) {
    noise_args <- c(noise_args, noise_wps_args_vector_from_list(noise_wps_args))
  }
  logs$noise <- run_wps_script(
    wps,
    work_dir,
    db_name,
    noise_script,
    noise_args,
    java = java,
    quiet = quiet
  )

  logs$isophones <- run_wps_script(
    wps,
    work_dir,
    db_name,
    isosurface_script,
    c(
      "-resultTable", "RECEIVERS_LEVEL",
      "-isoClass", paste(as.numeric(iso_levels), collapse = ","),
      "-resultTableField", as.character(iso_field[1]),
      "-smoothCoefficient", as.character(iso_smooth)
    ),
    java = java,
    quiet = quiet
  )

  output_file <- file.path(output_dir, "receivers_level.geojson")
  receivers_file <- file.path(output_dir, "receivers.geojson")
  emissions_file <- file.path(output_dir, "lw_roads.geojson")
  roads_file <- file.path(output_dir, "roads.geojson")
  isophones_file <- file.path(output_dir, "contouring_noise_map.geojson")
  logs$export_receivers_level <- run_wps_script(
    wps,
    work_dir,
    db_name,
    export_script,
    c("-exportPath", output_file, "-tableToExport", "RECEIVERS_LEVEL"),
    java = java,
    quiet = quiet
  )
  logs$export_receivers <- run_wps_script(
    wps,
    work_dir,
    db_name,
    export_script,
    c("-exportPath", receivers_file, "-tableToExport", "RECEIVERS"),
    java = java,
    quiet = quiet
  )
  logs$export_emissions <- run_wps_script(
    wps,
    work_dir,
    db_name,
    export_script,
    c("-exportPath", emissions_file, "-tableToExport", "LW_ROADS"),
    java = java,
    quiet = quiet
  )
  logs$export_roads <- run_wps_script(
    wps,
    work_dir,
    db_name,
    export_script,
    c("-exportPath", roads_file, "-tableToExport", "ROADS"),
    java = java,
    quiet = quiet
  )
  logs$export_isophones <- run_wps_script(
    wps,
    work_dir,
    db_name,
    export_script,
    c("-exportPath", isophones_file, "-tableToExport", "CONTOURING_NOISE_MAP"),
    java = java,
    quiet = quiet
  )
  if (!file.exists(output_file)) {
    stop("NoiseModelling finished without producing `RECEIVERS_LEVEL` output.", call. = FALSE)
  }
  receivers_level <- sf::st_read(output_file, quiet = TRUE)
  receivers <- if (file.exists(receivers_file)) sf::st_read(receivers_file, quiet = TRUE) else inputs$receivers
  emissions <- if (file.exists(emissions_file)) sf::st_read(emissions_file, quiet = TRUE) else NULL
  roads <- if (file.exists(roads_file)) sf::st_read(roads_file, quiet = TRUE) else inputs$roads
  isophones <- if (file.exists(isophones_file)) sf::st_read(isophones_file, quiet = TRUE) else NULL
  run_inputs <- inputs
  run_inputs$receivers <- receivers
  run_inputs$roads <- roads
  noise_map <- receivers_level_to_noise_map(receivers_level, receivers)
  noise_surface <- noise_map_to_noise_surface(
    noise_map,
    inputs = run_inputs,
    period = "DEN",
    field = "LAEQ",
    resolution = surface_resolution
  )
  result <- list(
    inputs = run_inputs,
    receivers_level = receivers_level,
    noise_map = noise_map,
    noise_surface = noise_surface,
    emissions = emissions,
    isophones = isophones,
    output_file = output_file,
    receivers_file = receivers_file,
    emissions_file = emissions_file,
    roads_file = roads_file,
    isophones_file = isophones_file,
    logs = logs,
    noisemodelling_path = nm_path
  )
  if (isTRUE(keep_files)) {
    result$work_dir <- work_dir
  }
  result
}

#' @noRd
resolve_noise_raster_inputs <- function(buildings,
                                        roads,
                                        greenspace,
                                        canopy_height,
                                        dem,
                                        datasource_canopy_height,
                                        datasource_greenspace,
                                        greenspace_year,
                                        greenspace_zoom,
                                        opentopo_key,
                                        min_tree_height,
                                        quiet = TRUE) {
  datasource_canopy_height <- normalize_shadow_source(
    datasource_canopy_height,
    choices = c("metachm", "ethchm")
  )
  datasource_greenspace <- normalize_shadow_source(
    datasource_greenspace,
    choices = c("esri", "sentinel2")
  )
  needs_canopy <- is.null(canopy_height) && !is.null(datasource_canopy_height)
  needs_greenspace <- is.null(greenspace) && !is.null(datasource_greenspace)
  needs_dem <- is.null(dem) && !is.null(opentopo_key)
  if (!needs_canopy && !needs_greenspace && !needs_dem) {
    return(list(canopy_height = canopy_height, greenspace = greenspace, dem = dem))
  }
  analysis_extent <- if (is.null(roads)) {
    sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(buildings)))
  } else {
    make_noise_analysis_extent(buildings, sf::st_transform(roads, sf::st_crs(buildings)))
  }
  bbox_vector <- as_wgs84_bbox_vector(analysis_extent)
  fetched <- list(canopy_height = canopy_height, greenspace = greenspace, dem = dem)
  if (needs_dem) {
    if (!quiet) cli::cli_alert_info("Downloading DEM for noise analysis ...")
    fetched$dem <- get_dem(bbox_vector, opentopo_key)
  }
  if (needs_canopy) {
    if (!quiet) cli::cli_alert_info("Downloading canopy height map for noise analysis ...")
    chm_layers <- suppressMessages(get_chm(
      bbox_vector,
      min_tree_height,
      datasource = datasource_canopy_height
    ))
    fetched$canopy_height <- chm_layers[[1]]
  }
  if (needs_greenspace) {
    if (!quiet) cli::cli_alert_info("Downloading greenspace map for noise analysis ...")
    fetched$greenspace <- get_greenspace(
      bbox = bbox_vector,
      buffer = NULL,
      type = datasource_greenspace,
      zoom = greenspace_zoom,
      year = greenspace_year,
      min_tree_height = min_tree_height
    )
  }
  fetched
}

#' @noRd
receivers_level_to_noise_map <- function(receivers_level, receivers) {
  if (is.null(receivers_level) || nrow(receivers_level) == 0 || is.null(receivers)) {
    return(NULL)
  }
  level <- if (inherits(receivers_level, "sf")) {
    sf::st_drop_geometry(receivers_level)
  } else {
    as.data.frame(receivers_level)
  }
  if (!all(c("IDRECEIVER", "PERIOD") %in% names(level))) {
    return(receivers)
  }
  value_cols <- intersect(
    c("LAEQ", "LEQ", "HZ63", "HZ125", "HZ250", "HZ500", "HZ1000", "HZ2000", "HZ4000", "HZ8000"),
    names(level)
  )
  if (length(value_cols) == 0) {
    return(receivers)
  }
  level$IDRECEIVER <- as.integer(level$IDRECEIVER)
  level$PERIOD <- toupper(as.character(level$PERIOD))
  level <- mark_no_result_receiver_levels(level, value_cols)
  level <- level[!is.na(level$IDRECEIVER) & !is.na(level$PERIOD) & nzchar(level$PERIOD), , drop = FALSE]
  if (nrow(level) == 0) {
    return(receivers)
  }

  wide <- unique(level["IDRECEIVER"])
  for (field in value_cols) {
    field_wide <- stats::reshape(
      level[, c("IDRECEIVER", "PERIOD", field), drop = FALSE],
      idvar = "IDRECEIVER",
      timevar = "PERIOD",
      direction = "wide"
    )
    names(field_wide) <- sub(paste0("^", field, "\\."), paste0(field, "_"), names(field_wide))
    wide <- merge(wide, field_wide, by = "IDRECEIVER", all = TRUE, sort = FALSE)
  }

  out <- receivers
  if (!"PK" %in% names(out)) {
    out$PK <- seq_len(nrow(out))
  }
  merge(out, wide, by.x = "PK", by.y = "IDRECEIVER", all.x = TRUE, sort = FALSE)
}

#' @noRd
mark_no_result_receiver_levels <- function(level, value_cols) {
  octave_cols <- intersect(c("HZ63", "HZ125", "HZ250", "HZ500", "HZ1000", "HZ2000", "HZ4000", "HZ8000"), names(level))
  if (length(octave_cols) == 0) {
    return(level)
  }
  octave_values <- as.data.frame(lapply(level[octave_cols], as.numeric))
  no_result <- Reduce(`&`, lapply(octave_values, function(x) is.na(x) | x <= -98))
  if (any(no_result, na.rm = TRUE)) {
    level[no_result, value_cols] <- NA_real_
  }
  level
}

#' @noRd
noise_map_to_noise_surface <- function(noise_map,
                                       inputs = NULL,
                                       period = "DEN",
                                       field = "LAEQ",
                                       resolution = NULL,
                                       breaks = c(35, 40, 45, 50, 55, 60, 65, 70, 75),
                                       max_distance = NULL,
                                       nodata = -99,
                                       idw_power = 2) {
  if (is.null(noise_map) || nrow(noise_map) == 0 || !inherits(noise_map, "sf")) {
    return(NULL)
  }
  period <- toupper(period[1])
  field <- toupper(field[1])
  value_col <- paste(field, period, sep = "_")
  if (!value_col %in% names(noise_map)) {
    return(NULL)
  }
  values <- as.numeric(noise_map[[value_col]])
  values[!is.finite(values)] <- NA_real_
  if (!is.null(nodata)) {
    values[values <= nodata] <- NA_real_
  }
  keep <- which(!is.na(values))
  if (length(keep) < 3) {
    return(NULL)
  }
  points <- noise_map[keep, ]
  coords <- sf::st_coordinates(points)[, c("X", "Y"), drop = FALSE]
  values <- values[keep]

  if (is.null(resolution) || !is.finite(resolution) || resolution <= 0) {
    resolution <- estimate_noise_surface_resolution(coords)
  }
  bbox_source <- if (!is.null(inputs) && !is.null(inputs$receivers)) inputs$receivers else points
  bbox <- sf::st_bbox(bbox_source)
  raster <- terra::rast(
    xmin = bbox[["xmin"]],
    xmax = bbox[["xmax"]],
    ymin = bbox[["ymin"]],
    ymax = bbox[["ymax"]],
    resolution = resolution,
    crs = sf::st_crs(points)$wkt
  )
  xy <- terra::xyFromCell(raster, seq_len(terra::ncell(raster)))
  raster_values <- idw_noise_values(
    xy = xy,
    points = coords,
    values = values,
    power = idw_power,
    max_distance = max_distance
  )
  terra::values(raster) <- raster_values
  names(raster) <- value_col

  bands <- noise_surface_bands(raster, breaks = breaks)
  structure(
    list(
      raster = raster,
      bands = bands,
      breaks = breaks,
      field = field,
      period = period,
      value_col = value_col,
      resolution = resolution
    ),
    class = "gloBFPr_noise_surface"
  )
}

#' @noRd
estimate_noise_surface_resolution <- function(coords) {
  dx <- sort(unique(coords[, 1]))
  dy <- sort(unique(coords[, 2]))
  steps <- c(diff(dx), diff(dy))
  steps <- steps[is.finite(steps) & steps > 0]
  if (length(steps) == 0) {
    return(10)
  }
  as.numeric(stats::median(steps))
}

#' @noRd
idw_noise_values <- function(xy,
                             points,
                             values,
                             power = 2,
                             max_distance = NULL,
                             chunk_size = 5000) {
  out <- rep(NA_real_, nrow(xy))
  chunks <- split(seq_len(nrow(xy)), ceiling(seq_len(nrow(xy)) / chunk_size))
  for (chunk in chunks) {
    dx <- outer(xy[chunk, 1], points[, 1], "-")
    dy <- outer(xy[chunk, 2], points[, 2], "-")
    dist <- sqrt(dx^2 + dy^2)
    exact <- dist == 0
    if (!is.null(max_distance) && is.finite(max_distance)) {
      dist[dist > max_distance] <- NA_real_
    }
    weights <- 1 / (dist^power)
    weights[!is.finite(weights)] <- NA_real_
    weighted_sum <- rowSums(weights * rep(values, each = nrow(weights)), na.rm = TRUE)
    weight_sum <- rowSums(weights, na.rm = TRUE)
    vals <- weighted_sum / weight_sum
    vals[weight_sum == 0] <- NA_real_
    exact_rows <- which(rowSums(exact) > 0)
    if (length(exact_rows) > 0) {
      exact_cols <- max.col(exact[exact_rows, , drop = FALSE], ties.method = "first")
      vals[exact_rows] <- values[exact_cols]
    }
    out[chunk] <- vals
  }
  out
}

#' @noRd
noise_surface_bands <- function(raster, breaks) {
  breaks <- sort(unique(as.numeric(breaks)))
  rcl <- rbind(
    c(-Inf, breaks[1], 1),
    cbind(breaks[-length(breaks)], breaks[-1], seq_along(breaks[-1]) + 1),
    c(breaks[length(breaks)], Inf, length(breaks) + 1)
  )
  classified <- terra::classify(raster, rcl = rcl, include.lowest = TRUE)
  names(classified) <- "noise_class"
  bands <- suppressWarnings(sf::st_as_sf(terra::as.polygons(
    classified,
    dissolve = TRUE,
    na.rm = TRUE
  )))
  if (nrow(bands) == 0) {
    return(bands)
  }
  bands$noise_class <- as.integer(bands$noise_class)
  lower <- c(-Inf, breaks)
  upper <- c(breaks, Inf)
  bands$level_min <- lower[bands$noise_class]
  bands$level_max <- upper[bands$noise_class]
  bands$label <- noise_band_labels(bands$level_min, bands$level_max)
  bands
}

#' @noRd
noise_band_labels <- function(lower, upper) {
  vapply(seq_along(lower), function(i) {
    if (!is.finite(lower[i])) return(paste0("< ", upper[i], " dB(A)"))
    if (!is.finite(upper[i])) return(paste0(">= ", lower[i], " dB(A)"))
    paste0(lower[i], "-", upper[i], " dB(A)")
  }, character(1))
}

#' Plot a NoiseModelling-style road-noise map
#'
#' @description
#' Draws NoiseModelling isosurface polygons underneath roads and buildings when
#' available. If `x` does not contain official isosurfaces, it falls back to an
#' interpolated receiver surface.
#'
#' @param x A result from `get_noise_map(run = TRUE)`, a `gloBFPr_noise_surface`
#' object, or an `sf` receiver noise map.
#' @param period Noise period to plot. Defaults to `"DEN"`.
#' @param field Noise field to plot. Defaults to `"LAEQ"`.
#' @param resolution Optional raster resolution in map units when `x` is an
#' `sf` receiver map.
#' @param breaks Noise class breakpoints in dB.
#' @param nodata Values at or below this threshold are treated as no-data.
#' Defaults to `-99`, NoiseModelling's no-result sentinel value.
#' @param palette Fill colors from quiet to loud.
#' @param road_col Road overlay color. Defaults to white with transparency.
#' Use `NA` to omit roads.
#' @param road_alpha Road overlay alpha from `0` fully transparent to `1`
#' opaque.
#' @param road_lwd Road overlay line width.
#' @param building_col Building fill color. Use `NA` to omit buildings.
#' @param legend Logical. Draw the dB(A) class legend in a dedicated panel to
#' the right of the map. Defaults to `TRUE`.
#' @param legend_width Width of the legend panel relative to the map panel.
#' Smaller values give the map more room. Defaults to `0.26`.
#' @param legend_cex Legend text size. Defaults to `0.85`.
#' @param scalebar Logical. Draw a distance scale bar just below the legend
#' (or in the bottom-left of the map when `legend = FALSE`). Defaults to
#' `TRUE`. Distances assume a projected CRS in metres; for geographic
#' coordinates they are approximated at the map's mid-latitude.
#' @param scalebar_unit Scale bar unit: `"auto"` (default), `"km"`, or `"m"`
#' to pick whichever keeps the label readable.
#' @param scalebar_cex Scale bar label size. Defaults to `0.7`.
#' @param north_arrow Logical. Draw a north arrow in the lower-left map area.
#'   Defaults to `TRUE`.
#' @param mar Margins (in lines) around the map panel. Defaults to a tight
#' margin so the map fills the device.
#' @param add Logical. If `TRUE`, add to the current plot.
#' @param ... Additional arguments passed to `plot()`.
#'
#' @return Invisibly returns the isosurface `sf` object or fallback
#' `gloBFPr_noise_surface` object used for plotting.
#' @export
plot_noise_map <- function(x,
                           period = "DEN",
                           field = "LAEQ",
                           resolution = NULL,
                           breaks = c(35, 40, 45, 50, 55, 60, 65, 70, 75),
                           nodata = -99,
                           palette = noise_map_palette(),
                           road_col = "white",
                           road_alpha = 0.35,
                           road_lwd = 1.4,
                           building_col = "black",
                           legend = TRUE,
                           legend_width = 0.26,
                           legend_cex = 0.85,
                           scalebar = TRUE,
                           scalebar_unit = c("auto", "km", "m"),
                           scalebar_cex = 0.7,
                           north_arrow = TRUE,
                           mar = c(0.2, 0.2, 0.2, 0.2),
                           add = FALSE,
                           ...) {
  scalebar_unit <- match.arg(scalebar_unit)
  inputs <- if (is.list(x) && !is.null(x$inputs)) x$inputs else NULL
  if (is.list(x) && !is.null(x$isophones) && inherits(x$isophones, "sf") && nrow(x$isophones) > 0) {
    plot_noise_isophones_layers(
      x$isophones,
      inputs = inputs,
      palette = palette,
      road_col = road_col,
      road_alpha = road_alpha,
      road_lwd = road_lwd,
      building_col = building_col,
      legend = legend,
      legend_width = legend_width,
      legend_cex = legend_cex,
      scalebar = scalebar,
      scalebar_unit = scalebar_unit,
      scalebar_cex = scalebar_cex,
      north_arrow = north_arrow,
      mar = mar,
      add = add,
      ...
    )
    return(invisible(x$isophones))
  }
  surface <- if (inherits(x, "gloBFPr_noise_surface")) {
    x
  } else if (is.list(x) && !is.null(x$noise_surface) && inherits(x$noise_surface, "gloBFPr_noise_surface") &&
             identical(toupper(period[1]), x$noise_surface$period) &&
             identical(toupper(field[1]), x$noise_surface$field)) {
    x$noise_surface
  } else {
    noise_map <- if (inherits(x, "sf")) x else x$noise_map
    noise_map_to_noise_surface(
      noise_map,
      inputs = inputs,
      period = period,
      field = field,
      resolution = resolution,
      breaks = breaks,
      nodata = nodata
    )
  }
  if (is.null(surface)) {
    stop("No plottable noise surface was available. Check that receiver levels contain finite values above `nodata`.", call. = FALSE)
  }
  plot_noise_surface_layers(
    surface,
    inputs = inputs,
    palette = palette,
    road_col = road_col,
    road_alpha = road_alpha,
    road_lwd = road_lwd,
    building_col = building_col,
    legend = legend,
    legend_width = legend_width,
    legend_cex = legend_cex,
    scalebar = scalebar,
    scalebar_unit = scalebar_unit,
    scalebar_cex = scalebar_cex,
    north_arrow = north_arrow,
    mar = mar,
    add = add,
    ...
  )
  invisible(surface)
}

#' @noRd
noise_map_palette <- function() {
  c(
    "#b8d6d3", "#c9e2d8", "#e8f3c4", "#fff1b8", "#f9cf88",
    "#f28b55", "#df4c4b", "#c51f63", "#8d006f", "#4b1f59",
    "#35164d"
  )
}

#' @noRd
noise_legend_labels <- function(labels) {
  labels <- sub("\\s*dB\\s*\\(A\\)\\s*$", "", labels)
  trimws(labels)
}

#' @noRd
noise_bbox_aspect <- function(geometry) {
  bb <- try(sf::st_bbox(geometry), silent = TRUE)
  if (inherits(bb, "try-error")) return(NULL)
  w <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
  h <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0) return(NULL)
  w / h
}

#' @noRd
noise_map_open_layout <- function(legend = TRUE,
                                  legend_width = 0.26,
                                  mar = c(0.2, 0.2, 0.2, 0.2),
                                  asp_ratio = NULL) {
  old_par <- graphics::par(no.readonly = TRUE)
  if (isTRUE(legend)) {
    din <- graphics::par("din")
    line_in <- graphics::par("csi")
    legend_in <- min(max(legend_width, 0.05) * din[1], din[1] * 0.6)
    map_in <- din[1] - legend_in
    if (!is.null(asp_ratio) && is.finite(asp_ratio) && asp_ratio > 0) {
      # With asp = 1 the map is usually height-limited, so give the map panel
      # only the width it can actually fill and hand the slack to the legend.
      map_h <- din[2] - (mar[1] + mar[3]) * line_in
      needed <- map_h * asp_ratio + (mar[2] + mar[4]) * line_in
      map_in <- min(map_in, max(needed, din[1] * 0.3))
    }
    graphics::layout(
      matrix(c(1L, 2L), nrow = 1L),
      widths = c(map_in, max(din[1] - map_in, din[1] * 0.05))
    )
  }
  graphics::par(mar = mar)
  old_par
}

#' @noRd
noise_map_close_layout <- function(old_par, legend = TRUE) {
  if (isTRUE(legend)) {
    graphics::layout(1L)
  }
  suppressWarnings(graphics::par(old_par))
  invisible(NULL)
}

#' @noRd
noise_map_scale <- function(geometry = NULL) {
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  span <- usr[2] - usr[1]
  if (!is.finite(span) || span <= 0 || !is.finite(pin[1]) || pin[1] <= 0) {
    return(NULL)
  }
  metres_per_inch <- span / pin[1]
  longlat <- FALSE
  if (!is.null(geometry)) {
    ll <- suppressWarnings(try(sf::st_is_longlat(geometry), silent = TRUE))
    longlat <- isTRUE(!inherits(ll, "try-error") && isTRUE(ll))
  }
  if (longlat) {
    # Degrees of longitude shrink with latitude; scale at the map's mid-latitude.
    mid_lat <- (usr[3] + usr[4]) / 2
    metres_per_inch <- metres_per_inch * 111320 * cos(mid_lat * pi / 180)
  }
  if (!is.finite(metres_per_inch) || metres_per_inch <= 0) {
    return(NULL)
  }
  metres_per_inch
}

#' @noRd
noise_scalebar_nice <- function(x) {
  if (!is.finite(x) || x <= 0) return(NA_real_)
  pow <- 10^floor(log10(x))
  frac <- x / pow
  nice <- if (frac >= 5) 5 else if (frac >= 2.5) 2.5 else if (frac >= 2) 2 else 1
  nice * pow
}

#' @noRd
noise_map_draw_scalebar <- function(metres_per_inch,
                                    x_in,
                                    y_in,
                                    max_width_in,
                                    unit = "km",
                                    cex = 0.7,
                                    col = "#333333") {
  if (is.null(metres_per_inch) || !is.finite(metres_per_inch) || metres_per_inch <= 0) {
    return(invisible(NULL))
  }
  if (!is.finite(max_width_in) || max_width_in <= 0.2) {
    return(invisible(NULL))
  }
  metres <- noise_scalebar_nice(max_width_in * 0.9 * metres_per_inch)
  if (!is.finite(metres) || metres <= 0) return(invisible(NULL))
  if (identical(unit, "auto")) unit <- if (metres >= 1000) "km" else "m"
  divisor <- if (identical(unit, "km")) 1000 else 1
  value <- metres / divisor
  bar_in <- metres / metres_per_inch
  height_in <- 0.07
  x0 <- graphics::grconvertX(x_in, "inches", "user")
  xm <- graphics::grconvertX(x_in + bar_in / 2, "inches", "user")
  x1 <- graphics::grconvertX(x_in + bar_in, "inches", "user")
  y0 <- graphics::grconvertY(y_in, "inches", "user")
  y1 <- graphics::grconvertY(y_in + height_in, "inches", "user")
  graphics::rect(x0, y0, xm, y1, col = col, border = col, xpd = NA)
  graphics::rect(xm, y0, x1, y1, col = "white", border = col, xpd = NA)
  label <- paste(format(value, trim = TRUE, scientific = FALSE), unit)
  graphics::text(x0, y0, "0", adj = c(0.5, 1.35), cex = cex, col = col, xpd = NA)
  graphics::text(x1, y0, label, adj = c(0.5, 1.35), cex = cex, col = col, xpd = NA)
  invisible(list(metres = metres, unit = unit, width_in = bar_in, height_in = height_in))
}

#' @noRd
noise_map_draw_north_arrow <- function(x_in,
                                       y_in,
                                       height_in = 0.28,
                                       cex = 0.7,
                                       col = "#333333") {
  x <- graphics::grconvertX(x_in, "inches", "user")
  y0 <- graphics::grconvertY(y_in, "inches", "user")
  y1 <- graphics::grconvertY(y_in + height_in, "inches", "user")
  graphics::arrows(
    x0 = x, y0 = y0, x1 = x, y1 = y1,
    length = 0.07, lwd = 1.1, col = col, xpd = NA
  )
  graphics::text(
    x, y1, "N",
    font = 2, cex = cex, col = col, adj = c(0.5, -0.25), xpd = NA
  )
  invisible(NULL)
}

#' @noRd
noise_map_scalebar_below_legend <- function(metres_per_inch,
                                            legend_info,
                                            unit = "km",
                                            cex = 0.7,
                                            north_arrow = FALSE) {
  if (is.null(metres_per_inch)) return(invisible(NULL))
  rect <- if (is.list(legend_info)) legend_info$rect else NULL
  if (is.null(rect)) return(invisible(NULL))
  left_in <- graphics::grconvertX(rect$left, "user", "inches")
  bottom_in <- graphics::grconvertY(rect$top - rect$h, "user", "inches") - 0.3
  right_in <- graphics::grconvertX(1, "npc", "inches")
  legend_w_in <- graphics::grconvertX(rect$left + rect$w, "user", "inches") - left_in
  arrow_w_in <- if (isTRUE(north_arrow)) 0.28 else 0
  gap_in <- if (isTRUE(north_arrow)) 0.12 else 0
  # Keep the bar roughly as wide as the legend box rather than the whole panel.
  max_width_in <- min(right_in - left_in - arrow_w_in - gap_in, max(legend_w_in, 1))
  if (isTRUE(north_arrow)) {
    noise_map_draw_north_arrow(
      x_in = left_in + arrow_w_in / 2,
      y_in = bottom_in,
      height_in = 0.26,
      cex = cex
    )
  }
  noise_map_draw_scalebar(
    metres_per_inch,
    x_in = left_in + arrow_w_in + gap_in,
    y_in = bottom_in,
    max_width_in = max_width_in,
    unit = unit,
    cex = cex
  )
}

#' @noRd
noise_map_scalebar_in_map <- function(metres_per_inch,
                                      unit = "km",
                                      cex = 0.7,
                                      north_arrow = FALSE) {
  if (is.null(metres_per_inch)) return(invisible(NULL))
  pin <- graphics::par("pin")
  left_in <- graphics::grconvertX(0.04, "npc", "inches")
  bottom_in <- graphics::grconvertY(0.08, "npc", "inches")
  arrow_w_in <- if (isTRUE(north_arrow)) 0.28 else 0
  gap_in <- if (isTRUE(north_arrow)) 0.12 else 0
  if (isTRUE(north_arrow)) {
    noise_map_draw_north_arrow(
      x_in = left_in + arrow_w_in / 2,
      y_in = bottom_in,
      height_in = 0.26,
      cex = cex
    )
  }
  noise_map_draw_scalebar(
    metres_per_inch,
    x_in = left_in + arrow_w_in + gap_in,
    y_in = bottom_in,
    max_width_in = pin[1] * 0.3 - arrow_w_in - gap_in,
    unit = unit,
    cex = cex
  )
}

#' @noRd
noise_map_north_arrow_in_map <- function(cex = 0.7, col = "#333333") {
  usr <- graphics::par("usr")
  if (!all(is.finite(usr))) return(invisible(NULL))
  w <- usr[2] - usr[1]
  h <- usr[4] - usr[3]
  if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0) {
    return(invisible(NULL))
  }

  x <- usr[1] + 0.075 * w
  y0 <- usr[3] + 0.120 * h
  y1 <- usr[3] + 0.205 * h
  pad_x <- 0.035 * w
  graphics::rect(
    x - pad_x, y0 - 0.035 * h,
    x + pad_x, y1 + 0.055 * h,
    col = grDevices::adjustcolor("white", alpha.f = 0.82),
    border = NA,
    xpd = NA
  )
  graphics::arrows(
    x0 = x, y0 = y0, x1 = x, y1 = y1,
    length = 0.1, lwd = 1.3, col = col, xpd = NA
  )
  graphics::text(
    x, y1 + 0.035 * h, "N",
    font = 2, cex = cex * 1.15, col = col, xpd = NA
  )
  invisible(NULL)
}

#' @noRd
noise_map_draw_legend <- function(labels,
                                  fills,
                                  title = "LAeq dB(A)",
                                  cex = 0.85,
                                  mar = c(0.2, 0.2, 0.2, 0.2),
                                  reserve_in = 0) {
  graphics::par(mar = c(mar[1], max(mar[2], 0.8), mar[3], 0.2))
  graphics::plot.new()
  pin <- graphics::par("pin")
  # Nudge the legend up so legend + scale bar stay centred on the map.
  y_centre <- if (reserve_in > 0 && pin[2] > 0) 0.5 + (reserve_in / 2) / pin[2] else 0.5
  info <- graphics::legend(
    x = 0,
    y = y_centre,
    yjust = 0.5,
    xjust = 0,
    legend = noise_legend_labels(labels),
    fill = fills,
    border = NA,
    bty = "n",
    title = title,
    title.adj = 0,
    cex = cex,
    y.intersp = 1.15,
    xpd = NA
  )
  invisible(info)
}

#' @noRd
plot_noise_surface_layers <- function(surface,
                                      inputs = NULL,
                                      palette = noise_map_palette(),
                                      road_col = "white",
                                      road_alpha = 0.35,
                                      road_lwd = 1.4,
                                      building_col = "black",
                                      legend = TRUE,
                                      legend_width = 0.26,
                                      legend_cex = 0.85,
                                      scalebar = FALSE,
                                      scalebar_unit = "km",
                                      scalebar_cex = 0.7,
                                      north_arrow = TRUE,
                                      mar = c(0.2, 0.2, 0.2, 0.2),
                                      add = FALSE,
                                      ...) {
  bands <- surface$bands
  if (is.null(bands) || nrow(bands) == 0) {
    stop("The noise surface has no contour bands to plot.", call. = FALSE)
  }
  class_values <- bands$noise_class
  fill <- palette[pmin(pmax(class_values, 1), length(palette))]
  if (!isTRUE(add)) {
    old_par <- noise_map_open_layout(
      legend = legend,
      legend_width = legend_width,
      mar = mar,
      asp_ratio = noise_bbox_aspect(bands)
    )
    on.exit(noise_map_close_layout(old_par, legend = legend), add = TRUE)
    graphics::plot(sf::st_geometry(bands), col = fill, border = NA, asp = 1, ...)
  } else {
    graphics::plot(sf::st_geometry(bands), col = fill, border = NA, add = TRUE, ...)
  }
  if (!is.null(inputs) && !is.null(inputs$roads) && !is.na(road_col) && road_alpha > 0) {
    graphics::plot(
      sf::st_geometry(inputs$roads),
      col = grDevices::adjustcolor(road_col, alpha.f = road_alpha),
      lwd = road_lwd,
      add = TRUE
    )
  }
  if (!is.null(inputs) && !is.null(inputs$buildings) && !is.na(building_col)) {
    graphics::plot(
      sf::st_geometry(inputs$buildings),
      col = building_col,
      border = NA,
      add = TRUE
    )
  }
  map_scale <- if (isTRUE(scalebar)) noise_map_scale(bands) else NULL
  if (!isTRUE(legend) || isTRUE(add)) {
    if (isTRUE(scalebar)) {
      noise_map_scalebar_in_map(
        map_scale,
        unit = scalebar_unit,
        cex = scalebar_cex,
        north_arrow = north_arrow
      )
    } else if (isTRUE(north_arrow)) {
      noise_map_north_arrow_in_map(cex = scalebar_cex)
    }
    return(invisible(NULL))
  }
  legend_labels <- unique(bands$label[order(bands$noise_class)])
  legend_classes <- unique(bands$noise_class[order(bands$noise_class)])
  legend_info <- noise_map_draw_legend(
    labels = legend_labels,
    fills = palette[pmin(pmax(legend_classes, 1), length(palette))],
    cex = legend_cex,
    mar = mar,
    reserve_in = if (isTRUE(scalebar)) 0.5 else 0
  )
  if (isTRUE(scalebar)) {
    noise_map_scalebar_below_legend(
      map_scale,
      legend_info,
      unit = scalebar_unit,
      cex = scalebar_cex,
      north_arrow = north_arrow
    )
  } else if (isTRUE(north_arrow)) {
    noise_map_north_arrow_in_map(cex = scalebar_cex)
  }
  invisible(NULL)
}

#' @noRd
plot_noise_isophones_layers <- function(isophones,
                                        inputs = NULL,
                                        palette = noise_map_palette(),
                                        road_col = "white",
                                        road_alpha = 0.35,
                                        road_lwd = 1.4,
                                        building_col = "black",
                                        legend = TRUE,
                                        legend_width = 0.26,
                                        legend_cex = 0.85,
                                        scalebar = FALSE,
                                        scalebar_unit = "km",
                                        scalebar_cex = 0.7,
                                        north_arrow = TRUE,
                                        mar = c(0.2, 0.2, 0.2, 0.2),
                                        add = FALSE,
                                        ...) {
  class_values <- isophone_class_values(isophones)
  fill <- palette[pmin(pmax(class_values, 1), length(palette))]
  if (!isTRUE(add)) {
    old_par <- noise_map_open_layout(
      legend = legend,
      legend_width = legend_width,
      mar = mar,
      asp_ratio = noise_bbox_aspect(isophones)
    )
    on.exit(noise_map_close_layout(old_par, legend = legend), add = TRUE)
    graphics::plot(sf::st_geometry(isophones), col = fill, border = NA, asp = 1, ...)
  } else {
    graphics::plot(sf::st_geometry(isophones), col = fill, border = NA, add = TRUE, ...)
  }
  if (!is.null(inputs) && !is.null(inputs$roads) && !is.na(road_col) && road_alpha > 0) {
    graphics::plot(
      sf::st_geometry(inputs$roads),
      col = grDevices::adjustcolor(road_col, alpha.f = road_alpha),
      lwd = road_lwd,
      add = TRUE
    )
  }
  if (!is.null(inputs) && !is.null(inputs$buildings) && !is.na(building_col)) {
    graphics::plot(
      sf::st_geometry(inputs$buildings),
      col = building_col,
      border = NA,
      add = TRUE
    )
  }
  map_scale <- if (isTRUE(scalebar)) noise_map_scale(isophones) else NULL
  if (!isTRUE(legend) || isTRUE(add)) {
    if (isTRUE(scalebar)) {
      noise_map_scalebar_in_map(
        map_scale,
        unit = scalebar_unit,
        cex = scalebar_cex,
        north_arrow = north_arrow
      )
    } else if (isTRUE(north_arrow)) {
      noise_map_north_arrow_in_map(cex = scalebar_cex)
    }
    return(invisible(NULL))
  }
  legend_info <- isophone_legend_info(isophones, class_values)
  legend_classes <- legend_info$class
  legend_labels <- legend_info$label
  legend_box <- noise_map_draw_legend(
    labels = legend_labels,
    fills = palette[pmin(pmax(legend_classes, 1), length(palette))],
    cex = legend_cex,
    mar = mar,
    reserve_in = if (isTRUE(scalebar)) 0.5 else 0
  )
  if (isTRUE(scalebar)) {
    noise_map_scalebar_below_legend(
      map_scale,
      legend_box,
      unit = scalebar_unit,
      cex = scalebar_cex,
      north_arrow = north_arrow
    )
  } else if (isTRUE(north_arrow)) {
    noise_map_north_arrow_in_map(cex = scalebar_cex)
  }
  invisible(NULL)
}

#' @noRd
isophone_class_values <- function(isophones) {
  candidates <- c("ISOLVL", "ISO_LVL", "ISOLINE", "IDISO", "LEVEL", "DB", "LAEQ")
  field <- candidates[candidates %in% toupper(names(isophones))]
  if (length(field) > 0) {
    actual <- names(isophones)[match(field[1], toupper(names(isophones)))]
    values <- parse_noise_isophone_values(isophones[[actual]])
    if (any(is.finite(values))) {
      unique_values <- sort(unique(values[is.finite(values)]))
      class <- match(values, unique_values)
      class[!is.finite(values)] <- 1L
      return(as.integer(class))
    }
  }
  seq_len(nrow(isophones))
}

#' @noRd
isophone_legend_info <- function(isophones, class_values) {
  candidates <- c("ISOLVL", "ISO_LVL", "ISOLINE", "IDISO", "LEVEL", "DB", "LAEQ")
  field <- candidates[candidates %in% toupper(names(isophones))]
  if (length(field) > 0) {
    actual <- names(isophones)[match(field[1], toupper(names(isophones)))]
    values <- parse_noise_isophone_values(isophones[[actual]])
    if (any(is.finite(values))) {
      unique_values <- sort(unique(values[is.finite(values)]))
      classes <- seq_along(unique_values)
      return(data.frame(
        class = classes,
        value = unique_values,
        label = isophone_labels(unique_values),
        stringsAsFactors = FALSE
      ))
    }
  }
  classes <- sort(unique(class_values))
  data.frame(
    class = classes,
    value = classes,
    label = paste("Class", classes, "dB(A)"),
    stringsAsFactors = FALSE
  )
}

#' @noRd
isophone_labels <- function(values) {
  if (length(values) == 0) return(character())
  values <- sort(values)
  if (all(values == floor(values)) && min(values) >= 1 && max(values) <= 11) {
    breaks <- c(35, 40, 45, 50, 55, 60, 65, 70, 75, 80)
    lower <- c(-Inf, breaks)
    upper <- c(breaks, Inf)
    return(noise_band_labels(lower[values], upper[values]))
  }
  if (all(values == floor(values)) && min(values) >= 0 && max(values) <= 10) {
    breaks <- c(35, 40, 45, 50, 55, 60, 65, 70, 75, 80)
    lower <- c(-Inf, breaks)
    upper <- c(breaks, Inf)
    return(noise_band_labels(lower[values + 1L], upper[values + 1L]))
  }
  noise_band_labels(values, c(values[-1], Inf))
}

#' @noRd
parse_noise_isophone_values <- function(values) {
  numeric_values <- suppressWarnings(as.numeric(values))
  if (any(is.finite(numeric_values))) {
    return(numeric_values)
  }
  text_values <- as.character(values)
  suppressWarnings(as.numeric(sub(".*?(-?[0-9]+(?:\\.[0-9]+)?).*", "\\1", text_values)))
}

#' @noRd
get_osm_roads_for_buildings <- function(x,
                                        overpass_url = "https://overpass-api.de/api/interpreter",
                                        timeout = 180) {
  bbox <- as_wgs84_bbox_vector(sf::st_as_sfc(sf::st_bbox(x)))
  query <- sprintf(
    paste0(
      "[out:json][timeout:%d];",
      "(",
      "way[\"highway\"][\"highway\"!~\"^(footway|path|cycleway|steps|bridleway|corridor|elevator|escape|pedestrian|track)$\"](%f,%f,%f,%f);",
      ");",
      "out tags geom;"
    ),
    timeout,
    bbox[2], bbox[1], bbox[4], bbox[3]
  )
  response <- httr2::request(overpass_url) |>
    httr2::req_body_form(data = query) |>
    httr2::req_perform()
  body <- httr2::resp_body_json(response, simplifyVector = FALSE)
  osm_roads_json_to_sf(body, sf::st_crs(x))
}

#' @noRd
osm_roads_json_to_sf <- function(body, target_crs) {
  elements <- body$elements
  if (is.null(elements) || length(elements) == 0) {
    stop("No OSM roads were found in the building bounding box.", call. = FALSE)
  }
  roads <- lapply(elements, function(element) {
    if (!identical(element$type, "way") || is.null(element$geometry) || length(element$geometry) < 2) {
      return(NULL)
    }
    coords <- do.call(rbind, lapply(element$geometry, function(point) {
      c(point$lon, point$lat)
    }))
    tags <- element$tags
    if (is.null(tags$highway)) return(NULL)
    list(
      osm_id = element$id,
      highway = tags$highway,
      name = if (!is.null(tags$name)) tags$name else NA_character_,
      maxspeed = if (!is.null(tags$maxspeed)) tags$maxspeed else NA_character_,
      lanes = if (!is.null(tags$lanes)) tags$lanes else NA_character_,
      geometry = sf::st_linestring(coords)
    )
  })
  roads <- Filter(Negate(is.null), roads)
  if (length(roads) == 0) {
    stop("No usable OSM road geometries were found in the building bounding box.", call. = FALSE)
  }
  out <- sf::st_sf(
    osm_id = vapply(roads, `[[`, numeric(1), "osm_id"),
    highway = vapply(roads, `[[`, character(1), "highway"),
    name = vapply(roads, `[[`, character(1), "name"),
    maxspeed = vapply(roads, `[[`, character(1), "maxspeed"),
    lanes = vapply(roads, `[[`, character(1), "lanes"),
    geometry = sf::st_sfc(lapply(roads, `[[`, "geometry"), crs = 4326)
  )
  sf::st_transform(out, target_crs)
}

#' @noRd
write_noisemodelling_shapefiles <- function(inputs, input_dir, quiet = TRUE) {
  layers <- list(
    BUILDINGS = inputs$buildings,
    ROADS = inputs$roads,
    GROUND = inputs$ground,
    RECEIVERS = inputs$receivers
  )
  if (!is.null(inputs$dem)) {
    layers$DEM <- dem_to_noisemodelling_points(inputs$dem)
  }
  layers <- Filter(Negate(is.null), layers)
  vapply(names(layers), function(name) {
    path <- file.path(input_dir, paste0(tolower(name), ".geojson"))
    write_noisemodelling_shapefile(layers[[name]], path, quiet = quiet)
    path
  }, character(1))
}

#' @noRd
write_noisemodelling_shapefile <- function(x, path, quiet = TRUE) {
  stem <- tools::file_path_sans_ext(basename(path))
  unlink(list.files(dirname(path), paste0("^", stem, "\\."), full.names = TRUE))
  sf::st_write(x, path, quiet = quiet)
  invisible(path)
}

#' @noRd
dem_to_noisemodelling_points <- function(dem) {
  dem <- dem[[1]]
  values <- terra::values(dem, mat = FALSE)
  keep <- which(!is.na(values))
  xy <- terra::xyFromCell(dem, keep)
  geom <- sf::st_sfc(
    lapply(seq_along(keep), function(i) sf::st_point(c(xy[i, 1], xy[i, 2], values[keep[i]]))),
    crs = terra::crs(dem)
  )
  sf::st_sf(PK = seq_along(keep), geometry = geom)
}

#' @noRd
resolve_noisemodelling_path <- function(nm_path, version, download_nm, quiet = TRUE) {
  if (!is.null(nm_path) && length(nm_path) > 0 && nzchar(nm_path[1])) {
    if (file.exists(noisemodelling_wps_script(nm_path[1]))) {
      return(normalizePath(nm_path[1], mustWork = TRUE))
    }
    stop(
      "NoiseModelling headless runner was not found at `nm_path`: ",
      nm_path[1],
      call. = FALSE
    )
  }
  candidates <- c(
    getOption("gloBFPr.noisemodelling.path", NULL),
    Sys.getenv("NOISEMODELLING_HOME", unset = NA_character_),
    file.path(noisemodelling_cache_dir(), paste0("NoiseModelling_without_gui-", version))
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (candidate in candidates) {
    if (file.exists(noisemodelling_wps_script(candidate))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  if (isTRUE(download_nm)) {
    return(install_noisemodelling(version = version, quiet = quiet))
  }
  stop(
    "NoiseModelling headless runner was not found. ",
    "Supply `nm_path`, set `NOISEMODELLING_HOME`, or run `install_noisemodelling()`.",
    call. = FALSE
  )
}

#' @noRd
noisemodelling_cache_dir <- function() {
  tools::R_user_dir("gloBFPr", which = "cache")
}

#' @noRd
noisemodelling_wps_script <- function(nm_path) {
  file.path(nm_path, "bin", if (.Platform$OS.type == "windows") "wps_scripts.bat" else "wps_scripts")
}

#' @noRd
nm_wps_file <- function(nm_path, folder, file) {
  path <- file.path(nm_path, "noisemodelling", "wps", folder, file)
  if (!file.exists(path)) {
    stop("NoiseModelling WPS script was not found: ", path, call. = FALSE)
  }
  path
}

#' @noRd
noisemodelling_srid <- function(x) {
  epsg <- sf::st_crs(x)$epsg
  if (is.na(epsg)) {
    stop("NoiseModelling execution requires an EPSG-coded projected CRS.", call. = FALSE)
  }
  if (epsg %in% c(4326, 3785)) {
    stop("NoiseModelling execution requires a metric projected CRS, not EPSG:", epsg, ".", call. = FALSE)
  }
  epsg
}

#' @noRd
noise_wps_args_vector <- function(...) {
  noise_wps_args_vector_from_list(list(...))
}

#' @noRd
noise_wps_args_vector_from_list <- function(args) {
  if (is.null(args) || length(args) == 0) return(character())
  if (is.null(names(args)) || any(!nzchar(names(args)))) {
    stop("`noise_wps_args` must be a named list.", call. = FALSE)
  }
  args <- args[!vapply(args, is.null, logical(1))]
  if (length(args) == 0) return(character())
  unlist(Map(function(name, value) {
    value <- noise_wps_value(value)
    if (!nzchar(value)) return(character())
    c(paste0("-", name), value)
  }, names(args), args), use.names = FALSE)
}

#' @noRd
noise_wps_value <- function(value) {
  if (length(value) == 0 || is.null(value) || all(is.na(value))) return("")
  if (is.logical(value)) return(if (isTRUE(value[1])) "true" else "false")
  if (length(value) > 1) return(paste(value, collapse = ","))
  as.character(value[1])
}

#' @noRd
format_favourable_occurrences <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.numeric(x)
  if (length(x) == 1) x <- rep(x, 16)
  if (length(x) != 16 || any(!is.finite(x))) {
    stop("`favourable_occurrences` must be one value or 16 finite values between 0 and 1.", call. = FALSE)
  }
  x <- pmin(pmax(x, 0), 1)
  paste(x, collapse = ",")
}

#' @noRd
check_noisemodelling_java <- function(wps, java = NULL) {
  java_bin <- resolve_java_bin(java)
  output <- suppressWarnings(system2(java_bin, "-version", stdout = TRUE, stderr = TRUE))
  version <- parse_java_major_version(output)
  if (is.na(version) || version < 11) {
    stop(
      "NoiseModelling 5.0.1 requires Java 11-21. ",
      "The detected Java runtime is: ", paste(output, collapse = " "), ". ",
      "Install a supported Java version, set `JAVA_HOME`, pass a Java executable path, ",
      "or use a version request such as `java = 17`.",
      call. = FALSE
    )
  }
  if (version > 21) {
    stop(
      "NoiseModelling 5.0.1 is not compatible with the detected Java runtime (Java ",
      version,
      "). ",
      "Its bundled Groovy runner can fail with newer class-file versions, for example ",
      "`Unsupported class file major version ", java_class_file_major(version), "`. ",
      "Use Java 11-21 instead; Java 17 is recommended. ",
      "Set `JAVA_HOME`, pass a Java executable path, or use `java = 17`.",
      call. = FALSE
    )
  }
  if (!file.exists(wps)) {
    stop("NoiseModelling runner was not found: ", wps, call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
resolve_java_bin <- function(java = NULL) {
  if (!is.null(java) && length(java) > 0 && !is.na(java[1]) && nzchar(as.character(java[1]))) {
    if (is_java_version_request(java[1])) {
      return(find_java_by_version(as.integer(java[1])))
    }
    java <- as.character(java[1])
    if (dir.exists(java)) {
      candidate <- file.path(java, "bin", if (.Platform$OS.type == "windows") "java.exe" else "java")
      if (file.exists(candidate)) return(candidate)
    }
    return(java)
  }
  java_home <- Sys.getenv("JAVA_HOME", unset = NA_character_)
  if (!is.na(java_home) && nzchar(java_home)) {
    candidate <- file.path(java_home, "bin", if (.Platform$OS.type == "windows") "java.exe" else "java")
    if (file.exists(candidate)) return(candidate)
  }
  "java"
}

#' @noRd
is_java_version_request <- function(java) {
  if (is.numeric(java) || is.integer(java)) return(is.finite(java))
  grepl("^[0-9]+$", as.character(java))
}

#' @noRd
find_java_by_version <- function(version) {
  if (is.na(version)) {
    stop("`java` must be a Java executable, Java home directory, or major version such as `java = 17`.", call. = FALSE)
  }
  candidates <- java_home_candidates(version)
  for (home in candidates) {
    java_bin <- file.path(home, "bin", if (.Platform$OS.type == "windows") "java.exe" else "java")
    if (file.exists(java_bin)) {
      output <- suppressWarnings(system2(java_bin, "-version", stdout = TRUE, stderr = TRUE))
      if (identical(parse_java_major_version(output), version)) {
        return(java_bin)
      }
    }
  }
  path_java <- Sys.which(if (.Platform$OS.type == "windows") "java.exe" else "java")
  if (nzchar(path_java)) {
    output <- suppressWarnings(system2(path_java, "-version", stdout = TRUE, stderr = TRUE))
    if (identical(parse_java_major_version(output), version)) {
      return(path_java)
    }
  }
  stop(
    "Java ", version, " was requested, but no matching Java runtime was found. ",
    "Install Java ", version, " or pass the Java executable path. ",
    "Supported Java versions for NoiseModelling 5.0.1 are 11-21; Java 17 is recommended.",
    call. = FALSE
  )
}

#' @noRd
java_home_candidates <- function(version) {
  candidates <- character()
  if (identical(Sys.info()[["sysname"]], "Darwin") && file.exists("/usr/libexec/java_home")) {
    output <- suppressWarnings(system2(
      "/usr/libexec/java_home",
      args = c("-v", as.character(version)),
      stdout = TRUE,
      stderr = TRUE
    ))
    status <- attr(output, "status")
    if (is.null(status) || identical(status, 0L)) {
      candidates <- c(candidates, output[nzchar(output)])
    }
  }
  java_home <- Sys.getenv("JAVA_HOME", unset = NA_character_)
  if (!is.na(java_home) && nzchar(java_home)) {
    candidates <- c(candidates, java_home)
  }
  candidates <- c(
    candidates,
    Sys.getenv(paste0("JAVA", version, "_HOME"), unset = NA_character_),
    Sys.getenv(paste0("JAVA_HOME_", version), unset = NA_character_),
    Sys.glob("/Library/Java/JavaVirtualMachines/*/Contents/Home"),
    Sys.glob("/usr/lib/jvm/*")
  )
  unique(normalizePath(candidates[!is.na(candidates) & nzchar(candidates)], mustWork = FALSE))
}

#' @noRd
parse_java_major_version <- function(output) {
  text <- paste(output, collapse = " ")
  match <- regmatches(text, regexpr('version "([^"]+)"', text))
  if (length(match) == 0 || !nzchar(match)) return(NA_integer_)
  version <- sub('version "([^"]+)".*', "\\1", match)
  parts <- strsplit(version, "\\.")[[1]]
  major <- suppressWarnings(as.integer(parts[1]))
  if (!is.na(major) && major == 1 && length(parts) > 1) {
    major <- suppressWarnings(as.integer(parts[2]))
  }
  major
}

#' @noRd
java_class_file_major <- function(java_major) {
  java_major + 44L
}

#' @noRd
patch_road_emission_wps_script <- function(script, work_dir) {
  lines <- readLines(script, warn = FALSE)
  patched <- lines
  patched <- sub(
    "createTableQuery\\.setLength\\(createTableQuery\\.length\\(\\) - 2\\)",
    "createTableQuery.delete(createTableQuery.length() - 2, createTableQuery.length())",
    patched
  )
  patched <- sub(
    "preparedInsertQuery\\.setLength\\(preparedInsertQuery\\.length\\(\\) - 2\\)",
    "preparedInsertQuery.delete(preparedInsertQuery.length() - 2, preparedInsertQuery.length())",
    patched
  )
  if (identical(lines, patched)) {
    return(script)
  }
  patched_dir <- file.path(work_dir, "patched_wps")
  dir.create(patched_dir, showWarnings = FALSE, recursive = TRUE)
  patched_script <- file.path(patched_dir, basename(script))
  writeLines(patched, patched_script, useBytes = TRUE)
  patched_script
}

#' @noRd
run_wps_script <- function(wps, work_dir, db_name, script, script_args, java = NULL, quiet = TRUE) {
  env <- character()
  java_bin <- resolve_java_bin(java)
  if (!identical(java_bin, "java")) {
    env <- paste0("JAVA_HOME=", dirname(dirname(java_bin)))
  }
  args <- c("-w", work_dir, "-d", db_name, "-s", script, script_args)
  if (!quiet) cli::cli_alert_info("Running NoiseModelling WPS script: {basename(script)}")
  output <- system2(wps, args = args, stdout = TRUE, stderr = TRUE, env = env)
  status <- attr(output, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop(
      "NoiseModelling WPS script failed: ", basename(script), "\n",
      paste(output, collapse = "\n"),
      call. = FALSE
    )
  }
  output
}
#' @noRd
fetch_greenspace_tile <- function(...) {
  if (!requireNamespace("greenSD", quietly = TRUE)) {
    stop("Package 'greenSD' is required for this function. Install it with: install.packages('greenSD')", call. = FALSE)
  }
  tryCatch(
    greenSD::get_tile_green(...),
    error = function(e) {
      message <- conditionMessage(e)
      if (grepl("lazy-load database|\\.rdb|R_decompress1|libdeflate|corrupt", message, ignore.case = TRUE)) {
        stop(
          "The installed `greenSD` package appears to be corrupt and cannot load its lazy-load database. ",
          "Restart R, reinstall `greenSD`, and then rerun the noise workflow. ",
          "For example: remove.packages('greenSD'); install.packages('greenSD'). ",
          "Original error: ", message,
          call. = FALSE
        )
      }
      stop(e)
    }
  )
}
