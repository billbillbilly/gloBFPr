#' search_3dglobdf
#' @description
#' Search and retrieve 3D building footprint data from 3D-GloBFP or
#' GlobalBuildingAtlas that intersect a given bounding box or area of interest
#' (a city), with options to return vector or raster outputs including building
#' polygons, binary presence rasters, and height-coded rasters.
#'
#' @param bbox `sf`, `sfc`, or a numeric vector (xmin, ymin, xmax, ymax)
#' defining the area of interest. This can be ignored if `place` is specified.
#' @param place vector (optional). A single line address,
#' e.g. ("1600 Pennsylvania Ave NW, Washington") or a vector of addresses
#' (c("Madrid", "Barcelona")).
#' @param crop logical. If `TRUE`, the resulting building footprint geometries
#' will be cropped to the input `bbox`. Default is `FALSE`.
#' @param data_source character. Building data source to query. Use `"GBF"` for
#' 3D-GloBFP (default) or `"GBA"`/`"gba"` for GlobalBuildingAtlas.
#' @param keep_source_id logical. If `TRUE`, keep the original source feature
#' identifier as `source_id` when it is available. Default is `FALSE`.
#' @param out_type character. Default is `'poly'`.
#' Output type(s) to return. Options include:
#'   \itemize{
#'     \item `"poly"`: building footprints as an `sf` polygon object.
#'     \item `"binary_rast"`: binary `terra` raster where buildings = 1.
#'     \item `"graduated_rast"`: `terra` raster encoding building height values.
#'     \item `"rast"`: a named list with both binary and graduated rasters.
#'     \item `"all"`: a named list including the polygon layer and both raster layers.
#'   }
#' @param mask logical (optional). Default is `FALSE`. If `TRUE`, masks the
#' graduated raster using the building footprint layer. Only used when `out_type`
#' is `"graduated_rast"`, `"rast"`, or `"all"`.
#' @param cell_size numeric (optional). Default is 1. Only used when `out_type`
#' is `"graduated_rast"`, `"rast"`, or `"all"`.
#' @param quiet logical. If `TRUE`, suppress cli messages and progress output.
#' Default is `TRUE`.
#'
#' @return Varies based on `out_type`:
#' \itemize{
#'   \item If `"poly"`: an `sf` object of building footprints. `MULTIPOLYGON`
#'   geometries are converted to `POLYGON` geometries while preserving one row
#'   per source feature. Polygons that touch or intersect share a `group_id`,
#'   which can be used to treat fragmented rows as one building group.
#'   \item If `"binary_rast"`: a binary `SpatRaster` (`terra`) indicating building presence.
#'   \item If `"graduated_rast"`: a quantitative `SpatRaster` of building heights.
#'   \item If `"rast"`: a named list with two `SpatRaster` objects: `binary` and `graduated`.
#'   \item If `"all"`: a named list with `poly` (sf), `binary`, and `graduated` rasters.
#' }
#'
#' @note
#' The downloading process may take some time, depending on the number and size
#' of building footprint tiles.
#'
#' This implementation for gloBFP-3D relies on the current structure of the dataset as hosted on Figshare.
#' It may break if the dataset owner changes the file organization or metadata format.
#'
#' The server of GlobalBuildingAtlas may have issues sometimes, so users may need to
#' switch over to gloBFP-3D, which means using `data_source="GBF"`
#'
#' When using `data_source = "GBA"`, the GlobalBuildingAtlas dataset does not
#' provide unique identifiers for individual building parcels. As a result, the
#' `group_id` assigned to overlapping or fragmented polygons may not accurately
#' reflect true building boundaries in all cases. This limitation may affect
#' morphological analyses at the individual building level (e.g., footprint area,
#' perimeter, building-level height statistics). However, the data remains suitable
#' for environmental simulation purposes such as noise mapping, solar radiation
#' analysis, and wind flow modelling, where parcel-level identity is less critical.
#'
#' @references
#' Che Yangzi, Li Xuecao, Liu Xiaoping, Wang Yuhao, Liao Weilin, Zheng Xianwei,
#' Zhang Xucai, Xu Xiaocong, Shi Qian, Zhu Jiajun, Zhang Honghui, Yuan Hua, &
#' Dai Yongjiu (2024). 3D-GloBFP: the first global three-dimensional building
#' footprint dataset. Earth Syst. Sci. Data, 16, 5357-5374
#'
#' Zhu X. X., Chen S., Zhang F., Shi Y., & Wang Y. (2025).
#' GlobalBuildingAtlas: an open global and complete dataset of building
#' polygons, heights and LoD1 3D models. Earth Syst. Sci. Data, 17, 6647-6668.
#'
#' @examples
#' \dontrun{
#' buildings <- gloBFPr::search_3dglobdf(bbox=c(-84.485519,45.636118,-84.462774,45.650639))
#' }
#'
#' @importFrom sf st_bbox
#' @importFrom sf st_read
#' @importFrom sf st_crop
#' @importFrom sf st_intersects
#' @importFrom utils download.file
#' @importFrom utils unzip
#' @importFrom utils URLencode
#' @importFrom dplyr bind_rows
#' @importFrom cli cli_alert_info
#' @importFrom cli cli_alert_success
#' @importFrom nominatimlite geo_lite_sf
#'
#' @export

search_3dglobdf <- function(bbox=NULL,
                            place=NULL,
                            crop=FALSE,
                            data_source="GBF",
                            keep_source_id=FALSE,
                            out_type='poly',
                            mask=FALSE,
                            cell_size=1,
                            quiet=TRUE) {
  if (inherits(bbox, 'NULL') && inherits(place, 'NULL')) {
    base::warning('Area of interest is missing: boox or place')
    return(NULL)
  }

  start_time <- Sys.time()

  if ((inherits(bbox, 'NULL') && !inherits(place, 'NULL')) || (!inherits(bbox, 'NULL') && !inherits(place, 'NULL'))) {
    city <- nominatimlite::geo_lite_sf(place, points_only = FALSE)
    city <- sf::st_transform(city, crs = 4326)
    bbox <- city$geometry
  } else if (!inherits(bbox, 'NULL') && inherits(place, 'NULL') ) {
    # check type of bbox
    if (is.numeric(bbox) && length(bbox) == 4) {
      # convert bbox to sf when it is not a sf polygon
      bbox <- sf::st_as_sfc(
        sf::st_bbox(
          c(xmin = unname(bbox[1]),
            ymin = unname(bbox[2]),
            xmax = unname(bbox[3]),
            ymax = unname(bbox[4])),
          crs = 4326
        )
      )
    }
  }

  # Ensure input is in WGS84
  bbox <- sf::st_transform(bbox, 4326)

  data_source <- normalize_building_data_source(data_source)
  if (data_source == "GBF") {
    all_data <- read_gbf_buildings(bbox, quiet = quiet)
  } else {
    all_data <- read_gba_buildings(bbox, quiet = quiet)
  }

  if (is.null(all_data) || nrow(all_data) == 0) {
    base::warning("No buildings were found for the provided bbox.")
    return(NULL)
  }

  utm_crs <- get_utm_crs(bbox)
  bbox_proj <- sf::st_transform(bbox, crs = utm_crs)
  all_data <- sf::st_transform(all_data, crs = utm_crs)
  # sf_data <- sf_data[!sf::st_is_empty(sf_data), ]
  all_data <- suppressWarnings(all_data[sf::st_intersects(all_data, bbox_proj, sparse = FALSE), ])
  # crop the data if 'crop' == true
  if (crop) {
    #all_data <- suppressWarnings(all_data[sf::st_intersects(all_data, bbox, sparse = FALSE), ])
    all_data <- suppressWarnings(sf::st_crop(all_data, bbox_proj))
  }
  all_data <- normalize_building_geometries(all_data)
  if (is.null(all_data) || nrow(all_data) == 0) {
    base::warning("No buildings were found for the provided bbox.")
    return(NULL)
  }
  all_data <- drop_internal_building_columns(all_data, keep_source_id = keep_source_id)
  all_data <- assign_building_group_id(all_data)

  # assign an id to each building
  all_data$id <- seq_len(nrow(all_data))

  # export data as sf polygons
  if(out_type == 'poly') {
    end_time <- Sys.time()
    process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    if (!quiet) time_taken(process_time)
    return(all_data)
  }

  # auto-generate raster outputs
  if (out_type %in% c("binary_rast", "graduated_rast", "rast", "all")) {
    if (isTRUE(mask) || out_type == 'binary_rast' || out_type == 'all') {
      binary <- rasterize_binary(all_data, bbox, res=cell_size)
    }

    if (isTRUE(mask) && out_type != "binary_rast") {
      graduated <- rasterize_height(all_data, bbox, res=cell_size, mask=binary)
    } else {
      graduated <- rasterize_height(all_data, bbox, res=cell_size)
    }

    if (out_type == "binary_rast") {
      end_time <- Sys.time()
      process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      if (!quiet) time_taken(process_time)
      return(binary)
    }
    if (out_type == "graduated_rast") {
      end_time <- Sys.time()
      process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      if (!quiet) time_taken(process_time)
      return(graduated)
    }
    if (out_type == "rast") {
      end_time <- Sys.time()
      process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      if (!quiet) time_taken(process_time)
      return(list(binary = binary, graduated = graduated))
    }
    if (out_type == "all") {
      end_time <- Sys.time()
      process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      if (!quiet) time_taken(process_time)
      return(list(poly = all_data, binary = binary, graduated = graduated))
    }
  }

  stop("Invalid out_type specified.")
}

#' get_fused_dsm
#' @description
#' Generate digital surface model using multiple datasets, including building height map,
#' canopy height map, and terrain model. Each building is given a single flat
#' roof elevation (its own base ground elevation, sampled at its centroid,
#' plus its height) rather than following the terrain slope beneath it pixel
#' by pixel, so buildings on sloped ground do not come out tilted or warped.
#' @param x sf. building footprint polygon, typically output from [search_3dglobdf()]
#' @param datasource_canopy_height character or `NULL`. Canopy height source.
#' Currently supports `"metachm"`, `"ethCHM"`, or `NULL`.
#' @param min_tree_height numeric. Minimum canopy height threshold in meters.
#' @param resolution numeric or `NULL`. Output raster resolution in meters. If
#' `NULL` (default), the finest native resolution among the downloaded DEM
#' and canopy height model is used, so the output is never silently degraded
#' to the coarser of the two source rasters. Set explicitly (e.g. `1`) to
#' force a finer grid than the native source data (e.g. when the DEM falls
#' back to coarse SRTM data), or a coarser one to speed up large areas.
#' @param opentopo_key character. OpenTopography API key used to download DEM data.
#' @param key Deprecated alias for `opentopo_key`.
#' @param quiet logical. If `TRUE`, suppress cli messages and progress output.
#' Default is `TRUE`.
#' @examples
#' \donttest{
#'  example <- gloBFPr::globfp_example
#'  dsm <- get_fused_dsm(x= example, opentopo_key = 'key')
#' }
#'
#' @export
get_fused_dsm <- function(x = NULL,
                          datasource_canopy_height = "metachm",
                          min_tree_height = 2,
                          resolution = NULL,
                          opentopo_key = NULL,
                          key = NULL,
                          quiet = TRUE) {
  if (inherits(x, 'NULL')) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygon generated by `search_3dglobdf()`.")
    return(x)
  }
  if (inherits(opentopo_key, 'NULL') && !inherits(key, 'NULL')) {
    opentopo_key <- key
  }
  if (inherits(opentopo_key, 'NULL')) {
    stop("API key for OpenTopography is missing.")
  }

  normalize_source <- function(value, choices) {
    if (inherits(value, "NULL")) return(NULL)
    value <- tolower(as.character(value[1]))
    if (value %in% c("none", "null", "na")) return(NULL)
    match.arg(value, choices)
  }
  datasource_canopy_height <- normalize_source(datasource_canopy_height, c("metachm", "ethchm"))

  old_terra_progress <- terra::terraOptions(print = FALSE)$progress
  terra::terraOptions(progress = 0)
  on.exit(terra::terraOptions(progress = old_terra_progress), add = TRUE)

  projected_poly <- x
  bbox <- get_bbox(x)
  bbox_vector <- bbox_poly_to_list(bbox)
  utm_crs <- get_utm_crs(bbox)

  if (!quiet) cli::cli_alert_info('Start downloading DSM raster inputs ...')
  dem <- get_dem(bbox_vector, opentopo_key)
  # Bilinear, not nearest-neighbor: DEM is a continuous surface, and "near"
  # produces a visibly blocky/stair-stepped terrain once resampled to a
  # finer grid below.
  dem_projected <- terra::project(dem, paste0("EPSG:", utm_crs), method = "bilinear")
  dem_res <- terra::res(dem_projected)[1]

  chm_layers <- NULL
  if (!inherits(datasource_canopy_height, "NULL")) {
    chm_args <- list(
      bbox = bbox_vector,
      min_height = min_tree_height,
      datasource = datasource_canopy_height
    )
    if ("quiet" %in% names(formals(get_chm))) {
      chm_args$quiet <- quiet
    }
    chm_layers <- suppressMessages(do.call(get_chm, chm_args))
  }

  # Resolve the output resolution: an explicit override if supplied,
  # otherwise the *finer* of the DEM's and CHM's native resolutions. This
  # fixes two prior issues: (1) when no canopy layer was requested, the
  # fused DSM inherited the DEM's native resolution verbatim, which can be
  # as coarse as ~30 m when USGS 1 m/10 m coverage isn't available and the
  # SRTM fallback is used; (2) the DEM was always used as the alignment
  # reference regardless of resolution, so even a finer CHM/building raster
  # got silently downsampled to match a coarser DEM.
  if (!inherits(resolution, "NULL")) {
    raster_res <- as.numeric(resolution[1])
    if (is.na(raster_res) || raster_res <= 0) {
      stop("`resolution` must be a positive number (meters).", call. = FALSE)
    }
  } else if (!inherits(chm_layers, "NULL")) {
    raster_res <- min(dem_res, terra::res(chm_layers[[1]])[1])
  } else {
    raster_res <- dem_res
  }

  # Build the building-height raster at the target resolution and use it as
  # the alignment grid, resampling the continuous DEM/CHM surfaces onto it
  # (bilinear) rather than the other way around.
  bh_all <- rasterize_height(projected_poly, bbox, raster_res)
  dem <- terra::resample(dem_projected, bh_all, method = "bilinear")
  chm <- if (!inherits(chm_layers, "NULL")) {
    terra::resample(chm_layers[[1]], bh_all, method = "bilinear")
  } else {
    NULL
  }

  bh_all <- terra::ifel(is.na(bh_all), 0, bh_all)
  chm <- if (!is.null(chm)) {
    terra::ifel(is.na(chm), 0, chm)
  } else {
    terra::ifel(is.na(dem), NA, 0)
  }
  surface <- terra::ifel(chm > bh_all, chm, bh_all)
  ground_dsm <- dem + surface

  # -- Flatten building roofs ------------------------------------------------
  # `dem` varies pixel-by-pixel with underlying terrain, so `dem + bh_all`
  # alone would tilt/warp every roof to follow the ground slope beneath it.
  # Buildings are rigid structures with a single roof elevation, so instead
  # each building gets one base ground elevation (DEM sampled at its
  # centroid, matching the convention used by get_3d_world()) added to its
  # height, producing a flat roof across its whole footprint.
  analysis_poly  <- prepare_group_analysis_buildings(projected_poly)
  if (!"id" %in% names(analysis_poly)) {
    analysis_poly$id <- seq_len(nrow(analysis_poly))
  }
  bldg_id_rast   <- rasterize_height(analysis_poly, bbox, raster_res, height_field = "id")
  building_mask  <- bldg_id_rast > 0

  centroids  <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(analysis_poly))))
  base_elev  <- terra::extract(dem, centroids[, 1:2, drop = FALSE], method = "bilinear")[, 1]
  dem_median <- stats::median(terra::values(dem, mat = FALSE), na.rm = TRUE)
  base_elev[!is.finite(base_elev)] <- dem_median

  base_elev_rast <- terra::subst(bldg_id_rast, from = analysis_poly$id, to = base_elev)
  flat_roof      <- base_elev_rast + bh_all

  # An overhanging tree taller than the roof at that pixel still wins.
  dsm <- terra::ifel(
    building_mask,
    terra::ifel(flat_roof > ground_dsm, flat_roof, ground_dsm),
    ground_dsm
  )

  return(dsm)
}
