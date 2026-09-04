#' Prepare gloBFPr spatial outputs for an OpenFOAM Docker workflow
#'
#' @description
#' Collects building, terrain, canopy, and ground-cover data from gloBFPr and
#' writes them into a structured OpenFOAM case folder ready for external CFD
#' simulation (via Docker).
#'
#' Output layers and their intended role in OpenFOAM:
#' \describe{
#'   \item{building STL}{Solid geometry for snappyHexMesh}
#'   \item{building height raster}{Auxiliary mesh-generation reference}
#'   \item{binary building raster}{Building mask / validation}
#'   \item{fused DSM}{Optional ground STL source (terrain + buildings +
#'     canopy surface)}
#'   \item{canopy height raster}{Standalone CHM for defining porous-zone
#'     extents and drag coefficients in fvOptions / topoSetDict}
#'   \item{ground roughness raster (z0)}{Per-cell aerodynamic roughness length
#'     derived from ESA WorldCover land cover, for nutURoughWallFunction.
#'     Building footprint cells and tree-cover cells are set to NA because they
#'     are handled by solid geometry and porous zones respectively.}
#' }
#'
#' @param case_dir Character. OpenFOAM case directory.
#' @param bbox Numeric vector c(xmin, ymin, xmax, ymax) in WGS-84 lon/lat.
#' @param place Optional place name passed to search_3dglobdf().
#' @param buildings_list Optional precomputed output from
#'   search_3dglobdf(..., out_type = "all").
#' @param data_source Character. Building source passed to search_3dglobdf().
#'   Default "GBF".
#' @param cell_size Numeric. Raster resolution in metres for building rasters.
#' @param crop Logical. Whether to crop buildings to bbox.
#' @param mask Logical. Whether to mask height raster by building footprints.
#' @param include_buildings Logical. Prepare building vector/raster/STL.
#' @param include_fused_dsm Logical. Prepare fused DSM.
#' @param include_tree_canopy Logical. Extract standalone canopy height raster
#'   (CHM) and include canopy in fused DSM. The CHM is written as a separate
#'   file so OpenFOAM porous-zone definitions can reference it directly.
#' @param canopy_source Character or NULL. "metachm", "ethCHM", or NULL.
#' @param min_tree_height Numeric. Minimum tree canopy height in metres.
#' @param include_morphology Logical. Add morphology metrics.
#' @param include_neighbors Logical. Add neighbour metrics.
#' @param include_greenspace Logical. Produce two greenspace outputs:
#'   (1) a ground roughness raster (z0, metres) from land-cover data
#'   for nutURoughWallFunction, and (2) distance-to-greenspace as a
#'   building-level attribute for post-processing / context.
#' @param mask_tree_cover Logical. When building the roughness raster, set
#'   tree-cover cells to NA so porous-zone drag is not counted twice.
#'   Default TRUE.
#' @param landcover_source Character. Land-cover dataset for the ground
#'   roughness raster. `"esa"` (default) uses ESA WorldCover (years 2020-2021);
#'   `"esri"` uses the Sentinel-2 10 m ESRI LULC Time Series (years 2017-2025),
#'   which provides more recent and historically consistent annual maps.
#' @param landcover_year Integer. Year of the land-cover product.
#'   For `"esa"`: 2020 or 2021. For `"esri"`: 2017-2025. Default 2021.
#' @param include_bgvi Logical. Add building green visibility index.
#' @param include_shadow Logical. Add shadow outputs.
#' @param include_radiation Logical. Add radiation outputs.
#' @param opentopo_key Character. OpenTopography API key for get_fused_dsm().
#' @param target_crs Optional projected target CRS. Usually leave as `NULL`:
#'   when the building data are in longitude/latitude, the local UTM zone is
#'   derived automatically from the data centroid, whichever of `bbox`,
#'   `place`, or `buildings_list` was supplied. Set this only to force a
#'   specific projection (e.g. a national grid). OpenFOAM domains must be
#'   metric, so a geographic CRS is rejected.
#' @param height_col Optional height column name. If NULL, guessed
#'   automatically.
#' @param default_height Numeric. Default building height if height column is
#'   missing (metres).
#' @param min_height Numeric. Minimum building height (metres).
#' @param domain_buffer Numeric. Buffer in metres around data extent.
#' @param zmax_buffer Numeric. Buffer in metres above maximum surface height.
#' @param stl_name Character. STL file name.
#' @param overwrite Logical. Whether to overwrite prepared files.
#' @param quiet Logical. Suppress messages.
#'
#' @return A list with:
#' \describe{
#'   \item{case_dir}{Absolute path to the OpenFOAM case directory.}
#'   \item{data_dir}{Absolute path to the gloBFPr data sub-directory.}
#'   \item{files}{Named list of output file paths (NULL when not produced).}
#'   \item{data}{Named list of in-memory spatial objects.}
#'   \item{domain}{Named list with xmin/xmax/ymin/ymax/zmin/zmax for
#'     blockMeshDict.}
#'   \item{origin}{Named numeric vector (x, y, z) of the local coordinate
#'     system origin in the input CRS.}
#'   \item{crs}{CRS of the building data.}
#'   \item{height_col}{Name of the height column used.}
#'   \item{n_buildings}{Number of building polygons.}
#'   \item{max_building_height}{Maximum building height in metres.}
#'   \item{max_surface_height}{Maximum surface height (buildings + DSM).}
#'   \item{include_tree_canopy}{Logical flag as supplied.}
#'   \item{canopy_source}{Canopy data source as supplied.}
#'   \item{min_tree_height}{Minimum tree height as supplied.}
#'   \item{landcover_source}{Land-cover dataset used (`"esa"` or `"esri"`).}
#'   \item{landcover_year}{Land-cover year used.}
#' }
#'
#' @export
prepare_openfoam_inputs <- function(
    case_dir,
    bbox                = NULL,
    place               = NULL,
    buildings_list      = NULL,
    data_source         = "GBF",
    cell_size           = 1,
    crop                = FALSE,
    mask                = TRUE,
    include_buildings   = TRUE,
    include_fused_dsm   = TRUE,
    include_tree_canopy = TRUE,
    canopy_source       = "metachm",
    min_tree_height     = 2,
    include_morphology  = TRUE,
    include_neighbors   = TRUE,
    include_greenspace  = TRUE,
    mask_tree_cover     = TRUE,
    landcover_source    = c("esa", "esri"),
    landcover_year      = 2021,
    include_bgvi        = FALSE,
    include_shadow      = FALSE,
    include_radiation   = FALSE,
    opentopo_key        = NULL,
    target_crs          = NULL,
    height_col          = NULL,
    default_height      = 10,
    min_height          = 2,
    domain_buffer       = 100,
    zmax_buffer         = 50,
    stl_name            = "buildings.stl",
    overwrite           = FALSE,
    quiet               = FALSE
) {

  # Internal helpers
  require_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package `", pkg, "` is required but not installed.", call. = FALSE)
    }
  }
  require_pkg("sf")
  require_pkg("terra")

  msg <- function(...) if (!isTRUE(quiet)) message(...)

  make_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE,
                                      showWarnings = FALSE)
  }

  safe_call <- function(label, expr) {
    msg("Preparing ", label, " ...")
    tryCatch(expr, error = function(e) {
      warning("Failed to prepare ", label, ": ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
  }

  guess_height_column <- function(x) {
    candidates <- c(
      "height", "Height", "HEIGHT",
      "bldg_h", "bldg_ht", "building_h", "building_height",
      "mean_height", "median_height", "height_m", "elevation"
    )
    matched <- candidates[candidates %in% names(x)]
    if (length(matched) > 0) return(matched[1])
    numeric_cols <- names(x)[vapply(x, is.numeric, logical(1))]
    height_like  <- numeric_cols[
      grepl("height|hgt|bldg|elev", numeric_cols, ignore.case = TRUE)
    ]
    if (length(height_like) > 0) return(height_like[1])
    NULL
  }

  clean_numeric_height <- function(x, default_height, min_height) {
    x <- suppressWarnings(as.numeric(x))
    x[!is.finite(x)] <- default_height
    x[x < min_height] <- min_height
    x
  }

  triangle_normal <- function(p1, p2, p3) {
    u <- p2 - p1; v <- p3 - p1
    n <- c(
      u[2] * v[3] - u[3] * v[2],
      u[3] * v[1] - u[1] * v[3],
      u[1] * v[2] - u[2] * v[1]
    )
    len <- sqrt(sum(n ^ 2))
    if (!is.finite(len) || len == 0) return(c(0, 0, 0))
    n / len
  }

  facet_to_stl <- function(p1, p2, p3) {
    n <- triangle_normal(p1, p2, p3)
    paste0(
      "  facet normal ", paste(format(n, scientific = FALSE), collapse = " "), "\n",
      "    outer loop\n",
      "      vertex ", paste(format(p1, scientific = FALSE), collapse = " "), "\n",
      "      vertex ", paste(format(p2, scientific = FALSE), collapse = " "), "\n",
      "      vertex ", paste(format(p3, scientific = FALSE), collapse = " "), "\n",
      "    endloop\n",
      "  endfacet\n"
    )
  }

  polygon_to_stl_facets <- function(coords, height) {
    if (nrow(coords) < 3) return(character(0))
    z0 <- 0; z1 <- height
    facets <- character(0)
    n <- nrow(coords)
    for (i in seq_len(n)) {
      j  <- if (i == n) 1L else i + 1L
      p1 <- c(coords[i, 1], coords[i, 2], z0)
      p2 <- c(coords[j, 1], coords[j, 2], z0)
      p3 <- c(coords[j, 1], coords[j, 2], z1)
      p4 <- c(coords[i, 1], coords[i, 2], z1)
      facets <- c(facets, facet_to_stl(p1, p2, p3), facet_to_stl(p1, p3, p4))
    }
    # Roof cap (fan triangulation from vertex 0)
    p0_top <- c(coords[1, 1], coords[1, 2], z1)
    for (i in 2:(n - 1)) {
      facets <- c(facets, facet_to_stl(
        p0_top,
        c(coords[i,     1], coords[i,     2], z1),
        c(coords[i + 1, 1], coords[i + 1, 2], z1)
      ))
    }
    # Bottom cap
    p0_bot <- c(coords[1, 1], coords[1, 2], z0)
    for (i in 2:(n - 1)) {
      facets <- c(facets, facet_to_stl(
        p0_bot,
        c(coords[i + 1, 1], coords[i + 1, 2], z0),
        c(coords[i,     1], coords[i,     2], z0)
      ))
    }
    facets
  }

  write_buildings_stl <- function(buildings, height_col, stl_path) {
    buildings <- sf::st_make_valid(buildings)
    buildings <- suppressWarnings(sf::st_collection_extract(buildings, "POLYGON"))
    buildings <- suppressWarnings(sf::st_cast(buildings, "POLYGON"))
    if (nrow(buildings) == 0)
      stop("No valid building polygons available for STL export.", call. = FALSE)
    con <- file(stl_path, open = "w", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines("solid buildings", con)
    for (i in seq_len(nrow(buildings))) {
      geom    <- sf::st_geometry(buildings[i, ])[[1]]
      exterior <- geom[[1]]
      coords  <- as.matrix(exterior[, 1:2, drop = FALSE])
      if (nrow(coords) > 1 && all(coords[1, ] == coords[nrow(coords), ]))
        coords <- coords[-nrow(coords), , drop = FALSE]
      h      <- buildings[[height_col]][i]
      facets <- polygon_to_stl_facets(coords, h)
      if (length(facets) > 0) writeLines(facets, con)
    }
    writeLines("endsolid buildings", con)
    invisible(stl_path)
  }

  shift_sf_to_local <- function(x, origin) {
    original_crs <- sf::st_crs(x)
    geom <- sf::st_geometry(x) - c(origin["x"], origin["y"])
    sf::st_geometry(x) <- geom
    sf::st_crs(x) <- original_crs   # the `-` operator drops CRS; restore it
    x
  }

  shift_raster_to_local <- function(r, origin) {
    if (is.null(r)) return(NULL)
    e <- terra::ext(r)
    terra::ext(r) <- terra::ext(
      e[1] - origin["x"], e[2] - origin["x"],
      e[3] - origin["y"], e[4] - origin["y"]
    )
    r
  }

  write_raster_if_not_null <- function(r, path) {
    if (is.null(r)) return(NULL)
    terra::writeRaster(r, path, overwrite = TRUE)
    normalizePath(path, mustWork = FALSE)
  }

  # Directories
  if (missing(case_dir) || !nzchar(case_dir))
    stop("`case_dir` must be provided.", call. = FALSE)
  if (!is.numeric(cell_size) || length(cell_size) != 1L ||
      !is.finite(cell_size) || cell_size <= 0)
    stop("`cell_size` must be one positive number.", call. = FALSE)
  if (!is.numeric(domain_buffer) || length(domain_buffer) != 1L ||
      !is.finite(domain_buffer) || domain_buffer < 0)
    stop("`domain_buffer` must be one non-negative number.", call. = FALSE)
  if (!is.numeric(zmax_buffer) || length(zmax_buffer) != 1L ||
      !is.finite(zmax_buffer) || zmax_buffer < 0)
    stop("`zmax_buffer` must be one non-negative number.", call. = FALSE)
  landcover_source <- match.arg(landcover_source)
  if (landcover_source == "esa" && !landcover_year %in% c(2020L, 2021L))
    stop("`landcover_year` must be 2020 or 2021 when `landcover_source = 'esa'`. ",
         "Use `landcover_source = 'esri'` for years 2017-2025.", call. = FALSE)
  if (landcover_source == "esri" && !landcover_year %in% 2017:2025)
    stop("`landcover_year` must be between 2017 and 2025 when `landcover_source = 'esri'`.",
         call. = FALSE)

  constant_dir  <- file.path(case_dir,    "constant")
  tri_dir       <- file.path(constant_dir, "triSurface")
  data_dir      <- file.path(constant_dir, "gloBFPr")
  raster_dir    <- file.path(data_dir,    "rasters")
  vector_dir    <- file.path(data_dir,    "vectors")
  metrics_dir   <- file.path(data_dir,    "metrics")
  metadata_dir  <- file.path(data_dir,    "metadata")

  if (dir.exists(data_dir) && !overwrite)
    stop(
      "OpenFOAM input data already exists at: ", data_dir, "\n",
      "Use `overwrite = TRUE` to replace it.",
      call. = FALSE
    )
  if (overwrite && dir.exists(data_dir))
    unlink(data_dir, recursive = TRUE, force = TRUE)

  for (d in c(case_dir, constant_dir, tri_dir, data_dir,
              raster_dir, vector_dir, metrics_dir, metadata_dir))
    make_dir(d)

  # Building data
  if (is.null(buildings_list)) {
    if (is.null(bbox) && is.null(place))
      stop("Provide either `bbox`, `place`, or `buildings_list`.", call. = FALSE)
    msg("Searching building data with search_3dglobdf() ...")
    buildings_list <- search_3dglobdf(
      bbox        = bbox,
      place       = place,
      crop        = crop,
      data_source = data_source,
      out_type    = "all",
      mask        = mask,
      cell_size   = cell_size,
      quiet       = quiet
    )
  }
  if (is.null(buildings_list))
    stop("No building data returned by search_3dglobdf().", call. = FALSE)

  buildings     <- buildings_list$poly
  binary        <- buildings_list$binary
  height_raster <- buildings_list$graduated

  if (is.null(buildings) || !inherits(buildings, "sf") || nrow(buildings) == 0)
    stop("No valid building polygon layer found.", call. = FALSE)
  if (is.null(sf::st_crs(buildings)))
    stop("Building polygons have no CRS.", call. = FALSE)

  # OpenFOAM coordinates and all buffer arguments are metric. Automatically use
  # the local UTM zone when the caller supplies longitude/latitude data. This
  # works regardless of whether the buildings came from `bbox`, `place`, or
  # `buildings_list`, because all three converge on `buildings` by this point.
  if (is.null(target_crs) && isTRUE(sf::st_is_longlat(buildings))) {
    target_crs <- get_utm_crs(buildings)
    msg("Auto-selected projected CRS: EPSG:", target_crs,
        " (UTM zone from data centroid)")
  }
  if (!is.null(target_crs))
    buildings <- sf::st_transform(buildings, target_crs)

  if (isTRUE(sf::st_is_longlat(buildings)))
    stop("Building data are still in longitude/latitude after CRS handling. ",
         "OpenFOAM domains must be metric - supply a projected `target_crs`.",
         call. = FALSE)

  # Supplied search results can contain rasters in a CRS different from the
  # polygon layer. Keep every layer in the same metric coordinate system.
  building_crs <- sf::st_crs(buildings)$wkt
  project_to_building_crs <- function(r, method) {
    if (is.null(r)) return(NULL)
    if (!terra::same.crs(r, building_crs))
      r <- terra::project(r, building_crs, method = method)
    r
  }
  binary        <- project_to_building_crs(binary, "near")
  height_raster <- project_to_building_crs(height_raster, "bilinear")

  if (is.null(height_col)) height_col <- guess_height_column(buildings)
  if (is.null(height_col)) {
    height_col           <- ".openfoam_height_m"
    buildings[[height_col]] <- default_height
  } else {
    buildings[[height_col]] <- clean_numeric_height(
      buildings[[height_col]],
      default_height = default_height,
      min_height     = min_height
    )
  }

  # Derive a WGS-84 bbox for functions that need it (canopy, roughness).
  # Prefer the explicit bbox argument; fall back to the building extent.
  bbox_wgs84 <- if (!is.null(bbox)) {
    as_wgs84_bbox_vector(bbox)
  } else {
    as_wgs84_bbox_vector(buildings)
  }


  # Fused DSM


  fused_dsm <- NULL
  if (isTRUE(include_fused_dsm)) {
    if (is.null(opentopo_key) || !nzchar(opentopo_key)) {
      warning(
        "`include_fused_dsm = TRUE`, but `opentopo_key` is missing. ",
        "Skipping fused DSM.",
        call. = FALSE
      )
    } else {
      canopy_arg <- if (isTRUE(include_tree_canopy)) canopy_source else NULL
      fused_dsm  <- safe_call(
        "fused DSM",
        get_fused_dsm(
          x                        = buildings,
          datasource_canopy_height = canopy_arg,
          min_tree_height          = min_tree_height,
          opentopo_key             = opentopo_key,
          quiet                    = quiet
        )
      )
    }
  }


  # Standalone canopy height raster

  # get_chm() returns list(filteredCHM, binaryCHM).  Index [[1]] is the
  # continuous height raster needed for porous-zone / drag-coefficient
  # definitions in OpenFOAM.  get_greenspace() and get_fused_dsm() only ever
  # use the binary [[2]] layer, so the height raster was silently discarded
  # before this fix.

  canopy_height_raster <- NULL
  if (isTRUE(include_tree_canopy)) {
    chm_result <- safe_call(
      "canopy height raster",
      get_chm(bbox_wgs84, min_height = min_tree_height,
              datasource = canopy_source, quiet = quiet)
    )
    if (!is.null(chm_result)) {
      canopy_height_raster <- chm_result[[1]]   # continuous CHM, not binary
      if (!is.null(target_crs))
        canopy_height_raster <- terra::project(canopy_height_raster, target_crs)
    }
  }


  # Optional vector / metric outputs (building-level)


  metrics          <- list()
  buildings_metrics <- buildings

  if (isTRUE(include_morphology)) {
    tmp <- safe_call("morphology metrics", get_morphology(buildings_metrics))
    if (!is.null(tmp)) { buildings_metrics <- tmp; metrics$morphology <- tmp }
  }

  if (isTRUE(include_neighbors)) {
    tmp <- safe_call("neighbor metrics", get_neighbors(buildings_metrics))
    if (!is.null(tmp)) { buildings_metrics <- tmp; metrics$neighbors <- tmp }
  }

  # Ground roughness raster (OpenFOAM nutURoughWallFunction input) AND
  # building-level distance-to-greenspace attribute (post-processing context).
  roughness_raster <- NULL
  if (isTRUE(include_greenspace)) {

    # 1. Ground roughness raster from ESA WorldCover land cover
    roughness_raster <- safe_call(
      "ground roughness raster (z0)",
      get_roughness_raster(
        bbox            = bbox_wgs84,
        source          = landcover_source,
        year            = landcover_year,
        crs             = sf::st_crs(buildings),
        building_mask   = buildings,
        mask_tree_cover = mask_tree_cover
      )
    )

    # 2. Building-level distance to nearest greenspace (context / post-processing)
    dng_source <- if (identical(canopy_source, "metachm")) "metachm" else "esri"
    tmp <- safe_call("distance to nearest greenspace",
                     get_dng(buildings_metrics, datasource = dng_source,
                             min_tree_height = min_tree_height, quiet = quiet))
    if (!is.null(tmp)) {
      buildings_metrics          <- tmp
      metrics$greenspace_distance <- tmp
    }
  }

  if (isTRUE(include_bgvi)) {
    tmp <- safe_call("building green view index", get_bgvi(buildings_metrics))
    if (!is.null(tmp)) { buildings_metrics <- tmp; metrics$bgvi <- tmp }
  }

  if (isTRUE(include_shadow)) {
    metrics$shadow_footprint <- safe_call("shadow footprint",
                                          get_shadow_footprint(buildings_metrics))
    metrics$shadow_height    <- safe_call("shadow height",
                                          get_shadow_height(buildings_metrics))
  }

  if (isTRUE(include_radiation)) {
    metrics$radiation <- safe_call("radiation", get_radiation(buildings_metrics))
  }


  # Local OpenFOAM coordinate system
  bb     <- sf::st_bbox(buildings_metrics)
  origin <- c(
    x = unname(bb["xmin"] - domain_buffer),
    y = unname(bb["ymin"] - domain_buffer),
    z = 0
  )

  buildings_local          <- shift_sf_to_local(buildings_metrics, origin)
  binary_local             <- shift_raster_to_local(binary,               origin)
  height_raster_local      <- shift_raster_to_local(height_raster,        origin)
  fused_dsm_local          <- shift_raster_to_local(fused_dsm,            origin)
  canopy_height_raster_local <- shift_raster_to_local(canopy_height_raster, origin)
  roughness_raster_local   <- shift_raster_to_local(roughness_raster,     origin)

  # Domain extents for blockMeshDict
  bb_local        <- sf::st_bbox(buildings_local)
  max_building_h  <- max(buildings_local[[height_col]], na.rm = TRUE)
  max_surface_h   <- max_building_h
  if (!is.null(fused_dsm_local)) {
    dsm_max <- suppressWarnings(
      max(terra::global(fused_dsm_local, "max", na.rm = TRUE)[[1]], na.rm = TRUE)
    )
    if (is.finite(dsm_max))
      max_surface_h <- max(max_surface_h, dsm_max, na.rm = TRUE)
  }

  domain <- list(
    xmin = 0,
    xmax = unname(bb_local["xmax"] + domain_buffer),
    ymin = 0,
    ymax = unname(bb_local["ymax"] + domain_buffer),
    zmin = 0,
    zmax = max_surface_h + zmax_buffer
  )


  # Write files


  files <- list(
    building_gpkg          = file.path(vector_dir,   "buildings_openfoam.gpkg"),
    building_stl           = file.path(tri_dir,       stl_name),
    building_binary_raster = file.path(raster_dir,   "building_binary.tif"),
    building_height_raster = file.path(raster_dir,   "building_height.tif"),
    fused_dsm              = file.path(raster_dir,   "fused_dsm.tif"),
    canopy_height_raster   = file.path(raster_dir,   "canopy_height.tif"),
    roughness_raster       = file.path(raster_dir,   "ground_roughness_z0.tif"),
    metadata_rds           = file.path(metadata_dir, "openfoam_inputs_metadata.rds"),
    buildings_rds          = file.path(metadata_dir, "buildings_openfoam.rds"),
    metrics_rds            = file.path(metadata_dir, "metrics_list.rds")
  )

  msg("Writing building GeoPackage ...")
  sf::st_write(buildings_local, files$building_gpkg,
               delete_dsn = TRUE, quiet = TRUE)

  if (isTRUE(include_buildings)) {
    msg("Writing building STL ...")
    write_buildings_stl(
      buildings  = buildings_local,
      height_col = height_col,
      stl_path   = files$building_stl
    )
  } else {
    files$building_stl <- NULL
  }

  files$building_binary_raster <- write_raster_if_not_null(
    binary_local, files$building_binary_raster)
  files$building_height_raster <- write_raster_if_not_null(
    height_raster_local, files$building_height_raster)
  files$fused_dsm              <- write_raster_if_not_null(
    fused_dsm_local, files$fused_dsm)
  files$canopy_height_raster   <- write_raster_if_not_null(
    canopy_height_raster_local, files$canopy_height_raster)
  files$roughness_raster       <- write_raster_if_not_null(
    roughness_raster_local, files$roughness_raster)

  saveRDS(buildings_local, files$buildings_rds)
  saveRDS(metrics,         files$metrics_rds)


  # Return value


  info <- list(
    case_dir           = normalizePath(case_dir, mustWork = FALSE),
    data_dir           = normalizePath(data_dir, mustWork = FALSE),
    files              = lapply(files, function(x)
      if (is.null(x)) NULL else normalizePath(x, mustWork = FALSE)),
    data               = list(
      buildings            = buildings_local,
      binary               = binary_local,
      building_height      = height_raster_local,
      fused_dsm            = fused_dsm_local,
      canopy_height        = canopy_height_raster_local,
      ground_roughness_z0  = roughness_raster_local,
      metrics              = metrics
    ),
    domain             = domain,
    origin             = origin,
    crs                = sf::st_crs(buildings),
    height_col         = height_col,
    n_buildings        = nrow(buildings_local),
    max_building_height = max_building_h,
    max_surface_height  = max_surface_h,
    include_tree_canopy = include_tree_canopy,
    canopy_source       = canopy_source,
    min_tree_height     = min_tree_height,
    landcover_source    = landcover_source,
    landcover_year      = landcover_year
  )

  saveRDS(info, files$metadata_rds)

  msg("OpenFOAM input preparation finished.")
  msg("Building count           : ", info$n_buildings)
  msg("Maximum building height  : ", round(info$max_building_height, 2), " m")
  msg("Maximum surface height   : ", round(info$max_surface_height,  2), " m")
  msg("Canopy height raster     : ",
      if (!is.null(files$canopy_height_raster)) files$canopy_height_raster
      else "(not produced)")
  msg("Ground roughness raster  : ",
      if (!is.null(files$roughness_raster)) files$roughness_raster
      else "(not produced)")
  msg("Suggested OpenFOAM domain:")
  msg(
    "  x = [", round(domain$xmin, 2), ", ", round(domain$xmax, 2), "],",
    "  y = [", round(domain$ymin, 2), ", ", round(domain$ymax, 2), "],",
    "  z = [", round(domain$zmin, 2), ", ", round(domain$zmax, 2), "]"
  )

  invisible(info)
}
