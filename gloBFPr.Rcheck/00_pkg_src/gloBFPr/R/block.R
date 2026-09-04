#' aggregate_block
#' @description Aggregate building-level metrics to block level using the output
#'   of [generate_block()]. Additive quantities (areas, volumes, population) are
#'   summed by default; per-building indices (shape metrics, elongation ratios,
#'   etc.) are averaged. Both defaults can be overridden per column via `.fns`.
#'   Two derived metrics are always added: `n_buildings` (count of buildings per
#'   block) and `coverage_ratio` (total building footprint area / block area).
#'
#' @param block_output list. The named list returned by [generate_block()],
#'   containing `$blocks` (sf polygons) and `$buildings` (sf with `block_id`).
#' @param .fns named list. Optional overrides mapping column names to aggregation
#'   functions, e.g. `list(vol = max, Height = median)`. Any column not named
#'   here uses the built-in default (sum or mean as described above).
#' @param population logical. If `TRUE`, fetch GHSL population at block level
#'   and add a `pop_total` column. Default `FALSE`.
#' @param population_year integer. GHSL population year to use when
#'   `population = TRUE`. One of 1975, 1980, ..., 2025, 2030. Default `2025`.
#' @param residential logical. If `TRUE`, fetch GHS built-up surface rasters at
#'   block level and compute `res_prop` (residential built-up surface fraction
#'   per block). Default `FALSE`. Note: per-building `res` flags are excluded
#'   from block-level aggregation - use this parameter for block-level
#'   residential proportion instead.
#' @param residential_year integer. GHS built-up surface year to use when
#'   `residential = TRUE`. One of 1975, 1980, ..., 2020, 2025, 2030. Default
#'   `2020`.
#' @param quiet logical. If `TRUE`, suppress cli messages. Default `FALSE`.
#' @returns The `blocks` sf object from `block_output` with one column per
#'   aggregated metric, plus `n_buildings` and `coverage_ratio`. When
#'   `population = TRUE`, also includes `pop_total`. When `residential = TRUE`,
#'   also includes `res_prop`.
#' @importFrom sf st_drop_geometry st_area st_transform st_crs
#' @importFrom cli cli_alert_info cli_alert_warning cli_alert_success
#' @export

aggregate_block <- function(block_output, .fns = NULL, population = FALSE,
                            population_year = 2025, residential = FALSE,
                            residential_year = 2020, quiet = FALSE) {
  if (!is.list(block_output) || inherits(block_output, "sf") ||
      !all(c("blocks", "buildings") %in% names(block_output))) {
    stop("`block_output` must be the list returned by `generate_block()`.", call. = FALSE)
  }

  blocks    <- block_output$blocks
  buildings <- block_output$buildings

  if (!"block_id" %in% names(buildings))
    stop("`buildings` must contain a `block_id` column.", call. = FALSE)

  # Drop buildings not assigned to any block
  b <- buildings[!is.na(buildings$block_id), ]
  if (nrow(b) == 0) {
    if (!quiet) cli::cli_alert_warning("No buildings have a block assignment.")
    return(blocks)
  }
  b_df <- sf::st_drop_geometry(b)

  # -- Determine numeric columns to aggregate --------------------------------
  # res / GHS sub-columns are per-building centroid extractions; exclude from
  # block aggregation - use residential = TRUE for block-level res_prop instead
  skip_cols    <- c("block_id", "group_id", "id", "source_id",
                    "res", "res_pct", "res_vals", "nres_vals", "total_built_vals")
  numeric_cols <- names(b_df)[vapply(b_df, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, skip_cols)

  # Logical columns (e.g. `res`) treated as numeric 0/1 so mean = proportion
  logical_cols <- names(b_df)[vapply(b_df, is.logical, logical(1))]
  logical_cols <- setdiff(logical_cols, skip_cols)
  for (lc in logical_cols) b_df[[lc]] <- as.integer(b_df[[lc]])

  all_cols <- c(numeric_cols, logical_cols)

  # -- Default aggregation rules ----------------------------------------------
  # Physical additive quantities -> sum; everything else -> mean
  sum_defaults <- c("g_area", "pmeter", "v_surf", "t_surf", "vol", "obb_vol", "pop_total")

  col_fns <- stats::setNames(
    lapply(all_cols, function(col) if (col %in% sum_defaults) sum else mean),
    all_cols
  )

  # Apply user overrides
  if (!is.null(.fns)) {
    unknown <- setdiff(names(.fns), all_cols)
    if (length(unknown) > 0 && !quiet)
      cli::cli_alert_warning(
        "Column(s) in `.fns` not found in buildings and will be ignored: {paste(unknown, collapse = ', ')}"
      )
    for (nm in intersect(names(.fns), all_cols)) col_fns[[nm]] <- .fns[[nm]]
  }

  # -- Aggregate -------------------------------------------------------------
  if (!quiet) cli::cli_alert_info("Aggregating {length(all_cols)} metric(s) to block level ...")

  agg <- lapply(all_cols, function(col) {
    tapply(b_df[[col]], b_df$block_id, col_fns[[col]], na.rm = TRUE)
  })
  names(agg) <- all_cols

  # -- Attach aggregated columns to blocks sf --------------------------------
  result <- blocks
  bid_chr <- as.character(result$block_id)

  for (col in all_cols) {
    vals <- agg[[col]]
    result[[col]] <- as.numeric(vals[bid_chr])
  }

  # -- Derived block metrics -------------------------------------------------
  # n_buildings: count of assigned buildings per block
  n_tbl <- table(b_df$block_id)
  result$n_buildings <- as.integer(n_tbl[bid_chr])

  # coverage_ratio: total building footprint / block polygon area
  # Use UTM for area accuracy regardless of input CRS
  utm_crs        <- get_utm_crs(get_bbox(result))
  blocks_utm     <- sf::st_transform(result, utm_crs)
  buildings_utm  <- sf::st_transform(b, utm_crs)

  bldg_area_sum  <- tapply(
    as.numeric(sf::st_area(buildings_utm)),
    sf::st_drop_geometry(buildings_utm)$block_id,
    sum, na.rm = TRUE
  )
  block_area     <- as.numeric(sf::st_area(blocks_utm))
  result$coverage_ratio <- as.numeric(bldg_area_sum[bid_chr]) / block_area

  # -- Block-level population ------------------------------------------------
  if (isTRUE(population)) {
    if (!quiet) cli::cli_alert_info("Fetching GHSL population at block level (year = {population_year}) ...")
    block_bbox <- get_bbox(blocks_utm)
    pop_result <- tryCatch(
      get_GHSpop(bbox = block_bbox, year = population_year,
                 polygons = blocks_utm, quiet = quiet),
      error = function(e) {
        if (!quiet) cli::cli_alert_warning("Population fetch failed: {conditionMessage(e)}")
        NULL
      }
    )
    if (!is.null(pop_result)) {
      result$pop_total <- pop_result$pop_total
    }
  }

  # -- Block-level GHS residential surface proportion ------------------------
  if (isTRUE(residential)) {
    if (!quiet) cli::cli_alert_info("Fetching GHS residential surface at block level (year = {residential_year}) ...")
    block_bbox <- get_bbox(blocks_utm)
    ghs <- tryCatch(
      get_GHSres(bbox = block_bbox, year = residential_year),
      error = function(e) {
        if (!quiet) cli::cli_alert_warning("GHS residential fetch failed: {conditionMessage(e)}")
        NULL
      }
    )
    if (!is.null(ghs)) {
      blocks_vect <- terra::vect(blocks_utm)
      ext_total <- suppressWarnings(terra::extract(ghs$total, blocks_vect, weights = TRUE))
      ext_res   <- suppressWarnings(terra::extract(ghs$res,   blocks_vect, weights = TRUE))

      val_total <- setdiff(names(ext_total), c("ID", "weight"))[1]
      val_res   <- setdiff(names(ext_res),   c("ID", "weight"))[1]

      # weighted sum of residential and total built-up surface per block
      wres <- stats::aggregate(
        ext_res[[val_res]] * ext_res$weight,
        by = list(ID = ext_res$ID), FUN = sum, na.rm = TRUE
      )
      wtotal <- stats::aggregate(
        ext_total[[val_total]] * ext_total$weight,
        by = list(ID = ext_total$ID), FUN = sum, na.rm = TRUE
      )

      res_prop_vals <- ifelse(wtotal$x > 0, wres$x / wtotal$x, NA_real_)

      # map back: terra IDs are 1-indexed rows of blocks_utm (same order as result)
      result$res_prop <- NA_real_
      result$res_prop[wres$ID] <- res_prop_vals
    }
  }

  if (!quiet) cli::cli_alert_success("Done. Returning blocks with {ncol(result) - 2} aggregated metric(s).")

  result
}

#' generate_block
#' @description Cluster given buildings into blocks based on street network.
#'   Uses a two-stage approach: (1) vector polygonization of the road network
#'   for well-formed areas, then (2) a raster fallback for buildings that fall
#'   in network gaps or dead-end pockets. Blocks with no buildings are dropped.
#'
#'   Before polygonization, dual carriageways (motorways and trunk roads
#'   represented as parallel lines) are simplified to single centrelines to
#'   prevent artificially narrow slivers between them from being misidentified
#'   as blocks. The simplification approach is conceptually adapted from
#'   UrbanWaterBlocks (Yin et al., 2025).
#'
#' @param x sf. Building footprint polygons, typically output from [search_3dglobdf()].
#' @param network sf or character. Optional road network or path to one. When
#'   supplied, `network_source` is ignored.
#' @param network_source character. Source for automatic network fetching when
#'   `network = NULL`. Either `"overture"` (Overture Maps via DuckDB parquet
#'   query, default, requires the `duckdb` and `DBI` packages) or `"osm"`
#'   (OpenStreetMap via the Overpass API, slower for large areas).
#' @param overture_release character or `NULL`. Overture Maps release string
#'   used when `network_source = "overture"`, e.g. `"2025-03-19.0"`. The
#'   default `NULL` (or `"auto"`) queries the public bucket for the newest
#'   release and caches it for the session. See
#'   <https://github.com/OvertureMaps/data/releases> for available releases.
#' @param res numeric. Raster resolution in metres for the fallback stage. Default `2`.
#' @param min_block_area numeric. Minimum block area in m^2 below which polygons
#'   are treated as slivers and merged into neighbours (raster stage) or dropped
#'   (polygonize stage). Default `500`.
#' @param dc_highway_types character vector. Highway class values treated as
#'   dual carriageway candidates (OSM `highway` tag or Overture `class` column).
#'   Default `c("motorway", "trunk")`.
#' @param dc_overlap_threshold numeric. Minimum overlap fraction (0-1) for a
#'   line to be considered a duplicate carriageway and removed. Default `0.7`.
#' @param quiet logical. If `TRUE`, suppress cli messages. Default `FALSE`.
#' @returns A named list with two elements:
#'   \describe{
#'     \item{`blocks`}{An `sf` polygon object, one row per block, with a
#'       `block_id` column. CRS matches the input `x`.}
#'     \item{`buildings`}{The input `x` with an added `block_id` integer column
#'       linking each building to its block. Buildings that cannot be assigned
#'       to any block receive `NA`.}
#'   }
#' @references
#'   Yin, H., et al. (2025). UrbanWaterBlocks: A python tool for block-based
#'   urban water management. Sustainable Cities and Society.
#' @importFrom sf st_transform st_centroid st_join st_within st_as_sf st_geometry
#' @importFrom sf st_union st_polygonize st_collection_extract st_area st_buffer
#' @importFrom sf st_crs st_sf st_read st_length st_drop_geometry
#' @importFrom terra rast rasterize vect patches zonal cellSize ifel focal cover
#' @importFrom terra as.polygons extract align ext
#' @importFrom cli cli_alert_info cli_alert_warning cli_alert_success
#' @export

generate_block <- function(x,
                           network              = NULL,
                           network_source       = c("overture", "osm"),
                           overture_release     = NULL,
                           res                  = 2,
                           min_block_area       = 500,
                           dc_highway_types     = c("motorway", "trunk"),
                           dc_overlap_threshold = 0.7,
                           quiet                = FALSE) {
  if (is.null(x)) {
    if (!quiet) cli::cli_alert_info("Please input building footprint polygon generated by `search_3dglobdf()`.")
    return(x)
  }

  building   <- x
  input_crs  <- sf::st_crs(x)
  bbox       <- get_bbox(building)
  utm_crs    <- get_utm_crs(bbox)
  building_proj <- sf::st_transform(building, utm_crs)

  # -- 1. Get / load street network ------------------------------------------
  if (!is.null(network)) {
    net_sf <- if (is.character(network)) sf::st_read(network, quiet = TRUE) else network
  } else {
    network_source <- match.arg(network_source)
    bb <- as_wgs84_bbox_vector(bbox)
    if (network_source == "overture") {
      if (!quiet) cli::cli_alert_info("Fetching road network from Overture Maps ...")
      net_sf <- fetch_overture_network(bb, release = overture_release)
    } else {
      if (!requireNamespace("osmdata", quietly = TRUE))
        stop("Install 'osmdata' or supply a network sf object via the `network` argument.", call. = FALSE)
      if (!quiet) cli::cli_alert_info("Fetching road network from OSM ...")
      net_sf <- fetch_osm_network(bb, quiet = quiet)
    }
  }

  net_proj <- sf::st_transform(net_sf, utm_crs)

  # -- 2. Dual carriageway simplification ------------------------------------
  # Adapted from UrbanWaterBlocks (Yin et al., 2025):
  # For each DC highway type, buffer the longest line in a spatial cluster and
  # drop neighbours whose length overlaps the buffer by > dc_overlap_threshold.
  if (!quiet) cli::cli_alert_info("Simplifying dual carriageways ...")
  net_proj <- simplify_dual_carriageways(
    net_proj,
    highway_types     = dc_highway_types,
    overlap_threshold = dc_overlap_threshold
  )

  net_geom <- net_proj[, "geometry"]

  # -- 3. Polygonize: fast-path vector approach -------------------------------
  if (!quiet) cli::cli_alert_info("Polygonizing road network ...")
  enclosures_geom <- sf::st_polygonize(sf::st_union(sf::st_geometry(net_geom)))
  enclosures_geom <- sf::st_collection_extract(enclosures_geom, "POLYGON")
  enclosures      <- sf::st_as_sf(enclosures_geom)
  names(enclosures)[1] <- "geometry"
  sf::st_geometry(enclosures) <- "geometry"

  # Drop slivers
  enclosures <- enclosures[as.numeric(sf::st_area(enclosures)) >= min_block_area, ]

  # -- 4. Assign buildings to enclosures via centroid-within join -------------
  centroids <- suppressWarnings(sf::st_centroid(building_proj))
  join <- sf::st_join(
    sf::st_sf(bldg_idx = seq_len(nrow(building_proj)),
              geometry = sf::st_geometry(centroids)),
    sf::st_sf(enc_id   = seq_len(nrow(enclosures)),
              geometry = sf::st_geometry(enclosures)),
    join = sf::st_within,
    left = TRUE
  )

  covered_mask  <- !is.na(join$enc_id)
  uncovered_idx <- which(!covered_mask)

  # -- 5. Raster fallback for unassigned buildings ----------------------------
  raster_blocks_sf <- NULL
  patch_vals_unc   <- NULL

  if (length(uncovered_idx) > 0) {
    if (!quiet)
      cli::cli_alert_info(
        "{length(uncovered_idx)} building(s) not covered by polygonization - applying raster fallback ..."
      )

    template <- terra::rast(
      ext        = terra::align(terra::ext(sf::st_bbox(building_proj)), res),
      resolution = res,
      crs        = sf::st_crs(building_proj)$wkt
    )

    # Optionally prune dead-end branches before rasterizing
    net_for_raster <- net_geom
    if (requireNamespace("sfnetworks", quietly = TRUE) &&
        requireNamespace("tidygraph", quietly = TRUE)) {
      net_for_raster <- tryCatch({
        sfnetworks::as_sfnetwork(net_geom, directed = FALSE) |>
          tidygraph::convert(sfnetworks::to_spatial_smooth) |>
          sfnetworks::activate("edges") |>
          sf::st_as_sf()
      }, error = function(e) net_geom)
    }

    net_buf   <- sf::st_buffer(net_for_raster, dist = res)
    road_rast <- terra::rasterize(terra::vect(net_buf), template,
                                  field = 1, background = NA)

    non_road   <- terra::ifel(is.na(road_rast), 1, NA)
    patch_rast <- terra::patches(non_road, directions = 8)

    cell_area  <- terra::cellSize(patch_rast, unit = "m")
    patch_area <- terra::zonal(cell_area, patch_rast, fun = "sum", na.rm = TRUE)
    small_ids  <- patch_area[[1]][patch_area[[2]] < min_block_area]

    if (length(small_ids) > 0) {
      patch_clean <- terra::ifel(!is.na(terra::match(patch_rast, small_ids)), NA, patch_rast)
      filled      <- terra::focal(patch_clean, w = 3, fun = "modal",
                                  na.policy = "only", na.rm = TRUE)
      patch_rast  <- terra::cover(patch_clean, filled)
    }

    unc_centroids  <- suppressWarnings(sf::st_centroid(building_proj[uncovered_idx, ]))
    patch_vals_unc <- terra::extract(patch_rast, terra::vect(unc_centroids))[, 2]

    building_patch_ids <- unique(patch_vals_unc[!is.na(patch_vals_unc)])
    if (length(building_patch_ids) > 0) {
      patch_with_bldg  <- terra::ifel(!is.na(terra::match(patch_rast, building_patch_ids)), patch_rast, NA)
      raster_blocks_sf <- terra::as.polygons(patch_with_bldg) |>
        sf::st_as_sf() |>
        sf::st_transform(utm_crs)
      names(raster_blocks_sf)[1] <- "patch_val"
    }
  }

  # -- 6. Combine enclosure + raster blocks; drop those without buildings -----
  enc_with_bldg <- unique(join$enc_id[covered_mask])

  enc_blocks <- if (length(enc_with_bldg) > 0) {
    b <- enclosures[enc_with_bldg, ]
    b$block_id <- seq_len(nrow(b))
    b[, c("block_id", "geometry")]
  } else {
    NULL
  }

  n_enc <- if (!is.null(enc_blocks)) nrow(enc_blocks) else 0L

  if (!is.null(raster_blocks_sf)) {
    raster_blocks_sf$block_id <- n_enc + seq_len(nrow(raster_blocks_sf))
    raster_blocks_sf <- raster_blocks_sf[, c("block_id", "geometry")]
  }

  all_blocks <- if (!is.null(enc_blocks) && !is.null(raster_blocks_sf)) {
    rbind(enc_blocks, raster_blocks_sf)
  } else if (!is.null(enc_blocks)) {
    enc_blocks
  } else {
    raster_blocks_sf
  }

  if (is.null(all_blocks) || nrow(all_blocks) == 0) {
    if (!quiet) cli::cli_alert_warning("No blocks could be generated.")
    return(NULL)
  }

  # -- 7. Assign block_id back to each building ------------------------------
  block_ids <- rep(NA_integer_, nrow(building_proj))

  if (!is.null(enc_blocks) && any(covered_mask)) {
    enc_id_to_block         <- stats::setNames(enc_blocks$block_id, enc_with_bldg)
    block_ids[covered_mask] <- enc_id_to_block[as.character(join$enc_id[covered_mask])]
  }

  if (!is.null(raster_blocks_sf) && length(uncovered_idx) > 0) {
    patch_to_block           <- stats::setNames(raster_blocks_sf$block_id,
                                                raster_blocks_sf$patch_val)
    block_ids[uncovered_idx] <- patch_to_block[as.character(patch_vals_unc)]
  }

  building$block_id <- block_ids

  # -- 8. Reproject outputs to input CRS -------------------------------------
  all_blocks <- sf::st_transform(all_blocks, input_crs)

  if (!quiet) cli::cli_alert_success("Generated {nrow(all_blocks)} block(s).")

  list(blocks = all_blocks, buildings = building)
}

#' simplify_dual_carriageways
#' @description Remove redundant parallel lines representing dual carriageways
#'   from a road network sf object. For each group of spatially proximate lines
#'   of the target highway types, the longest line is kept and shorter neighbours
#'   are removed if more than `overlap_threshold` of their length overlaps the
#'   buffer of the longest line.
#'
#'   The approach is conceptually adapted from UrbanWaterBlocks (Yin et al., 2025).
#'
#' @param net sf. Line geometries of the road network, in a projected (metric) CRS.
#' @param highway_types character vector. Values of the `highway` column to
#'   treat as dual carriageway candidates. Default `c("motorway", "trunk")`.
#' @param overlap_threshold numeric. Fraction of a line's length that must fall
#'   inside the buffer of the longest neighbour to be considered a duplicate.
#'   Default `0.7`.
#' @returns The input `net` sf with duplicate carriageway lines removed.
#' @references
#'   Yin, H., et al. (2025). UrbanWaterBlocks: A python tool for block-based
#'   urban water management. Sustainable Cities and Society.
#' @noRd

fetch_osm_network <- function(bbox, max_retries = 6, quiet = FALSE) {
  highway_values <- c("motorway", "trunk", "primary", "secondary", "tertiary",
                      "residential", "unclassified", "living_street", "road")
  last_error <- NULL

  for (attempt in seq_len(max_retries)) {
    result <- rlang::try_fetch(
      osmdata::opq(bbox = bbox) |>
        osmdata::add_osm_feature(key = "highway", value = highway_values) |>
        osmdata::osmdata_sf() |>
        (`[[`)("osm_lines"),
      error = function(e) e
    )

    if (!inherits(result, "error")) {
      if (is.null(result) || nrow(result) == 0)
        stop("No road network retrieved for the bounding box.", call. = FALSE)
      return(result)
    }

    last_error <- result
    wait <- 2^(attempt - 1L)   # 1, 2, 4, 8, 16, 32 seconds

    if (!quiet)
      cli::cli_alert_warning(
        "OSM request failed (attempt {attempt}/{max_retries}): {conditionMessage(result)}. Retrying in {wait}s ..."
      )

    Sys.sleep(wait)
  }

  stop(
    "Failed to fetch road network from OSM after ", max_retries, " attempts. ",
    "Last error: ", conditionMessage(last_error), "\n",
    "Consider using network_source = 'overture' or supplying a pre-downloaded network.",
    call. = FALSE
  )
}

#' @noRd

fetch_overture_network <- function(bbox, release = NULL) {
  if (!requireNamespace("duckdb", quietly = TRUE))
    stop("Install the 'duckdb' package to use network_source = 'overture'.", call. = FALSE)

  xmin <- bbox[["xmin"]]; ymin <- bbox[["ymin"]]
  xmax <- bbox[["xmax"]]; ymax <- bbox[["ymax"]]

  if (!requireNamespace("DBI", quietly = TRUE))
    stop("Install the 'DBI' package to use network_source = 'overture'.", call. = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  DBI::dbExecute(con, "INSTALL httpfs;  LOAD httpfs;")
  DBI::dbExecute(con, "INSTALL spatial; LOAD spatial;")
  DBI::dbExecute(con, "SET s3_region = 'us-west-2';")

  # Drivable road classes - mirrors the OSM highway tags used by the OSM path
  keep_classes <- paste0("'",
    paste(c("motorway", "motorway_link",
            "trunk",    "trunk_link",
            "primary",  "primary_link",
            "secondary","secondary_link",
            "tertiary", "tertiary_link",
            "residential", "unclassified",
            "living_street", "road"),
          collapse = "','"),
    "'")

  release <- resolve_overture_release(con, release)

  s3_path <- sprintf(
    "s3://overturemaps-us-west-2/release/%s/theme=transportation/type=segment/*",
    release
  )

  query <- sprintf("
    SELECT
      ST_AsWKB(geometry) AS wkb,
      class               AS highway
    FROM read_parquet('%s', hive_partitioning = 1)
    WHERE subtype = 'road'
      AND class IN (%s)
      AND bbox.xmax >= %f AND bbox.xmin <= %f
      AND bbox.ymax >= %f AND bbox.ymin <= %f
  ", s3_path, keep_classes, xmin, xmax, ymin, ymax)

  res <- DBI::dbGetQuery(con, query)

  if (nrow(res) == 0)
    stop("No road segments returned from Overture Maps for the bounding box.", call. = FALSE)

  # Convert WKB -> sfc, then build sf
  geom <- sf::st_as_sfc(
    structure(as.list(res$wkb), class = "WKB"),
    crs = 4326
  )
  sf::st_sf(highway = res$highway, geometry = geom)
}

#' @importFrom sf st_length st_buffer st_intersection st_geometry st_drop_geometry
#' @importFrom sf st_intersects
#' @noRd

simplify_dual_carriageways <- function(net,
                                       highway_types     = c("motorway", "trunk"),
                                       overlap_threshold = 0.7) {
  if (!"highway" %in% names(net)) return(net)

  # Separate DC candidates from the rest
  is_dc  <- net$highway %in% highway_types
  dc     <- net[is_dc, ]
  other  <- net[!is_dc, ]

  if (nrow(dc) == 0) return(net)

  lengths_dc <- as.numeric(sf::st_length(dc))

  # Buffer width by highway type (metres); wider for motorways
  buffer_width <- ifelse(dc$highway == "motorway", 25, 15)

  # Find spatial clusters: lines that intersect each other's buffers
  dc_buf     <- sf::st_buffer(dc, dist = buffer_width)
  neighbours <- sf::st_intersects(dc, dc_buf, sparse = TRUE)

  keep <- rep(TRUE, nrow(dc))

  for (i in seq_len(nrow(dc))) {
    if (!keep[i]) next

    nb_idx <- neighbours[[i]]
    nb_idx <- nb_idx[nb_idx != i & keep[nb_idx]]
    if (length(nb_idx) == 0) next

    # Among this line and its neighbours, find the longest
    group_idx    <- c(i, nb_idx)
    group_len    <- lengths_dc[group_idx]
    longest_idx  <- group_idx[which.max(group_len)]

    if (longest_idx != i) next  # will be handled when we reach the longest

    # Buffer the longest line and measure overlap for each shorter neighbour
    buf_longest <- sf::st_buffer(dc[longest_idx, ], dist = buffer_width[longest_idx])

    for (j in nb_idx) {
      seg_len <- lengths_dc[j]
      if (seg_len == 0) next
      overlap_geom <- suppressWarnings(
        sf::st_intersection(sf::st_geometry(dc[j, ]), sf::st_geometry(buf_longest))
      )
      if (length(overlap_geom) == 0 || all(sf::st_is_empty(overlap_geom))) next
      overlap_len  <- tryCatch(
        as.numeric(sf::st_length(overlap_geom)),
        error = function(e) 0
      )
      if (sum(overlap_len, na.rm = TRUE) / seg_len >= overlap_threshold) {
        keep[j] <- FALSE
      }
    }
  }

  rbind(dc[keep, ], other)
}
