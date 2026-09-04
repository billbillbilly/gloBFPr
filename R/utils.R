#' @importFrom igraph make_empty_graph add_edges components
#' @importFrom terra ext
#' @importFrom terra vect
#' @importFrom terra rast
#' @importFrom terra rasterize
#' @importFrom terra align
#' @importFrom terra mask project
#' @importFrom sf st_centroid st_within
#' @importFrom sf st_coordinates
#' @importFrom sf st_transform
#' @importFrom sf st_crs st_convex_hull
#' @importFrom sf st_as_sfc st_multipoint

#### Footprint data processing ####
#' @noRd
normalize_building_data_source <- function(data_source) {
  if (inherits(data_source, "NULL") || length(data_source) == 0 || is.na(data_source[1])) {
    stop('Invalid data_source specified. Use "GBF" or "GBA".', call. = FALSE)
  }
  data_source <- toupper(as.character(data_source[1]))
  if (data_source %in% c("GBA", "GBF")) {
    return(data_source)
  }
  stop('Invalid data_source specified. Use "GBF" or "GBA".', call. = FALSE)
}

#' @noRd
normalize_building_geometries <- function(x) {
  if (is.null(x) || nrow(x) == 0) {
    return(NULL)
  }
  x <- x[!sf::st_is_empty(x), ]
  if (nrow(x) == 0) {
    return(NULL)
  }
  geometry <- sf::st_geometry(x)
  geometry <- lapply(geometry, largest_polygon_part)
  geometry <- sf::st_sfc(geometry, crs = sf::st_crs(x))
  sf::st_geometry(x) <- geometry
  x
}

#' @noRd
drop_internal_building_columns <- function(x, keep_source_id = FALSE) {
  if (isTRUE(keep_source_id) && ".source_id" %in% names(x)) {
    x$source_id <- x$.source_id
  }
  internal_cols <- intersect(names(x), c(".source_id"))
  if (length(internal_cols) > 0) {
    x <- x[, setdiff(names(x), internal_cols)]
  }
  x
}

#' @noRd
assign_building_group_id <- function(x) {
  if (is.null(x) || nrow(x) == 0) {
    return(x)
  }

  # Buildings that touch or intersect are treated as one connected group.
  neighbors <- sf::st_intersects(x, x, sparse = TRUE)

  # Build upper-triangle edge list and find connected components via igraph,
  # which is much faster than an R-level BFS loop for large datasets.
  edges <- do.call(rbind, lapply(seq_along(neighbors), function(i) {
    j <- neighbors[[i]]
    j <- j[j > i]
    if (length(j) == 0) return(NULL)
    cbind(i, j)
  }))

  g <- igraph::make_empty_graph(n = nrow(x), directed = FALSE)
  if (!is.null(edges)) {
    g <- igraph::add_edges(g, t(edges))
  }
  x$group_id <- igraph::components(g)$membership
  x
}

#' @noRd
largest_polygon_part <- function(geometry) {
  polygons <- suppressWarnings(sf::st_cast(sf::st_sfc(geometry), "POLYGON"))
  if (length(polygons) == 0) {
    return(sf::st_polygon())
  }
  if (length(polygons) == 1) {
    return(polygons[[1]])
  }
  areas <- as.numeric(sf::st_area(polygons))
  polygons[[which.max(areas)]]
}

#' @noRd
read_gbf_buildings <- function(bbox, quiet = TRUE) {
  # find all areas of spatial grid that intersect with bbox
  intersecting <- globfp3d_metadata[sf::st_intersects(globfp3d_metadata, bbox, sparse = FALSE), ]

  if (nrow(intersecting) == 0) {
    return(NULL)
  }

  # Store the original 'timeout' option and ensure it's reset upon function exit
  original_timeout <- getOption('timeout')
  on.exit(options(timeout = original_timeout), add = TRUE)
  options(timeout=9999)

  workers <- min(max(1, as.numeric(future::availableCores()) - 1), nrow(intersecting))
  os <- Sys.info()[["sysname"]]

  bbox_wkt <- sf::st_as_text(sf::st_union(sf::st_geometry(bbox)))

  if (workers == 1) {
    # Single tile: run directly instead of through a (sequential) future. This
    # also avoids future's "added, removed, or modified connections" warning,
    # which is triggered by connections opened inside download.file()/st_read().
    result_list <- lapply(
      intersecting$download_url,
      read_gbf_tile,
      bbox_wkt = bbox_wkt,
      quiet = quiet
    )
  } else {
    # future (>= 1.40) warns when a future expression leaves a connection
    # behind. Here the connection is created by download.file()/sf::st_read(),
    # i.e. outside our control, so the check is silenced locally.
    old_misuse <- getOption("future.connections.onMisuse")
    options(future.connections.onMisuse = "ignore")
    on.exit(options(future.connections.onMisuse = old_misuse), add = TRUE)

    # Restore whatever plan the user had, rather than forcing sequential.
    oplan <- future::plan()
    on.exit(future::plan(oplan), add = TRUE)
    if (os == "Windows" || isFALSE(parallelly::supportsMulticore())) {
      future::plan(future::multisession, workers = workers)
    } else {
      future::plan(future::multicore, workers = workers)
    }

    result_list <- furrr::future_map(
      intersecting$download_url,
      read_gbf_tile,
      bbox_wkt = bbox_wkt,
      quiet = quiet,
      .options = furrr::furrr_options(seed = TRUE),
      .progress = !quiet
    )
  }
  result_list <- Filter(Negate(is.null), result_list)

  if (length(result_list) == 0) {
    return(NULL)
  }

  result_list <- lapply(result_list, function(x) {
    #x <- sf::st_cast(x, "POLYGON")  # ensure same geometry type
    x <- x[, intersect(names(x), names(result_list[[1]]))]  # keep common columns only
    return(x)
  })
  # Combine all into one sf object
  all_data <- dplyr::bind_rows(result_list)
  all_data <- all_data[,c('Height','geometry')]
  return(all_data)
}

#' @noRd
read_gbf_tile <- function(download_url, bbox_wkt, quiet = TRUE) {
  d_mode <- 'auto'
  if (Sys.info()[["sysname"]] == "Windows") {
    d_mode <- 'wininet'
  }

  temp_zip <- tempfile(fileext = ".zip")
  unzip_dir <- tempfile()
  on.exit(unlink(c(temp_zip, unzip_dir), recursive = TRUE), add = TRUE)

  utils::download.file(download_url,
                       destfile = temp_zip,
                       method = d_mode,
                       quiet = TRUE)

  utils::unzip(temp_zip, exdir = unzip_dir)

  shp_files <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE)
  if (length(shp_files) == 0) {
    return(NULL)
  }

  tryCatch({
    sf::st_read(shp_files[1], quiet = TRUE, wkt_filter = bbox_wkt)
  }, error = function(e) {
    if (!quiet) base::message("Failed to read shapefile: ", shp_files[1])
    return(NULL)
  })
}

#' @noRd
read_gba_buildings <- function(bbox, quiet = TRUE) {
  if (!quiet) cli::cli_alert_info("Downloading GlobalBuildingAtlas WFS features ...")
  all_data <- read_gba_bbox_recursive(bbox, quiet = quiet)
  if (is.null(all_data) || nrow(all_data) == 0) {
    return(NULL)
  }
  all_data <- normalize_gba_columns(all_data)
  all_data <- dedupe_gba_features(all_data)
  keep_cols <- intersect(c("Height", ".source_id", "geometry"), names(all_data))
  all_data <- all_data[, keep_cols]
  return(all_data)
}

#' @noRd
read_gba_bbox_recursive <- function(bbox, quiet = TRUE, depth = 0L, max_depth = 5L) {
  if (gba_bbox_needs_split(bbox) && depth < max_depth) {
    pieces <- lapply(split_bbox_quadrants(bbox), function(tile) {
      read_gba_bbox_recursive(tile, quiet = quiet, depth = depth + 1L, max_depth = max_depth)
    })
    pieces <- Filter(Negate(is.null), pieces)
    if (length(pieces) == 0) {
      return(NULL)
    }
    return(dplyr::bind_rows(pieces))
  }

  result <- tryCatch(
    read_gba_bbox(bbox),
    error = function(e) e
  )
  if (!inherits(result, "error")) {
    return(result)
  }

  if (depth >= max_depth) {
    stop("Failed to read GlobalBuildingAtlas WFS data: ", conditionMessage(result), call. = FALSE)
  }
  if (grepl("404|not found", conditionMessage(result), ignore.case = TRUE)) {
    return(NULL)
  }
  if (!quiet) {
    cli::cli_alert_info("GlobalBuildingAtlas WFS request failed; retrying with smaller bbox tiles ...")
  }

  pieces <- lapply(split_bbox_quadrants(bbox), function(tile) {
    tryCatch(
      read_gba_bbox_recursive(tile, quiet = quiet, depth = depth + 1L, max_depth = max_depth),
      error = function(e) {
        if (!quiet) base::message("Failed to read GBA tile: ", conditionMessage(e))
        NULL
      }
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0) {
    return(NULL)
  }
  dplyr::bind_rows(pieces)
}

#' @noRd
read_gba_bbox <- function(bbox) {
  bbox_values <- sf::st_bbox(bbox)
  bbox_param <- paste(
    bbox_values[["xmin"]],
    bbox_values[["ymin"]],
    bbox_values[["xmax"]],
    bbox_values[["ymax"]],
    "EPSG:4326",
    sep = ","
  )
  params <- c(
    service = "WFS",
    version = "2.0.0",
    request = "GetFeature",
    typeNames = "global3D:lod1_global",
    outputFormat = "application/json",
    srsName = "EPSG:4326",
    BBOX = bbox_param
  )
  query <- paste(
    names(params),
    vapply(params, utils::URLencode, character(1), reserved = TRUE),
    sep = "=",
    collapse = "&"
  )
  url <- paste0("https://tubvsig-so2sat-vm1.srv.mwn.de/geoserver/ows?", query)

  response <- httr2::request(url) |>
    httr2::req_timeout(20) |>
    httr2::req_perform()
  body <- httr2::resp_body_string(response)
  if (grepl('"features"[[:space:]]*:[[:space:]]*\\[[[:space:]]*\\]', body)) {
    return(NULL)
  }

  geojson_file <- tempfile(fileext = ".geojson")
  on.exit(unlink(geojson_file), add = TRUE)
  writeLines(body, geojson_file, useBytes = TRUE)
  sf::st_read(geojson_file, quiet = TRUE)
}

#' @noRd
gba_bbox_needs_split <- function(bbox, max_span = 0.01) {
  bbox_values <- sf::st_bbox(bbox)
  width <- bbox_values[["xmax"]] - bbox_values[["xmin"]]
  height <- bbox_values[["ymax"]] - bbox_values[["ymin"]]
  width > max_span || height > max_span
}

#' @noRd
split_bbox_quadrants <- function(bbox) {
  bbox_values <- sf::st_bbox(bbox)
  xmid <- mean(c(bbox_values[["xmin"]], bbox_values[["xmax"]]))
  ymid <- mean(c(bbox_values[["ymin"]], bbox_values[["ymax"]]))
  boxes <- list(
    c(bbox_values[["xmin"]], bbox_values[["ymin"]], xmid, ymid),
    c(xmid, bbox_values[["ymin"]], bbox_values[["xmax"]], ymid),
    c(bbox_values[["xmin"]], ymid, xmid, bbox_values[["ymax"]]),
    c(xmid, ymid, bbox_values[["xmax"]], bbox_values[["ymax"]])
  )
  lapply(boxes, function(values) {
    sf::st_as_sfc(sf::st_bbox(
      c(xmin = values[1], ymin = values[2], xmax = values[3], ymax = values[4]),
      crs = sf::st_crs(bbox)
    ))
  })
}

#' @noRd
normalize_gba_columns <- function(all_data) {
  if (nrow(all_data) == 0) {
    return(NULL)
  }
  if ("id" %in% names(all_data)) {
    all_data$.source_id <- all_data$id
  } else if ("ogc_fid" %in% names(all_data)) {
    all_data$.source_id <- all_data$ogc_fid
  }
  if ("height" %in% names(all_data) && !"Height" %in% names(all_data)) {
    names(all_data)[names(all_data) == "height"] <- "Height"
  }
  if (!"Height" %in% names(all_data)) {
    stop("GlobalBuildingAtlas response is missing the height field.", call. = FALSE)
  }
  all_data
}

#' @noRd
dedupe_gba_features <- function(all_data) {
  if (is.null(all_data) || nrow(all_data) <= 1) {
    return(all_data)
  }
  if (".source_id" %in% names(all_data)) {
    source_id <- as.character(all_data$.source_id)
    has_source_id <- !is.na(source_id) & nzchar(source_id)
    keep <- rep(TRUE, nrow(all_data))
    keep[has_source_id] <- !duplicated(source_id[has_source_id])
    all_data <- all_data[keep, ]
  }
  if (nrow(all_data) <= 1) {
    return(all_data)
  }
  # Use digits = 6 (~0.1 m precision) so near-identical polygons from adjacent
  # GBA tile boundaries are treated as duplicates and removed.
  exact_key <- paste(
    as.character(all_data$Height),
    sf::st_as_text(sf::st_geometry(all_data), digits = 6),
    sep = "|"
  )
  all_data[!duplicated(exact_key), ]
}

#' @noRd
rasterize_binary <- function(poly, bbox, res) {
  utm_crs <- get_utm_crs(bbox)
  proj_ <- reproj(bbox, poly, utm_crs, res)
  bbox_proj <- proj_[[1]]
  poly_proj <- proj_[[2]]
  bbox_raster <- proj_[[3]]

  template <- terra::rast(ext = bbox_raster,
                          resolution = res,
                          crs = sf::st_crs(bbox_proj)$wkt)
  binary <- terra::rasterize(terra::vect(poly_proj),
                             template,
                             field = 1,
                             background = 0)
  return(binary)
}

#' @noRd
rasterize_height <- function(poly, bbox, res, mask=NULL, height_field = "Height") {
  if (!height_field %in% names(poly)) {stop("Missing height field in polygon.")}
  utm_crs <- get_utm_crs(bbox)
  proj_ <- reproj(bbox, poly, utm_crs, res)
  bbox_proj <- proj_[[1]]
  poly_proj <- terra::vect(proj_[[2]])
  bbox_raster <- proj_[[3]]

  template <- terra::rast(ext = bbox_raster,
                          resolution = res,
                          crs = sf::st_crs(bbox_proj)$wkt)
  graduated <- terra::rasterize(poly_proj,
                                template,
                                field = height_field,
                                fun = "max",
                                background = 0)
  if (inherits(mask, "SpatRaster")) {
    mask <- terra::ifel(mask == 1, 1, NA)
    graduated <- terra::mask(graduated, mask)
  }
  return(graduated)
}

#### Metrics calculation ####

#' @noRd
prepare_group_analysis_buildings <- function(x) {
  if (!"group_id" %in% names(x)) {
    return(x)
  }

  group_values <- unique(x$group_id)
  group_indices <- lapply(group_values, function(g) which(x$group_id == g))
  group_geometries <- lapply(group_indices, function(idx) {
    sf::st_union(sf::st_geometry(x[idx, ]))[[1]]
  })
  group_geometry <- sf::st_sfc(group_geometries, crs = sf::st_crs(x))
  group_height <- vapply(group_indices, function(idx) max(x$Height[idx], na.rm = TRUE), numeric(1))

  units <- sf::st_sf(
    id = seq_along(group_values),
    group_id = group_values,
    Height = group_height,
    geometry = group_geometry
  )
  units
}

#' @noRd
copy_group_results <- function(target, grouped, cols) {
  cols <- intersect(cols, names(grouped))
  if (length(cols) == 0) {
    return(target)
  }
  if ("group_id" %in% names(target) && "group_id" %in% names(grouped)) {
    match_idx <- match(target$group_id, grouped$group_id)
    for (col in cols) {
      target[[col]] <- grouped[[col]][match_idx]
    }
  } else {
    for (col in cols) {
      target[[col]] <- grouped[[col]]
    }
  }
  target
}

#' @noRd
prepare_morphology_units <- function(x) {
  if (!"group_id" %in% names(x)) {
    return(x)
  }

  group_values <- unique(x$group_id)
  group_indices <- lapply(group_values, function(g) which(x$group_id == g))
  group_geometries <- lapply(group_indices, function(idx) {
    sf::st_union(sf::st_geometry(x[idx, ]))[[1]]
  })
  group_geometry <- sf::st_sfc(group_geometries, crs = sf::st_crs(x))
  group_height <- vapply(group_indices, function(idx) max(x$Height[idx], na.rm = TRUE), numeric(1))

  units <- sf::st_sf(
    group_id = group_values,
    Height = group_height,
    geometry = group_geometry
  )
  units$.group_parts <- I(lapply(group_indices, function(idx) x[idx, ]))

  detail <- lapply(group_indices, morphology_group_properties, x = x)
  units$.g_area <- vapply(detail, function(d) d$g_area, numeric(1))
  units$.pmeter <- vapply(detail, function(d) d$pmeter, numeric(1))
  units$.v_surf <- vapply(detail, function(d) d$v_surf, numeric(1))
  units$.t_surf <- vapply(detail, function(d) d$t_surf, numeric(1))
  units$.vol <- vapply(detail, function(d) d$vol, numeric(1))
  units$.obb_vol <- vapply(detail, function(d) d$obb_vol, numeric(1))
  units
}

#' @noRd
morphology_group_properties <- function(idx, x) {
  parts <- x[idx, ]
  part_area <- as.numeric(sf::st_area(parts))
  part_perimeter <- as.numeric(sf::st_perimeter(parts))
  part_height <- as.numeric(parts$Height)
  union_geometry <- sf::st_union(sf::st_geometry(parts))
  g_area <- sum(part_area, na.rm = TRUE)
  pmeter <- as.numeric(sf::st_perimeter(union_geometry))
  v_surf <- group_vertical_surface(parts, part_perimeter, part_height)
  vol <- sum(part_area * part_height, na.rm = TRUE)
  obb_vol <- as.numeric(sf::st_area(sf::st_minimum_rotated_rectangle(union_geometry)) *
                          max(part_height, na.rm = TRUE))
  list(
    g_area = g_area,
    pmeter = pmeter,
    v_surf = v_surf,
    t_surf = v_surf + g_area,
    vol = vol,
    obb_vol = obb_vol
  )
}

#' @noRd
group_vertical_surface <- function(parts, part_perimeter, part_height) {
  v_surf <- sum(part_perimeter * part_height, na.rm = TRUE)
  if (nrow(parts) <= 1) {
    return(v_surf)
  }

  geometries <- sf::st_geometry(parts)
  for (i in seq_len(nrow(parts) - 1)) {
    for (j in (i + 1):nrow(parts)) {
      shared <- suppressWarnings(sf::st_intersection(
        sf::st_boundary(geometries[i]),
        sf::st_boundary(geometries[j])
      ))
      if (length(shared) == 0 || all(sf::st_is_empty(shared))) {
        next
      }
      shared_length <- sum(as.numeric(sf::st_length(shared)), na.rm = TRUE)
      if (!is.finite(shared_length) || shared_length <= 0) {
        next
      }
      v_surf <- v_surf - 2 * shared_length * min(part_height[i], part_height[j], na.rm = TRUE)
    }
  }
  v_surf
}

#' @noRd
get_morphology_property <- function(x, name) {
  detail_name <- paste0(".", name)
  if (detail_name %in% names(x)) {
    return(x[[detail_name]])
  }
  switch(
    name,
    g_area = as.numeric(sf::st_area(x)),
    pmeter = as.numeric(sf::st_perimeter(x)),
    v_surf = as.numeric(sf::st_perimeter(x)) * x$Height,
    t_surf = as.numeric(sf::st_perimeter(x)) * x$Height + as.numeric(sf::st_area(x)),
    vol = as.numeric(sf::st_area(x)) * x$Height,
    obb_vol = as.numeric(sf::st_area(sf::st_minimum_rotated_rectangle(x)) * x$Height),
    stop("Unknown morphology property: ", name, call. = FALSE)
  )
}

#' @noRd
cal_hemisphericality <- function(poly_) {
  # minimum_area_rectangle <- sf::st_minimum_rotated_rectangle(poly)
  # coords_2d <- sf::st_coordinates(minimum_area_rectangle)
  # coords_2d <- coords_2d[1:4, c("X", "Y")]
  # pt_1 <- as.vector(c(coords_2d[1,1], coords_2d[1,2]), 0)
  # pt_2 <- as.vector(c(coords_2d[1,1], coords_2d[1,2]), poly$Height)
  # pt_3 <- as.vector(c(coords_2d[3,1], coords_2d[3,2]), 0)
  # delta_1_2 <- pt_1 - pt_2
  # delta_1_3 <- pt_1 - pt_3
  # a2 <- delta_1_2[3]^2
  # b2 <- delta_1_3[1]^2 + delta_1_3[2]^2
  # radius <- sqrt(a2 + (sqrt(b2)/2)^2)
  # volume_hemisphere <- pi * radius^3 / 2

  base_area <- poly_$g_area
  building_volume <- poly_$vol
  # Compute radius of hemisphere with same base area
  radius <- sqrt(base_area / pi)
  # Volume of hemisphere with this radius
  vol_hemisphere <- (2/3) * pi * radius^3
  # deviation from ideal hemisphere volume
  hemisphericality <- (building_volume - vol_hemisphere) / vol_hemisphere
  return(hemisphericality)
}

#' @noRd
cal_convexity <- function(poly_) {
  vertices <- sf::st_coordinates(poly_)[, 1:2]
  multipoint <- sf::st_multipoint(vertices)
  multipoint <- sf::st_sfc(multipoint, crs = sf::st_crs(poly_))
  convex_hull <- sf::st_convex_hull(multipoint)
  return(as.numeric(poly_$g_area / sf::st_area(convex_hull)))
}

#' @noRd
cal_accessibility <- function(poly_) {
  # create a 3D voxel grid
  filtered_grid <- ploy2grid(poly_)
  # compute centroid of the volume
  centroid_x <- mean(filtered_grid$X)
  centroid_y <- mean(filtered_grid$Y)
  centroid_z <- mean(filtered_grid$Z)

  # compute mean Euclidean distance to the centroid
  dists <- sqrt((filtered_grid$X - centroid_x)^2 +
                  (filtered_grid$Y - centroid_y)^2 +
                  (filtered_grid$Z - centroid_z)^2)
  return(mean(dists))
}

#' @noRd
cal_accessibility_group <- function(parts) {
  filtered_grid <- group_ploy2grid(parts)
  if (nrow(filtered_grid) == 0) {
    return(NA_real_)
  }
  centroid_x <- mean(filtered_grid$X)
  centroid_y <- mean(filtered_grid$Y)
  centroid_z <- mean(filtered_grid$Z)
  dists <- sqrt((filtered_grid$X - centroid_x)^2 +
                  (filtered_grid$Y - centroid_y)^2 +
                  (filtered_grid$Z - centroid_z)^2)
  mean(dists)
}

#' @noRd
cal_mean_pairwise_distance <- function(poly_) {
  # create a 3D voxel grid
  filtered_grid <- ploy2grid(poly_)
 # Compute all pairwise distances
  coords <- as.matrix(filtered_grid)
  dist <- mean_pairwise_distance(coords)
  return(dist)
}

#' @noRd
cal_mean_pairwise_distance_group <- function(parts) {
  filtered_grid <- group_ploy2grid(parts)
  if (nrow(filtered_grid) == 0) {
    return(NA_real_)
  }
  coords <- as.matrix(filtered_grid)
  mean_pairwise_distance(coords)
}

#' @noRd
cal_volume_exchange_ratio <- function(poly_) {
  # Get volume of minimum enclosing sphere
  coords <- sf::st_coordinates(sf::st_convex_hull(poly_))[, 1:2]
  center <- base::colMeans(coords)
  radii <- sqrt(base::rowSums((coords - matrix(center, nrow(coords), 2, byrow = TRUE))^2))
  radius <- max(radii)
  vol_sphere <- (4/3) * pi * radius^3
  # Compute deviation
  vol_exch <- (vol_sphere - poly_$vol) / vol_sphere
  return(vol_exch)
}

#' @noRd
ploy2grid <- function(poly_) {
  height <- as.numeric(poly_$Height[1])
  bbox_ <- sf::st_bbox(poly_)
  dx <- dy <- 5  # horizontal spacing
  dz <- 2        # vertical spacing
  x_seq <- seq(bbox_["xmin"], bbox_["xmax"], by = dx)
  y_seq <- seq(bbox_["ymin"], bbox_["ymax"], by = dy)
  z_seq <- seq(0, height, by = dz)
  grid_3d <- expand.grid(X = x_seq, Y = y_seq, Z = z_seq)
  # keep only points that fall inside the footprint (x, y only)
  xy_points <- sf::st_as_sf(grid_3d, coords = c("X", "Y"), crs = st_crs(poly_))
  inside <- sf::st_within(xy_points, poly_, sparse = FALSE)[, 1]
  filtered_grid <- grid_3d[inside, ]
  return(filtered_grid)
}

#' @noRd
group_ploy2grid <- function(parts) {
  grids <- lapply(seq_len(nrow(parts)), function(i) ploy2grid(parts[i, ]))
  grids <- grids[vapply(grids, nrow, integer(1)) > 0]
  if (length(grids) == 0) {
    return(data.frame(X = numeric(), Y = numeric(), Z = numeric()))
  }
  unique(do.call(rbind, grids))
}

#' @noRd
cal_elongation_ratios <- function(poly_) {
  # Compute minimum bounding rectangle
  min_rect <- sf::st_minimum_rotated_rectangle(poly_)
  coords <- sf::st_coordinates(min_rect)[, 1:2]

  # Get edge lengths
  edges <- sqrt(rowSums((coords - coords[c(2:5, 1), ])^2))
  side_lengths <- sort(round(edges[1:2], 4))

  x_length <- side_lengths[1]  # shorter side (width)
  y_length <- side_lengths[2]  # longer side (length)
  z_length <- as.numeric(poly_$Height[1])  # height

  ratio_x <- x_length / z_length
  ratio_y <- y_length / z_length
  ratio_z <- z_length / max(x_length, y_length)

  return(list(ratio_x, ratio_y, ratio_z))
}

#' @noRd
directional_green_feature <- function(binary_green, p, direction, field_of_view) {
  direction_bearings <- c(
    north = 0,
    northeast = 45,
    east = 90,
    southeast = 135,
    south = 180,
    southwest = 225,
    west = 270,
    northwest = 315
  )

  center <- direction_bearings[[direction]]
  if (is.null(center)) {
    stop(sprintf("Unsupported BGVI direction: %s", direction), call. = FALSE)
  }
  directional_green_feature_bearing(binary_green, p, center, field_of_view)
}

#' @noRd
directional_green_feature_bearing <- function(binary_green, p, center, field_of_view) {
  center <- as.numeric(center[1]) %% 360
  field_of_view <- as.numeric(field_of_view[1])
  if (field_of_view >= 360) return(binary_green)

  cells <- seq_len(terra::ncell(binary_green))
  xy <- terra::xyFromCell(binary_green, cells)
  dx <- xy[, 1] - p[1]
  dy <- xy[, 2] - p[2]
  bearing <- (atan2(dx, dy) * 180 / pi + 360) %% 360
  angular_distance <- abs(((bearing - center + 180) %% 360) - 180)

  green_values <- terra::values(binary_green, mat = FALSE)
  directional_values <- ifelse(
    !is.na(green_values) & green_values == 1 & angular_distance <= field_of_view / 2,
    1,
    0
  )
  out <- binary_green
  terra::values(out) <- directional_values
  out
}

#' @noRd
bgvi_sector_mask <- function(template, p, orientation, field_of_view) {
  out <- template
  terra::values(out) <- 0
  if (inherits(orientation, "NULL") || length(orientation) == 0L || is.na(orientation[1])) {
    terra::values(out) <- 1
    return(out)
  }
  center <- bgvi_orientation_to_bearing(orientation)
  if (field_of_view >= 360) {
    terra::values(out) <- 1
    return(out)
  }
  cells <- seq_len(terra::ncell(template))
  xy <- terra::xyFromCell(template, cells)
  dx <- xy[, 1] - p[1]
  dy <- xy[, 2] - p[2]
  bearing <- (atan2(dx, dy) * 180 / pi + 360) %% 360
  angular_distance <- abs(((bearing - center + 180) %% 360) - 180)
  values <- ifelse(!is.na(terra::values(template, mat = FALSE)) &
                     angular_distance <= field_of_view / 2, 1, 0)
  terra::values(out) <- values
  out
}

#' @noRd
bgvi_orientation_to_bearing <- function(orientation) {
  direction_bearings <- c(
    north = 0,
    northeast = 45,
    east = 90,
    southeast = 135,
    south = 180,
    southwest = 225,
    west = 270,
    northwest = 315
  )
  if (is.character(orientation)) {
    orientation <- tolower(orientation[1])
    bearing <- direction_bearings[[orientation]]
    if (is.null(bearing)) {
      stop(sprintf(
        "`orientation` must be a number of degrees clockwise from north or one of: %s.",
        paste(names(direction_bearings), collapse = ", ")
      ), call. = FALSE)
    }
    return(bearing)
  }
  bearing <- as.numeric(orientation[1])
  if (is.na(bearing)) {
    stop("`orientation` must be a valid bearing or direction name.", call. = FALSE)
  }
  bearing %% 360
}

#' @noRd
resolve_bgvi_building_index <- function(x, building) {
  if (length(building) != 1) {
    stop("`building` must identify exactly one building.", call. = FALSE)
  }
  idx <- NA_integer_
  if ("id" %in% names(x)) {
    idx <- match(as.character(building), as.character(x$id))
  }
  if (is.na(idx)) {
    idx <- suppressWarnings(as.integer(building))
  }
  if (is.na(idx) || idx < 1L || idx > nrow(x)) {
    stop("`building` must be a valid row number or `id` value.", call. = FALSE)
  }
  idx
}

#' @noRd
resolve_bgvi_observer_height <- function(building,
                                         level = c("bottom", "top"),
                                         floor = NULL,
                                         height = NULL) {
  level <- match.arg(level)
  if (!is.null(height)) {
    height <- as.numeric(height[1])
    if (is.na(height) || height < 0) {
      stop("`height` must be a non-negative observer offset in metres.", call. = FALSE)
    }
    return(height)
  }
  if (!is.null(floor)) {
    floor <- as.integer(floor[1])
    if (is.na(floor) || floor < 1L) {
      stop("`floor` must be a positive integer.", call. = FALSE)
    }
    max_floor <- max(1L, round(as.numeric(building$Height[1]) / 3))
    if (floor > max_floor) {
      stop(sprintf("`floor` is above the estimated top floor (%s).", max_floor), call. = FALSE)
    }
    return(1.7 + (floor - 1L) * 3)
  }
  if (identical(level, "bottom")) {
    return(1.7)
  }
  max(1.7, 1.7 + as.numeric(building$Height[1]) - 3)
}

#' @noRd
viewshed_to_spatraster <- function(viewshed) {
  terra::rast(
    viewshed@visible,
    extent = terra::ext(viewshed@extent, xy = TRUE),
    crs = terra::crs(viewshed@crs)
  )
}

#' @noRd
align_bgvi_viewshed_raster <- function(viewshed_raster, template) {
  if (!terra::same.crs(viewshed_raster, template)) {
    viewshed_raster <- terra::project(viewshed_raster, template, method = "near")
  }
  terra::resample(viewshed_raster, template, method = "near")
}

#' @noRd
bgvi_raster_has_values <- function(r) {
  values <- terra::values(r, mat = FALSE)
  any(!is.na(values) & is.finite(values))
}

#' @noRd
bgvi_raster_bbox <- function(r) {
  e <- terra::ext(r)
  c(
    xmin = terra::xmin(e),
    xmax = terra::xmax(e),
    ymin = terra::ymin(e),
    ymax = terra::ymax(e)
  )
}

#' @noRd
bgvi_raster_geometry <- function(r) {
  bb <- bgvi_raster_bbox(r)
  sf::st_as_sfc(sf::st_bbox(bb, crs = terra::crs(r)))
}

#' @noRd
draw_bgvi_raster_layer <- function(r, col) {
  if (!bgvi_raster_has_values(r)) return(invisible(FALSE))
  terra::plot(
    r,
    col = col,
    legend = FALSE,
    axes = FALSE,
    add = TRUE
  )
  invisible(TRUE)
}

#' @noRd
bgvi_plot_raster <- function(result, buildings = NULL) {
  r <- result$viewshed_raster
  terra::values(r) <- 1

  sector <- if (!is.null(result$orientation)) {
    result$sector_mask == 1
  } else {
    r == 1
  }
  visible <- result$viewshed_raster > 0 & sector
  green <- result$binary_green == 1
  visible_green <- result$visible_green == 1
  radius_mask <- if (!is.null(result$radius) && !is.null(result$viewpoint)) {
    xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
    dist <- sqrt((xy[, 1] - result$viewpoint[1])^2 + (xy[, 2] - result$viewpoint[2])^2)
    mask <- r
    terra::values(mask) <- dist <= result$radius
    mask == 1
  } else {
    r == r
  }

  r <- terra::ifel(green, 3, r)
  r <- terra::ifel(visible, 2, r)
  r <- terra::ifel(visible_green, 4, r)

  if (!is.null(buildings) && inherits(buildings, "sf")) {
    building_mask <- terra::rasterize(terra::vect(buildings), r, field = 1, background = NA)
    r <- terra::ifel(!is.na(building_mask) & radius_mask, 5, r)
  }
  target_mask <- terra::rasterize(terra::vect(result$building), r, field = 1, background = NA)
  r <- terra::ifel(!is.na(target_mask) & radius_mask, 6, r)
  r
}

#' @noRd
plot_bgvi_viewshed_result <- function(result,
                                      buildings,
                                      scalebar = TRUE,
                                      scalebar_unit = c("auto", "km", "m"),
                                      scalebar_cex = 0.7,
                                      north_arrow = TRUE,
                                      ...) {
  scalebar_unit <- match.arg(scalebar_unit)
  plot_raster <- bgvi_plot_raster(result, buildings = buildings)
  bb <- bgvi_raster_bbox(result$viewshed_raster)
  map_geom <- bgvi_raster_geometry(result$viewshed_raster)
  cols <- c(
    "#F4F2EC",
    "#D2D6D3",
    "#B9D8A8",
    "#159A63",
    "#8F8F8F",
    "#305CDE"
  )
  old_par <- open_shared_map_layout(map_geom, legend = TRUE)
  on.exit(close_shared_map_layout(old_par, legend = TRUE), add = TRUE)

  terra::plot(
    plot_raster,
    col = cols,
    breaks = seq(0.5, 6.5, by = 1),
    legend = FALSE,
    axes = FALSE,
    main = "",
    ...
  )
  draw_bgvi_viewshed_caption(bb, result)
  graphics::points(result$viewpoint[1], result$viewpoint[2],
                   pch = 21, bg = "white", col = "black", cex = 1.1)

  map_scale <- if (isTRUE(scalebar)) noise_map_scale(map_geom) else NULL
  legend_info <- noise_map_draw_legend(
    labels = c("Base map", "Viewshed area", "Green base", "Visible green",
               "Building", "Target building"),
    fills = cols,
    title = "Layer",
    reserve_in = if (isTRUE(scalebar)) 0.5 else 0
  )
  draw_shared_map_ornaments(
    map_geom,
    scalebar = scalebar,
    scalebar_unit = scalebar_unit,
    scalebar_cex = scalebar_cex,
    north_arrow = north_arrow,
    legend_info = legend_info,
    map_scale = map_scale
  )
  invisible(plot_raster)
}

#' @noRd
format_bgvi_direction <- function(orientation) {
  if (is.null(orientation)) {
    return("all")
  }
  if (is.numeric(orientation)) {
    return(sprintf("%g deg", orientation %% 360))
  }
  as.character(orientation[1])
}

#' @noRd
draw_bgvi_viewshed_caption <- function(bb, result, cex = 0.72) {
  direction <- format_bgvi_direction(result$orientation)
  fov <- if (is.null(result$field_of_view)) 360 else result$field_of_view
  viewshed_area <- if (!is.null(result$viewshed_area)) result$viewshed_area else NA_real_
  visible_green_area <- if (!is.null(result$green_area)) result$green_area else NA_real_
  caption <- sprintf(
    "Direction: %s FOV: %g deg viewshed area: %.3f km2 visible green area: %.3f km2",
    direction,
    fov,
    viewshed_area / 1e6,
    visible_green_area / 1e6
  )

  width <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
  height <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  if (!is.finite(height) || height <= 0) {
    usr <- graphics::par("usr")
    height <- usr[4L] - usr[3L]
  }
  caption_cex <- cex
  while (caption_cex > 0.48 &&
         graphics::strwidth(caption, units = "user", cex = caption_cex) > width) {
    caption_cex <- caption_cex - 0.04
  }
  lines <- caption
  if (graphics::strwidth(caption, units = "user", cex = caption_cex) > width) {
    avg_char_width <- max(graphics::strwidth("M", units = "user", cex = caption_cex), 1e-6)
    wrap_width <- max(20L, floor(width / avg_char_width))
    lines <- strwrap(caption, width = wrap_width)
  }
  line_height <- 0.045 * height
  for (i in seq_along(lines)) {
    graphics::text(
      x = bb[["xmin"]],
      y = bb[["ymin"]] - 0.055 * height - (i - 1L) * line_height,
      labels = lines[i],
      adj = c(0, 0.5),
      font = 1,
      cex = caption_cex,
      xpd = NA
    )
  }
  invisible(NULL)
}

#' @noRd
gvi_from_viewshed <- function(viewshed, building, binary_green) {
  if (!requireNamespace("viewscape", quietly = TRUE)) {
    stop("Package 'viewscape' is required for this function. Install it with: install.packages('viewscape')", call. = FALSE)
  }
  v_area <- length(as.vector(viewshed@visible[viewshed@visible > 0])) *
    viewshed@resolution[1]^2
  green_proportion <- viewscape::calculate_feature(
    viewshed = viewshed,
    feature = binary_green,
    type = 2,
    exclude_value = 0
  )
  green_area <- v_area * green_proportion
  gvi <- green_area / max(v_area - building$g_area, 1e-6)
  list(gvi = min(gvi, 1), green_area = green_area, viewshed_area = v_area)
}

#' @noRd
get_gvi <- function(dsm, p, height, r, building, binary_chm,
                    directions = NULL, field_of_view = 45) {
  if (!requireNamespace("viewscape", quietly = TRUE)) {
    stop("Package 'viewscape' is required for this function. Install it with: install.packages('viewscape')", call. = FALSE)
  }
  tryCatch({
    # Viewshed
    v <- viewscape::compute_viewshed(dsm = dsm, viewpoints = p,
                                     offset_viewpoint = height,
                                     r = r
    )
    result <- gvi_from_viewshed(v, building, binary_chm)

    if (inherits(directions, "NULL")) {
      return(result)
    }

    directional_results <- lapply(directions, function(direction) {
      directional_feature <- directional_green_feature(
        binary_green = binary_chm,
        p = p,
        direction = direction,
        field_of_view = field_of_view
      )
      gvi_from_viewshed(v, building, directional_feature)
    })
    names(directional_results) <- directions

    return(c(list(overall = result), directional_results))
  }, error = function(e) {
    stop(sprintf("GVI calculation failed: %s", e$message))
  })
}
# Note:
# In fact, if the building itself occupies nearly the entire viewshed,
# or the canopy overlaps very closely, v_area - building$g_area could be <= 0,
# yet still surrounded by trees - in which case GVI = 1 is conceptually valid.

# max(v_area - building$g_area, 1e-6) prevents division by zero while allowing close overlap
# min(gvi, 1) caps the value at 1 if needed for interpretation as a proportion
# Still returns 0 if the viewshed fails or has no values

#### Data collection and processing ####
#' @noRd
download_to_file <- function(url, destfile, quiet = FALSE) {
  tryCatch({
    httr2::request(url) |>
      httr2::req_timeout(9999) |>
      httr2::req_perform(path = destfile)
    invisible(destfile)
  }, error = function(e) {
    stop(sprintf("Failed to download %s: %s", url, e$message), call. = FALSE)
  })
}

#' @noRd
get_GHSpop <- function(bbox = NULL, year = NULL, points = NULL, polygons = NULL, quiet = FALSE) {
  original_timeout <- getOption('timeout')
  on.exit(options(timeout = original_timeout), add = TRUE)
  options(timeout=9999)

  years <- c(2030, 2025, 2020, 2015, 2010, 2005, 2000, 1995, 1990, 1985, 1980, 1975)
  result_list <- list()
  temp_paths <- c()

  if (!year %in% years) {
    stop(sprintf("Input year %d is not in allowed range. Skipping.", year))
  }

  intersected_tiles <- ghsl_tiles[sf::st_intersects(ghsl_tiles, bbox, sparse = FALSE), ]
  if (nrow(intersected_tiles) == 0) {
    stop("No GHSL population tiles intersect the input extent")
  }

  on.exit(unlink(temp_paths, recursive = TRUE), add = TRUE)

  if (!inherits(points, "NULL") || !inherits(polygons, "NULL")) {
    target <- if (!inherits(points, "NULL")) points else polygons
    pop_total <- rep(NA_real_, nrow(target))
    cell_id <- rep(NA_character_, nrow(target))

    for (i in seq_len(nrow(intersected_tiles))) {
      temp_zip <- tempfile(fileext = ".zip")
      unzip_dir <- tempfile()
      temp_paths <- c(temp_paths, temp_zip, unzip_dir)

      url_ <- get_GHSurl(year, intersected_tiles$tile_id[i], 'pop')
      download_to_file(url_, temp_zip, quiet = quiet)
      utils::unzip(temp_zip, exdir = unzip_dir)
      tif_files <- list.files(unzip_dir, pattern = "\\.tif$", full.names = TRUE)
      if (length(tif_files) == 0) next

      rast_data <- terra::rast(tif_files[1])
      if (!inherits(points, "NULL")) {
        points_src <- sf::st_transform(points, terra::crs(rast_data))
        point_xy <- sf::st_coordinates(points_src)
        point_cells <- terra::cellFromXY(rast_data, point_xy)
        extracted <- suppressWarnings(
          terra::extract(rast_data, terra::vect(points_src), raw = TRUE, ID = FALSE)[, 1]
        )
        matched <- is.na(pop_total) & !is.na(extracted)
        pop_total[matched] <- extracted[matched]
        cell_id[matched] <- paste(intersected_tiles$tile_id[i], point_cells[matched], sep = ":")
      } else {
        polygons_src <- sf::st_transform(polygons, terra::crs(rast_data))
        extracted <- suppressWarnings(
          terra::extract(rast_data, terra::vect(polygons_src), weights = TRUE)
        )
        if (nrow(extracted) == 0) next
        value_col <- setdiff(names(extracted), c("ID", "weight"))[1]
        weighted_pop <- extracted[[value_col]] * extracted$weight
        tile_sum <- stats::aggregate(
          weighted_pop,
          by = list(ID = extracted$ID),
          FUN = function(x) sum(x, na.rm = TRUE)
        )
        has_values <- stats::aggregate(
          !is.na(extracted[[value_col]]),
          by = list(ID = extracted$ID),
          FUN = any
        )
        tile_sum$x[!has_values$x] <- NA_real_
        matched <- !is.na(tile_sum$x)
        pop_total[tile_sum$ID[matched]] <- rowSums(
          cbind(
            ifelse(is.na(pop_total[tile_sum$ID[matched]]), 0, pop_total[tile_sum$ID[matched]]),
            tile_sum$x[matched]
          )
        )
      }
    }

    if (all(is.na(pop_total))) {
      stop("No population values extracted from downloaded GHSL rasters")
    }
    if (!quiet) cli::cli_alert_success('Finished downloading population data')
    if (!inherits(points, "NULL")) {
      return(list(pop_total = pop_total, cell_id = cell_id))
    }
    return(list(pop_total = pop_total))
  }

  for (i in seq_len(nrow(intersected_tiles))) {
    temp_zip <- tempfile(fileext = ".zip")
    unzip_dir <- tempfile()
    temp_paths <- c(temp_paths, temp_zip, unzip_dir)

    url_ <- get_GHSurl(year, intersected_tiles$tile_id[i], 'pop')
    download_to_file(url_, temp_zip, quiet = quiet)
    utils::unzip(temp_zip, exdir = unzip_dir)
    tif_files <- list.files(unzip_dir, pattern = "\\.tif$", full.names = TRUE)
    if (length(tif_files) == 0) next

    rast_data <- terra::rast(tif_files[1])
    bbox_src <- sf::st_transform(bbox, terra::crs(rast_data))
    rast_data <- suppressWarnings(terra::crop(rast_data, terra::vect(bbox_src), snap = "out"))
    if (terra::ncell(rast_data) == 0) next
    result_list[[length(result_list) + 1]] <- rast_data
  }

  if (length(result_list) == 0) {
    stop("No population rasters downloaded")
  }
  if (!quiet) cli::cli_alert_success('Finished downloading population data')

  r <- if (length(result_list) == 1) result_list[[1]] else do.call(terra::merge, result_list)
  r <- r / (100 * 100)
  utm_crs <- get_utm_crs(bbox)
  terra::project(r, paste0('EPSG:', utm_crs), method = 'near')
}

#' @noRd
get_GHSres <- function(bbox = NULL, year = NULL, quiet = TRUE) {
  # Store the original 'timeout' option and ensure it's reset upon function exit
  original_timeout <- getOption('timeout')
  on.exit(options(timeout = original_timeout), add = TRUE)
  options(timeout=9999)

  # GHS population grid
  years <- c(2030, 2025, 2020, 2018, 2015, 2010, 2005, 2000, 1995, 1990, 1985, 1980, 1975)
  result_list_total <- list()
  result_list_nres <- list()
  temp_paths <- c()  # store paths for later cleanup

  if (year %in% years) {
    intersected_tiles <- ghsl_tiles[sf::st_intersects(ghsl_tiles, bbox, sparse = FALSE), ]
    for (i in seq_len(nrow(intersected_tiles))) {
      temp_total_zip <- tempfile(fileext = ".zip")
      temp_nres_zip <- tempfile(fileext = ".zip")
      urls <- get_GHSurl(year, intersected_tiles$tile_id[i], type = 'b_surf')
      download_to_file(urls[[1]], temp_total_zip, quiet = quiet)
      download_to_file(urls[[2]], temp_nres_zip, quiet = quiet)
      unzip_total_dir <- tempfile()
      unzip_nres_dir <- tempfile()
      utils::unzip(temp_total_zip, exdir = unzip_total_dir)
      utils::unzip(temp_nres_zip, exdir = unzip_nres_dir)
      total_tif_files <- list.files(unzip_total_dir, pattern = "\\.tif$", full.names = TRUE)
      nres_tif_files <- list.files(unzip_nres_dir, pattern = "\\.tif$", full.names = TRUE)
      if (length(total_tif_files) == 0) next
      total_rast_data <- terra::rast(total_tif_files[1])
      nres_rast_data <- terra::rast(nres_tif_files[1])
      result_list_total[[length(result_list_total) + 1]] <- total_rast_data
      result_list_nres[[length(result_list_nres) + 1]] <- nres_rast_data
      temp_paths <- c(temp_paths,
                      temp_total_zip, temp_nres_zip,
                      unzip_total_dir, unzip_nres_dir)
    }

    if (length(result_list_total) == 0) {
      stop("No building surface raster downloaded")
    }

    # Combine all into one terra raster object
    r_total <- if (length(result_list_total) == 1) result_list_total[[1]] else do.call(terra::merge, result_list_total)
    r_nres <- if (length(result_list_nres) == 1) result_list_nres[[1]] else do.call(terra::merge, result_list_nres)
    # reproject rasters
    utm_crs <- get_utm_crs(bbox)
    r_total <- terra::project(r_total, paste0('EPSG:', utm_crs), method = 'near')
    r_nres <- terra::project(r_nres, paste0('EPSG:', utm_crs), method = 'near')
    # calculate residential
    r_res <- r_total - r_nres
    # crop
    crop_ext <- terra::ext(terra::project(terra::vect(bbox), terra::crs(r_res)))
    r_total <- terra::crop(r_total, crop_ext)
    r_nres <- terra::crop(r_nres, crop_ext)
    r_res <- terra::crop(r_res, crop_ext)

    # Ensure cleanup
    on.exit(unlink(temp_paths, recursive = TRUE), add = TRUE)
    return(list(total = r_total, nres = r_nres, res = r_res))
  } else {
    stop(sprintf("Input year %d is not in allowed range. Skipping.", year))
  }
}

#' @noRd
get_GHSurl <- function(year, id, type) {
  if (type == 'pop') {
    # source: https://human-settlement.emergency.copernicus.eu/download.php?ds=pop
    return(
      paste0(
        'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2025_GLOBE_R2023A_54009_100/V1-0/tiles/GHS_POP_E',
        year,
        '_GLOBE_R2023A_54009_100_V1_0_',
        id,
        '.zip'
      )
    )
  } else if (type == 'b_surf') {
    return(
      list(
        paste0(
          'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_E',
          year,
          '_GLOBE_R2023A_54009_100/V1-0/tiles/GHS_BUILT_S_E',
          year,
          '_GLOBE_R2023A_54009_100_V1_0_',
          id,
          '.zip'
        ),
        paste0(
          'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_NRES_E',
          year,
          '_GLOBE_R2023A_54009_100/V1-0/tiles/GHS_BUILT_S_NRES_E',
          year,
          '_GLOBE_R2023A_54009_100_V1_0_',
          id,
          '.zip'
        )
      )
    )
  }
}



#' @noRd
get_dem <- function(bbox, key) {
  if (!requireNamespace("dsmSearch", quietly = TRUE)) {
    stop("Package 'dsmSearch' is required for this function. Install it with: install.packages('dsmSearch')", call. = FALSE)
  }
  if (missing(key)) stop("missing api key")
  bbox <- as_wgs84_bbox_vector(bbox)

  us_poly <- suppressMessages(tigris::nation(progress_bar = FALSE))
  bbox_poly <- sf::st_transform(
    sf::st_as_sfc(
      sf::st_bbox(
        c(xmin = unname(bbox[1]),
          ymin = unname(bbox[2]),
          xmax = unname(bbox[3]),
          ymax = unname(bbox[4])),
        crs = 4326
      )
    ), sf::st_crs(us_poly)
  )

  dem <- NULL

  if (any(sf::st_intersects(bbox_poly, us_poly, sparse = FALSE))) {
    # USGS sources
    dem <- tryCatch({
      dsmSearch::get_dsm_30(bbox = bbox, key = key, datatype = 'USGS1m')
    }, error = function(e1) {
      tryCatch({
        dsmSearch::get_dsm_30(bbox = bbox, key = key, datatype = 'USGS10m')
      }, error = function(e2) {
        dsmSearch::get_dsm_30(bbox = bbox, key = key, datatype = 'SRTMGL1')
      })
    })
  } else {
    # Global SRTM fallback
    dem <- tryCatch({
      dsmSearch::get_dsm_30(bbox = bbox, key = key, datatype = 'SRTMGL1')
    }, error = function(e) {
      NULL
    })
  }

  if (is.null(dem)) {
    stop("Failed to get elevation data within the area of interest.")
  }

  return(dem)
}


#' Retry a raster download after purging corrupt cached tiles
#'
#' @description
#' Tile downloads (Meta CHM, ETH CHM, DEM) are cached to `tempdir()` by the
#' upstream package. A truncated download stays cached for the whole session
#' and fails on every subsequent read, surfacing as GDAL *warnings*
#' ("TIFFReadEncodedStrip failed", "IReadBlock failed") followed by empty data
#' rather than a clean error - so a plain `tryCatch()` never sees it.
#'
#' This helper watches for those warnings, deletes the specific files named in
#' them, and re-evaluates the expression.
#'
#' @param expr Expression that downloads and returns a SpatRaster.
#' @param tries Integer. Total attempts. Default 3.
#' @param quiet Logical. Suppress retry messages.
#' @noRd
with_tile_retry <- function(expr, tries = 3L, quiet = TRUE) {
  ex <- substitute(expr)
  pf <- parent.frame()
  bad_pattern <- "TIFFReadEncodedStrip|TIFFFillStrip|IReadBlock|too few values"

  res <- NULL
  for (i in seq_len(tries)) {
    corrupt   <- FALSE
    bad_files <- character(0)

    res <- withCallingHandlers(
      tryCatch(
        eval(ex, pf),
        error = function(e) {
          if (grepl(bad_pattern, conditionMessage(e))) corrupt <<- TRUE
          NULL
        }
      ),
      warning = function(w) {
        m <- conditionMessage(w)
        if (grepl(bad_pattern, m)) {
          corrupt <<- TRUE
          hit <- regmatches(m, regexpr("/[^,[:space:]]+\\.tif", m))
          if (length(hit)) bad_files <<- c(bad_files, hit)
          invokeRestart("muffleWarning")
        }
      }
    )

    if (!corrupt && !is.null(res)) return(res)
    if (i >= tries) break

    # Purge the exact files named in the warnings; fall back to clearing
    # cached .tif files in tempdir() if none could be parsed out.
    bad_files <- unique(bad_files[file.exists(bad_files)])
    if (!length(bad_files))
      bad_files <- list.files(tempdir(), pattern = "\\.tif$",
                              full.names = TRUE, recursive = TRUE)
    if (length(bad_files)) unlink(bad_files, force = TRUE)

    if (!isTRUE(quiet))
      message("  Corrupt tile cache detected - purged ", length(bad_files),
              " file(s), retrying (attempt ", i + 1L, " of ", tries, ") ...")
    Sys.sleep(2)
  }

  if (is.null(res))
    stop("Tile download failed after ", tries, " attempts. The remote server ",
         "may be serving truncated tiles; try again later.", call. = FALSE)
  res
}

#' @noRd
get_chm <- function(bbox, min_height, datasource = "metachm", quiet = TRUE) {
  if (!requireNamespace("dsmSearch", quietly = TRUE)) {
    stop("Package 'dsmSearch' is required for this function. Install it with: install.packages('dsmSearch')", call. = FALSE)
  }
  bbox <- as_wgs84_bbox_vector(bbox)
  datasource <- tolower(as.character(datasource[1]))
  datasource <- match.arg(datasource, c("metachm", "ethchm"))
  datatype <- switch(
    datasource,
    metachm = "metaCHM",
    ethchm = "ethCHM"
  )
  # get CHM - retried with a cache purge if the tile comes back truncated
  chm <- with_tile_retry(
    suppressMessages(dsmSearch::get_dsm_30(bbox = bbox, datatype = datatype)),
    tries = 3L,
    quiet = quiet
  )
  # reporject chm
  bbox <- sf::st_as_sfc(
    sf::st_bbox(
      c(xmin = unname(bbox[1]),
        ymin = unname(bbox[2]),
        xmax = unname(bbox[3]),
        ymax = unname(bbox[4])),
      crs = 4326
    )
  )
  utm_crs <- get_utm_crs(bbox)
  chm <- terra::project(chm, paste0('EPSG:', utm_crs), method = 'near')
  # filtered CHM based on the minimum tree height
  filteredCHM <- terra::ifel(chm < min_height, 0, chm)
  binaryCHM <- terra::ifel(chm < min_height, 0, 1)
  return(list(filteredCHM, binaryCHM))
}

#' @noRd
get_greenspace <- function(bbox = NULL, buffer = NULL,
                           type = NULL, zoom = 17, year = NULL,
                           min_tree_height = 2) {
  if (inherits(type, "NULL")) {
    stop("Please input greenspace datasource.")
  }
  type <- match.arg(type, c("metachm", "esri", "sentinel2"))
  bbox_vector <- if (is.numeric(bbox) && length(bbox) == 4) {
    as_wgs84_bbox_vector(bbox)
  } else {
    as_wgs84_bbox_vector(bbox)
  }
  bbox_query <- unname(bbox_vector)

  if (type == "metachm") {
    g <- get_chm(bbox_vector, min_tree_height)[[2]]
  } else if (type == "esri") {
    g <- fetch_greenspace_tile(bbox = bbox_query, zoom = zoom,
                               provider = "esri")
 	utm_crs <- get_utm_crs(bbox)
    g <- terra::project(g$green, paste0('EPSG:', utm_crs), method = 'near')
  } else if (type == "sentinel2") {
    g <- fetch_greenspace_tile(bbox = bbox_query, zoom = zoom,
                               provider = "eox", year = year)
 	utm_crs <- get_utm_crs(bbox)
    g <- terra::project(g$green, paste0('EPSG:', utm_crs), method = 'near')
  }
  if (is.null(buffer)) {
    return(g)
  } else {
    return(terra::crop(g, terra::vect(buffer), mask = TRUE))
  }

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

#' Fetch ESA WorldCover land cover and reclassify to aerodynamic roughness
#' length (z0) values suitable for OpenFOAM nutURoughWallFunction.
#'
#' ESA WorldCover class codes and their z0 assignments (metres):
#'   10  Tree cover          -> 1.00  (set to NA if mask_tree_cover = TRUE,
#'                                     so porous-zone treatment is not doubled)
#'   20  Shrubland           -> 0.20
#'   30  Grassland           -> 0.05
#'   40  Cropland            -> 0.05
#'   50  Built-up            -> 0.50  (set to NA over building footprints so
#'                                     solid STL geometry is not doubled)
#'   60  Bare / sparse veg   -> 0.01
#'   70  Snow and ice        -> 0.001
#'   80  Permanent water     -> 0.0002
#'   90  Herbaceous wetland  -> 0.05
#'   95  Mangroves           -> 0.50
#'  100  Moss and lichen     -> 0.01
#'
#' @param bbox  Bounding box accepted by as_wgs84_bbox_vector().
#' @param source Character. Land-cover dataset: `"esa"` (ESA WorldCover,
#'   years 2020-2021) or `"esri"` (Sentinel-2 10 m ESRI LULC Time Series,
#'   years 2017-2025). Default `"esa"`.
#' @param year  Land-cover year. For `source = "esa"`: 2020 or 2021.
#'   For `source = "esri"`: 2017-2025. Default 2021.
#' @param crs   Optional target CRS for the output raster.
#' @param building_mask  Optional sf polygon layer of building footprints.
#'   Cells that fall inside buildings are set to NA (solid geometry handles
#'   them; roughness should not be applied twice).
#' @param mask_tree_cover Logical. If TRUE (default), cells classified as tree
#'   cover are set to NA so that porous-zone drag is not accumulated on top of
#'   a roughness penalty. For `"esa"`: class 10. For `"esri"`: class 2.
#'
#' @return A single-layer SpatRaster named "z0_roughness_m".
#' @noRd
get_roughness_raster <- function(bbox,
                                 source = c("esa", "esri"),
                                 year = 2021,
                                 crs = NULL,
                                 building_mask = NULL,
                                 mask_tree_cover = TRUE) {
  source <- match.arg(source)

  if (!requireNamespace("greenSD", quietly = TRUE))
    stop("Package 'greenSD' is required. Install with: install.packages('greenSD')",
         call. = FALSE)

  bbox_wgs84 <- as_wgs84_bbox_vector(bbox)

  # greenSD::get_esa_wc() fails on a *named* bbox vector (its internal lapply
  # over tiles errors with "!anyNA(x) is not TRUE"), so strip names here.
  bbox_gs <- unname(as.numeric(bbox_wgs84))

  # -- Fetch land-cover raster -----------------------------------------------
  if (source == "esa") {
    if (!year %in% c(2020L, 2021L))
      stop("ESA WorldCover is only available for years 2020 and 2021. ",
           "Use `landcover_source = 'esri'` for years 2017-2025.", call. = FALSE)
    lc <- greenSD::get_esa_wc(bbox_gs, datatype = "landcover", year = year)
  } else {
    # Sentinel-2 10 m ESRI LULC Time Series (2017-2025)
    if (!year %in% 2017:2025)
      stop("ESRI LULC is available for years 2017-2025.", call. = FALSE)
    lc <- greenSD::get_esa_wc(bbox_gs, datatype = "lulc", year = year)
  }

  if (is.null(lc))
    stop(
      "No land-cover tiles found for the requested area and year ",
      "(source = '", source, "', year = ", year, ").\n",
      if (source == "esri")
        paste0("  The most recent years of the Sentinel-2 LULC Time Series are ",
               "often not yet published, and coverage varies by region.\n",
               "  Try an earlier year, e.g. `landcover_year = 2023`.")
      else
        "  ESA WorldCover covers 2020 and 2021 only.",
      call. = FALSE
    )

  # greenSD may return a list or a multi-layer time-series stack; take layer 1
  if (is.list(lc) && !inherits(lc, "SpatRaster")) lc <- lc[[1L]]
  if (terra::nlyr(lc) > 1L) lc <- lc[[1L]]

  # Drop category and colour tables. Categorical rasters make `classify()`
  # and `lc == code` operate on labels rather than the raw class codes, and a
  # colour table carried through classify() produces a mis-rendered z0 raster.
  if (terra::is.factor(lc)) levels(lc) <- NULL
  try(terra::coltab(lc) <- NULL, silent = TRUE)
  lc <- terra::as.int(lc)

  # -- Reproject to UTM for metric accuracy ----------------------------------
  utm_crs <- get_utm_crs(bbox_wgs84)
  lc <- terra::project(lc, paste0("EPSG:", utm_crs), method = "near")

  # -- Reclassify: class code -> z0 (m) --------------------------------------
  if (source == "esa") {
    # ESA WorldCover class codes
    rcl <- matrix(
      c(
         10,  1.0000,   # Tree cover
         20,  0.2000,   # Shrubland
         30,  0.0500,   # Grassland
         40,  0.0500,   # Cropland
         50,  0.0300,   # Built-up -> paved ground between resolved buildings
         60,  0.0100,   # Bare / sparse vegetation
         70,  0.0010,   # Snow and ice
         80,  0.0002,   # Permanent water bodies
         90,  0.0500,   # Herbaceous wetland
         95,  0.5000,   # Mangroves (tall vegetation, not building-resolved)
        100,  0.0100    # Moss and lichen
      ),
      ncol = 2, byrow = TRUE
    )
    tree_class <- 10L
  } else {
    # ESRI LULC class codes (Sentinel-2 10 m Time Series)
    # Class 3 and 6 are not defined in the ESRI schema.
    rcl <- matrix(
      c(
         1,  0.0002,  # Water - rivers, ponds, lakes, oceans
         2,  1.0000,  # Trees - tall dense vegetation (>= ~4.5 m); masked as porous zone
         4,  0.0500,  # Flooded vegetation - rice paddies, flooded mangroves, emergent veg
         5,  0.0500,  # Crops - cereals, grasses, soy, fallow (not at tree height)
         7,  0.0300,  # Built area -> paved ground between resolved buildings
         8,  0.0100,  # Bare ground - rock, sand, desert, dry lake beds, mines
         9,  0.0010,  # Snow / ice - glaciers, permanent snowpack
        10,  0.0300,  # Clouds - no land-cover information; fallback value
        11,  0.1000   # Rangeland - parks, lawns, pastures, savannas, sparse shrub/grass
      ),
      ncol = 2, byrow = TRUE
    )
    tree_class <- 2L
  }

  z0 <- terra::classify(lc, rcl, others = 0.03)
  names(z0) <- "z0_roughness_m"

  # -- Mask tree-cover cells (handled as porous zones, not roughness) --------
  if (isTRUE(mask_tree_cover))
    z0 <- terra::ifel(lc == tree_class, NA, z0)

  # -- Mask building footprint cells (solid STL geometry) -------------------
  if (!is.null(building_mask)) {
    bm_proj <- sf::st_transform(building_mask, paste0("EPSG:", utm_crs))
    building_r <- terra::rasterize(
      terra::vect(bm_proj), z0, field = 1L, background = 0L
    )
    z0 <- terra::ifel(building_r == 1L, NA, z0)
  }

  # -- Reproject to caller-supplied CRS if requested -------------------------
  # `crs` may arrive as an sf crs object, a sp CRS, an EPSG integer, or a
  # string; terra::project() only accepts a character CRS or a SpatRaster.
  if (!is.null(crs)) {
    crs_chr <- if (inherits(crs, "crs")) {
      if (!is.na(crs$wkt)) crs$wkt else as.character(crs$input)
    } else if (is.numeric(crs)) {
      paste0("EPSG:", as.integer(crs))
    } else {
      as.character(crs)
    }

    # "near" (not bilinear): z0 is a lookup from discrete land-cover classes.
    # Interpolating across class boundaries invents roughness lengths that
    # correspond to no class, and smears NA tree/building cells into
    # neighbouring cells.
    if (!is.na(crs_chr) && nzchar(crs_chr))
      z0 <- terra::project(z0, crs_chr, method = "near")
  }

  z0
}

#' @noRd
filter_patch_area <- function(r, min_area, unit = "m2", directions = 8) {
  stopifnot(inherits(r, "SpatRaster"))
  unit <- match.arg(unit, c("m2", "ha", "km2"))

  # patch
  greens <- terra::ifel(r == 1, 1, NA)
  cl <- terra::patches(greens, directions = directions)

  # Cell area in m^2 (works for lon/lat and projected)
  cell_area_m2 <- terra::cellSize(cl, unit = "m")
  # Sum area per patch (zonal)
  z <- terra::zonal(cell_area_m2, cl, fun = "sum", na.rm = TRUE)  # columns: zone, sum
  if (is.null(z) || nrow(z) == 0) {
    out <- terra::ifel(is.na(cl), 0, 0)
    names(out) <- "greenspace_filtered"
    return(out)
  }

  # Map patch area back to each cell
  area_r <- terra::subst(cl, from = z[[1]], to = z[[2]])

  # Threshold (convert min_area to m^2)
  thr_m2 <- switch(unit,
                   m2  = min_area,
                   ha  = min_area * 1e4,
                   km2 = min_area * 1e6
                   )
  keep_mask <- !is.na(cl) & (area_r >= thr_m2)

  out <- terra::ifel(keep_mask, 1, 0)
  names(out) <- "greenspace_filtered"
  return(out)
}

#' @noRd
merge_elev <- function(building, dem, chm = NULL) {
  # Prioritize layers: canopy/building surface above terrain.
  out <- dem
  if (!is.null(building)) {
    building <- terra::resample(building, out, method = "bilinear")
    out <- terra::ifel(!is.na(building) & building != 0, building, out)
  }
  if (!is.null(chm)) {
    chm <- terra::resample(chm, out, method = "bilinear")
    out <- terra::ifel(!is.na(chm) & chm != 0, chm, out)
  }
  out
}

#' @importFrom sf st_buffer st_centroid
#' @noRd
get_buffer <- function(x = NULL, radius = NULL) {
  bbox <- get_bbox(x)
  utm_crs <- get_utm_crs(bbox)
  x <- sf::st_transform(x, utm_crs)
  ct <- sf::st_centroid(x)
  buffer_ <- sf::st_buffer(ct, dist = radius) # utm
  bbox <- get_bbox(buffer_) # WGS 84
  return(list(buffer=buffer_, bbox=bbox, centroid=ct))
}

#### Projection tools ####
#' @noRd
get_utm_crs <- function(bbox) {
  if (is.numeric(bbox) && length(bbox) == 4) {
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
  centroid <- sf::st_centroid(sf::st_union(bbox))
  coords <- sf::st_coordinates(centroid)
  lon <- coords[1]
  lat <- coords[2]
  zone <- floor((lon + 180) / 6) + 1
  epsg <- if (lat >= 0) {
    32600 + zone
  } else {
    32700 + zone
  }
  return(epsg)
}

#' @noRd
reproj <- function(bbox, poly, utm_crs, res) {
  bbox_proj <- sf::st_transform(bbox, crs = utm_crs)
  poly_proj <- sf::st_transform(poly, crs = utm_crs)
  # Force bbox to be axis-aligned in projected CRS
  bbox_aligned <- sf::st_as_sfc(sf::st_bbox(bbox_proj), crs = utm_crs)

  # bbox_raster <- terra::ext(sf::st_bbox(bbox_proj))

  bbox_raster <- terra::ext(sf::st_bbox(bbox_aligned))
  # Snap the extent to nearest multiple of resolution
  bbox_raster <- terra::align(bbox_raster, res)
  return(list(bbox_aligned, poly_proj, bbox_raster))
}

#' @noRd
get_bbox <- function(x) {
  bbox <- sf::st_as_sfc(sf::st_bbox(x), crs = sf::st_crs(x))
  bbox <- sf::st_transform(bbox, crs = 4326)
  return(bbox)
}

#' @noRd
bbox_poly_to_list <- function(bbox) {
  coor <- sf::st_coordinates(bbox)
  return(
    c(min(coor[,1]), min(coor[,2]), max(coor[,1]), max(coor[,2]))
  )
}

#' @noRd
as_wgs84_bbox_vector <- function(bbox) {
  if (inherits(bbox, c("sf", "sfc"))) {
    if (is.na(sf::st_crs(bbox))) {
      stop("`bbox` must have a coordinate reference system.", call. = FALSE)
    }
    bbox <- get_bbox(bbox)
    bbox <- bbox_poly_to_list(bbox)
  } else if (inherits(bbox, "bbox")) {
    bbox <- sf::st_as_sfc(bbox)
    if (is.na(sf::st_crs(bbox))) {
      stop("`bbox` must have a coordinate reference system.", call. = FALSE)
    }
    bbox <- sf::st_transform(bbox, crs = 4326)
    bbox <- bbox_poly_to_list(bbox)
  } else {
    bbox <- as.numeric(bbox)
  }

  if (length(bbox) != 4 || anyNA(bbox)) {
    stop("`bbox` must be a four-value bounding box.", call. = FALSE)
  }
  bbox <- c(
    xmin = unname(bbox[1]),
    ymin = unname(bbox[2]),
    xmax = unname(bbox[3]),
    ymax = unname(bbox[4])
  )
  if (bbox["xmin"] > bbox["xmax"] || bbox["ymin"] > bbox["ymax"] ||
      bbox["xmin"] < -180 || bbox["xmax"] > 180 ||
      bbox["ymin"] < -90 || bbox["ymax"] > 90) {
    stop(
      "`bbox` must be in EPSG:4326 longitude/latitude coordinates before calling dsmSearch.",
      call. = FALSE
    )
  }
  bbox
}

#' @noRd
unify_layers <- function(bbox, ...) {
  utm_crs <- get_utm_crs(bbox)
  input_layers <- list(...)

  # Reproject to UTM
  input_layers <- lapply(input_layers, function(x) terra::project(x, paste0("epsg:", utm_crs)))

  # Use the first layer as reference
  ref <- input_layers[[1]]

  # Align all layers to match reference (extent, resolution, projection)
  aligned_layers <- lapply(input_layers, function(x) {
    terra::resample(x, ref, method = "near")
  })

  return(aligned_layers)
}

#' @noRd
normalize_bgvi_source <- function(value, choices) {
  if (inherits(value, "NULL")) return(NULL)
  value <- tolower(as.character(value[1]))
  if (value %in% c("none", "null", "na")) return(NULL)
  match.arg(value, choices)
}

#' @noRd
validate_bgvi_directions <- function(directions) {
  valid_directions <- c(
    "southwest", "southeast", "northeast", "northwest",
    "north", "east", "west", "south"
  )
  if (inherits(directions, "NULL")) return(NULL)
  directions <- unique(tolower(as.character(directions)))
  invalid_directions <- setdiff(directions, valid_directions)
  if (length(invalid_directions) > 0) {
    stop(
      sprintf(
        "`directions` contains invalid value(s): %s. Valid values are: %s.",
        paste(invalid_directions, collapse = ", "),
        paste(valid_directions, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  directions
}

#' @noRd
validate_bgvi_field_of_view <- function(field_of_view) {
  field_of_view <- as.numeric(field_of_view[1])
  if (is.na(field_of_view) || field_of_view <= 0 || field_of_view > 360) {
    stop("`field_of_view` must be a numeric value greater than 0 and less than or equal to 360.", call. = FALSE)
  }
  field_of_view
}

#' @noRd
prepare_bgvi_inputs <- function(x,
                                datasource_canopy_height = "metachm",
                                datasource_greenspace = NULL,
                                min_tree_height = 2,
                                zoom = 17,
                                radius = 800,
                                year = NULL,
                                resolution = NULL,
                                key = NULL,
                                quiet = FALSE) {
  datasource_canopy_height <- normalize_bgvi_source(datasource_canopy_height, c("metachm", "ethchm"))
  datasource_greenspace <- normalize_bgvi_source(datasource_greenspace, c("esri", "sentinel2"))
  if (inherits(datasource_canopy_height, "NULL") && inherits(datasource_greenspace, "NULL")) {
    stop("At least one of `datasource_canopy_height` or `datasource_greenspace` must be supplied.", call. = FALSE)
  }

  input_bbox <- get_bbox(x)
  utm_crs <- get_utm_crs(input_bbox)
  projected_poly <- sf::st_transform(x, utm_crs)
  analysis_poly <- prepare_group_analysis_buildings(projected_poly)
  analysis_bbox <- get_bbox(sf::st_buffer(projected_poly, dist = radius))
  bbox_vector <- bbox_poly_to_list(analysis_bbox)

  if (!quiet) cli::cli_alert_info('Start downloading BGVI raster inputs ...')
  dem <- get_dem(bbox_vector, key)
  chm_layers <- NULL
  if (!inherits(datasource_canopy_height, "NULL")) {
    chm_layers <- suppressMessages(get_chm(
      bbox_vector,
      min_tree_height,
      datasource = datasource_canopy_height
    ))
  }
  greenspace <- NULL
  if (!inherits(datasource_greenspace, "NULL")) {
    greenspace <- get_greenspace(
      bbox = analysis_bbox,
      buffer = NULL,
      type = datasource_greenspace,
      zoom = zoom,
      year = year,
      min_tree_height = min_tree_height
    )
  }

  analysis_poly$g_area <- as.numeric(sf::st_area(analysis_poly))
  analysis_utm_crs <- get_utm_crs(analysis_bbox)
  dem_projected <- terra::project(dem, paste0("EPSG:", analysis_utm_crs), method = "bilinear")
  dem_res <- terra::res(dem_projected)[1]

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

  bh_all <- rasterize_height(projected_poly, analysis_bbox, raster_res)
  dem <- terra::resample(dem_projected, bh_all, method = "bilinear")
  chm <- if (!inherits(chm_layers, "NULL")) {
    terra::resample(chm_layers[[1]], bh_all, method = "bilinear")
  } else {
    terra::ifel(is.na(dem), NA, 0)
  }
  canopy_green_aligned <- if (!inherits(chm_layers, "NULL")) {
    terra::resample(chm_layers[[2]], bh_all, method = "near")
  } else {
    NULL
  }
  tile_green_aligned <- if (!inherits(greenspace, "NULL")) {
    terra::resample(greenspace, bh_all, method = "near")
  } else {
    NULL
  }

  binary_green <- terra::ifel(is.na(dem), NA, 0)
  if (!is.null(canopy_green_aligned)) {
    binary_green <- terra::ifel(canopy_green_aligned == 1, 1, binary_green)
  }
  if (!is.null(tile_green_aligned)) {
    binary_green <- terra::ifel(tile_green_aligned == 1, 1, binary_green)
  }

  flat_id_poly <- sf::st_transform(analysis_poly, analysis_utm_crs)
  if (!"id" %in% names(flat_id_poly)) {
    flat_id_poly$id <- seq_len(nrow(flat_id_poly))
  }
  bldg_id_rast <- rasterize_height(flat_id_poly, analysis_bbox, raster_res, height_field = "id")

  flat_centroids <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(flat_id_poly))))
  flat_base_elev <- terra::extract(dem, flat_centroids[, 1:2, drop = FALSE], method = "bilinear")[, 1]
  dem_median <- stats::median(terra::values(dem, mat = FALSE), na.rm = TRUE)
  flat_base_elev[!is.finite(flat_base_elev)] <- dem_median

  base_elev_rast <- terra::subst(bldg_id_rast, from = flat_id_poly$id, to = flat_base_elev)
  binary_green <- terra::ifel(binary_green == 1, 1, 0)

  list(
    projected_poly = projected_poly,
    analysis_poly = analysis_poly,
    analysis_bbox = analysis_bbox,
    dem = dem,
    chm = chm,
    binary_green = binary_green,
    bh_all = bh_all,
    base_elev = base_elev_rast,
    bldg_id = bldg_id_rast,
    radius = radius,
    raster_res = raster_res
  )
}

#' @noRd
compute_gvi_per_building <- function(building,
                                     dem_path,
                                     chm_path,
                                     binary_green_path,
                                     bh_all_path,
                                     base_elev_path = NULL,
                                     bldg_id_path = NULL,
                                     radius,
                                     floor,
                                     floor_step,
                                     short_building_threshold,
                                     directions = NULL,
                                     field_of_view = 45) {
  # Read raster files inside each worker.
  dem <- terra::rast(dem_path)
  chm <- terra::rast(chm_path)
  binary_green <- terra::rast(binary_green_path)
  bh_all <- terra::rast(bh_all_path)
  base_elev <- if (!is.null(base_elev_path)) terra::rast(base_elev_path) else NULL
  bldg_id   <- if (!is.null(bldg_id_path)) terra::rast(bldg_id_path) else NULL

  scene <- prepare_bgvi_building_scene(
    building = building,
    dem = dem,
    chm = chm,
    binary_green = binary_green,
    bh_all = bh_all,
    base_elev = base_elev,
    bldg_id = bldg_id,
    radius = radius
  )
  p <- scene$viewpoint
  dsm_ <- scene$dsm
  binary_green <- scene$binary_green

  get_gvi_values <- function(height) {
    if (inherits(directions, "NULL")) {
      return(get_gvi(dsm_, p, height, radius, building, binary_green))
    }
    get_gvi(
      dsm_,
      p,
      height,
      radius,
      building,
      binary_green,
      directions = directions,
      field_of_view = field_of_view
    )
  }
  overall_gvi <- function(values) {
    v <- if ("overall" %in% names(values)) values$overall else values
    as.numeric(v$gvi)
  }
  overall_green_area <- function(values) {
    v <- if ("overall" %in% names(values)) values$overall else values
    as.numeric(v$green_area)
  }

  height_bottom <- 1.7
  bottom_values <- get_gvi_values(height_bottom)
  bottom_gvi <- overall_gvi(bottom_values)
  bottom_green_area <- overall_green_area(bottom_values)
  top_values <- if (building$Height < short_building_threshold) {
    bottom_values
  } else {
    height_top <- 1.7 + building$Height - 3
    get_gvi_values(height_top)
  }
  top_gvi <- overall_gvi(top_values)
  top_green_area <- overall_green_area(top_values)
  direction_results <- list()

  if (isTRUE(floor)) {
    floor_ids <- unique(c(seq.int(1L, building$estimated_floors, by = floor_step),
                          building$estimated_floors))
    GVIs <- numeric(length(floor_ids))
    GreenAreas <- numeric(length(floor_ids))
    directional_GVIs <- if (!inherits(directions, "NULL")) {
      stats::setNames(vector("list", length(directions)), directions)
    } else {
      NULL
    }
    for (j in seq_along(floor_ids)) {
      height <- 1.7 + (floor_ids[j] - 1) * 3
      floor_values <- get_gvi_values(height)
      GVIs[j] <- overall_gvi(floor_values)
      GreenAreas[j] <- overall_green_area(floor_values)
      if (!inherits(directions, "NULL")) {
        for (direction in directions) {
          directional_GVIs[[direction]][j] <- as.numeric(floor_values[[direction]]$gvi)
        }
      }
    }
    if (!inherits(directions, "NULL")) {
      for (direction in directions) {
        direction_results[[direction]] <- list(
          mean = mean(directional_GVIs[[direction]]),
          bottom = as.numeric(bottom_values[[direction]]$gvi),
          top = as.numeric(top_values[[direction]]$gvi),
          min = min(directional_GVIs[[direction]]),
          max = max(directional_GVIs[[direction]]),
          sd = if (length(directional_GVIs[[direction]]) > 1) {
            stats::sd(directional_GVIs[[direction]])
          } else {
            NA_real_
          }
        )
      }
    }
    return(list(
      mean = mean(GVIs),
      bottom = bottom_gvi,
      top = top_gvi,
      min = min(GVIs),
      max = max(GVIs),
      sd = if (length(GVIs) > 1) stats::sd(GVIs) else NA_real_,
      bottom_green_area = bottom_green_area,
      top_green_area = top_green_area,
      mean_green_area = mean(GreenAreas),
      directions = direction_results
    ))
  } else {
    if (building$Height < short_building_threshold) {
      mean_gvi <- bottom_gvi
      mean_green_area <- bottom_green_area
    } else {
      mean_gvi <- mean(c(top_gvi, bottom_gvi))
      mean_green_area <- mean(c(top_green_area, bottom_green_area))
    }
    if (!inherits(directions, "NULL")) {
      for (direction in directions) {
        direction_bottom <- as.numeric(bottom_values[[direction]]$gvi)
        direction_top <- as.numeric(top_values[[direction]]$gvi)
        direction_mean <- if (building$Height < short_building_threshold) {
          direction_bottom
        } else {
          mean(c(direction_top, direction_bottom))
        }
        direction_values <- if (building$Height < short_building_threshold) {
          direction_bottom
        } else {
          c(direction_bottom, direction_top)
        }
        direction_results[[direction]] <- list(
          mean = direction_mean,
          bottom = direction_bottom,
          top = direction_top,
          min = min(direction_values),
          max = max(direction_values),
          sd = if (length(direction_values) > 1) stats::sd(direction_values) else NA_real_
        )
      }
    }
    return(list(
      mean = mean_gvi,
      bottom = bottom_gvi,
      top = top_gvi,
      min = NA_real_,
      max = NA_real_,
      sd = NA_real_,
      bottom_green_area = bottom_green_area,
      top_green_area = top_green_area,
      mean_green_area = mean_green_area,
      directions = direction_results
    ))
  }
}

#' @noRd
prepare_bgvi_building_scene <- function(building,
                                        dem,
                                        chm,
                                        binary_green,
                                        bh_all,
                                        base_elev = NULL,
                                        bldg_id = NULL,
                                        radius) {
  # Compute centroid
  centroid <- suppressWarnings(sf::st_centroid(building$geometry))
  p <- as.vector(sf::st_coordinates(centroid))

  # Crop all rasters to the local viewshed radius.
  buffer <- sf::st_buffer(centroid, dist = radius)
  dem <- terra::crop(dem, terra::vect(buffer), mask = TRUE)
  chm <- terra::crop(chm, terra::vect(buffer), mask = TRUE)
  binary_green <- terra::crop(binary_green, terra::vect(buffer), mask = TRUE)
  bh_all <- terra::crop(bh_all, terra::vect(buffer), mask = TRUE)
  if (!is.null(base_elev)) base_elev <- terra::crop(base_elev, terra::vect(buffer), mask = TRUE)
  if (!is.null(bldg_id))   bldg_id   <- terra::crop(bldg_id, terra::vect(buffer), mask = TRUE)

  # Flatten the target building footprint: it is the observer, not an obstacle.
  target_mask <- terra::rasterize(terra::vect(building), bh_all, field = 1, background = 0)
  bh_without_target <- terra::ifel(target_mask == 1, 0, bh_all)
  chm_without_target <- terra::ifel(target_mask == 1, 0, chm)
  binary_green <- terra::ifel(target_mask == 1, 0, binary_green)

  surface <- terra::ifel(chm_without_target > bh_without_target,
                         chm_without_target,
                         bh_without_target)
  ground_dsm <- dem + surface

  # Flatten neighboring buildings' roofs: each building keeps a single flat
  # roof elevation (its own base ground elevation + height) instead of
  # following the terrain slope pixel by pixel, matching get_fused_dsm().
  # The target building's own footprint was already zeroed out above (its
  # `bh_without_target` is 0 there), so it correctly reduces to ground level
  # rather than getting a flat roof of its own.
  if (!is.null(base_elev) && !is.null(bldg_id)) {
    building_mask <- bldg_id > 0
    flat_roof <- base_elev + bh_without_target
    dsm_ <- terra::ifel(
      building_mask,
      terra::ifel(flat_roof > ground_dsm, flat_roof, ground_dsm),
      ground_dsm
    )
  } else {
    dsm_ <- ground_dsm
  }
  list(
    dsm = dsm_,
    binary_green = binary_green,
    target_mask = target_mask,
    viewpoint = p,
    centroid = centroid
  )
}

#### utils ####
time_taken <- function(process_time) {
  if (process_time >= 60) {
    cli::cli_alert_success(paste0("Completed. Time taken: ", base::round(process_time/60), " minutes."))
  } else {
    cli::cli_alert_success(paste0("Completed. Time taken: ", base::round(process_time), " seconds."))
  }
}

#### Shadow and radiation helpers ####
prepare_shadow_buildings <- function(x, height_field) {
  if (!inherits(x, "sf")) {
    stop("`x` must be an sf polygon object.", call. = FALSE)
  }
  if (!height_field %in% names(x)) {
    stop("Missing height field `", height_field, "` in `x`.", call. = FALSE)
  }
  if (is.na(sf::st_crs(x))) {
    stop("`x` must have a valid CRS.", call. = FALSE)
  }
  geom_type <- unique(as.character(sf::st_geometry_type(x)))
  if (!all(geom_type %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("`x` must contain polygon geometries.", call. = FALSE)
  }
  out <- x
  if (isTRUE(sf::st_is_longlat(out))) {
    out <- sf::st_transform(out, get_utm_crs(out))
  }
  out[[height_field]] <- as.numeric(out[[height_field]])
  if (anyNA(out[[height_field]])) {
    stop("`", height_field, "` must contain numeric building heights.", call. = FALSE)
  }
  out
}

resolve_solar_inputs <- function(x,
                                 azimuth = NULL,
                                 elevation = NULL,
                                 solar_time = NULL,
                                 time_zone = NULL) {
  if (!is.null(solar_time) || !is.null(time_zone)) {
    if (is.null(solar_time) || is.null(time_zone)) {
      stop("Both `solar_time` and `time_zone` must be supplied together.", call. = FALSE)
    }
    return(solar_position_from_time(x, solar_time, time_zone))
  }
  validate_azimuth_elevation(azimuth, elevation)
}

validate_azimuth_elevation <- function(azimuth, elevation) {
  if (is.null(azimuth) || is.null(elevation)) {
    stop("Provide either `azimuth` and `elevation`, or `solar_time` and `time_zone`.", call. = FALSE)
  }
  azimuth <- as.numeric(unlist(azimuth, use.names = FALSE))
  elevation <- as.numeric(unlist(elevation, use.names = FALSE))
  if (length(azimuth) != length(elevation)) {
    stop("`azimuth` and `elevation` must have the same length.", call. = FALSE)
  }
  if (length(azimuth) == 0 || anyNA(azimuth) || anyNA(elevation)) {
    stop("`azimuth` and `elevation` must contain numeric values.", call. = FALSE)
  }
  cbind(azimuth = azimuth %% 360, elevation = elevation)
}

normalize_solar_time <- function(solar_time, time_zone) {
  if (!is.character(time_zone) || length(time_zone) != 1 || is.na(time_zone)) {
    stop("`time_zone` must be a single character time zone.", call. = FALSE)
  }
  if (is.list(solar_time)) {
    solar_time <- unlist(solar_time, use.names = FALSE)
  }
  if (!is.character(solar_time) || length(solar_time) == 0 || anyNA(solar_time)) {
    stop(
      "`solar_time` must be a character string or character vector, such as ",
      "\"2026-06-21 15:00:00\".",
      call. = FALSE
    )
  }
  as.POSIXct(solar_time, tz = time_zone)
}

make_shadow_sun_ids <- function(solar_time, n) {
  if (is.null(solar_time)) {
    return(paste0("solar_", seq_len(n)))
  }
  if (is.list(solar_time)) {
    solar_time <- unlist(solar_time, use.names = FALSE)
  }
  solar_time <- as.character(solar_time)
  if (length(solar_time) == n) {
    return(solar_time)
  }
  paste0("solar_", seq_len(n))
}

solar_position_from_time <- function(x, solar_time, time_zone) {
  solar_time <- normalize_solar_time(solar_time, time_zone)
  if (anyNA(solar_time)) {
    stop("`solar_time` could not be converted to POSIXct with `time_zone`.", call. = FALSE)
  }
  centroid <- sf::st_transform(sf::st_centroid(sf::st_union(sf::st_geometry(x))), 4326)
  xy <- sf::st_coordinates(centroid)[1, ]
  calc_solar_position(lon = xy[["X"]], lat = xy[["Y"]], time = solar_time)
}

calc_solar_position <- function(lon, lat, time) {
  utc <- as.POSIXct(format(time, tz = "UTC", usetz = TRUE), tz = "UTC")
  lt <- as.POSIXlt(utc, tz = "UTC")
  day <- lt$yday + 1
  hour <- lt$hour + lt$min / 60 + lt$sec / 3600
  gamma <- 2 * pi / 365 * (day - 1 + (hour - 12) / 24)
  eqtime <- 229.18 * (
    0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -
      0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma)
  )
  decl <- 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) -
    0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) -
    0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)
  time_offset <- eqtime + 4 * lon
  tst <- (hour * 60 + time_offset) %% 1440
  ha <- deg2rad(tst / 4 - 180)
  lat_rad <- deg2rad(lat)
  cos_zenith <- sin(lat_rad) * sin(decl) + cos(lat_rad) * cos(decl) * cos(ha)
  cos_zenith <- pmin(pmax(cos_zenith, -1), 1)
  zenith <- acos(cos_zenith)
  elevation <- 90 - rad2deg(zenith)
  azimuth <- rad2deg(atan2(
    sin(ha),
    cos(ha) * sin(lat_rad) - tan(decl) * cos(lat_rad)
  )) + 180
  cbind(azimuth = azimuth %% 360, elevation = elevation)
}

compute_shadow_footprints <- function(buildings, solar_pos, height_field, b) {
  if (solar_pos[[2]] <= 0) {
    stop("Solar elevation must be above the horizon for shadow footprints.", call. = FALSE)
  }
  geoms <- lapply(seq_len(nrow(buildings)), function(i) {
    cast_shadow_geometry(
      sf::st_geometry(buildings[i, ]),
      buildings[[height_field]][i],
      solar_pos,
      b
    )
  })
  out <- buildings
  sf::st_geometry(out) <- do.call(c, geoms)
  out
}

compute_canopy_shadow_footprints <- function(canopy, buildings, solar_pos, b) {
  if (is.null(canopy) || nrow(canopy$xy) == 0) {
    return(empty_shadow_sf(sf::st_crs(buildings)))
  }
  if (solar_pos[[2]] <= 0) {
    stop("Solar elevation must be above the horizon for canopy shadow footprints.", call. = FALSE)
  }
  shift <- canopy_shadow_shift(canopy$height, solar_pos)
  half_cell <- canopy$cell_size / 2
  geoms <- lapply(seq_len(nrow(canopy$xy)), function(i) {
    x <- canopy$xy[i, 1]
    y <- canopy$xy[i, 2]
    source <- sf::st_polygon(list(rbind(
      c(x - half_cell, y - half_cell),
      c(x - half_cell, y + half_cell),
      c(x + half_cell, y + half_cell),
      c(x + half_cell, y - half_cell),
      c(x - half_cell, y - half_cell)
    )))
    shifted <- source + c(shift$dx[i], shift$dy[i])
    side_geoms <- shadow_side_polygons(
      sf::st_sfc(source, crs = sf::st_crs(buildings)),
      c(dx = shift$dx[i], dy = shift$dy[i])
    )
    geom <- sf::st_union(c(
      sf::st_sfc(source, crs = sf::st_crs(buildings)),
      sf::st_sfc(shifted, crs = sf::st_crs(buildings)),
      side_geoms
    ))
    if (is.finite(b) && b > 0) {
      geom <- sf::st_buffer(sf::st_buffer(geom, b), -b)
    }
    sf::st_make_valid(geom)
  })
  sf::st_sf(
    shadow_source = rep("canopy", length(geoms)),
    canopy_height = canopy$height,
    ground = canopy$ground,
    geometry = do.call(c, geoms)
  )
}

empty_shadow_sf <- function(crs) {
  sf::st_sf(
    shadow_source = character(),
    geometry = sf::st_sfc(crs = crs)
  )
}

combine_shadow_sources <- function(building_shadows, canopy_shadows) {
  crs <- sf::st_crs(building_shadows)
  building_geom <- sf::st_geometry(building_shadows)
  canopy_geom <- sf::st_geometry(canopy_shadows)
  sf::st_sf(
    shadow_source = c(building_shadows$shadow_source, canopy_shadows$shadow_source),
    solar_index = c(building_shadows$solar_index, canopy_shadows$solar_index),
    sun_id = c(building_shadows$sun_id, canopy_shadows$sun_id),
    azimuth = c(building_shadows$azimuth, canopy_shadows$azimuth),
    elevation = c(building_shadows$elevation, canopy_shadows$elevation),
    canopy_height = c(rep(NA_real_, nrow(building_shadows)), canopy_shadows$canopy_height),
    ground = c(rep(NA_real_, nrow(building_shadows)), canopy_shadows$ground),
    geometry = sf::st_sfc(c(as.list(building_geom), as.list(canopy_geom)), crs = crs)
  )
}

combine_shadow_footprint_overlaps <- function(shadows) {
  if (nrow(shadows) == 0) {
    return(shadows)
  }
  source_values <- unique(shadows$shadow_source)
  out <- lapply(source_values, function(source) {
    subset <- shadows[shadows$shadow_source == source, ]
    geom <- sf::st_make_valid(sf::st_union(sf::st_geometry(subset)))
    sf::st_sf(
      shadow_source = source,
      solar_count = length(unique(subset$solar_index)),
      solar_indices = paste(sort(unique(subset$solar_index)), collapse = ","),
      sun_ids = paste(sort(unique(subset$sun_id)), collapse = ","),
      azimuth_min = min(subset$azimuth, na.rm = TRUE),
      azimuth_max = max(subset$azimuth, na.rm = TRUE),
      elevation_min = min(subset$elevation, na.rm = TRUE),
      elevation_max = max(subset$elevation, na.rm = TRUE),
      geometry = sf::st_sfc(geom, crs = sf::st_crs(shadows))
    )
  })
  do.call(rbind, out)
}

open_shared_map_layout <- function(geometry,
                                   legend = TRUE,
                                   legend_width = 0.26,
                                   mar = c(2.0, 0.2, 0.2, 0.2)) {
  noise_map_open_layout(
    legend = legend,
    legend_width = legend_width,
    mar = mar,
    asp_ratio = noise_bbox_aspect(geometry)
  )
}

close_shared_map_layout <- function(old_par, legend = TRUE) {
  noise_map_close_layout(old_par, legend = legend)
}

draw_shared_map_frame <- function(geometry, main = NULL) {
  bb <- sf::st_bbox(geometry)
  graphics::plot(
    NA,
    type = "n",
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "",
    xlim = c(bb[["xmin"]], bb[["xmax"]]),
    ylim = c(bb[["ymin"]], bb[["ymax"]])
  )
  if (!is.null(main) && nzchar(main)) {
    draw_shared_map_title_below(bb, main)
  }
  invisible(bb)
}

draw_shared_map_title_below <- function(bb, main, cex = 1.1) {
  if (is.null(main) || !nzchar(main)) return(invisible(NULL))
  h <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  if (!is.finite(h) || h <= 0) h <- graphics::par("usr")[4L] - graphics::par("usr")[3L]
  graphics::text(
    x = mean(c(bb[["xmin"]], bb[["xmax"]])),
    y = bb[["ymin"]] - 0.08 * h,
    labels = main,
    font = 2,
    cex = cex,
    xpd = NA
  )
  invisible(NULL)
}

draw_shared_map_ornaments <- function(geometry,
                                      scalebar = TRUE,
                                      scalebar_unit = "auto",
                                      scalebar_cex = 0.7,
                                      north_arrow = TRUE,
                                      legend_info = NULL,
                                      map_scale = NULL) {
  if (!isTRUE(scalebar)) {
    if (isTRUE(north_arrow)) {
      noise_map_north_arrow_in_map(cex = scalebar_cex)
    }
    return(invisible(NULL))
  }
  if (is.null(map_scale)) {
    map_scale <- noise_map_scale(geometry)
  }
  if (!is.null(legend_info)) {
    noise_map_scalebar_below_legend(
      map_scale,
      legend_info,
      unit = scalebar_unit,
      cex = scalebar_cex,
      north_arrow = north_arrow
    )
  } else {
    noise_map_scalebar_in_map(
      map_scale,
      unit = scalebar_unit,
      cex = scalebar_cex,
      north_arrow = north_arrow
    )
  }
  invisible(NULL)
}

draw_shared_colorbar_panel <- function(cols,
                                       lo,
                                       hi,
                                       title = "",
                                       cex = 0.85,
                                       mar = c(0.2, 0.2, 0.2, 0.2),
                                       reserve_in = 0) {
  graphics::par(mar = c(mar[1], max(mar[2], 0.8), mar[3], 0.2))
  graphics::plot.new()
  pin <- graphics::par("pin")
  y_centre <- if (reserve_in > 0 && pin[2] > 0) 0.5 + (reserve_in / 2) / pin[2] else 0.5
  bar_h <- 0.58
  bar_w <- 0.16
  left <- 0.05
  bottom <- y_centre - bar_h / 2
  n <- length(cols)
  ys <- seq(bottom, bottom + bar_h, length.out = n + 1L)
  graphics::rect(left, ys[-(n + 1L)], left + bar_w, ys[-1L],
                 col = cols, border = NA, xpd = NA)
  graphics::rect(left, bottom, left + bar_w, bottom + bar_h,
                 col = NA, border = "grey30", lwd = 0.6, xpd = NA)
  if (is.finite(lo) && is.finite(hi)) {
    ticks <- seq(lo, hi, length.out = 5L)
    tick_y <- seq(bottom, bottom + bar_h, length.out = 5L)
    graphics::segments(left + bar_w, tick_y, left + bar_w + 0.035, tick_y,
                       col = "grey30", lwd = 0.6, xpd = NA)
    graphics::text(left + bar_w + 0.055, tick_y,
                   labels = format(signif(ticks, 3L), trim = TRUE),
                   adj = c(0, 0.5), cex = cex * 0.75, xpd = NA)
  }
  if (nzchar(title)) {
    graphics::text(left, bottom + bar_h + 0.08, title,
                   adj = c(0, 0), cex = cex, font = 2, xpd = NA)
  }
  invisible(list(rect = list(left = left, top = bottom + bar_h, w = bar_w, h = bar_h)))
}

plot_shadow_footprints <- function(buildings,
                                   shadows,
                                   canopy = NULL,
                                   plot_overlap_gradient = FALSE,
                                   scalebar = TRUE,
                                   scalebar_unit = c("auto", "km", "m"),
                                   scalebar_cex = 0.7,
                                   north_arrow = TRUE) {
  scalebar_unit <- match.arg(scalebar_unit)
  map_geom <- sf::st_geometry(buildings)
  if (!is.null(shadows) && nrow(shadows) > 0) {
    map_geom <- c(map_geom, sf::st_geometry(shadows))
  }
  old_par <- open_shared_map_layout(map_geom, legend = TRUE)
  on.exit(close_shared_map_layout(old_par, legend = TRUE), add = TRUE)
  draw_shared_map_frame(map_geom)
  has_canopy <- !is.null(canopy) || ("shadow_source" %in% names(shadows) &&
    any(shadows$shadow_source == "canopy", na.rm = TRUE))

  if (nrow(shadows) == 0) {
    if (!is.null(canopy)) {
      plot_canopy_cells(canopy, crs = sf::st_crs(buildings), add = TRUE)
    }
    graphics::plot(sf::st_geometry(buildings), col = "black", border = NA, add = TRUE)
    canopy_legend <- canopy_legend_items(canopy)
    map_scale <- if (isTRUE(scalebar)) noise_map_scale(map_geom) else NULL
    legend_info <- noise_map_draw_legend(
      labels = if (is.null(canopy)) "Building" else c("Building", canopy_legend$labels),
      fills = if (is.null(canopy)) "black" else c("black", canopy_legend$fill),
      title = "Layer",
      reserve_in = if (isTRUE(scalebar)) 0.5 else 0
    )
    draw_shared_map_ornaments(map_geom, scalebar, scalebar_unit,
                              scalebar_cex, north_arrow, legend_info, map_scale)
    return(invisible(NULL))
  }

  if (has_canopy) {
    source_style <- shadow_source_plot_style(shadows, plot_overlap_gradient)
    for (source in names(source_style$fill)) {
      subset <- shadows[shadows$shadow_source == source, ]
      if (nrow(subset) == 0) next
      if (source == "canopy") {
        subset <- order_canopy_shadows(subset)
      }
      graphics::plot(
        sf::st_geometry(subset),
        col = source_style$fill[[source]],
        border = NA,
        add = TRUE
      )
    }

    if (!is.null(canopy)) {
      plot_canopy_cells(canopy, crs = sf::st_crs(buildings), add = TRUE)
    }
    graphics::plot(sf::st_geometry(buildings), col = "black", border = NA, add = TRUE)
    canopy_legend <- canopy_legend_items(canopy)
    map_scale <- if (isTRUE(scalebar)) noise_map_scale(map_geom) else NULL
    legend_info <- noise_map_draw_legend(
      labels = c("Building", canopy_legend$labels, "Building shadow", "Canopy shadow"),
      fills = c("black", canopy_legend$fill, source_style$fill[["building"]], source_style$fill[["canopy"]]),
      title = "Layer",
      reserve_in = if (isTRUE(scalebar)) 0.5 else 0
    )
    draw_shared_map_ornaments(map_geom, scalebar, scalebar_unit,
                              scalebar_cex, north_arrow, legend_info, map_scale)
    return(invisible(NULL))
  }

  shadow_cols <- shadow_plot_cols(shadows, plot_overlap_gradient)

  if ("sun_id" %in% names(shadows)) {
    for (sun_id in names(shadow_cols)) {
      subset <- shadows[shadows$sun_id == sun_id, ]
      if (nrow(subset) == 0) next
      graphics::plot(
        sf::st_geometry(subset),
        col = shadow_cols[[sun_id]],
        border = NA,
        add = TRUE
      )
    }
  } else {
    graphics::plot(
      sf::st_geometry(shadows),
      col = unname(shadow_cols[[1]]),
      border = NA,
      add = TRUE
    )
  }

  graphics::plot(sf::st_geometry(buildings), col = "black", border = NA, add = TRUE)
  map_scale <- if (isTRUE(scalebar)) noise_map_scale(map_geom) else NULL
  legend_labels <- "Shadow"
  legend_cols <- unname(shadow_cols[[1]])
  if ("sun_id" %in% names(shadows) &&
      length(shadow_cols) > 1 &&
      !isTRUE(plot_overlap_gradient)) {
    legend_labels <- paste("Shadow", names(shadow_cols), sep = ": ")
    legend_cols <- unname(shadow_cols)
  }
  legend_info <- noise_map_draw_legend(
    labels = c("Building", legend_labels),
    fills = c("black", legend_cols),
    title = "Layer",
    reserve_in = if (isTRUE(scalebar)) 0.5 else 0
  )
  draw_shared_map_ornaments(map_geom, scalebar, scalebar_unit,
                            scalebar_cex, north_arrow, legend_info, map_scale)
  invisible(NULL)
}

order_canopy_shadows <- function(shadows) {
  if (!"canopy_height" %in% names(shadows)) {
    return(shadows)
  }
  height <- as.numeric(shadows$canopy_height)
  shadows[order(height, na.last = TRUE), ]
}

plot_canopy_cells <- function(canopy, crs = NA, add = TRUE) {
  if (is.null(canopy) || is.null(canopy$xy) || nrow(canopy$xy) == 0) {
    return(invisible(NULL))
  }
  canopy_cells <- canopy_cells_sf(canopy, crs = crs)
  graphics::plot(
    sf::st_geometry(canopy_cells),
    col = canopy_height_cols(canopy_cells$canopy_height),
    border = NA,
    add = add
  )
  invisible(canopy_cells)
}

canopy_height_cols <- function(height) {
  if (length(height) == 0) {
    return(character())
  }
  height <- as.numeric(height)
  cols <- grDevices::colorRampPalette(c("#BBF7D0", "#16A34A"))(100)
  if (all(is.na(height))) {
    return(rep(cols[50], length(height)))
  }
  range_height <- range(height, na.rm = TRUE)
  if (!is.finite(diff(range_height)) || diff(range_height) == 0) {
    return(rep(cols[70], length(height)))
  }
  idx <- round((height - range_height[1]) / diff(range_height) * 99) + 1
  idx[is.na(idx)] <- 50
  cols[pmin(pmax(idx, 1), 100)]
}

canopy_legend_items <- function(canopy) {
  if (is.null(canopy) || is.null(canopy$height) || length(canopy$height) == 0) {
    return(list(labels = "Tree canopy", fill = "#16A34A"))
  }
  height <- as.numeric(canopy$height)
  if (all(is.na(height)) || diff(range(height, na.rm = TRUE)) == 0) {
    return(list(labels = "Tree canopy", fill = canopy_height_cols(height[1])))
  }
  range_height <- range(height, na.rm = TRUE)
  list(
    labels = c(
      sprintf("Tree canopy low (%.1f m)", range_height[1]),
      sprintf("Tree canopy high (%.1f m)", range_height[2])
    ),
    fill = canopy_height_cols(range_height)
  )
}

canopy_cells_sf <- function(canopy, crs = NA) {
  half_cell <- canopy$cell_size / 2
  geoms <- lapply(seq_len(nrow(canopy$xy)), function(i) {
    x <- canopy$xy[i, 1]
    y <- canopy$xy[i, 2]
    sf::st_polygon(list(rbind(
      c(x - half_cell, y - half_cell),
      c(x - half_cell, y + half_cell),
      c(x + half_cell, y + half_cell),
      c(x + half_cell, y - half_cell),
      c(x - half_cell, y - half_cell)
    )))
  })
  sf::st_sf(
    canopy_height = canopy$height,
    geometry = sf::st_sfc(geoms, crs = crs)
  )
}

shadow_source_plot_style <- function(shadows, plot_overlap_gradient = FALSE) {
  multiple_suns <- "sun_id" %in% names(shadows) && length(unique(shadows$sun_id)) > 1
  if (isTRUE(plot_overlap_gradient) && multiple_suns) {
    shadow_fill <- grDevices::adjustcolor("grey20", alpha.f = 0.22)
    fill <- c(
      building = shadow_fill,
      canopy = shadow_fill
    )
  } else {
    shadow_fill <- grDevices::adjustcolor("grey45", alpha.f = 0.38)
    fill <- c(
      building = shadow_fill,
      canopy = shadow_fill
    )
  }
  border <- c(
    building = grDevices::adjustcolor("grey15", alpha.f = 0.75),
    canopy = grDevices::adjustcolor("#146C43", alpha.f = 0.85)
  )
  list(fill = fill, border = border)
}

shadow_plot_cols <- function(shadows, plot_overlap_gradient = FALSE) {
  if (isTRUE(plot_overlap_gradient) &&
      "sun_id" %in% names(shadows) &&
      length(unique(shadows$sun_id)) > 1) {
    return(stats::setNames(
      rep(grDevices::adjustcolor("grey20", alpha.f = 0.22), length(unique(shadows$sun_id))),
      unique(shadows$sun_id)
    ))
  }
  if ("sun_id" %in% names(shadows)) {
    sun_ids <- unique(shadows$sun_id)
    palette <- grDevices::hcl.colors(max(length(sun_ids), 3), "Dark 3")
    return(stats::setNames(
      grDevices::adjustcolor(palette[seq_along(sun_ids)], alpha.f = 0.35),
      sun_ids
    ))
  }
  c(shadow = grDevices::adjustcolor("grey20", alpha.f = 0.35))
}

cast_shadow_geometry <- function(geometry, height, solar_pos, b) {
  shift <- shadow_shift(height, solar_pos)
  original <- sf::st_geometry(sf::st_cast(sf::st_sf(geometry = geometry), "POLYGON"))
  parts <- lapply(seq_along(original), function(i) {
    poly <- original[i]
    shifted <- shift_geometry(poly, shift)
    sides <- shadow_side_polygons(poly, shift)
    sf::st_union(c(poly, shifted, sides))
  })
  geom <- sf::st_union(do.call(c, parts))
  if (is.finite(b) && b > 0) {
    geom <- sf::st_buffer(sf::st_buffer(geom, b), -b)
  }
  sf::st_make_valid(geom)
}

shadow_shift <- function(height, solar_pos) {
  elev <- solar_pos[[2]]
  if (elev <= 0) return(c(dx = NA_real_, dy = NA_real_, length = Inf))
  length <- height / tan(deg2rad(elev))
  az <- deg2rad(solar_pos[[1]])
  c(dx = -sin(az) * length, dy = -cos(az) * length, length = length)
}

shift_geometry <- function(geometry, shift) {
  shifted <- geometry + c(shift[["dx"]], shift[["dy"]])
  sf::st_crs(shifted) <- sf::st_crs(geometry)
  shifted
}

shadow_side_polygons <- function(poly, shift) {
  coords <- sf::st_coordinates(poly)
  rings <- split(as.data.frame(coords), coords[, "L1"])
  side_geoms <- list()
  for (ring in rings) {
    xy <- as.matrix(ring[, c("X", "Y")])
    if (nrow(xy) < 2) next
    for (i in seq_len(nrow(xy) - 1)) {
      p1 <- xy[i, ]
      p2 <- xy[i + 1, ]
      side <- rbind(
        p1,
        p2,
        p2 + c(shift[["dx"]], shift[["dy"]]),
        p1 + c(shift[["dx"]], shift[["dy"]]),
        p1
      )
      side_geoms[[length(side_geoms) + 1]] <- sf::st_polygon(list(side))
    }
  }
  sf::st_sfc(side_geoms, crs = sf::st_crs(poly))
}

prepare_shadow_location <- function(shadow_locations,
                                    buildings,
                                    height_field,
                                    solar_pos,
                                    cell_size,
                                    extent_buffer) {
  if (is.null(shadow_locations)) {
    return(make_shadow_template(buildings, height_field, solar_pos, cell_size, extent_buffer))
  }
  if (inherits(shadow_locations, "SpatRaster")) {
    return(align_shadow_location_spatraster(shadow_locations, buildings))
  }
  if (inherits(shadow_locations, "sf")) {
    return(align_location_points(shadow_locations, buildings, "shadow_locations"))
  }
  stop("`shadow_locations` must be an sf point layer, a terra SpatRaster, or NULL.", call. = FALSE)
}

align_shadow_location_spatraster <- function(shadow_locations, buildings) {
  if (is.na(sf::st_crs(terra::crs(shadow_locations)))) {
    stop("`shadow_locations` must have a valid CRS.", call. = FALSE)
  }
  if (sf::st_crs(terra::crs(shadow_locations)) != sf::st_crs(buildings)) {
    shadow_locations <- terra::project(shadow_locations, sf::st_crs(buildings)$wkt)
  }
  shadow_locations
}

align_location_points <- function(x, buildings, arg_name) {
  if (is.na(sf::st_crs(x))) {
    stop("`", arg_name, "` must have a valid CRS.", call. = FALSE)
  }
  geom_type <- unique(as.character(sf::st_geometry_type(x)))
  if (!all(geom_type %in% c("POINT", "MULTIPOINT"))) {
    stop("`", arg_name, "` must be an sf point layer.", call. = FALSE)
  }
  if (sf::st_crs(x) != sf::st_crs(buildings)) {
    x <- sf::st_transform(x, sf::st_crs(buildings))
  }
  x
}

make_shadow_template <- function(buildings, height_field, solar_pos, cell_size, extent_buffer) {
  if (is.null(extent_buffer)) {
    heights <- buildings[[height_field]]
    extent_buffer <- max(heights, na.rm = TRUE) * 3
    positive_elev <- solar_pos[, 2][solar_pos[, 2] > 0]
    if (length(positive_elev) > 0) {
      extent_buffer <- max(heights, na.rm = TRUE) / tan(min(positive_elev) * pi / 180)
    }
    extent_buffer <- max(extent_buffer, cell_size)
  }
  bbox <- sf::st_bbox(buildings)
  terra::rast(
    xmin = bbox[["xmin"]] - extent_buffer,
    xmax = bbox[["xmax"]] + extent_buffer,
    ymin = bbox[["ymin"]] - extent_buffer,
    ymax = bbox[["ymax"]] + extent_buffer,
    resolution = cell_size,
    crs = sf::st_crs(buildings)$wkt
  )
}

compute_shadow_height_points <- function(location, buildings, solar_pos, height_field, b, canopy = NULL) {
  coords <- sf::st_coordinates(location)[, c("X", "Y"), drop = FALSE]
  out <- matrix(NA_real_, nrow = nrow(coords), ncol = nrow(solar_pos))
  for (j in seq_len(nrow(solar_pos))) {
    out[, j] <- shadow_height_at_xy(coords, buildings, solar_pos[j, ], height_field, b, canopy)
  }
  colnames(out) <- paste0("solar_", seq_len(nrow(solar_pos)))
  out
}

compute_shadow_height_spatraster <- function(location, buildings, solar_pos, height_field, b, canopy = NULL) {
  xy <- terra::xyFromCell(location, seq_len(terra::ncell(location)))
  out <- vector("list", nrow(solar_pos))
  for (j in seq_len(nrow(solar_pos))) {
    r <- location[[1]]
    terra::values(r) <- shadow_height_at_xy(xy, buildings, solar_pos[j, ], height_field, b, canopy)
    names(r) <- paste0("solar_", j)
    out[[j]] <- r
  }
  if (length(out) == 1) out[[1]] else do.call(c, out)
}

shadow_height_at_xy <- function(xy, buildings, solar_pos, height_field, b, canopy = NULL) {
  if (solar_pos[[2]] <= 0) {
    return(rep(Inf, nrow(xy)))
  }
  building_rings <- shadow_building_ring_data(buildings, height_field)
  shift <- shadow_shift(max(buildings[[height_field]], na.rm = TRUE), solar_pos)
  unit_shift <- shift[1:2] / max(buildings[[height_field]], na.rm = TRUE)
  shadowed <- building_shadow_height_cpp(
    xy = xy,
    rings = building_rings$rings,
    ring_building = building_rings$ring_building,
    ring_part = building_rings$ring_part,
    ring_is_hole = building_rings$ring_is_hole,
    heights = building_rings$heights,
    dx_unit = unname(unit_shift[["dx"]]),
    dy_unit = unname(unit_shift[["dy"]]),
    tan_elev = tan(deg2rad(solar_pos[[2]]))
  )
  if (!is.null(canopy)) {
    canopy_height <- canopy_shadow_height_at_xy(xy, canopy, solar_pos)
    shadowed <- pmax(shadowed, canopy_height, na.rm = TRUE)
  }
  shadowed[is.infinite(shadowed)] <- NA_real_
  shadowed
}

shadow_building_ring_data <- function(buildings, height_field) {
  rings <- list()
  ring_building <- integer()
  ring_part <- integer()
  ring_is_hole <- logical()
  part_id <- 0L
  for (i in seq_len(nrow(buildings))) {
    polygons <- sf::st_cast(sf::st_geometry(buildings[i, ]), "POLYGON", warn = FALSE)
    for (poly_i in seq_along(polygons)) {
      part_id <- part_id + 1L
      coords <- sf::st_coordinates(polygons[poly_i])
      ring_ids <- unique(coords[, "L1"])
      for (ring_id in ring_ids) {
        ring_coords <- coords[coords[, "L1"] == ring_id, c("X", "Y"), drop = FALSE]
        if (nrow(ring_coords) < 4) next
        rings[[length(rings) + 1L]] <- unname(ring_coords)
        ring_building <- c(ring_building, i)
        ring_part <- c(ring_part, part_id)
        ring_is_hole <- c(ring_is_hole, ring_id != 1)
      }
    }
  }
  list(
    rings = rings,
    ring_building = ring_building,
    ring_part = ring_part,
    ring_is_hole = ring_is_hole,
    heights = buildings[[height_field]]
  )
}

prepare_radiation_grid <- function(grid, buildings, height_field, grid_res, offset,
                                   ground = FALSE, ground_res = NULL, dem = NULL) {
  if (!is.null(grid)) {
    return(align_location_points(grid, buildings, "grid"))
  }
  surface <- surface_grid_sf(buildings, height_field, grid_res, offset)
  if (isTRUE(ground)) {
    gres <- if (is.null(ground_res)) grid_res else ground_res
    gnd  <- ground_grid_sf(buildings, gres, offset, dem)
    if (nrow(gnd) > 0L) surface <- rbind(surface, gnd)
  }
  surface
}

surface_grid_sf <- function(buildings, height_field, grid_res, offset) {
  roof <- roof_grid_sf(buildings, height_field, grid_res, offset)
  facade <- facade_grid_sf(buildings, height_field, grid_res, offset)
  rbind(roof, facade)
}

roof_grid_sf <- function(buildings, height_field, grid_res, offset) {
  out <- lapply(seq_len(nrow(buildings)), function(i) {
    pts <- sf::st_make_grid(buildings[i, ], cellsize = grid_res, what = "centers")
    pts <- pts[lengths(sf::st_within(pts, sf::st_geometry(buildings[i, ]))) > 0]
    if (length(pts) == 0) {
      pts <- sf::st_centroid(sf::st_geometry(buildings[i, ]))
    }
    xy <- sf::st_coordinates(pts)
    n <- nrow(xy)
    sf::st_sf(
      building_id = rep(i, n),
      surface = rep("roof", n),
      z = rep(unname(buildings[[height_field]][i]) + offset, n),
      nx = rep(0, n),
      ny = rep(0, n),
      nz = rep(1, n),
      svf = rep(1, n),
      geometry = sf::st_sfc(
        lapply(seq_len(n), function(k) sf::st_point(c(xy[k, 1], xy[k, 2]))),
        crs = sf::st_crs(buildings)
      )
    )
  })
  do.call(rbind, out)
}

facade_grid_sf <- function(buildings, height_field, grid_res, offset) {
  out <- list()
  for (i in seq_len(nrow(buildings))) {
    rings <- polygon_rings(sf::st_geometry(buildings[i, ]))
    for (ring in rings) {
      xy <- ring
      # Determine ring winding order via signed area (shoelace formula).
      # Positive = CCW (y-up), negative = CW.  `c(edge[2], -edge[1])` gives the
      # outward normal for CCW rings; flip sign for CW rings so normals always
      # point outward regardless of how the source data wound the polygon.
      n_ring <- nrow(xy)
      signed_area <- sum(
        xy[seq_len(n_ring - 1), 1] * xy[seq_len(n_ring - 1) + 1, 2] -
        xy[seq_len(n_ring - 1) + 1, 1] * xy[seq_len(n_ring - 1), 2]
      ) / 2
      orient <- if (signed_area >= 0) 1L else -1L
      for (e in seq_len(nrow(xy) - 1)) {
        p1 <- xy[e, ]
        p2 <- xy[e + 1, ]
        edge <- p2 - p1
        edge_len <- sqrt(sum(edge^2))
        if (!is.finite(edge_len) || edge_len == 0) next
        n_h <- max(1, ceiling(edge_len / grid_res))
        n_v <- max(1, ceiling(buildings[[height_field]][i] / grid_res))
        mids <- (seq_len(n_h) - 0.5) / n_h
        zs <- (seq_len(n_v) - 0.5) / n_v * buildings[[height_field]][i] + offset
        pts <- do.call(rbind, lapply(mids, function(m) p1 + m * edge))
        pts <- pts[rep(seq_len(nrow(pts)), each = length(zs)), , drop = FALSE]
        z <- rep(zs, times = length(mids))
        normal <- orient * c(edge[2], -edge[1]) / edge_len
        n <- nrow(pts)
        out[[length(out) + 1]] <- sf::st_sf(
          building_id = rep(i, n),
          surface = rep("facade", n),
          z = z,
          nx = rep(unname(normal[1]), n),
          ny = rep(unname(normal[2]), n),
          nz = rep(0, n),
          svf = rep(0.5, n),
          geometry = sf::st_sfc(
            lapply(seq_len(n), function(k) sf::st_point(c(pts[k, 1], pts[k, 2]))),
            crs = sf::st_crs(buildings)
          )
        )
      }
    }
  }
  do.call(rbind, out)
}

# `dem` is accepted for API symmetry but intentionally unused: ground `z` is a
# height above local ground, and terrain enters through surface_ground_elevation().
ground_grid_sf <- function(buildings, grid_res, offset, dem = NULL) {
  bbox_poly <- sf::st_as_sfc(sf::st_bbox(buildings))
  pts <- sf::st_make_grid(bbox_poly, cellsize = grid_res, what = "centers")
  # Drop points that fall inside any building footprint
  in_building <- lengths(sf::st_within(pts, sf::st_union(sf::st_geometry(buildings)))) > 0L
  pts <- pts[!in_building]
  if (length(pts) == 0L) {
    return(sf::st_sf(
      building_id = integer(0), surface = character(0), z = numeric(0),
      nx = numeric(0), ny = numeric(0), nz = numeric(0), svf = numeric(0),
      geometry = sf::st_sfc(crs = sf::st_crs(buildings))
    ))
  }
  pts_sf <- sf::st_sf(geometry = pts, crs = sf::st_crs(buildings))
  xy <- sf::st_coordinates(pts_sf)
  n  <- nrow(xy)
  # `z` must be HEIGHT ABOVE LOCAL GROUND, matching roof and facade points.
  # Shadow heights from shadow_height_at_xy() / canopy_shadow_height_at_xy()
  # are also expressed above local ground, so adding absolute DEM elevation
  # here would make every shadow test fail and leave the ground always sunlit.
  # Terrain is accounted for separately via surface_ground_elevation().
  sf::st_sf(
    building_id = rep(NA_integer_, n),
    surface     = rep("ground", n),
    z           = rep(offset, n),
    nx          = rep(0, n),
    ny          = rep(0, n),
    nz          = rep(1, n),
    svf         = rep(1, n),
    geometry    = sf::st_sfc(
      lapply(seq_len(n), function(k) sf::st_point(c(xy[k, 1L], xy[k, 2L]))),
      crs = sf::st_crs(buildings)
    )
  )
}

polygon_rings <- function(geometry) {
  poly <- sf::st_cast(sf::st_sf(geometry = geometry), "POLYGON")
  rings <- list()
  for (i in seq_along(sf::st_geometry(poly))) {
    coords <- sf::st_coordinates(sf::st_geometry(poly)[i])
    split_coords <- split(as.data.frame(coords), coords[, "L1"])
    rings <- c(rings, lapply(split_coords, function(x) as.matrix(x[, c("X", "Y")])))
  }
  rings
}

check_radiation_vectors <- function(solar_pos, solar_normal, solar_diffuse) {
  if (length(solar_normal) != nrow(solar_pos) || length(solar_diffuse) != nrow(solar_pos)) {
    stop("`solar_normal` and `solar_diffuse` must match the number of solar positions.", call. = FALSE)
  }
}

check_canopy_transmissivity <- function(canopy_transmissivity) {
  if (!is.numeric(canopy_transmissivity) || length(canopy_transmissivity) != 1 ||
      is.na(canopy_transmissivity) || canopy_transmissivity < 0 ||
      canopy_transmissivity > 1) {
    stop("`canopy_transmissivity` must be a numeric value from 0 to 1.", call. = FALSE)
  }
}

resolve_shadow_raster_inputs <- function(buildings,
                                         solar_pos,
                                         height_field,
                                         canopy_height,
                                         dem,
                                         datasource_canopy_height,
                                         key,
                                         min_tree_height,
                                         raster_buffer,
                                         quiet) {
  datasource_canopy_height <- normalize_shadow_source(
    datasource_canopy_height,
    choices = c("metachm", "ethchm")
  )
  needs_canopy <- is.null(canopy_height) && !is.null(datasource_canopy_height)
  uses_canopy <- needs_canopy || !is.null(canopy_height)
  needs_dem <- is.null(dem) && !is.null(key) && uses_canopy
  if (!needs_canopy && !needs_dem) {
    return(list(canopy_height = canopy_height, dem = dem))
  }
  if (needs_dem && is.null(key)) {
    stop("`key` is required to retrieve DEM data internally.", call. = FALSE)
  }
  analysis_bbox <- shadow_analysis_bbox(
    buildings = buildings,
    solar_pos = solar_pos,
    height_field = height_field,
    raster_buffer = raster_buffer
  )
  bbox_vector <- as_wgs84_bbox_vector(analysis_bbox)
  fetched <- list(canopy_height = canopy_height, dem = dem)
  if (needs_dem) {
    if (!quiet) cli::cli_alert_info("Downloading DEM for shadow/radiation analysis ...")
    fetched$dem <- get_dem(bbox_vector, key)
  }
  if (needs_canopy) {
    if (!quiet) cli::cli_alert_info("Downloading canopy height map for shadow/radiation analysis ...")
    chm_layers <- suppressMessages(get_chm(
      bbox_vector,
      min_tree_height,
      datasource = datasource_canopy_height
    ))
    fetched$canopy_height <- chm_layers[[1]]
  }
  fetched
}

normalize_shadow_source <- function(value, choices) {
  if (inherits(value, "NULL")) return(NULL)
  value <- tolower(as.character(value[1]))
  if (value %in% c("none", "null", "na")) return(NULL)
  match.arg(value, choices)
}

shadow_analysis_bbox <- function(buildings, solar_pos, height_field, raster_buffer) {
  if (is.null(raster_buffer)) {
    raster_buffer <- estimate_shadow_buffer(buildings, solar_pos, height_field)
  }
  sf::st_buffer(sf::st_as_sfc(sf::st_bbox(buildings)), dist = raster_buffer)
}

estimate_shadow_buffer <- function(buildings, solar_pos, height_field) {
  positive_elev <- solar_pos[, 2][solar_pos[, 2] > 0]
  if (length(positive_elev) == 0) {
    return(max(buildings[[height_field]], na.rm = TRUE) * 3)
  }
  buffer <- max(buildings[[height_field]], na.rm = TRUE) /
    tan(min(positive_elev) * pi / 180)
  max(buffer, max(buildings[[height_field]], na.rm = TRUE), 100)
}

estimate_svf_buffer <- function(buildings, height_field, canopy_height = NULL, max_distance = NULL) {
  if (!is.null(max_distance) && is.finite(max_distance)) {
    return(max_distance)
  }
  max_height <- max(buildings[[height_field]], na.rm = TRUE)
  if (inherits(canopy_height, "SpatRaster")) {
    canopy_max <- suppressWarnings(terra::global(canopy_height[[1]], "max", na.rm = TRUE)[1, 1])
    if (is.finite(canopy_max)) {
      max_height <- max(max_height, canopy_max, na.rm = TRUE)
    }
  }
  max(max_height * 5, 100)
}

make_svf_template <- function(buildings, grid_res, extent_buffer) {
  bbox <- sf::st_bbox(buildings)
  terra::rast(
    xmin = bbox[["xmin"]] - extent_buffer,
    xmax = bbox[["xmax"]] + extent_buffer,
    ymin = bbox[["ymin"]] - extent_buffer,
    ymax = bbox[["ymax"]] + extent_buffer,
    resolution = grid_res,
    crs = sf::st_crs(buildings)$wkt
  )
}

compute_svf_raster <- function(template,
                               buildings,
                               height_field,
                               canopy,
                               dem,
                               res_angle,
                               observer_height,
                               max_distance) {
  xy <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
  building_rings <- shadow_building_ring_data(buildings, height_field)
  observer_ground <- rep(0, nrow(xy))
  if (!is.null(dem)) {
    dem <- align_spatraster_to_template(dem, template, "dem")
    observer_ground <- terra::extract(dem, xy)[, 1]
    observer_ground[is.na(observer_ground)] <- 0
  }
  building_ground <- rep(0, nrow(buildings))
  if (!is.null(dem)) {
    building_xy <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(buildings)))[, c("X", "Y"), drop = FALSE]
    building_ground <- terra::extract(dem, building_xy)[, 1]
    building_ground[is.na(building_ground)] <- 0
  }
  canopy_xy <- matrix(numeric(), ncol = 2)
  canopy_h <- numeric()
  canopy_ground <- numeric()
  canopy_cell_size <- 0
  if (!is.null(canopy)) {
    canopy_xy <- canopy$xy
    canopy_h <- canopy$height
    canopy_ground <- canopy$ground
    canopy_cell_size <- canopy$cell_size
  }
  azimuth <- seq(0, 360 - res_angle, by = res_angle)
  values <- svf_cpp(
    xy = xy,
    azimuth = azimuth,
    rings = building_rings$rings,
    ring_building = building_rings$ring_building,
    building_height = building_rings$heights,
    building_ground = building_ground,
    canopy_xy = canopy_xy,
    canopy_height = canopy_h,
    canopy_ground = canopy_ground,
    observer_ground = observer_ground,
    observer_height = rep(observer_height, nrow(xy)),
    canopy_cell_size = canopy_cell_size,
    max_distance = if (is.null(max_distance)) Inf else max_distance
  )
  building_mask <- terra::rasterize(terra::vect(buildings), template, field = 1, background = NA)
  values[!is.na(terra::values(building_mask, mat = FALSE))] <- NA_real_
  out <- template
  terra::values(out) <- values
  names(out) <- "svf"
  out
}

surface_ground_elevation <- function(surface, dem) {
  ground <- rep(0, nrow(surface))
  if (!is.null(dem)) {
    xy <- sf::st_coordinates(surface)[, c("X", "Y"), drop = FALSE]
    ground <- terra::extract(dem, xy)[, 1]
    ground[is.na(ground)] <- 0
  }
  ground
}

building_ground_elevation <- function(buildings, dem) {
  ground <- rep(0, nrow(buildings))
  if (!is.null(dem)) {
    building_xy <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(buildings)))[, c("X", "Y"), drop = FALSE]
    ground <- terra::extract(dem, building_xy)[, 1]
    ground[is.na(ground)] <- 0
  }
  ground
}

surface_svf <- function(surface,
                        buildings,
                        height_field,
                        canopy,
                        dem,
                        res_angle,
                        max_distance) {
  if (!is.null(dem)) {
    dem <- align_spatraster_to_buildings(dem, buildings, "dem")
  }
  xy <- sf::st_coordinates(surface)[, c("X", "Y"), drop = FALSE]
  building_rings <- shadow_building_ring_data(buildings, height_field)
  observer_ground <- surface_ground_elevation(surface, dem)
  building_ground <- building_ground_elevation(buildings, dem)
  canopy_xy <- matrix(numeric(), ncol = 2)
  canopy_h <- numeric()
  canopy_ground <- numeric()
  canopy_cell_size <- 0
  if (!is.null(canopy)) {
    canopy_xy <- canopy$xy
    canopy_h <- canopy$height
    canopy_ground <- canopy$ground
    canopy_cell_size <- canopy$cell_size
  }
  azimuth <- seq(0, 360 - res_angle, by = res_angle)
  pmin(pmax(svf_cpp(
    xy = xy,
    azimuth = azimuth,
    rings = building_rings$rings,
    ring_building = building_rings$ring_building,
    building_height = building_rings$heights,
    building_ground = building_ground,
    canopy_xy = canopy_xy,
    canopy_height = canopy_h,
    canopy_ground = canopy_ground,
    observer_ground = observer_ground,
    observer_height = surface$z,
    canopy_cell_size = canopy_cell_size,
    max_distance = if (is.null(max_distance)) Inf else max_distance
  ), 0), 1)
}

plot_svf_raster <- function(buildings,
                            svf_raster,
                            scalebar = TRUE,
                            scalebar_unit = c("auto", "km", "m"),
                            scalebar_cex = 0.7,
                            north_arrow = TRUE) {
  scalebar_unit <- match.arg(scalebar_unit)
  cols <- grDevices::hcl.colors(100, "Cividis", rev = TRUE)
  old_par <- open_shared_map_layout(sf::st_geometry(buildings), legend = TRUE)
  on.exit(close_shared_map_layout(old_par, legend = TRUE), add = TRUE)

  graphics::plot(
    svf_raster,
    col = cols,
    legend = FALSE,
    axes = FALSE,
    main = ""
  )
  draw_shared_map_title_below(sf::st_bbox(buildings), "Sky View Factor")
  graphics::plot(sf::st_geometry(buildings), col = "black", border = NA, add = TRUE)
  map_scale <- if (isTRUE(scalebar)) noise_map_scale(sf::st_geometry(buildings)) else NULL
  rng <- range(terra::values(svf_raster), na.rm = TRUE)
  legend_info <- draw_shared_colorbar_panel(
    cols,
    rng[1L],
    rng[2L],
    title = "SVF",
    reserve_in = if (isTRUE(scalebar)) 0.5 else 0
  )
  draw_shared_map_ornaments(sf::st_geometry(buildings), scalebar, scalebar_unit,
                            scalebar_cex, north_arrow, legend_info, map_scale)
  invisible(NULL)
}

prepare_canopy_obstacles <- function(canopy_height, dem, buildings, min_tree_height) {
  if (is.null(canopy_height)) {
    return(NULL)
  }
  if (!inherits(canopy_height, "SpatRaster")) {
    stop("`canopy_height` must be a terra SpatRaster.", call. = FALSE)
  }
  canopy_height <- align_spatraster_to_buildings(canopy_height, buildings, "canopy_height")
  if (!is.null(dem)) {
    if (!inherits(dem, "SpatRaster")) {
      stop("`dem` must be a terra SpatRaster.", call. = FALSE)
    }
    dem <- align_spatraster_to_template(dem, canopy_height, "dem")
  }
  canopy_vals <- terra::values(canopy_height[[1]], mat = FALSE)
  keep <- !is.na(canopy_vals) & canopy_vals >= min_tree_height
  if (!any(keep)) {
    return(NULL)
  }
  cells <- which(keep)
  xy <- terra::xyFromCell(canopy_height, cells)
  ground <- rep(0, length(cells))
  if (!is.null(dem)) {
    ground <- terra::values(dem[[1]], mat = FALSE)[cells]
    ground[is.na(ground)] <- 0
  }
  list(
    xy = xy,
    height = canopy_vals[cells],
    ground = ground,
    dem = dem,
    cell_size = max(terra::res(canopy_height))
  )
}

align_spatraster_to_buildings <- function(r, buildings, arg_name) {
  if (is.na(sf::st_crs(terra::crs(r)))) {
    stop("`", arg_name, "` must have a valid CRS.", call. = FALSE)
  }
  if (sf::st_crs(terra::crs(r)) != sf::st_crs(buildings)) {
    r <- terra::project(r, sf::st_crs(buildings)$wkt)
  }
  r
}

align_spatraster_to_template <- function(r, template, arg_name) {
  if (is.na(sf::st_crs(terra::crs(r)))) {
    stop("`", arg_name, "` must have a valid CRS.", call. = FALSE)
  }
  if (terra::crs(r) != terra::crs(template)) {
    r <- terra::project(r, terra::crs(template))
  }
  if (!terra::compareGeom(r, template, stopOnError = FALSE, crs = TRUE,
                          ext = TRUE, rowcol = TRUE, res = TRUE)) {
    r <- terra::resample(r, template, method = "bilinear")
  }
  r
}

canopy_shadow_height_at_xy <- function(xy, canopy, solar_pos) {
  if (is.null(canopy) || nrow(canopy$xy) == 0) {
    return(rep(NA_real_, nrow(xy)))
  }
  if (solar_pos[[2]] <= 0) {
    return(rep(Inf, nrow(xy)))
  }
  target_ground <- rep(0, nrow(xy))
  if (!is.null(canopy$dem)) {
    target_ground <- terra::extract(canopy$dem, xy)[, 1]
    target_ground[is.na(target_ground)] <- 0
  }
  shift_per_height <- canopy_shadow_unit_shift(solar_pos)
  out <- canopy_shadow_height_cpp(
    xy = xy,
    canopy_xy = canopy$xy,
    canopy_height = canopy$height,
    canopy_ground = canopy$ground,
    target_ground = target_ground,
    dx_unit = shift_per_height[["dx"]],
    dy_unit = shift_per_height[["dy"]],
    tan_elev = tan(deg2rad(solar_pos[[2]])),
    half_cell = canopy$cell_size / sqrt(2)
  )
  out[is.infinite(out)] <- NA_real_
  out
}

canopy_shadow_unit_shift <- function(solar_pos) {
  elev <- solar_pos[[2]]
  if (elev <= 0) {
    return(c(dx = NA_real_, dy = NA_real_))
  }
  length_per_height <- 1 / tan(deg2rad(elev))
  az <- deg2rad(solar_pos[[1]])
  c(dx = -sin(az) * length_per_height, dy = -cos(az) * length_per_height)
}

canopy_shadow_shift <- function(height, solar_pos) {
  elev <- solar_pos[[2]]
  if (elev <= 0) {
    return(list(dx = rep(NA_real_, length(height)), dy = rep(NA_real_, length(height))))
  }
  length <- height / tan(deg2rad(elev))
  az <- deg2rad(solar_pos[[1]])
  list(dx = -sin(az) * length, dy = -cos(az) * length)
}

compute_surface_radiation <- function(surface, buildings, solar_pos, solar_normal,
                                      solar_diffuse, height_field, canopy,
                                      canopy_transmissivity,
                                      dem = NULL,
                                      svf_res_angle = 15,
                                      radius = Inf) {
  xy <- sf::st_coordinates(surface)[, c("X", "Y"), drop = FALSE]
  z <- surface$z
  svf_values <- surface_svf(
    surface = surface,
    buildings = buildings,
    height_field = height_field,
    canopy = canopy,
    dem = dem,
    res_angle = svf_res_angle,
    max_distance = radius
  )
  direct_by_time <- matrix(0, nrow = nrow(surface), ncol = nrow(solar_pos))
  diffuse_by_time <- matrix(0, nrow = nrow(surface), ncol = nrow(solar_pos))
  for (j in seq_len(nrow(solar_pos))) {
    elev <- solar_pos[j, 2]
    if (elev <= 0) next
    sun_vec <- sun_vector(solar_pos[j, ])
    incidence <- pmax(0, surface$nx * sun_vec[1] + surface$ny * sun_vec[2] + surface$nz * sun_vec[3])
    building_shadow_h <- shadow_height_at_xy(
      xy, buildings, solar_pos[j, ], height_field, b = 0
    )
    canopy_shadow_h <- canopy_shadow_height_at_xy(xy, canopy, solar_pos[j, ])
    shaded_by_building <- !is.na(building_shadow_h) & building_shadow_h > z
    shaded_by_canopy <- !is.na(canopy_shadow_h) & canopy_shadow_h > z
    transmission <- rep(1, nrow(surface))
    transmission[shaded_by_canopy] <- canopy_transmissivity
    transmission[shaded_by_building] <- 0
    direct_by_time[, j] <- solar_normal[j] * incidence * transmission
    diffuse_by_time[, j] <- solar_diffuse[j] * svf_values
  }
  list(
    svf = svf_values,
    direct = rowSums(direct_by_time),
    diffuse = rowSums(diffuse_by_time),
    total = rowSums(direct_by_time) + rowSums(diffuse_by_time),
    by_time = list(direct = direct_by_time, diffuse = diffuse_by_time)
  )
}

draw_colorbar <- function(cols, lo, hi, title = "", n_ticks = 5L) {
  # Vertical gradient colour bar drawn just outside the right edge of the
  # current panel.  All geometry is in *user* coordinates derived from
  # par("usr"), and xpd = NA allows drawing into the figure margin.
  usr <- graphics::par("usr")   # c(x1, x2, y1, y2), user coordinates
  xw  <- usr[2L] - usr[1L]
  yh  <- usr[4L] - usr[3L]
  if (!is.finite(xw) || !is.finite(yh) || xw <= 0 || yh <= 0) {
    return(invisible(NULL))
  }
  bar_x0 <- usr[2L] + xw * 0.03
  bar_x1 <- bar_x0  + xw * 0.035
  bar_y0 <- usr[3L] + yh * 0.20
  bar_y1 <- usr[3L] + yh * 0.80

  old_xpd <- graphics::par("xpd")
  on.exit(graphics::par(xpd = old_xpd), add = TRUE)
  graphics::par(xpd = NA)

  n  <- length(cols)
  ys <- seq(bar_y0, bar_y1, length.out = n + 1L)
  graphics::rect(bar_x0, ys[-(n + 1L)], bar_x1, ys[-1L],
                 col = cols, border = NA)
  graphics::rect(bar_x0, bar_y0, bar_x1, bar_y1,
                 col = NA, border = "grey30", lwd = 0.6)

  # Numeric tick labels along the right side of the bar
  if (is.finite(lo) && is.finite(hi)) {
    tick_v <- seq(lo, hi, length.out = n_ticks)
    tick_y <- seq(bar_y0, bar_y1, length.out = n_ticks)
    graphics::segments(bar_x1, tick_y, bar_x1 + xw * 0.008, tick_y,
                       col = "grey30", lwd = 0.6)
    graphics::text(bar_x1 + xw * 0.014, tick_y,
                   labels = format(signif(tick_v, 3L), trim = TRUE),
                   adj = c(0, 0.5), cex = 0.6)
  }
  # Title above the bar
  if (nzchar(title)) {
    graphics::text((bar_x0 + bar_x1) / 2, bar_y1 + yh * 0.04,
                   labels = title, adj = c(0.5, 0), cex = 0.65)
  }
  invisible(NULL)
}

plot_radiation_surface <- function(buildings,
                                   radiation,
                                   scalebar = TRUE,
                                   scalebar_unit = c("auto", "km", "m"),
                                   scalebar_cex = 0.7,
                                   north_arrow = TRUE) {
  scalebar_unit <- match.arg(scalebar_unit)
  if (nrow(radiation) == 0) {
    old_par <- open_shared_map_layout(sf::st_geometry(buildings), legend = FALSE)
    on.exit(close_shared_map_layout(old_par, legend = FALSE), add = TRUE)
    draw_shared_map_frame(sf::st_geometry(buildings), main = "Radiation")
    graphics::plot(sf::st_geometry(buildings), col = "grey95", border = "grey40", add = TRUE)
    draw_shared_map_ornaments(sf::st_geometry(buildings), scalebar, scalebar_unit,
                              scalebar_cex, north_arrow, NULL)
    return(invisible(NULL))
  }
  rad_cols <- radiation_palette(100L)
  make_idx <- function(vals, rng) {
    if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(50L, length(vals)))
    idx <- round((vals - rng[1L]) / diff(rng) * 99) + 1L
    idx[is.na(idx)] <- 1L
    pmin(pmax(idx, 1L), 100L)
  }
  ground     <- radiation$surface == "ground"
  facade     <- radiation$surface == "facade"
  roof       <- radiation$surface == "roof"
  has_ground <- any(ground)

  if (has_ground) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    on.exit(graphics::layout(1L), add = TRUE)
    graphics::layout(matrix(c(1L, 2L, 3L, 4L), nrow = 1L),
                     widths = c(1, 1, 1, 0.32))
    graphics::par(mar = c(2.0, 0.2, 0.2, 0.2))

    gnd_sf <- radiation[ground, ]
    bb     <- sf::st_bbox(buildings)
    n_gnd  <- sum(ground)
    gres   <- sqrt(
      (bb["xmax"] - bb["xmin"]) * (bb["ymax"] - bb["ymin"]) / max(n_gnd, 1L)
    )
    tmpl <- terra::rast(
      xmin = unname(bb["xmin"]), xmax = unname(bb["xmax"]),
      ymin = unname(bb["ymin"]), ymax = unname(bb["ymax"]),
      resolution = gres, crs = sf::st_crs(buildings)$wkt
    )
    bld_v <- terra::vect(sf::st_geometry(buildings))

    # Shared W/m^2 scale across all 2D map panels so colours are comparable.
    total_rng <- range(radiation$total, na.rm = TRUE)
    plot_rng <- total_rng
    if (!is.finite(diff(plot_rng)) || diff(plot_rng) == 0) {
      pad <- max(abs(plot_rng[1L]) * 0.01, 1)
      plot_rng <- c(plot_rng[1L] - pad, plot_rng[2L] + pad)
    }
    brks <- seq(plot_rng[1L], plot_rng[2L], length.out = 101L)

    # Both panels are set up with an identical empty base plot (same limits,
    # same aspect) so their plot regions match and the titles line up.
    new_panel <- function(main) {
      graphics::plot(NA, type = "n", asp = 1L, axes = FALSE,
                     xlab = "", ylab = "", main = "",
                     xlim = c(bb[["xmin"]], bb[["xmax"]]),
                     ylim = c(bb[["ymin"]], bb[["ymax"]]))
      draw_shared_map_title_below(bb, main)
    }

    # --- Left panel: ground total radiation ---
    gnd_r <- terra::rasterize(terra::vect(gnd_sf), tmpl,
                              field = "total", fun = "mean")
    gnd_r <- terra::mask(gnd_r, bld_v, inverse = TRUE)
    new_panel("Ground total radiation (W/m2)")
    # Drawn with base graphics::image rather than terra::plot, because
    # terra::plot mutates par("mar") and would shift the next panel's title.
    gm <- terra::as.matrix(gnd_r, wide = TRUE)
    graphics::image(
      x = seq(terra::xmin(gnd_r), terra::xmax(gnd_r), length.out = ncol(gm) + 1L),
      y = seq(terra::ymin(gnd_r), terra::ymax(gnd_r), length.out = nrow(gm) + 1L),
      z = t(gm[nrow(gm):1L, , drop = FALSE]),
      col = rad_cols, breaks = brks, add = TRUE, useRaster = TRUE
    )
    graphics::plot(sf::st_geometry(buildings),
                   col = "grey92", border = "grey35", lwd = 0.6, add = TRUE)
    map_scale <- if (isTRUE(scalebar)) noise_map_scale(sf::st_geometry(buildings)) else NULL

    # --- Middle panel: facade total radiation ---
    new_panel("Facade total radiation (W/m2)")
    graphics::plot(sf::st_geometry(buildings), col = "grey95",
                   border = "grey45", add = TRUE)
    if (any(facade)) {
      graphics::points(
        sf::st_coordinates(radiation[facade, ]),
        pch = 16, cex = 0.35,
        col = rad_cols[make_idx(radiation$total[facade], plot_rng)]
      )
    }
    graphics::plot(sf::st_geometry(buildings), col = NA, border = "grey30", add = TRUE)

    # --- Right panel: roof total radiation ---
    new_panel("Roof total radiation (W/m2)")
    graphics::plot(sf::st_geometry(buildings), col = "grey95",
                   border = "grey45", add = TRUE)
    if (any(roof)) {
      graphics::points(
        sf::st_coordinates(radiation[roof, ]),
        pch = 16, cex = 0.75,
        col = rad_cols[make_idx(radiation$total[roof], plot_rng)]
      )
    }
    graphics::plot(sf::st_geometry(buildings), col = NA, border = "grey30", add = TRUE)
    # Single shared colour bar for all panels.
    legend_info <- draw_shared_colorbar_panel(
      rad_cols,
      total_rng[1L],
      total_rng[2L],
      title = "W/m2",
      reserve_in = if (isTRUE(scalebar)) 0.5 else 0
    )
    draw_shared_map_ornaments(
      sf::st_geometry(buildings),
      scalebar = scalebar,
      scalebar_unit = scalebar_unit,
      scalebar_cex = scalebar_cex,
      north_arrow = north_arrow,
      legend_info = legend_info,
      map_scale = map_scale
    )

  } else {
    # Building-only layout
    old_par <- open_shared_map_layout(sf::st_geometry(buildings), legend = TRUE)
    on.exit(close_shared_map_layout(old_par, legend = TRUE), add = TRUE)
    total_rng <- range(radiation$total, na.rm = TRUE)
    draw_shared_map_frame(sf::st_geometry(buildings), main = "Total radiation (W/m2)")
    graphics::plot(sf::st_geometry(buildings), col = "grey95", border = "grey45", add = TRUE)
    if (any(facade)) {
      graphics::plot(
        sf::st_geometry(radiation[facade, ]),
        pch = 16, cex = 0.25,
        col = "grey60",
        add = TRUE
      )
    }
    if (any(roof)) {
      graphics::plot(
        sf::st_geometry(radiation[roof, ]),
        pch = 16, cex = 0.75,
        col = rad_cols[make_idx(radiation$total[roof], total_rng)],
        add = TRUE
      )
    }
    graphics::plot(sf::st_geometry(buildings), col = NA, border = "grey30", add = TRUE)
    map_scale <- if (isTRUE(scalebar)) noise_map_scale(sf::st_geometry(buildings)) else NULL
    legend_info <- draw_shared_colorbar_panel(
      rad_cols,
      total_rng[1L],
      total_rng[2L],
      title = "W/m2",
      reserve_in = if (isTRUE(scalebar)) 0.5 else 0
    )
    draw_shared_map_ornaments(sf::st_geometry(buildings), scalebar, scalebar_unit,
                              scalebar_cex, north_arrow, legend_info, map_scale)
  }
  invisible(NULL)
}

plot_radiation_canopy_impact <- function(buildings,
                                         radiation,
                                         no_canopy,
                                         scalebar = TRUE,
                                         scalebar_unit = c("auto", "km", "m"),
                                         scalebar_cex = 0.7,
                                         north_arrow = TRUE) {
  scalebar_unit <- match.arg(scalebar_unit)
  if (nrow(radiation) == 0 || is.null(no_canopy$total)) {
    return(invisible(NULL))
  }
  impact <- radiation
  impact$total_diff <- radiation$total - no_canopy$total
  cols <- grDevices::colorRampPalette(c("#2C7BB6", "#FFFFBF", "#D7191C"))(100L)
  make_idx <- function(vals, rng) {
    if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(50L, length(vals)))
    idx <- round((vals - rng[1L]) / diff(rng) * 99) + 1L
    idx[is.na(idx)] <- 50L
    pmin(pmax(idx, 1L), 100L)
  }
  ground <- impact$surface == "ground"
  facade <- impact$surface == "facade"
  roof <- impact$surface == "roof"
  max_abs <- max(abs(impact$total_diff), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  diff_rng <- c(-max_abs, max_abs)

  if (any(ground)) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    on.exit(graphics::layout(1L), add = TRUE)
    graphics::layout(matrix(c(1L, 2L, 3L, 4L), nrow = 1L),
                     widths = c(1, 1, 1, 0.32))
    graphics::par(mar = c(2.0, 0.2, 0.2, 0.2))

    bb <- sf::st_bbox(buildings)
    n_gnd <- sum(ground)
    gres <- sqrt(
      (bb["xmax"] - bb["xmin"]) * (bb["ymax"] - bb["ymin"]) / max(n_gnd, 1L)
    )
    tmpl <- terra::rast(
      xmin = unname(bb["xmin"]), xmax = unname(bb["xmax"]),
      ymin = unname(bb["ymin"]), ymax = unname(bb["ymax"]),
      resolution = gres, crs = sf::st_crs(buildings)$wkt
    )
    bld_v <- terra::vect(sf::st_geometry(buildings))
    brks <- seq(diff_rng[1L], diff_rng[2L], length.out = 101L)
    new_panel <- function(main) {
      graphics::plot(NA, type = "n", asp = 1L, axes = FALSE,
                     xlab = "", ylab = "", main = "",
                     xlim = c(bb[["xmin"]], bb[["xmax"]]),
                     ylim = c(bb[["ymin"]], bb[["ymax"]]))
      draw_shared_map_title_below(bb, main)
    }

    gnd_r <- terra::rasterize(terra::vect(impact[ground, ]), tmpl,
                              field = "total_diff", fun = "mean")
    gnd_r <- terra::mask(gnd_r, bld_v, inverse = TRUE)
    new_panel("Ground canopy impact (W/m2)")
    gm <- terra::as.matrix(gnd_r, wide = TRUE)
    graphics::image(
      x = seq(terra::xmin(gnd_r), terra::xmax(gnd_r), length.out = ncol(gm) + 1L),
      y = seq(terra::ymin(gnd_r), terra::ymax(gnd_r), length.out = nrow(gm) + 1L),
      z = t(gm[nrow(gm):1L, , drop = FALSE]),
      col = cols, breaks = brks, add = TRUE, useRaster = TRUE
    )
    graphics::plot(sf::st_geometry(buildings),
                   col = "grey92", border = "grey35", lwd = 0.6, add = TRUE)
    map_scale <- if (isTRUE(scalebar)) noise_map_scale(sf::st_geometry(buildings)) else NULL

    new_panel("Facade canopy impact (W/m2)")
    graphics::plot(sf::st_geometry(buildings), col = "grey95",
                   border = "grey45", add = TRUE)
    if (any(facade)) {
      graphics::points(
        sf::st_coordinates(impact[facade, ]),
        pch = 16, cex = 0.35,
        col = cols[make_idx(impact$total_diff[facade], diff_rng)]
      )
    }
    graphics::plot(sf::st_geometry(buildings), col = NA, border = "grey30", add = TRUE)

    new_panel("Roof canopy impact (W/m2)")
    graphics::plot(sf::st_geometry(buildings), col = "grey95",
                   border = "grey45", add = TRUE)
    if (any(roof)) {
      graphics::points(
        sf::st_coordinates(impact[roof, ]),
        pch = 16, cex = 0.75,
        col = cols[make_idx(impact$total_diff[roof], diff_rng)]
      )
    }
    graphics::plot(sf::st_geometry(buildings), col = NA, border = "grey30", add = TRUE)
    legend_info <- draw_shared_colorbar_panel(
      cols,
      diff_rng[1L],
      diff_rng[2L],
      title = "W/m2",
      reserve_in = if (isTRUE(scalebar)) 0.5 else 0
    )
    draw_shared_map_ornaments(
      sf::st_geometry(buildings),
      scalebar = scalebar,
      scalebar_unit = scalebar_unit,
      scalebar_cex = scalebar_cex,
      north_arrow = north_arrow,
      legend_info = legend_info,
      map_scale = map_scale
    )
  } else {
    old_par <- open_shared_map_layout(sf::st_geometry(buildings), legend = TRUE)
    on.exit(close_shared_map_layout(old_par, legend = TRUE), add = TRUE)
    draw_shared_map_frame(
      sf::st_geometry(buildings),
      main = "Canopy impact - radiation difference (W/m2)"
    )
    graphics::plot(sf::st_geometry(buildings), col = "grey95", border = "grey45", add = TRUE)
    if (any(facade)) {
      graphics::plot(
        sf::st_geometry(impact[facade, ]),
        pch = 16, cex = 0.25,
        col = cols[make_idx(impact$total_diff[facade], diff_rng)],
        add = TRUE
      )
    }
    if (any(roof)) {
      graphics::plot(
        sf::st_geometry(impact[roof, ]),
        pch = 16, cex = 0.75,
        col = cols[make_idx(impact$total_diff[roof], diff_rng)],
        add = TRUE
      )
    }
    graphics::plot(sf::st_geometry(buildings), col = NA, border = "grey30", add = TRUE)
    map_scale <- if (isTRUE(scalebar)) noise_map_scale(sf::st_geometry(buildings)) else NULL
    legend_info <- draw_shared_colorbar_panel(
      cols,
      diff_rng[1L],
      diff_rng[2L],
      title = "W/m2",
      reserve_in = if (isTRUE(scalebar)) 0.5 else 0
    )
    draw_shared_map_ornaments(sf::st_geometry(buildings), scalebar, scalebar_unit,
                              scalebar_cex, north_arrow, legend_info, map_scale)
  }
  invisible(NULL)
}

plot_radiation_surface_3d <- function(buildings, radiation, height_field = "Height") {
  if (nrow(radiation) == 0L) {
    graphics::plot.new()
    graphics::mtext("No radiation samples", side = 1, line = 0.6, font = 2)
    return(invisible(NULL))
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1L, 3L), mar = c(2, 1, 1, 1), oma = c(0, 0, 0, 4))

  fields <- c("direct", "diffuse", "total")
  titles <- c("Direct radiation", "Diffuse radiation", "Total radiation")
  bld_h  <- as.numeric(buildings[[height_field]])

  # ---- Consistent 3D projector fixed to all footprint corners + full height ----
  all_xy <- do.call(rbind, lapply(seq_len(nrow(buildings)), function(b) {
    sf::st_coordinates(sf::st_geometry(buildings[b, ]))[, c("X", "Y"), drop = FALSE]
  }))
  proj_fn <- make_3d_projector(
    all_xy[, "X"], all_xy[, "Y"],
    c(0, max(bld_h, na.rm = TRUE))
  )

  # Plot limits from projected radiation samples
  sxy  <- sf::st_coordinates(radiation)[, c("X", "Y"), drop = FALSE]
  sprj <- proj_fn(sxy[, "X"], sxy[, "Y"], radiation$z)
  xlim <- range(sprj$x, na.rm = TRUE)
  ylim <- range(sprj$y, na.rm = TRUE)

  # Building order: farthest first (painter's algorithm).
  # `depth` grows with distance, so the draw order is DECREASING depth.
  # The key is evaluated at z = 0 so that building height cannot contaminate
  # the planar ordering - a tall building must not sort as if it were nearer.
  bld_ctr <- suppressWarnings(
    sf::st_coordinates(sf::st_centroid(sf::st_geometry(buildings)))
  )
  bld_ord <- order(
    proj_fn(bld_ctr[, "X"], bld_ctr[, "Y"], rep(0, nrow(bld_ctr)))$depth,
    decreasing = TRUE
  )

  # Pre-extract exterior ring + winding orientation per building
  bld_rings <- lapply(seq_len(nrow(buildings)), function(b) {
    rings <- polygon_rings(sf::st_geometry(buildings[b, ]))
    ext   <- rings[[1L]]   # exterior ring is always first from polygon_rings()
    n     <- nrow(ext)
    sa    <- sum(ext[seq_len(n - 1L), 1L] * ext[2L:n, 2L] -
                 ext[2L:n, 1L] * ext[seq_len(n - 1L), 2L]) / 2
    list(ring = ext, orient = if (sa >= 0) 1L else -1L)
  })

  # One colour scale shared by all three panels so the same colour means the
  # same W/m^2 everywhere and the panels are directly comparable.
  vrng <- range(unlist(lapply(fields, function(f) radiation[[f]])), na.rm = TRUE)

  for (fi in seq_along(fields)) {
    fld  <- fields[fi]
    vals <- radiation[[fld]]

    graphics::plot(xlim, ylim, type = "n", axes = FALSE,
                   xlab = "", ylab = "", asp = 1L, main = "")
    graphics::mtext(titles[fi], side = 1, line = 0.6, font = 2, cex = 1.1)

    for (b in bld_ord) {
      h <- bld_h[b]
      if (!is.finite(h) || h <= 0) next

      ring   <- bld_rings[[b]]$ring
      orient <- bld_rings[[b]]$orient
      n_edge <- nrow(ring) - 1L

      fac_idx <- which(radiation$building_id == b & radiation$surface == "facade")
      rof_idx <- which(radiation$building_id == b & radiation$surface == "roof")

      # Sort walls farthest-first, evaluated at ground level for the same
      # reason as the building ordering above.
      w_ord <- if (n_edge > 0L) {
        mx <- (ring[seq_len(n_edge), 1L] + ring[seq_len(n_edge) + 1L, 1L]) / 2
        my <- (ring[seq_len(n_edge), 2L] + ring[seq_len(n_edge) + 1L, 2L]) / 2
        order(proj_fn(mx, my, rep(0, n_edge))$depth, decreasing = TRUE)
      } else integer(0L)

      # ---- Draw facade walls as filled horizontal strips ----
      for (e in w_ord) {
        p1  <- ring[e, ]
        p2  <- ring[e + 1L, ]
        edv <- p2 - p1
        elen <- sqrt(sum(edv^2))
        if (!is.finite(elen) || elen < 1e-9) next

        # Outward normal (consistent with facade_grid_sf after orientation fix)
        nrm <- orient * c(edv[2L], -edv[1L]) / elen

        # Backface culling: skip walls pointing away from the camera, otherwise
        # their base edges show through as a spurious footprint outline.
        if (!facing_camera(nrm[1L], nrm[2L], proj_fn)) next

        eidx <- if (length(fac_idx) > 0L) {
          fac_idx[abs(radiation$nx[fac_idx] - nrm[1L]) < 1e-4 &
                  abs(radiation$ny[fac_idx] - nrm[2L]) < 1e-4]
        } else integer(0L)

        # Full-wall quad, reused for the silhouette outline
        wall_q <- proj_fn(c(p1[1L], p2[1L], p2[1L], p1[1L]),
                          c(p1[2L], p2[2L], p2[2L], p1[2L]),
                          c(0, 0, h, h))

        if (length(eidx) == 0L) {
          # No radiation samples for this wall - draw neutral grey
          graphics::polygon(wall_q$x, wall_q$y, col = "grey75",
                            border = "grey40", lwd = 0.4)
          next
        }

        ez  <- radiation$z[eidx]
        ev  <- vals[eidx]
        z_u <- sort(unique(ez))
        n_z <- length(z_u)
        # Strip boundaries that cover exactly 0..h with no gaps
        bnd <- if (n_z == 1L) {
          c(0, h)
        } else {
          c(0, (z_u[-n_z] + z_u[-1L]) / 2, h)
        }

        for (zi in seq_len(n_z)) {
          sv   <- mean(ev[abs(ez - z_u[zi]) < 1e-9], na.rm = TRUE)
          scol <- radiation_value_cols_range(sv, vrng)
          qp   <- proj_fn(c(p1[1L], p2[1L], p2[1L], p1[1L]),
                          c(p1[2L], p2[2L], p2[2L], p1[2L]),
                          c(bnd[zi], bnd[zi], bnd[zi + 1L], bnd[zi + 1L]))
          # Border matched to the fill keeps strips seamless (no hairline gaps)
          graphics::polygon(qp$x, qp$y, col = scol, border = scol, lwd = 0.3)
        }
        # Outline the whole wall face so pale walls stay legible on white
        graphics::polygon(wall_q$x, wall_q$y, col = NA,
                          border = "grey40", lwd = 0.4)
      }

      # ---- Draw roof as filled projected footprint polygon ----
      rcol <- if (length(rof_idx) > 0L) {
        radiation_value_cols_range(mean(vals[rof_idx], na.rm = TRUE), vrng)
      } else {
        "grey85"
      }
      rp <- proj_fn(ring[, 1L], ring[, 2L], rep(h, nrow(ring)))
      graphics::polygon(rp$x, rp$y, col = rcol,
                        border = "grey40", lwd = 0.4)
    }

    # Single shared colour bar, drawn once beside the last panel
    if (fi == length(fields)) {
      draw_colorbar(
        radiation_palette(100L),
        lo = vrng[1L], hi = vrng[2L],
        title = "W/m^2"
      )
    }
  }
  invisible(NULL)
}

# Like project_radiation_3d but with normalization parameters fixed at construction
# time so that multiple calls use the same coordinate system.
make_3d_projector <- function(x_ref, y_ref, z_ref, azimuth = 35, elevation = 25) {
  x_ctr   <- mean(range(as.numeric(x_ref), na.rm = TRUE))
  y_ctr   <- mean(range(as.numeric(y_ref), na.rm = TRUE))
  z_min   <- min(as.numeric(z_ref), na.rm = TRUE)
  x0      <- as.numeric(x_ref) - x_ctr
  y0      <- as.numeric(y_ref) - y_ctr
  z0      <- as.numeric(z_ref) - z_min
  xy_span <- max(diff(range(x0, na.rm = TRUE)), diff(range(y0, na.rm = TRUE)), 1)
  z_span  <- diff(range(z0, na.rm = TRUE))
  z_scale <- if (is.finite(z_span) && z_span > 0) xy_span / z_span * 0.18 else 1
  az <- azimuth * pi / 180; el <- elevation * pi / 180
  ca <- cos(az); sa <- sin(az); ce <- cos(el); se <- sin(el)
  fn <- function(x, y, z) {
    x0 <- as.numeric(x) - x_ctr
    y0 <- as.numeric(y) - y_ctr
    z0 <- (as.numeric(z) - z_min) * z_scale
    xr <- x0 * ca - y0 * sa
    yr <- x0 * sa + y0 * ca
    list(x = xr, y = yr * ce + z0 * se, depth = yr * se - z0 * ce)
  }
  # Expose rotation terms so callers can perform backface culling: a wall with
  # outward normal (nx, ny) faces the camera when its rotated y-component is
  # negative, i.e. moving along the normal decreases projected depth.
  attr(fn, "sa") <- sa
  attr(fn, "ca") <- ca
  fn
}

# TRUE when a wall with outward normal (nx, ny) is visible to the camera.
facing_camera <- function(nx, ny, proj_fn) {
  sa <- attr(proj_fn, "sa"); ca <- attr(proj_fn, "ca")
  if (is.null(sa) || is.null(ca)) return(TRUE)
  (nx * sa + ny * ca) < 0
}

project_radiation_3d <- function(x, y, z, azimuth = 35, elevation = 25) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  z <- as.numeric(z)
  x0 <- x - mean(range(x, na.rm = TRUE))
  y0 <- y - mean(range(y, na.rm = TRUE))
  z0 <- z - min(z, na.rm = TRUE)
  xy_span <- max(diff(range(x0, na.rm = TRUE)), diff(range(y0, na.rm = TRUE)), 1)
  z_span <- diff(range(z0, na.rm = TRUE))
  if (is.finite(z_span) && z_span > 0) {
    z0 <- z0 * xy_span / z_span * 0.18
  }
  az <- deg2rad(azimuth)
  el <- deg2rad(elevation)
  xr <- x0 * cos(az) - y0 * sin(az)
  yr <- x0 * sin(az) + y0 * cos(az)
  list(
    x = xr,
    y = yr * cos(el) + z0 * sin(el),
    depth = yr * sin(el) - z0 * cos(el)
  )
}

# Single source of truth for the radiation colour ramp so that map surfaces and
# colour bars can never disagree.  `rev = TRUE` puts pale yellow at the LOW end
# and saturated dark red at the HIGH end; without it the ramp runs the other way
# and high-radiation surfaces render near-white, which reads as washed out.
radiation_palette <- function(n = 100L) {
  # Drop the palest ~20% of the ramp.  The raw YlOrRd low end is almost white,
  # which makes low-radiation surfaces blend into the page background and read
  # as "transparent".  Starting at a visible cream keeps every surface opaque.
  full <- grDevices::hcl.colors(round(n * 1.25), "YlOrRd", rev = TRUE)
  full[(length(full) - n + 1L):length(full)]
}

radiation_value_cols <- function(values) {
  values <- as.numeric(values)
  if (length(values) == 0) {
    return(character())
  }
  radiation_value_cols_range(values, range(values, na.rm = TRUE))
}

# Like radiation_value_cols but with a pre-supplied range so that the colour
# scale is consistent across multiple calls (e.g. per-building rendering).
radiation_value_cols_range <- function(values, value_range) {
  cols   <- radiation_palette(100L)
  values <- as.numeric(values)
  if (length(values) == 0L) return(character(0L))
  if (all(is.na(values))) return(rep(cols[1L], length(values)))
  if (!is.finite(diff(value_range)) || diff(value_range) == 0) {
    return(rep(cols[50L], length(values)))
  }
  idx <- round((values - value_range[1L]) / diff(value_range) * 99) + 1L
  idx[is.na(idx)] <- 1L
  cols[pmin(pmax(idx, 1L), 100L)]
}

sun_vector <- function(solar_pos) {
  az <- deg2rad(solar_pos[[1]])
  elev <- deg2rad(solar_pos[[2]])
  c(sin(az) * cos(elev), cos(az) * cos(elev), sin(elev))
}

deg2rad <- function(x) x * pi / 180
rad2deg <- function(x) x * 180 / pi

#### Noise mapping ####
#' @noRd
prepare_noise_buildings <- function(buildings,
                                    height_col,
                                    population = FALSE,
                                    population_field = NULL,
                                    population_year = 2025,
                                    quiet = TRUE) {
  out <- buildings
  if (isTRUE(population) && is.null(population_field) && !"POP" %in% names(out) && !"pop_total" %in% names(out)) {
    population_fun <- get0("get_pop", mode = "function", inherits = TRUE)
    if (is.null(population_fun)) {
      population_fun <- get0("get_pop_density", mode = "function", inherits = TRUE)
    }
    if (is.null(population_fun)) {
      stop("No building population function is available.", call. = FALSE)
    }
    out <- population_fun(out, year = population_year, quiet = quiet)
  }
  out$PK <- seq_len(nrow(out))
  out$HEIGHT <- as.numeric(out[[height_col]])
  pop_field <- resolve_noise_population_field(out, population_field)
  if (!is.null(pop_field)) {
    out$POP <- as.numeric(out[[pop_field]])
    out$POP[is.na(out$POP) | out$POP < 0] <- 0
    out$POP[!is.finite(out$POP)] <- 0
  }
  out <- out[!is.na(out$HEIGHT) & out$HEIGHT > 0, ]
  if (nrow(out) == 0) {
    stop("No buildings with positive finite height were found.", call. = FALSE)
  }
  keep <- unique(c("PK", "id", "group_id", "HEIGHT", "POP", attr(out, "sf_column")))
  out[, intersect(keep, names(out))]
}

#' @noRd
resolve_noise_population_field <- function(buildings, population_field = NULL) {
  if (!is.null(population_field)) {
    if (!population_field %in% names(buildings)) {
      stop("`population_field` was not found in `x`: ", population_field, call. = FALSE)
    }
    return(population_field)
  }
  candidates <- c("POP", "pop_total", "population", "pop")
  found <- candidates[candidates %in% names(buildings)]
  if (length(found) == 0) return(NULL)
  found[1]
}

#' @noRd
prepare_noise_roads <- function(roads, quiet = TRUE) {
  required_traffic <- c("speed_kmh", "light_veh_h", "heavy_veh_h")
  if (!all(required_traffic %in% names(roads))) {
    roads <- infer_osm_traffic(roads, quiet = quiet)
  }
  roads$PK <- seq_len(nrow(roads))
  keep <- unique(c(
    "PK", "highway", "highway_class", "name", "maxspeed", "lanes", "speed_kmh",
    "light_veh_h", "heavy_veh_h", "traffic_source",
    "LV_D", "LV_E", "LV_N", "HGV_D", "HGV_E", "HGV_N",
    "MV_D", "MV_E", "MV_N", "WAV_D", "WAV_E", "WAV_N", "WBV_D", "WBV_E", "WBV_N",
    "LV_SPD_D", "LV_SPD_E", "LV_SPD_N",
    "HGV_SPD_D", "HGV_SPD_E", "HGV_SPD_N",
    "MV_SPD_D", "MV_SPD_E", "MV_SPD_N",
    "WAV_SPD_D", "WAV_SPD_E", "WAV_SPD_N",
    "WBV_SPD_D", "WBV_SPD_E", "WBV_SPD_N",
    "PVMT", "TEMP_D", "TEMP_E", "TEMP_N", "WAY",
    attr(roads, "sf_column")
  ))
  roads[, intersect(keep, names(roads))]
}

#' @noRd
prepare_noise_ground <- function(buildings,
                                 greenspace = NULL,
                                 canopy_height = NULL,
                                 ground_default = 0,
                                 green_ground = 1,
                                 min_tree_height = 2) {
  bbox_poly <- sf::st_as_sfc(sf::st_bbox(buildings))
  bbox_sf <- sf::st_sf(G = ground_default, source = "default", geometry = bbox_poly)
  sf::st_crs(bbox_sf) <- sf::st_crs(buildings)
  green_layers <- list()

  green_sf <- green_input_to_ground(greenspace, buildings, green_ground, "greenspace")
  if (!is.null(green_sf)) green_layers[[length(green_layers) + 1L]] <- green_sf

  canopy_sf <- canopy_to_ground(canopy_height, buildings, min_tree_height, green_ground)
  if (!is.null(canopy_sf)) green_layers[[length(green_layers) + 1L]] <- canopy_sf

  if (length(green_layers) == 0) {
    return(bbox_sf)
  }
  green <- do.call(rbind, green_layers)
  green <- suppressWarnings(sf::st_intersection(green, bbox_sf[, "geometry"]))
  if (nrow(green) == 0) {
    return(bbox_sf)
  }
  green <- green[, c("G", "source", attr(green, "sf_column"))]
  rbind(bbox_sf, green)
}

#' @noRd
make_noise_analysis_extent <- function(buildings, roads) {
  bbox_values <- c(
    xmin = min(sf::st_bbox(buildings)[["xmin"]], sf::st_bbox(roads)[["xmin"]]),
    ymin = min(sf::st_bbox(buildings)[["ymin"]], sf::st_bbox(roads)[["ymin"]]),
    xmax = max(sf::st_bbox(buildings)[["xmax"]], sf::st_bbox(roads)[["xmax"]]),
    ymax = max(sf::st_bbox(buildings)[["ymax"]], sf::st_bbox(roads)[["ymax"]])
  )
  bbox <- structure(bbox_values, class = "bbox", crs = sf::st_crs(buildings))
  sf::st_sf(geometry = sf::st_as_sfc(bbox))
}

#' @noRd
green_input_to_ground <- function(greenspace, template, green_ground, source) {
  if (is.null(greenspace)) return(NULL)
  if (inherits(greenspace, "sf")) {
    out <- sf::st_transform(greenspace, sf::st_crs(template))
    out$G <- green_ground
    out$source <- source
    return(out[, c("G", "source", attr(out, "sf_column"))])
  }
  if (inherits(greenspace, "SpatRaster")) {
    r <- terra::project(greenspace[[1]], sf::st_crs(template)$wkt, method = "near")
    r <- terra::ifel(!is.na(r) & r > 0, 1, NA)
    return(raster_mask_to_ground(r, green_ground, source))
  }
  stop("`greenspace` must be an sf object or terra SpatRaster when supplied.", call. = FALSE)
}

#' @noRd
canopy_to_ground <- function(canopy_height, template, min_tree_height, green_ground) {
  if (is.null(canopy_height)) return(NULL)
  if (!inherits(canopy_height, "SpatRaster")) {
    stop("`canopy_height` must be a terra SpatRaster when supplied.", call. = FALSE)
  }
  r <- terra::project(canopy_height[[1]], sf::st_crs(template)$wkt, method = "near")
  r <- terra::ifel(!is.na(r) & r >= min_tree_height, 1, NA)
  raster_mask_to_ground(r, green_ground, "canopy")
}

#' @noRd
raster_mask_to_ground <- function(r, green_ground, source) {
  if (all(is.na(terra::values(r, mat = FALSE)))) return(NULL)
  poly <- suppressWarnings(sf::st_as_sf(terra::as.polygons(r, dissolve = TRUE, na.rm = TRUE)))
  if (nrow(poly) == 0) return(NULL)
  poly$G <- green_ground
  poly$source <- source
  poly[, c("G", "source", attr(poly, "sf_column"))]
}

#' @noRd
make_noise_receiver_grid <- function(buildings, roads, resolution) {
  bbox <- sf::st_bbox(buildings)
  if ((bbox[["xmax"]] - bbox[["xmin"]]) < resolution) {
    pad <- (resolution - (bbox[["xmax"]] - bbox[["xmin"]])) / 2
    bbox[["xmin"]] <- bbox[["xmin"]] - pad
    bbox[["xmax"]] <- bbox[["xmax"]] + pad
  }
  if ((bbox[["ymax"]] - bbox[["ymin"]]) < resolution) {
    pad <- (resolution - (bbox[["ymax"]] - bbox[["ymin"]])) / 2
    bbox[["ymin"]] <- bbox[["ymin"]] - pad
    bbox[["ymax"]] <- bbox[["ymax"]] + pad
  }
  bbox_values <- c(
    xmin = bbox[["xmin"]],
    ymin = bbox[["ymin"]],
    xmax = bbox[["xmax"]],
    ymax = bbox[["ymax"]]
  )
  combined_bbox <- structure(bbox_values, class = "bbox", crs = sf::st_crs(buildings))
  grid <- sf::st_make_grid(sf::st_as_sfc(combined_bbox), cellsize = resolution, what = "centers")
  coords <- sf::st_coordinates(grid)
  grid_z <- sf::st_sfc(
    lapply(seq_len(nrow(coords)), function(i) sf::st_point(c(coords[i, "X"], coords[i, "Y"], 4))),
    crs = sf::st_crs(buildings)
  )
  receivers <- sf::st_sf(PK = seq_along(grid_z), id = seq_along(grid_z), HEIGHT = 4, geometry = grid_z)
  outside <- receivers[lengths(sf::st_intersects(receivers, buildings)) == 0, ]
  if (nrow(outside) > 0) {
    return(outside)
  }
  centroid <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_union(sf::st_geometry(buildings)))))
  fallback_geom <- sf::st_sfc(sf::st_point(c(centroid[1, "X"], centroid[1, "Y"], 4)), crs = sf::st_crs(buildings))
  sf::st_sf(PK = 1L, id = 1L, HEIGHT = 4, geometry = fallback_geom)
}

#' @noRd
write_noise_gpkg <- function(gpkg, buildings, roads, ground, receivers = NULL, quiet = TRUE) {
  if (file.exists(gpkg)) unlink(gpkg)
  sf::st_write(buildings, gpkg, layer = "BUILDINGS", quiet = quiet)
  if (!is.null(roads)) {
    sf::st_write(roads, gpkg, layer = "ROADS", quiet = quiet, append = TRUE)
  }
  sf::st_write(ground, gpkg, layer = "GROUND", quiet = quiet, append = TRUE)
  if (!is.null(receivers)) {
    sf::st_write(receivers, gpkg, layer = "RECEIVERS", quiet = quiet, append = TRUE)
  }
  invisible(gpkg)
}

#' @noRd
normalize_osm_highway <- function(highway) {
  vapply(highway, function(x) {
    if (length(x) == 0 || all(is.na(x))) return(NA_character_)
    x <- as.character(x[[1]])
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    strsplit(x, ";|,")[[1]][1]
  }, character(1))
}

#' @noRd
parse_osm_speed_kmh <- function(maxspeed) {
  speed <- parse_osm_number(maxspeed)
  text <- tolower(vapply(maxspeed, function(x) {
    if (length(x) == 0 || all(is.na(x))) return("")
    as.character(x[[1]])
  }, character(1)))
  speed[grepl("mph", text) & !is.na(speed)] <- speed[grepl("mph", text) & !is.na(speed)] * 1.609344
  speed
}

#' @noRd
parse_osm_number <- function(x) {
  vapply(x, function(value) {
    if (length(value) == 0 || all(is.na(value))) return(NA_real_)
    text <- as.character(value[[1]])
    match <- regmatches(text, regexpr("[0-9]+(\\.[0-9]+)?", text))
    if (length(match) == 0 || !nzchar(match)) return(NA_real_)
    as.numeric(match)
  }, numeric(1))
}

#' @noRd
parse_osm_oneway <- function(x) {
  text <- tolower(vapply(x, function(value) {
    if (length(value) == 0 || all(is.na(value))) return("")
    as.character(value[[1]])
  }, character(1)))
  text %in% c("yes", "true", "1", "-1")
}
