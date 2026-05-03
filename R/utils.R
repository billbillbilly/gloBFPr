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
#' @importFrom viewscape compute_viewshed calculate_feature

#### Footprint data processing ####
#' @noMd
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

#' @noMd
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

#' @noMd
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

#' @noMd
cal_convexity <- function(poly_) {
  vertices <- sf::st_coordinates(poly_)[, 1:2]
  multipoint <- sf::st_multipoint(vertices)
  multipoint <- sf::st_sfc(multipoint, crs = sf::st_crs(poly_))
  convex_hull <- sf::st_convex_hull(multipoint)
  return(as.numeric(poly_$g_area / sf::st_area(convex_hull)))
}

#' @noMd
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

#' @noMd
cal_mean_pairwise_distance <- function(poly_) {
  # create a 3D voxel grid
  filtered_grid <- ploy2grid(poly_)
 # Compute all pairwise distances
  coords <- as.matrix(filtered_grid)
  dist <- mean_pairwise_distance(coords)
  return(dist)
}

#' @noMd
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

#' @noMd
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

#' @noMd
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

#' @noMd
get_gvi <- function(dsm, p, height, r, building, binary_chm) {
  tryCatch({
    # Viewshed
    v <- viewscape::compute_viewshed(dsm = dsm, viewpoints = p,
                                     offset_viewpoint = height,
                                     r = r
    )
    # Viewshed area
    v_area <- length(as.vector(v@visible[v@visible == 1])) * v@resolution[1]^2
    # Visible canopy area
    canopy_proportion <- viewscape::calculate_feature(viewshed = v,
                                                      feature = binary_chm,
                                                      type = 2,
                                                      exclude_value = 0)
    canopy_area <- v_area * canopy_proportion

    # Compute GVI (allow >1 if canopy exceeds viewshed minus building)
    # print(paste0("building area: ", as.numeric(building$g_area),
    #              "; viewshed area: ", v_area,
    #              "; visible green area: ", canopy_area,
    #              "; canopy proportion: ", canopy_proportion
    #              )
    #       )
    gvi <- canopy_area / max(v_area - building$g_area, 1e-6)  # use small constant to avoid division by 0
    gvi <- min(gvi, 1)  # cap at 1 if desired

    return(gvi)
  }, error = function(e) {
    stop(sprintf("GVI calculation failed: %s", e$message))
  })
}
# Note:
# In fact, if the building itself occupies nearly the entire viewshed,
# or the canopy overlaps very closely, v_area - building$g_area could be ≤ 0,
# yet still surrounded by trees — in which case GVI = 1 is conceptually valid.

# max(v_area - building$g_area, 1e-6) prevents division by zero while allowing close overlap
# min(gvi, 1) caps the value at 1 if needed for interpretation as a proportion
# Still returns 0 if the viewshed fails or has no values

#### Data collection and processing ####
#' @noMd
get_GHSpop <- function(bbox = NULL, year = NULL) {
  # Store the original 'timeout' option and ensure it's reset upon function exit
  original_timeout <- getOption('timeout')
  on.exit(options(timeout = original_timeout), add = TRUE)
  options(timeout=9999)

  d_mode <- 'auto'
  # check os
  os <- Sys.info()[["sysname"]]
  d_mode <- if (Sys.info()[["sysname"]] == "Windows") 'wb' else 'auto'

  # GHS population grid
  years <- c(2030, 2025, 2020, 2015, 2010, 2005, 2000, 1995, 1990, 1985, 1980, 1975)
  result_list <- list()
  temp_paths <- c()  # store paths for later cleanup

  if (year %in% years) {
    intersected_tiles <- ghsl_tiles[sf::st_intersects(ghsl_tiles, bbox, sparse = FALSE), ]
    for (i in seq_len(nrow(intersected_tiles))) {
      temp_zip <- tempfile(fileext = ".zip")
      url_ <- get_GHSurl(year, intersected_tiles$tile_id[i], 'pop')
      utils::download.file(url_,
                           destfile = temp_zip,
                           mode = d_mode,
                           quiet = TRUE)
      unzip_dir <- tempfile()
      utils::unzip(temp_zip, exdir = unzip_dir)
      tif_files <- list.files(unzip_dir, pattern = "\\.tif$", full.names = TRUE)
      if (length(tif_files) == 0) next
      rast_data <- terra::rast(tif_files[1])
      result_list[[length(result_list) + 1]] <- rast_data
      temp_paths <- c(temp_paths, temp_zip, unzip_dir)
      # unlink(c(temp_zip, unzip_dir), recursive = TRUE)
    }
    if (length(result_list) == 0) {
      stop("No population rasters downloaded")
    }
    cli::cli_alert_success('Finished downloading population data')

    # Combine all into one terra raster object
    r <- if (length(result_list) == 1) result_list[[1]] else do.call(terra::merge, result_list)

    # reproject raster
    utm_crs <- get_utm_crs(bbox)
    r <- terra::project(r, paste0('EPSG:', utm_crs), method = 'near')

    # calculate population density
    r <- r / (100*100)

    on.exit(unlink(temp_paths, recursive = TRUE), add = TRUE)
    return(r)
  } else {
    stop(sprintf("Input year %d is not in allowed range. Skipping.", year))
  }
}

#' @noMd
get_GHSres <- function(bbox = NULL, year = NULL) {
  # Store the original 'timeout' option and ensure it's reset upon function exit
  original_timeout <- getOption('timeout')
  on.exit(options(timeout = original_timeout), add = TRUE)
  options(timeout=9999)

  d_mode <- 'auto'
  # check os
  os <- Sys.info()[["sysname"]]
  d_mode <- if (Sys.info()[["sysname"]] == "Windows") 'wb' else 'auto'

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
      utils::download.file(urls[[1]],
                           destfile = temp_total_zip,
                           mode = d_mode,
                           quiet = TRUE)
      utils::download.file(urls[[2]],
                           destfile = temp_nres_zip,
                           mode = d_mode,
                           quiet = TRUE)
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

#' @noMd
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



#' @noMd
get_dem <- function(bbox, key) {
  if (missing(key)) stop("missing api key")

  us_poly <- suppressMessages(tigris::nation(progress_bar = FALSE))
  bbox_poly <- sf::st_transform(
    sf::st_as_sfc(
      sf::st_bbox(
        c(xmin = bbox[1],
          ymin = bbox[2],
          xmax = bbox[3],
          ymax = bbox[4]),
        crs = 4326
      )
    ), sf::st_crs(us_poly)
  )

  dem <- NULL

  if (any(sf::st_intersects(bbox_poly, us_poly, sparse = FALSE))) {
    # USGS sources
    dem <- tryCatch({
      dsmSearch::get_dsm_30(bbox = bbox, key = key, datatype = 'usgs1m')
    }, error = function(e1) {
      tryCatch({
        dsmSearch::get_dsm_30(bbox = bbox, key = key, datatype = 'usgs10m')
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


#' @noMd
get_chm <- function(bbox, min_height) {
  # get CHM
  chm <- suppressMessages(dsmSearch::get_dsm_30(bbox = bbox, datatype = 'metaCHM'))
  # reporject chm
  bbox <- sf::st_as_sfc(
    sf::st_bbox(
      c(xmin = bbox[1],
        ymin = bbox[2],
        xmax = bbox[3],
        ymax = bbox[4]),
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

#' @noMd
get_greenspace <- function(bbox = NULL, buffer = NULL,
                           type = NULL, zoom = 17, year = NULL,
                           min_tree_height = 2) {
  if (inherits(type, "NULL")) {
    stop("Please input greenspace datasource.")
  }
  type <- match.arg(type, c("metachm", "esri", "sentinel2"))
  bbox_vector <- if (is.numeric(bbox) && length(bbox) == 4) {
    bbox
  } else {
    bbox_poly_to_list(bbox)
  }

  if (type == "metachm") {
    g <- get_chm(bbox_vector, min_tree_height)[[2]]
  } else if (type == "esri") {
    g <- greenSD::get_tile_green(bbox = bbox_vector, zoom = zoom,
                                 provider = "esri")
    utm_crs <- get_utm_crs(bbox)
    g <- terra::project(g$green, paste0('EPSG:', utm_crs), method = 'near')
  } else if (type == "sentinel2") {
    g <- greenSD::get_tile_green(bbox = bbox_vector, zoom = zoom,
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

#' @noMd
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

#' @noMd
merge_elev <- function(building, dem, chm=NULL) {
  # prioritize layers:  (chm >) building > dem
  bc <- terra::overlay(r1, building, fun = function(x, y) {
    ifelse(x < y & x != 0, x, x + y)
  })
}

#' @importFrom sf st_buffer st_centroid
#' @noMd
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
#' @noMd
get_utm_crs <- function(bbox) {
  if (is.numeric(bbox) && length(bbox) == 4) {
    bbox <- sf::st_as_sfc(
      sf::st_bbox(
        c(xmin = bbox[1],
          ymin = bbox[2],
          xmax = bbox[3],
          ymax = bbox[4]),
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

#' @noMd
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

#' @noMd
get_bbox <- function(x) {
  bbox <- sf::st_as_sfc(sf::st_bbox(x), crs = sf::st_crs(x))
  bbox <- sf::st_transform(bbox, crs = 4326)
  return(bbox)
}

#' @noMd
bbox_poly_to_list <- function(bbox) {
  coor <- sf::st_coordinates(bbox)
  return(
    c(min(coor[,1]), min(coor[,2]), max(coor[,1]), max(coor[,2]))
  )
}

#' @noMd
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

#' @noMd
compute_gvi_per_building <- function(building,
                                     dem_path,
                                     chm_path,
                                     binary_chm_path,
                                     bh_all_path,
                                     radius,
                                     floor,
                                     floor_step) {
  # Read raster files inside each worker.
  dem <- terra::rast(dem_path)
  chm <- terra::rast(chm_path)
  binary_chm <- terra::rast(binary_chm_path)
  bh_all <- terra::rast(bh_all_path)

  # Compute centroid
  centroid <- suppressWarnings(sf::st_centroid(building$geometry))
  p <- as.vector(sf::st_coordinates(centroid))

  # Crop all rasters to the local viewshed radius.
  buffer <- sf::st_buffer(centroid, dist = radius)
  dem <- terra::crop(dem, terra::vect(buffer), mask = TRUE)
  chm <- terra::crop(chm, terra::vect(buffer), mask = TRUE)
  binary_chm <- terra::crop(binary_chm, terra::vect(buffer), mask = TRUE)
  bh_all <- terra::crop(bh_all, terra::vect(buffer), mask = TRUE)

  # Flatten the target building footprint: it is the observer, not an obstacle.
  target_mask <- terra::rasterize(terra::vect(building), bh_all, field = 1, background = 0)
  bh_without_target <- terra::ifel(target_mask == 1, 0, bh_all)
  chm_without_target <- terra::ifel(target_mask == 1, 0, chm)
  binary_chm <- terra::ifel(target_mask == 1, 0, binary_chm)

  surface <- terra::ifel(chm_without_target > bh_without_target,
                         chm_without_target,
                         bh_without_target)
  dsm_ <- dem + surface

  if (isTRUE(floor)) {
    floor_ids <- unique(c(seq.int(1L, building$estimated_floors, by = floor_step),
                          building$estimated_floors))
    GVIs <- numeric(length(floor_ids))
    for (j in seq_along(floor_ids)) {
      height <- 1.7 + (floor_ids[j] - 1) * 3
      GVIs[j] <- get_gvi(dsm_, p, height, radius, building, binary_chm)
    }
    return(list(
      mean = mean(GVIs),
      min = min(GVIs),
      max = max(GVIs),
      sd = if (length(GVIs) > 1) stats::sd(GVIs) else NA_real_
    ))
  } else {
    height_bottom <- 1.7
    if (building$Height < 6) {
      mean_gvi <- get_gvi(dsm_, p, height_bottom, radius, building, binary_chm)
    } else {
      height_top <- 1.7 + building$Height - 3
      gvi_top <- get_gvi(dsm_, p, height_top, radius, building, binary_chm)
      gvi_bottom <- get_gvi(dsm_, p, height_bottom, radius, building, binary_chm)
      mean_gvi <- mean(c(gvi_top, gvi_bottom))
    }
    return(list(mean = mean_gvi, min = NA_real_, max = NA_real_, sd = NA_real_))
  }
}

#### utils ####
time_taken <- function(process_time) {
  if (process_time >= 60) {
    cli::cli_alert_success(paste0("Completed. Time taken: ", base::round(process_time/60), " minutes."))
  } else {
    cli::cli_alert_success(paste0("Completed. Time taken: ", base::round(process_time), " seconds."))
  }
}
