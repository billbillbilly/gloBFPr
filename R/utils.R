#' @importFrom terra ext
#' @importFrom terra vect
#' @importFrom terra rast
#' @importFrom terra rasterize
#' @importFrom terra align
#' @importFrom terra mask
#' @importFrom sf st_centroid
#' @importFrom sf st_coordinates
#' @importFrom sf st_transform
#' @importFrom sf st_crs st_convex_hull
#' @importFrom sf st_as_sfc st_multipoint

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
add_mph <- function(projected_poly, type) {
  #### Stats ####
  projected_poly$ground_vertex_count <- nrow(sf::st_coordinates(projected_poly)) - 1

  #### Properties ####
  projected_poly$ground_area <- as.numeric(sf::st_area(projected_poly))
  projected_poly$perimeter <- as.numeric(sf::st_perimeter(projected_poly))
  projected_poly$vertical_surface <- projected_poly$perimeter * projected_poly$Height
  projected_poly$total_surface <- projected_poly$vertical_surface + projected_poly$ground_area
  projected_poly$actual_volume <- projected_poly$ground_area * projected_poly$Height
  projected_poly$oobb_volume <- sf::st_area(sf::st_minimum_rotated_rectangle(projected_poly)) * projected_poly$Height

  #### Indices ####
  projected_poly$squareness <- projected_poly$perimeter / projected_poly$ground_area # Perimeter-Area Ratio
  projected_poly$rectangularity <- projected_poly$ground_area / sf::st_area(sf::st_minimum_rotated_rectangle(projected_poly))
  projected_poly$fractality <- 1 - log(projected_poly) / (1.5 * log(projected_poly$total_surface))
  projected_poly$hemisphericality <- 0
  projected_poly$convexity <- 0
  projected_poly$cuboidness <- 0
  projected_poly$med <- 0
  projected_poly$mpd <- 0
  projected_poly$vol_exch <- 0
  projected_poly$elo_short <- 0
  projected_poly$elo_long <- 0
  for (i in 1:nrow(projected_poly)) {
    poly_ <- projected_poly[i,]
    projected_poly$hemisphericality[i] <- cal_hemisphericality(poly_)
    projected_poly$convexity[i] <- cal_convexity(poly_)
    projected_poly$cuboidness[i] <- cal_cuboidness(poly_)
    projected_poly$med[i] <- cal_accessibility(poly_)
    projected_poly$mpd[i] <- cal_mean_pairwise_distance(poly_)
    projected_poly$vol_exch[i] <- cal_volume_exchange_ratio(poly_)
    elongation_ratios <- cal_elongation_ratios(poly_)
    projected_poly$elo_short[i] <- elongation_ratios[[1]]
    projected_poly$elo_long[i] <- elongation_ratios[[2]]
  }
  return(projected_poly)
}

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

  vol_hemisphere <- (2/3) * pi^(-0.5) * poly_$ground_area^(3/2)
  hemisphericality <- (poly_$actual_volume - vol_hemisphere) / vol_hemisphere
  return(hemisphericality)
}

#' @noMd
cal_convexity <- function(poly_) {
  vertices <- sf::st_coordinates(poly_)[, 1:2]
  multipoint <- sf::st_multipoint(vertices)
  multipoint <- sf::st_sfc(multipoint, crs = sf::st_crs(poly_))
  convex_hull <- sf::st_convex_hull(multipoint)
  return(poly_$ground_area / sf::st_area(convex_hull))
}

#' @noMd
cal_mean_accessibility <- function(poly_) {
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
  n <- nrow(coords)

  total_dist <- 0
  count <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- sqrt(sum((coords[i, ] - coords[j, ])^2))
      total_dist <- total_dist + d
      count <- count + 1
    }
  }
  return(total_dist / count)
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
  vol_exch <- (vol_sphere - poly_$actual_volume) / vol_sphere
  return(vol_exch)
}

#' @noMd
ploy2grid <- function(poly_) {
  height <- poly_$Height
  bbox_ <- sf::st_bbox(poly_)
  dx <- dy <- 2  # horizontal spacing
  dz <- 2        # vertical spacing
  x_seq <- seq(bbox["xmin"], bbox["xmax"], by = dx)
  y_seq <- seq(bbox["ymin"], bbox["ymax"], by = dy)
  z_seq <- seq(0, height, by = dz)
  grid_3d <- expand.grid(X = x_seq, Y = y_seq, Z = z_seq)
  # keep only points that fall inside the footprint (x, y only)
  xy_points <- sf::st_as_sf(grid_3d, coords = c("X", "Y"), crs = st_crs(poly))
  inside <- sf::st_within(xy_points, poly_, sparse = FALSE)[, 1]
  filtered_grid <- grid_3d[inside, ]
  return(filtered_grid)
}

"  double mean_pairwise_distance(Rcpp::NumericMatrix coords) {
    int n = coords.nrow();
    int dim = coords.ncol();
    double total_dist = 0.0;
    long count = 0;

    for (int i = 0; i < n - 1; ++i) {
      for (int j = i + 1; j < n; ++j) {
        double d = 0.0;
        for (int k = 0; k < dim; ++k) {
          double diff = coords(i, k) - coords(j, k);
          d += diff * diff;
        }
        total_dist += std::sqrt(d);
        count += 1;
      }
    }

    return total_dist / count;
  }
"

cal_elongation_ratios <- function(poly_) {
  # Compute minimum bounding rectangle
  min_rect <- sf::st_minimum_rotated_rectangle(poly_)
  coords <- sf::st_coordinates(min_rect)[, 1:2]

  # Get edge lengths
  edges <- sqrt(base::rowSums((coords - coords[c(2:5, 1), ])^2))
  lengths <- sort(unique(round(edges, 8)))  # may have duplicate lengths due to rectangle

  # Assign shortest and longest edge
  S <- lengths[1]
  L <- lengths[length(lengths)]

  # Use height attribute
  H <- poly_$Height

  # Compute elongation ratios
  ratio_short <- S / H
  ratio_long <- L / H

  return(list(ratio_short, ratio_long))
}

add_neighbor  <- function(projected_poly) {
  centroids <- sf::st_centroid(projected_poly)
}

#' @noMd
add_pop <- function(bbox, projected_poly, year) {
  bbox <- sf::st_transform(bbox, crs = 4326)
  d_mode <- 'auto'
  # check os
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    d_mode <- 'wb'
  }
  # GHS population grid
  years <- c(2030, 2025, 2020, 2015, 2010, 2005, 2000, 1995, 1990, 1985, 1980, 1975)
  result_list <- list()
  if (year %in% years) {
    intersected_tiles <- ghsl_tiles[sf::st_intersects(ghsl_tiles, bbox, sparse = FALSE), ]
    for (i in seq_len(nrow(intersected_tiles))) {
      temp_zip <- tempfile(fileext = ".zip")
      utils::download.file(intersected_tiles$tile_id[i],
                           destfile = temp_zip,
                           mode = d_mode,
                           quiet = TRUE)
      unzip_dir <- tempfile()
      utils::unzip(temp_zip, exdir = unzip_dir)
      tif_files <- list.files(unzip_dir, pattern = "\\.tif$", full.names = TRUE)
      if (length(tif_files) == 0) next
      rast_data <- terra::rast(tif_files[1])
      result_list[[length(result_list) + 1]] <- rast_data
      unlink(c(temp_zip, unzip_dir), recursive = TRUE)
    }

    if (length(result_list) == 0) {
      base::warning("No population rasters downloaded. Returning original polygons.")
      return(projected_poly)
    }

    # Combine all into one terra raster object
    r <- if (length(result_list) == 1) result_list[[1]] else do.call(terra::merge, result_list)

    # reproject raster
    utm_crs <- get_utm_crs(bbox)
    r <- terra::project(r, paste0('EPSG:', utm_crs))

    # calculate population density
    r <- r / (100*100)

    # use all centroids of projected_poly to extract value from the raster
    pop_density_df <- terra::extract(r,
                                     terra::vect(sf::st_centroid(projected_poly)),
                                     raw = FALSE,
                                     ID = FALSE)
    # link pop_density_df back to projected_poly
    projected_poly$pop_density <- pop_density_df[[1]]

  } else {
    base::warning(sprintf("Input year %d is not in allowed range. Skipping.", year))
    return(projected_poly)
  }
  return(projected_poly)
}

#' @noMd
get_GHSurl <- function(year, id) {
  # source: https://human-settlement.emergency.copernicus.eu/download.php?ds=pop
  return(
    paste0(
      'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2030_GLOBE_R2023A_54009_100/V1-0/tiles/GHS_POP_E',
      year,
      '_GLOBE_R2023A_54009_100_V1_0_',
      id,
      '.zip'
    )
  )
}

#### Projection tools ####
#' @noMd
get_utm_crs <- function(bbox) {
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
