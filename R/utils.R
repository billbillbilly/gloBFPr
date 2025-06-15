#' @importFrom terra ext
#' @importFrom terra vect
#' @importFrom terra rast
#' @importFrom terra rasterize
#' @importFrom terra align
#' @importFrom terra mask
#' @importFrom sf st_centroid st_within
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

#' @noMd
get_dem <- function(bbox, key) {
  us_poly <- suppressMessages(tigris::nation(progress_bar = FALSE))
  bbox_poly <- sf::st_transform(bbox, sf::st_crs(us_poly))

  dem <- NULL

  if (any(sf::st_intersects(bbox_poly, us_poly, sparse = FALSE))) {
    # USGS sources
    dem <- tryCatch({
      dsmSearch::get_dsm_30(bbox, key, datatype = 'usgs1m')
    }, error = function(e1) {
      tryCatch({
        dsmSearch::get_dsm_30(bbox, key, datatype = 'usgs10m')
      }, error = function(e2) {
        dsmSearch::get_dsm_30(bbox, key, datatype = 'SRTMGL1')
      })
    })
  } else {
    # Global SRTM fallback
    dem <- tryCatch({
      dsmSearch::get_dsm_30(bbox, key, datatype = 'SRTMGL1')
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
  chm <- dsmSearch::get_dsm_30(bbox=bbox, datatype='metaCHM')
  # filtered CHM based on the minimum tree height
  filteredCHM <- terra::ifel(chm < min_height, 0, chm)
  isolatedCHM <- terra::ifel(chm < min_height, NA, 1)
  return(list(filteredCHM, isolatedCHM))
}

get_gvi <- function(dsm, p, height, building, binary_chm) {
  tryCatch({
    v <- terra::viewshed(dsm, p, height, 0)
    v <- terra::ifel(v == 1, v, NA)

    v_values <- terra::values(v)
    if (all(is.na(v_values))) return(0)

    # Viewshed area
    v_area <- sum(terra::cellSize(v, unit = "m") * v_values, na.rm = TRUE)
    if (v_area == 0 || is.na(v_area)) return(0)

    # Resample canopy to align
    canopy_aligned <- terra::resample(binary_chm, v, method = "near")
    visible_canopy <- terra::mask(canopy_aligned, v)
    visible_canopy <- terra::ifel(visible_canopy == 1, 1, NA)
    overlap_area <- sum(terra::cellSize(visible_canopy, unit = "m") * terra::values(visible_canopy), na.rm = TRUE)

    # Compute GVI (allow >1 if canopy exceeds viewshed minus building)
    gvi <- overlap_area / max(v_area - building$g_area, 1e-6)  # use small constant to avoid division by 0
    gvi <- min(gvi, 1)  # cap at 1 if desired

    return(gvi)
  }, error = function(e) {
    warning(sprintf("GVI calculation failed: %s", e$message))
    return(0)
  })
}
# Note:
# In fact, if the building itself occupies nearly the entire viewshed,
# or the canopy overlaps very closely, v_area - building$g_area could be ≤ 0,
# yet still surrounded by trees — in which case GVI = 1 is conceptually valid.

# max(v_area - building$g_area, 1e-6) prevents division by zero while allowing close overlap
# min(gvi, 1) caps the value at 1 if needed for interpretation as a proportion
# Still returns 0 if the viewshed fails or has no values

#' @noMd
merge_elev <- function(building, dem, chm=NULL) {
  # prioritize layers:  (chm >) building > dem
  bc <- terra::overlay(r1, building, fun = function(x, y) {
    ifelse(x < y & x != 0, x, x + y)
  })
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

#' @noMd
get_bbox <- function(x) {
  bbox <- sf::st_as_sfc(sf::st_bbox(x), crs = sf::st_crs(x))
  bbox <- sf::st_transform(bbox, crs = 4326)
  return(bbox)
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
