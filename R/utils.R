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

#' @noMd
add_2d_3d <- function(projected_poly, type) {
  if (type == '2d' || type == 'all') {
    #### 2D metrics ####
    projected_poly$area_2d <- as.numeric(sf::st_area(projected_poly))
    projected_poly$perimeter <- as.numeric(sf::st_perimeter(projected_poly))
    # Perimeter-Area Ratio
    projected_poly$PARA <- projected_poly$perimeter / projected_poly$area_2d
    projected_poly$rectangularity <- projected_poly$area_2d / sf::st_area(sf::st_minimum_rotated_rectangle(projected_poly))
  } else if (type == '3d' || type == 'all') {
    #### 3D metrics ####
    projected_poly$vertical_surface <- projected_poly$perimeter * projected_poly$Height
    projected_poly$total_surface <- projected_poly$vertical_surface + projected_poly$area_2d
    projected_poly$volume <- projected_poly$area_2d * projected_poly$Height
    projected_poly$hemisphericality <- 0
    projected_poly$convexity <- 0
    projected_poly$cuboidness <- 0
    for (i in 1:nrow(projected_poly)) {
      poly_ <- projected_poly[i,]
      projected_poly$hemisphericality[i] <- cal_hemisphericality(poly_)
      projected_poly$convexity[i] <- cal_convexity(poly_)
      projected_poly$cuboidness[i] <- cal_cuboidness(poly_)
    }
    projected_poly$fractality <- 1 - log(projected_poly) / (1.5 * log(projected_poly$total_surface))

  } else {
    return(projected_poly)
  }
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

  vol_hemisphere <- (2/3) * pi^(-0.5) * poly_$area_2d^(3/2)
  hemisphericality <- (poly_$volume - vol_hemisphere) / volume_hemisphere
  return(hemisphericality)
}

#' @noMd
cal_convexity <- function(poly_) {
  vertices <- sf::st_coordinates(poly_)[, 1:2]
  multipoint <- sf::st_multipoint(vertices)
  multipoint <- sf::st_sfc(multipoint, crs = sf::st_crs(poly_))
  return(sf::st_convex_hull(multipoint))
}

#' @noMd


#' @noMd
add_pop <- function(bbox, projected_poly, year) {
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

#' @noMd


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
