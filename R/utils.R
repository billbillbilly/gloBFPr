#' @importFrom terra ext
#' @importFrom terra vect
#' @importFrom terra rast
#' @importFrom terra rasterize
#' @importFrom terra align
#' @importFrom terra mask
#' @importFrom sf st_centroid
#' @importFrom sf st_coordinates
#' @importFrom sf st_transform
#' @importFrom sf st_crs
#' @importFrom sf st_as_sfc

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
