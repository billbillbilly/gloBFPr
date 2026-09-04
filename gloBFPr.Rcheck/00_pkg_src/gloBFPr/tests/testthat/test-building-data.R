testthat::test_that("get_fused_dsm combines DEM with max building and canopy heights", {
  make_square <- function(xmin, ymin, xmax, ymax) {
    sf::st_polygon(list(matrix(
      c(xmin, ymin,
        xmin, ymax,
        xmax, ymax,
        xmax, ymin,
        xmin, ymin),
      ncol = 2,
      byrow = TRUE
    )))
  }

  buildings <- sf::st_sf(
    id = seq_len(2),
    Height = c(30, 8),
    geometry = sf::st_sfc(
      make_square(3.0000, 0.0000, 3.0002, 0.0002),
      make_square(3.0004, 0.0000, 3.0006, 0.0002),
      crs = 4326
    )
  )

  make_template <- function(bbox) {
    bbox_sfc <- sf::st_as_sfc(sf::st_bbox(
      c(xmin = bbox[1], ymin = bbox[2], xmax = bbox[3], ymax = bbox[4]),
      crs = 4326
    ))
    utm_crs <- getFromNamespace("get_utm_crs", "gloBFPr")(bbox_sfc)
    bbox_proj <- sf::st_transform(bbox_sfc, utm_crs)
    terra::rast(
      ext = terra::ext(sf::st_bbox(bbox_proj)),
      resolution = 5,
      crs = sf::st_crs(bbox_proj)$wkt
    )
  }

  mocked_get_chm <- function(bbox, min_height, datasource = "metachm") {
    r <- make_template(bbox)
    filtered <- r
    binary <- r
    terra::values(filtered) <- 15
    terra::values(binary) <- 1
    list(filtered, binary)
  }

  mocked_get_dem <- function(bbox, key) {
    r <- make_template(bbox)
    terra::values(r) <- 100
    r
  }

  namespace <- asNamespace("gloBFPr")
  old_get_chm <- get("get_chm", envir = namespace)
  old_get_dem <- get("get_dem", envir = namespace)
  unlockBinding("get_chm", namespace)
  unlockBinding("get_dem", namespace)
  assign("get_chm", mocked_get_chm, envir = namespace)
  assign("get_dem", mocked_get_dem, envir = namespace)
  lockBinding("get_chm", namespace)
  lockBinding("get_dem", namespace)
  on.exit({
    unlockBinding("get_chm", namespace)
    unlockBinding("get_dem", namespace)
    assign("get_chm", old_get_chm, envir = namespace)
    assign("get_dem", old_get_dem, envir = namespace)
    lockBinding("get_chm", namespace)
    lockBinding("get_dem", namespace)
  }, add = TRUE)

  dsm <- gloBFPr::get_fused_dsm(buildings, key = "fake-key", quiet = TRUE)
  centroids <- sf::st_coordinates(suppressWarnings(
    sf::st_transform(sf::st_centroid(buildings), sf::st_crs(dsm))
  ))
  elevations <- terra::extract(dsm, centroids)[, 1]

  testthat::expect_s4_class(dsm, "SpatRaster")
  testthat::expect_equal(elevations, c(130, 115))
})

testthat::test_that("get_fused_dsm can build DSM without canopy height", {
  make_square <- function(xmin, ymin, xmax, ymax) {
    sf::st_polygon(list(matrix(
      c(xmin, ymin,
        xmin, ymax,
        xmax, ymax,
        xmax, ymin,
        xmin, ymin),
      ncol = 2,
      byrow = TRUE
    )))
  }

  buildings <- sf::st_sf(
    id = 1,
    Height = 20,
    geometry = sf::st_sfc(make_square(3.0000, 0.0000, 3.0002, 0.0002), crs = 4326)
  )

  mocked_get_dem <- function(bbox, key) {
    bbox_sfc <- sf::st_as_sfc(sf::st_bbox(
      c(xmin = bbox[1], ymin = bbox[2], xmax = bbox[3], ymax = bbox[4]),
      crs = 4326
    ))
    utm_crs <- getFromNamespace("get_utm_crs", "gloBFPr")(bbox_sfc)
    bbox_proj <- sf::st_transform(bbox_sfc, utm_crs)
    r <- terra::rast(
      ext = terra::ext(sf::st_bbox(bbox_proj)),
      resolution = 5,
      crs = sf::st_crs(bbox_proj)$wkt
    )
    terra::values(r) <- 100
    r
  }

  namespace <- asNamespace("gloBFPr")
  old_get_dem <- get("get_dem", envir = namespace)
  unlockBinding("get_dem", namespace)
  assign("get_dem", mocked_get_dem, envir = namespace)
  lockBinding("get_dem", namespace)
  on.exit({
    unlockBinding("get_dem", namespace)
    assign("get_dem", old_get_dem, envir = namespace)
    lockBinding("get_dem", namespace)
  }, add = TRUE)

  dsm <- gloBFPr::get_fused_dsm(
    buildings,
    datasource_canopy_height = NULL,
    key = "fake-key",
    quiet = TRUE
  )
  centroid <- sf::st_coordinates(suppressWarnings(
    sf::st_transform(sf::st_centroid(buildings), sf::st_crs(dsm))
  ))
  elevation <- terra::extract(dsm, centroid)[, 1]

  testthat::expect_equal(elevation, 120)
})
