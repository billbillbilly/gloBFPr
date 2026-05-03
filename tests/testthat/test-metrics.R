testthat::test_that("runs correctly", {
  res_pop <- gloBFPr::get_pop_density(year = 2025)
  res_morphology <- gloBFPr::get_morphology()
  res_neighbors <- gloBFPr::get_neighbors()
  res_dng <- gloBFPr::get_dng()
  res_residential <- gloBFPr::get_residential()
  res_residential_missing_year <- gloBFPr::get_residential(gloBFPr::globfp_example)

  testthat::expect_type(res_pop, "NULL")
  testthat::expect_type(res_morphology, "NULL")
  testthat::expect_type(res_neighbors, "NULL")
  testthat::expect_type(res_dng, "NULL")
  testthat::expect_type(res_residential, "NULL")
  testthat::expect_type(res_residential_missing_year, "NULL")
})

testthat::test_that("get_residential classifies by residential percentage", {
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
    id = seq_len(4),
    Height = rep(10, 4),
    geometry = sf::st_sfc(
      make_square(3.0000, 0.0000, 3.0002, 0.0002),
      make_square(3.0020, 0.0000, 3.0022, 0.0002),
      make_square(3.0040, 0.0000, 3.0042, 0.0002),
      make_square(3.0060, 0.0000, 3.0062, 0.0002),
      crs = 4326
    )
  )

  mocked_get_GHSres <- function(bbox, year) {
    utm_crs <- getFromNamespace("get_utm_crs", "gloBFPr")(bbox)
    projected <- sf::st_transform(buildings, utm_crs)
    centroids <- sf::st_coordinates(sf::st_centroid(projected))

    template <- terra::rast(
      xmin = min(centroids[, "X"]) - 50,
      xmax = max(centroids[, "X"]) + 50,
      ymin = min(centroids[, "Y"]) - 50,
      ymax = max(centroids[, "Y"]) + 50,
      resolution = 50,
      crs = paste0("EPSG:", utm_crs)
    )

    total <- template
    nres <- template
    terra::values(total) <- NA_real_
    terra::values(nres) <- NA_real_

    cells <- terra::cellFromXY(template, centroids)
    total_values <- c(100, 100, 0, NA)
    nres_values <- c(20, 30, 0, NA)
    total[cells] <- total_values
    nres[cells] <- nres_values

    list(total = total, nres = nres, res = total - nres)
  }

  namespace <- asNamespace("gloBFPr")
  old_get_GHSres <- get("get_GHSres", envir = namespace)
  unlockBinding("get_GHSres", namespace)
  assign("get_GHSres", mocked_get_GHSres, envir = namespace)
  lockBinding("get_GHSres", namespace)
  on.exit({
    unlockBinding("get_GHSres", namespace)
    assign("get_GHSres", old_get_GHSres, envir = namespace)
    lockBinding("get_GHSres", namespace)
  }, add = TRUE)

  result <- gloBFPr::get_residential(buildings, year = 2025, threshold = 70)

  testthat::expect_equal(result$total_built_vals, c(100, 100, 0, NA))
  testthat::expect_equal(result$nres_vals, c(20, 30, 0, NA))
  testthat::expect_equal(result$res_vals, c(80, 70, 0, NA))
  testthat::expect_equal(result$res_pct, c(80, 70, NA, NA))
  testthat::expect_equal(result$res, c(TRUE, TRUE, FALSE, FALSE))
})

testthat::test_that("get_dng uses filtered green candidates and handles empty results", {
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

  building <- sf::st_sf(
    id = 1,
    Height = 10,
    geometry = sf::st_sfc(make_square(3.0000, 0.0000, 3.0002, 0.0002), crs = 4326)
  )

  captured_min_tree_height <- NULL
  mocked_get_greenspace <- function(bbox = NULL, buffer = NULL, type = NULL,
                                    zoom = 17, year = NULL, min_tree_height = 2) {
    captured_min_tree_height <<- min_tree_height
    buffer_vect <- terra::vect(buffer)
    template <- terra::rast(
      ext = terra::ext(buffer_vect),
      resolution = 20,
      crs = sf::st_crs(buffer)$wkt
    )
    terra::values(template) <- 0

    centroid <- sf::st_coordinates(sf::st_centroid(buffer))
    green_xy <- matrix(c(centroid[1, "X"] + 40, centroid[1, "Y"]), ncol = 2)
    template[terra::cellFromXY(template, green_xy)] <- 1
    template
  }

  namespace <- asNamespace("gloBFPr")
  old_get_greenspace <- get("get_greenspace", envir = namespace)
  unlockBinding("get_greenspace", namespace)
  assign("get_greenspace", mocked_get_greenspace, envir = namespace)
  lockBinding("get_greenspace", namespace)
  on.exit({
    unlockBinding("get_greenspace", namespace)
    assign("get_greenspace", old_get_greenspace, envir = namespace)
    lockBinding("get_greenspace", namespace)
  }, add = TRUE)

  result <- gloBFPr::get_dng(
    building,
    datasource = "metachm",
    min_tree_height = 7,
    radius = 100,
    min_area = 1
  )

  testthat::expect_equal(captured_min_tree_height, 7)
  testthat::expect_true(is.numeric(result$dng))
  testthat::expect_true(is.finite(result$dng))
  testthat::expect_gt(result$dng, 0)

  no_candidate <- gloBFPr::get_dng(
    building,
    datasource = "metachm",
    radius = 100,
    min_area = 1e9
  )

  testthat::expect_true(is.na(no_candidate$dng))
})

testthat::test_that("get_greenspace metachm returns the binary CHM raster", {
  bbox <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = 3.0000, ymin = 0.0000, xmax = 3.0002, ymax = 0.0002),
    crs = 4326
  ))

  mocked_get_chm <- function(bbox, min_height) {
    testthat::expect_true(is.numeric(bbox))
    testthat::expect_length(bbox, 4)
    r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                     crs = "EPSG:3857")
    filtered <- r
    binary <- r
    terra::values(filtered) <- c(0, min_height + 1, 0, min_height + 2)
    terra::values(binary) <- c(0, 1, 0, 1)
    list(filtered, binary)
  }

  namespace <- asNamespace("gloBFPr")
  old_get_chm <- get("get_chm", envir = namespace)
  unlockBinding("get_chm", namespace)
  assign("get_chm", mocked_get_chm, envir = namespace)
  lockBinding("get_chm", namespace)
  on.exit({
    unlockBinding("get_chm", namespace)
    assign("get_chm", old_get_chm, envir = namespace)
    lockBinding("get_chm", namespace)
  }, add = TRUE)

  result <- getFromNamespace("get_greenspace", "gloBFPr")(
    bbox = bbox,
    type = "metachm",
    min_tree_height = 5
  )

  testthat::expect_s4_class(result, "SpatRaster")
  testthat::expect_equal(as.vector(terra::values(result)), c(0, 1, 0, 1))
})

testthat::test_that("get_bgvi flattens target building and samples floors", {
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

  building <- sf::st_sf(
    id = 1,
    Height = 30,
    geometry = sf::st_sfc(make_square(3.0000, 0.0000, 3.0002, 0.0002), crs = 4326)
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
      resolution = 10,
      crs = sf::st_crs(bbox_proj)$wkt
    )
  }

  mocked_get_chm <- function(bbox, min_height) {
    r <- make_template(bbox)
    filtered <- r
    binary <- r
    terra::values(filtered) <- 0
    terra::values(binary) <- 0
    list(filtered, binary)
  }

  mocked_get_dem <- function(bbox, key) {
    r <- make_template(bbox)
    terra::values(r) <- 0
    r
  }

  mocked_get_gvi <- function(dsm, p, height, r, building, binary_chm) {
    target_elevation <- terra::extract(dsm, matrix(p, ncol = 2))[, 1]
    testthat::expect_equal(target_elevation, 0)
    height / 100
  }

  namespace <- asNamespace("gloBFPr")
  old_get_chm <- get("get_chm", envir = namespace)
  old_get_dem <- get("get_dem", envir = namespace)
  old_get_gvi <- get("get_gvi", envir = namespace)
  unlockBinding("get_chm", namespace)
  unlockBinding("get_dem", namespace)
  unlockBinding("get_gvi", namespace)
  assign("get_chm", mocked_get_chm, envir = namespace)
  assign("get_dem", mocked_get_dem, envir = namespace)
  assign("get_gvi", mocked_get_gvi, envir = namespace)
  lockBinding("get_chm", namespace)
  lockBinding("get_dem", namespace)
  lockBinding("get_gvi", namespace)
  on.exit({
    unlockBinding("get_chm", namespace)
    unlockBinding("get_dem", namespace)
    unlockBinding("get_gvi", namespace)
    assign("get_chm", old_get_chm, envir = namespace)
    assign("get_dem", old_get_dem, envir = namespace)
    assign("get_gvi", old_get_gvi, envir = namespace)
    lockBinding("get_chm", namespace)
    lockBinding("get_dem", namespace)
    lockBinding("get_gvi", namespace)
  }, add = TRUE)

  result <- gloBFPr::get_bgvi(
    building,
    floor = TRUE,
    floor_step = 4,
    workers = 1,
    key = "fake-key"
  )

  expected_gvis <- c(1.7, 13.7, 25.7, 28.7) / 100
  testthat::expect_equal(result$estimated_floors, 10)
  testthat::expect_equal(result$mean_gvi, mean(expected_gvis))
  testthat::expect_equal(result$min_gvi, min(expected_gvis))
  testthat::expect_equal(result$max_gvi, max(expected_gvis))
  testthat::expect_equal(result$sd_gvi, stats::sd(expected_gvis))
  testthat::expect_error(
    gloBFPr::get_bgvi(building, datasource = "esri", key = "fake-key", workers = 1),
    "should be"
  )
})
