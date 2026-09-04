testthat::test_that("runs correctly", {
  res_pop <- gloBFPr:::get_pop(year = 2025)
  res_morphology <- gloBFPr::get_morphology()
  res_neighbors <- gloBFPr::get_neighbors()
  res_dng <- gloBFPr::get_dng()
  res_residential <- gloBFPr:::get_residential()
  res_residential_missing_year <- gloBFPr:::get_residential(gloBFPr::globfp_example)

  testthat::expect_type(res_pop, "NULL")
  testthat::expect_type(res_morphology, "NULL")
  testthat::expect_type(res_neighbors, "NULL")
  testthat::expect_type(res_dng, "NULL")
  testthat::expect_type(res_residential, "NULL")
  testthat::expect_type(res_residential_missing_year, "NULL")
})

testthat::test_that("get_neighbors reports min and max neighbor distances", {
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
      make_square(3.0004, 0.0000, 3.0006, 0.0002),
      make_square(3.0000, 0.0004, 3.0002, 0.0006),
      make_square(3.0004, 0.0004, 3.0006, 0.0006),
      crs = 4326
    )
  )

  result <- gloBFPr::get_neighbors(buildings, radius = 100, quiet = TRUE)

  testthat::expect_true(all(c(
    "n_count", "m_ndist", "min_ndist", "max_ndist", "sd_ndist"
  ) %in% names(result)))

  with_neighbors <- result$n_count > 0
  testthat::expect_true(any(with_neighbors))
  testthat::expect_true(all(is.finite(result$min_ndist[with_neighbors])))
  testthat::expect_true(all(is.finite(result$max_ndist[with_neighbors])))
  testthat::expect_true(all(result$min_ndist[with_neighbors] <= result$m_ndist[with_neighbors]))
  testthat::expect_true(all(result$m_ndist[with_neighbors] <= result$max_ndist[with_neighbors]))
})

testthat::test_that("get_morphology treats group_id features as one building", {
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
    id = seq_len(3),
    group_id = c(1, 1, 2),
    Height = c(10, 20, 5),
    geometry = sf::st_sfc(
      make_square(0, 0, 10, 10),
      make_square(10, 0, 20, 10),
      make_square(40, 0, 50, 10),
      crs = 3857
    )
  )

  result <- gloBFPr::get_morphology(
    buildings,
    metrics = c("g_area", "pmeter", "v_surf", "t_surf", "vol", "obb_vol",
                "pa_ratio", "rec", "fra", "cbn"),
    quiet = TRUE
  )

  testthat::expect_equal(nrow(result), 3)
  testthat::expect_equal(as.numeric(sf::st_area(result)), c(100, 100, 100))
  testthat::expect_equal(result$g_area, c(200, 200, 100))
  testthat::expect_equal(result$pmeter, c(60, 60, 40))
  testthat::expect_equal(result$v_surf, c(1000, 1000, 200))
  testthat::expect_equal(result$t_surf, c(1200, 1200, 300))
  testthat::expect_equal(result$vol, c(3000, 3000, 500))
  testthat::expect_equal(result$obb_vol, c(4000, 4000, 500))
  testthat::expect_equal(result$group_id, c(1, 1, 2))
})

testthat::test_that("get_pop allocates cell population by building volume", {
  make_square <- function(xmin, ymin, xmax, ymax) {
    sf::st_polygon(list(matrix(
      c(xmin, ymin, xmin, ymax, xmax, ymax, xmax, ymin, xmin, ymin),
      ncol = 2, byrow = TRUE
    )))
  }
  buildings <- sf::st_sf(
    id = seq_len(2), Height = c(10, 20),
    geometry = sf::st_sfc(
      make_square(3.0000, 0.0000, 3.0002, 0.0002),
      make_square(3.0004, 0.0000, 3.0006, 0.0002),
      crs = 4326
    )
  )
  captured_points <- NULL
  mocked_get_GHSpop <- function(bbox = NULL, year = NULL, points = NULL,
                                polygons = NULL, quiet = FALSE) {
    captured_points <<- points
    testthat::expect_equal(year, 2025)
    testthat::expect_true(quiet)
    testthat::expect_null(polygons)
    list(pop_total = c(300, 300), cell_id = c("tile:1", "tile:1"))
  }
  namespace <- asNamespace("gloBFPr")
  old_get_GHSpop <- get("get_GHSpop", envir = namespace)
  unlockBinding("get_GHSpop", namespace)
  assign("get_GHSpop", mocked_get_GHSpop, envir = namespace)
  lockBinding("get_GHSpop", namespace)
  on.exit({
    unlockBinding("get_GHSpop", namespace)
    assign("get_GHSpop", old_get_GHSpop, envir = namespace)
    lockBinding("get_GHSpop", namespace)
  }, add = TRUE)
  result <- gloBFPr:::get_pop(buildings, year = 2025, quiet = TRUE)
  testthat::expect_s3_class(captured_points, "sf")
  testthat::expect_equal(result$pop_total, c(100, 200), tolerance = 1e-6)
  testthat::expect_false("pop_den" %in% names(result))
})

testthat::test_that("get_pop treats group_id features as one building", {
  make_square <- function(xmin, ymin, xmax, ymax) {
    sf::st_polygon(list(matrix(
      c(xmin, ymin, xmin, ymax, xmax, ymax, xmax, ymin, xmin, ymin),
      ncol = 2, byrow = TRUE
    )))
  }
  buildings <- sf::st_sf(
    id = seq_len(3), group_id = c(1, 1, 2), Height = c(10, 20, 8),
    geometry = sf::st_sfc(
      make_square(0, 0, 10, 10), make_square(10, 0, 20, 10), make_square(40, 0, 50, 10),
      crs = 3857
    )
  )
  mocked_get_GHSpop <- function(bbox = NULL, year = NULL, points = NULL,
                                polygons = NULL, quiet = FALSE) {
    n <- if (!is.null(points)) nrow(points) else nrow(polygons)
    testthat::expect_equal(n, 2)
    list(pop_total = c(380, 380), cell_id = c("tile:1", "tile:1"))
  }
  namespace <- asNamespace("gloBFPr")
  old_get_GHSpop <- get("get_GHSpop", envir = namespace)
  unlockBinding("get_GHSpop", namespace)
  assign("get_GHSpop", mocked_get_GHSpop, envir = namespace)
  lockBinding("get_GHSpop", namespace)
  on.exit({
    unlockBinding("get_GHSpop", namespace)
    assign("get_GHSpop", old_get_GHSpop, envir = namespace)
    lockBinding("get_GHSpop", namespace)
  }, add = TRUE)
  result <- gloBFPr:::get_pop(buildings, year = 2025, quiet = TRUE)
  testthat::expect_equal(result$pop_total, c(300, 300, 80), tolerance = 1e-6)
  testthat::expect_false("pop_den" %in% names(result))
})

testthat::test_that("population allocation falls back from volume to area", {
  allocate_population_cpp <- gloBFPr:::allocate_population_cpp
  area_weighted <- allocate_population_cpp(
    pop_total = c(120, 120), cell_group = c(1L, 1L),
    volume = c(0, NA_real_), area = c(10, 30)
  )
  testthat::expect_equal(area_weighted, c(30, 90))
  equal_weighted <- allocate_population_cpp(
    pop_total = c(90, 90, 90), cell_group = c(1L, 1L, 1L),
    volume = c(0, NA_real_, 0), area = c(0, NA_real_, 0)
  )
  testthat::expect_equal(equal_weighted, c(30, 30, 30))
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
    centroids <- sf::st_coordinates(suppressWarnings(sf::st_centroid(projected)))

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

  result <- gloBFPr:::get_residential(buildings, year = 2025, threshold = 70)

  testthat::expect_equal(result$total_built_vals, c(100, 100, 0, NA))
  testthat::expect_equal(result$nres_vals, c(20, 30, 0, NA))
  testthat::expect_equal(result$res_vals, c(80, 70, 0, NA))
  testthat::expect_equal(result$res_pct, c(80, 70, NA, NA))
  testthat::expect_equal(result$res, c(TRUE, TRUE, FALSE, FALSE))
})

testthat::test_that("get_neighbors computes on group features", {
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
    id = seq_len(3),
    group_id = c(1, 1, 2),
    Height = c(10, 20, 10),
    geometry = sf::st_sfc(
      make_square(3.0000, 0.0000, 3.0002, 0.0002),
      make_square(3.0002, 0.0000, 3.0004, 0.0002),
      make_square(3.0008, 0.0000, 3.0010, 0.0002),
      crs = 4326
    )
  )

  result <- gloBFPr::get_neighbors(buildings, radius = 200, quiet = TRUE)

  testthat::expect_equal(result$n_count[1], result$n_count[2])
  testthat::expect_equal(result$m_ndist[1], result$m_ndist[2])
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
    if (is.null(buffer)) buffer <- bbox
    buffer <- sf::st_transform(buffer, getFromNamespace("get_utm_crs", "gloBFPr")(buffer))
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

testthat::test_that("get_dng copies group distance back to features", {
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
    group_id = c(1, 1),
    Height = c(10, 12),
    geometry = sf::st_sfc(
      make_square(3.0000, 0.0000, 3.0002, 0.0002),
      make_square(3.0002, 0.0000, 3.0004, 0.0002),
      crs = 4326
    )
  )

  mocked_get_greenspace <- function(bbox = NULL, buffer = NULL, type = NULL,
                                    zoom = 17, year = NULL, min_tree_height = 2) {
    if (is.null(buffer)) buffer <- bbox
    buffer <- sf::st_transform(buffer, getFromNamespace("get_utm_crs", "gloBFPr")(buffer))
    buffer_vect <- terra::vect(buffer)
    template <- terra::rast(
      ext = terra::ext(buffer_vect),
      resolution = 20,
      crs = sf::st_crs(buffer)$wkt
    )
    terra::values(template) <- 0
    centroid <- sf::st_coordinates(suppressWarnings(sf::st_centroid(buffer)))
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
    buildings,
    datasource = "metachm",
    radius = 100,
    min_area = 1
  )

  testthat::expect_equal(result$dng[1], result$dng[2])
  testthat::expect_true(is.finite(result$dng[1]))
})

testthat::test_that("build_network_graph nodes shared intersections and weights edges", {
  build_network_graph <- gloBFPr:::build_network_graph
  network_distance_to_targets <- gloBFPr:::network_distance_to_targets

  # Two ways meeting at (100, 0): an L-shaped street 200 m long end to end.
  roads <- sf::st_sf(
    id = seq_len(2),
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(100, 0, 100, 100), ncol = 2, byrow = TRUE)),
      crs = 32617
    )
  )

  # max_segment = Inf keeps the raw vertices so the node count is predictable.
  g <- build_network_graph(roads, max_segment = Inf)
  testthat::expect_equal(igraph::vcount(g$graph), 3L)
  testthat::expect_equal(igraph::ecount(g$graph), 2L)
  testthat::expect_equal(sum(igraph::E(g$graph)$weight), 200)

  # Origin sits 5 m off one end, target 5 m off the other: 5 + 200 + 5.
  routed <- network_distance_to_targets(
    g,
    matrix(c(0, -5), ncol = 2),
    matrix(c(105, 100), ncol = 2)
  )
  testthat::expect_equal(routed$distance, 210)
  testthat::expect_equal(routed$index, 1L)
})

testthat::test_that("network routing picks the network-nearest, not the euclid-nearest, patch", {
  build_network_graph <- gloBFPr:::build_network_graph
  network_distance_to_targets <- gloBFPr:::network_distance_to_targets

  roads <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0, 100, 100), ncol = 2, byrow = TRUE)),
      crs = 32617
    )
  )
  g <- build_network_graph(roads, max_segment = Inf)

  # Target 1 is a long way round the L; target 2 is further as the crow flies
  # but sits directly on the street.
  targets <- matrix(c(105, 100,
                      50, 0), ncol = 2, byrow = TRUE)
  routed <- network_distance_to_targets(g, matrix(c(0, -5), ncol = 2), targets)

  testthat::expect_equal(routed$index, 2L)
  testthat::expect_equal(routed$distance, 55)
})

testthat::test_that("get_dng routes along a supplied network and flags the method", {
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

  # Green pixel 40 m due east of the building centroid.
  mocked_get_greenspace <- function(bbox = NULL, buffer = NULL, type = NULL,
                                    zoom = 17, year = NULL, min_tree_height = 2) {
    if (is.null(buffer)) buffer <- bbox
    buffer <- sf::st_transform(buffer, getFromNamespace("get_utm_crs", "gloBFPr")(buffer))
    buffer_vect <- terra::vect(buffer)
    template <- terra::rast(
      ext = terra::ext(buffer_vect),
      resolution = 20,
      crs = sf::st_crs(buffer)$wkt
    )
    terra::values(template) <- 0
    centroid <- sf::st_coordinates(suppressWarnings(sf::st_centroid(buffer)))
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

  utm_crs <- gloBFPr:::get_utm_crs(gloBFPr:::get_bbox(building))
  centroid <- sf::st_coordinates(
    suppressWarnings(sf::st_centroid(sf::st_transform(building, utm_crs)))
  )

  # A detour street: north 60 m, east 40 m, back south 60 m. Any route to the
  # green pixel must go the long way round, so the network distance has to
  # exceed the 40 m straight-line distance.
  detour <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(
        c(centroid[1, "X"],       centroid[1, "Y"],
          centroid[1, "X"],       centroid[1, "Y"] + 60,
          centroid[1, "X"] + 40, centroid[1, "Y"] + 60,
          centroid[1, "X"] + 40, centroid[1, "Y"]),
        ncol = 2, byrow = TRUE
      )),
      crs = utm_crs
    )
  )

  euclid <- gloBFPr::get_dng(
    building, datasource = "metachm", radius = 100, min_area = 1
  )
  routed <- gloBFPr::get_dng(
    building, datasource = "metachm", radius = 100, min_area = 1,
    network = detour, quiet = TRUE
  )

  testthat::expect_equal(euclid$dng_method, "euclidean")
  testthat::expect_equal(routed$dng_method, "network")
  testthat::expect_gt(routed$dng, euclid$dng)
})

testthat::test_that("get_dng falls back to euclidean when the network is unusable", {
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

  mocked_get_greenspace <- function(bbox = NULL, buffer = NULL, type = NULL,
                                    zoom = 17, year = NULL, min_tree_height = 2) {
    if (is.null(buffer)) buffer <- bbox
    buffer <- sf::st_transform(buffer, getFromNamespace("get_utm_crs", "gloBFPr")(buffer))
    buffer_vect <- terra::vect(buffer)
    template <- terra::rast(
      ext = terra::ext(buffer_vect),
      resolution = 20,
      crs = sf::st_crs(buffer)$wkt
    )
    terra::values(template) <- 0
    centroid <- sf::st_coordinates(suppressWarnings(sf::st_centroid(buffer)))
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

  utm_crs <- gloBFPr:::get_utm_crs(gloBFPr:::get_bbox(building))
  empty_net <- sf::st_sf(
    id = integer(0),
    geometry = sf::st_sfc(crs = utm_crs)
  )

  result <- gloBFPr::get_dng(
    building, datasource = "metachm", radius = 100, min_area = 1,
    network = empty_net, quiet = TRUE
  )

  testthat::expect_equal(result$dng_method, "euclidean")
  testthat::expect_true(is.finite(result$dng))

  testthat::expect_error(
    gloBFPr::get_dng(
      building, datasource = "metachm", radius = 100, min_area = 1,
      network = 42, quiet = TRUE
    ),
    "must be NULL"
  )
})

testthat::test_that("plot_bgvi_viewshed returns single-building diagnostic layers", {
  testthat::skip_if_not_installed("viewscape")

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
    id = c("target", "neighbor"),
    Height = c(18, 12),
    geometry = sf::st_sfc(
      make_square(3.0000, 0.0000, 3.0002, 0.0002),
      make_square(3.0006, 0.0000, 3.0008, 0.0002),
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
      resolution = 10,
      crs = sf::st_crs(bbox_proj)$wkt
    )
  }

  mocked_get_chm <- function(bbox, min_height, datasource = "metachm") {
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
  mocked_get_greenspace <- function(bbox = NULL, buffer = NULL, type = NULL,
                                    zoom = 17, year = NULL, min_tree_height = 2) {
    bbox_vector <- if (is.numeric(bbox) && length(bbox) == 4) {
      bbox
    } else {
      getFromNamespace("bbox_poly_to_list", "gloBFPr")(bbox)
    }
    r <- make_template(bbox_vector)
    terra::values(r) <- 1
    r
  }

  namespace <- asNamespace("gloBFPr")
  old_get_chm <- get("get_chm", envir = namespace)
  old_get_dem <- get("get_dem", envir = namespace)
  old_get_greenspace <- get("get_greenspace", envir = namespace)
  unlockBinding("get_chm", namespace)
  unlockBinding("get_dem", namespace)
  unlockBinding("get_greenspace", namespace)
  assign("get_chm", mocked_get_chm, envir = namespace)
  assign("get_dem", mocked_get_dem, envir = namespace)
  assign("get_greenspace", mocked_get_greenspace, envir = namespace)
  lockBinding("get_chm", namespace)
  lockBinding("get_dem", namespace)
  lockBinding("get_greenspace", namespace)
  on.exit({
    unlockBinding("get_chm", namespace)
    unlockBinding("get_dem", namespace)
    unlockBinding("get_greenspace", namespace)
    assign("get_chm", old_get_chm, envir = namespace)
    assign("get_dem", old_get_dem, envir = namespace)
    assign("get_greenspace", old_get_greenspace, envir = namespace)
    lockBinding("get_chm", namespace)
    lockBinding("get_dem", namespace)
    lockBinding("get_greenspace", namespace)
  }, add = TRUE)

  result <- gloBFPr::plot_bgvi_viewshed(
    buildings,
    building = "target",
    floor = 3,
    orientation = 90,
    field_of_view = 90,
    datasource_greenspace = "esri",
    radius = 120,
    key = "fake-key",
    plot = FALSE,
    quiet = TRUE
  )

  testthat::expect_equal(result$height, 7.7)
  testthat::expect_equal(result$building$id, "target")
  testthat::expect_s4_class(result$viewshed_raster, "SpatRaster")
  testthat::expect_s4_class(result$visible_green, "SpatRaster")
  viewscape_raster <- viewscape:::filter_invisible(result$viewshed, TRUE)
  testthat::expect_equal(
    as.vector(terra::ext(result$viewshed_raster)),
    as.vector(terra::ext(viewscape_raster)),
    tolerance = 1e-6
  )
  testthat::expect_true(terra::compareGeom(
    result$viewshed_raster,
    result$binary_green,
    stopOnError = FALSE
  ))
  testthat::expect_true(is.finite(result$gvi))
  testthat::expect_true(sum(terra::values(result$sector_mask) == 1, na.rm = TRUE) <
	                          terra::ncell(result$sector_mask))
  testthat::expect_true(any(terra::values(result$plot_raster) == 4, na.rm = TRUE))
  p <- result$viewpoint
  corner_building <- sf::st_sf(
    id = "outside-radius",
    Height = 8,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(p[1] + 85, p[2] + 85,
        p[1] + 85, p[2] + 95,
        p[1] + 95, p[2] + 95,
        p[1] + 95, p[2] + 85,
        p[1] + 85, p[2] + 85),
      ncol = 2,
      byrow = TRUE
    ))), crs = sf::st_crs(result$building))
  )
  radius_limited <- getFromNamespace("bgvi_plot_raster", "gloBFPr")(
    result,
    buildings = rbind(result$building[, c("id", "Height", "geometry")], corner_building)
  )
  corner_cell <- terra::cellFromXY(radius_limited, matrix(c(p[1] + 90, p[2] + 90), ncol = 2))
  testthat::expect_false(terra::values(radius_limited)[corner_cell] %in% c(5, 6))

  empty_result <- result
  terra::values(empty_result$viewshed_raster) <- 0
  terra::values(empty_result$visible_green) <- NA_real_
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  testthat::expect_silent(
    getFromNamespace("plot_bgvi_viewshed_result", "gloBFPr")(
      empty_result,
      buildings = buildings
    )
  )
})

testthat::test_that("get_greenspace metachm returns the binary CHM raster", {
  bbox <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = 3.0000, ymin = 0.0000, xmax = 3.0002, ymax = 0.0002),
    crs = 4326
  ))

  mocked_get_chm <- function(bbox, min_height, datasource = "metachm") {
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

  mocked_get_chm <- function(bbox, min_height, datasource = "metachm") {
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

  mocked_get_greenspace <- function(bbox = NULL, buffer = NULL, type = NULL,
                                    zoom = 17, year = NULL, min_tree_height = 2) {
    bbox_vector <- if (is.numeric(bbox) && length(bbox) == 4) {
      bbox
    } else {
      getFromNamespace("bbox_poly_to_list", "gloBFPr")(bbox)
    }
    r <- make_template(bbox_vector)
    terra::values(r) <- 1
    r
  }

  mocked_get_gvi <- function(dsm, p, height, r, building, binary_chm,
                             directions = NULL, field_of_view = 45) {
    target_elevation <- terra::extract(dsm, matrix(p, ncol = 2))[, 1]
    if (!isTRUE(all.equal(target_elevation, 0))) {
      stop("Target building footprint was not flattened.")
    }
    value <- height / 100
    green_area <- value * 500
    if (is.null(directions)) return(list(gvi = value, green_area = green_area))
    directional_offsets <- seq_along(directions) / 1000
    directional <- lapply(seq_along(directions), function(i) {
      list(gvi = value + directional_offsets[i], green_area = green_area)
    })
    names(directional) <- directions
    c(list(overall = list(gvi = value, green_area = green_area)), directional)
  }

  namespace <- asNamespace("gloBFPr")
  old_get_chm <- get("get_chm", envir = namespace)
  old_get_dem <- get("get_dem", envir = namespace)
  old_get_greenspace <- get("get_greenspace", envir = namespace)
  old_get_gvi <- get("get_gvi", envir = namespace)
  unlockBinding("get_chm", namespace)
  unlockBinding("get_dem", namespace)
  unlockBinding("get_greenspace", namespace)
  unlockBinding("get_gvi", namespace)
  assign("get_chm", mocked_get_chm, envir = namespace)
  assign("get_dem", mocked_get_dem, envir = namespace)
  assign("get_greenspace", mocked_get_greenspace, envir = namespace)
  assign("get_gvi", mocked_get_gvi, envir = namespace)
  lockBinding("get_chm", namespace)
  lockBinding("get_dem", namespace)
  lockBinding("get_greenspace", namespace)
  lockBinding("get_gvi", namespace)
  on.exit({
    unlockBinding("get_chm", namespace)
    unlockBinding("get_dem", namespace)
    unlockBinding("get_greenspace", namespace)
    unlockBinding("get_gvi", namespace)
    assign("get_chm", old_get_chm, envir = namespace)
    assign("get_dem", old_get_dem, envir = namespace)
    assign("get_greenspace", old_get_greenspace, envir = namespace)
    assign("get_gvi", old_get_gvi, envir = namespace)
    lockBinding("get_chm", namespace)
    lockBinding("get_dem", namespace)
    lockBinding("get_greenspace", namespace)
    lockBinding("get_gvi", namespace)
  }, add = TRUE)

  result <- gloBFPr::get_bgvi(
    building,
    floor = TRUE,
    floor_step = 4,
    workers = 1,
    key = "fake-key",
    quiet = TRUE
  )

  expected_gvis <- c(1.7, 13.7, 25.7, 28.7) / 100
  testthat::expect_equal(result$estimated_floors, 10)
  testthat::expect_equal(result$mean_gvi, mean(expected_gvis))
  testthat::expect_equal(result$bottom_gvi, 1.7 / 100)
  testthat::expect_equal(result$top_gvi, 28.7 / 100)
  testthat::expect_equal(result$min_gvi, min(expected_gvis))
  testthat::expect_equal(result$max_gvi, max(expected_gvis))
  testthat::expect_equal(result$sd_gvi, stats::sd(expected_gvis))
  testthat::expect_equal(result$bottom_green_area, 1.7 / 100 * 500, tolerance = 1e-6)
  testthat::expect_equal(result$top_green_area, 28.7 / 100 * 500, tolerance = 1e-6)
  expected_green_areas <- c(1.7, 13.7, 25.7, 28.7) / 100 * 500
  testthat::expect_equal(result$mean_green_area, mean(expected_green_areas), tolerance = 1e-6)
  greenspace_only <- gloBFPr::get_bgvi(
    building,
    datasource_canopy_height = NULL,
    datasource_greenspace = "esri",
    workers = 1,
    key = "fake-key",
    quiet = TRUE
  )
  testthat::expect_equal(greenspace_only$mean_gvi, mean(c(1.7, 28.7) / 100))
  testthat::expect_equal(greenspace_only$bottom_gvi, 1.7 / 100)
  testthat::expect_equal(greenspace_only$top_gvi, 28.7 / 100)
  short_threshold <- gloBFPr::get_bgvi(
    building,
    short_building_threshold = 40,
    workers = 1,
    key = "fake-key",
    quiet = TRUE
  )
  testthat::expect_equal(short_threshold$mean_gvi, 1.7 / 100)
  testthat::expect_equal(short_threshold$bottom_gvi, 1.7 / 100)
  testthat::expect_equal(short_threshold$top_gvi, 1.7 / 100)
  directional <- gloBFPr::get_bgvi(
    building,
    directions = c("south", "east"),
    field_of_view = 45,
    workers = 1,
    key = "fake-key",
    quiet = TRUE
  )
  testthat::expect_true(all(c(
    "gvi_bottom_south", "gvi_top_south", "mean_gvi_south",
    "min_gvi_south", "max_gvi_south", "sd_gvi_south",
    "gvi_bottom_east", "gvi_top_east", "mean_gvi_east"
  ) %in% names(directional)))
  testthat::expect_equal(directional$gvi_bottom_south, 1.7 / 100 + 0.001)
  testthat::expect_equal(directional$gvi_top_south, 28.7 / 100 + 0.001)
  testthat::expect_equal(directional$mean_gvi_south, mean(c(1.7 / 100 + 0.001, 28.7 / 100 + 0.001)))
  testthat::expect_equal(directional$min_gvi_south, 1.7 / 100 + 0.001)
  testthat::expect_equal(directional$max_gvi_south, 28.7 / 100 + 0.001)
  testthat::expect_equal(directional$gvi_bottom_east, 1.7 / 100 + 0.002)
  testthat::expect_equal(directional$gvi_top_east, 28.7 / 100 + 0.002)
  grouped_building <- sf::st_sf(
    id = seq_len(2),
    group_id = c(1, 1),
    Height = c(10, 30),
    geometry = sf::st_sfc(
      make_square(3.0000, 0.0000, 3.0002, 0.0002),
      make_square(3.0002, 0.0000, 3.0004, 0.0002),
      crs = 4326
    )
  )
  grouped_gvi <- gloBFPr::get_bgvi(
    grouped_building,
    workers = 1,
    key = "fake-key",
    quiet = TRUE
  )
  testthat::expect_equal(nrow(grouped_gvi), 2)
  testthat::expect_equal(grouped_gvi$mean_gvi, rep(mean(c(1.7, 28.7) / 100), 2))
  testthat::expect_equal(grouped_gvi$top_gvi, rep(28.7 / 100, 2))
  testthat::expect_error(
    gloBFPr::get_bgvi(
      building,
      directions = "up",
      key = "fake-key",
      workers = 1,
      quiet = TRUE
    ),
    "invalid"
  )
  testthat::expect_error(
    gloBFPr::get_bgvi(
      building,
      datasource = "metachm",
      key = "fake-key",
      workers = 1,
      quiet = TRUE
    ),
    "matches multiple formal|unused argument"
  )
  testthat::expect_error(
    gloBFPr::get_bgvi(
      building,
      datasource_canopy_height = NULL,
      datasource_greenspace = NULL,
      key = "fake-key",
      workers = 1,
      quiet = TRUE
    ),
    "At least one"
  )
})

testthat::test_that("plot_bgvi_viewshed rejects removed datasource argument", {
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
    Height = 10,
    geometry = sf::st_sfc(make_square(3.0000, 0.0000, 3.0002, 0.0002), crs = 4326)
  )

  testthat::expect_error(
    gloBFPr::plot_bgvi_viewshed(
      building,
      datasource = "metachm",
      plot = FALSE,
      quiet = TRUE
    ),
    "matches multiple formal|no longer supported"
  )
})
