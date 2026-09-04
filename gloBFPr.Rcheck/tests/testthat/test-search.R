testthat::test_that("runs correctly", {
  buildings <- testthat::expect_warning(
    search_3dglobdf(bbox=c(-135.459881,81.000000,-133.526287,82.000000),
                    crop = FALSE),
    "No buildings were found",
    fixed = TRUE
  )

  testthat::expect_null(buildings)
})

testthat::test_that("search_3dglobdf accepts GlobalBuildingAtlas data source", {
  make_square_coords <- function(xmin, ymin, xmax, ymax) {
    list(matrix(
      c(xmin, ymin,
        xmin, ymax,
        xmax, ymax,
        xmax, ymin,
        xmin, ymin),
      ncol = 2,
      byrow = TRUE
    ))
  }

  mocked_read_gba_buildings <- function(bbox, quiet = TRUE) {
    sf::st_sf(
      Height = 7.5,
      .source_id = "gba-source-1",
      geometry = sf::st_sfc(
        sf::st_multipolygon(list(
          make_square_coords(-84.480, 45.640, -84.479, 45.641),
          make_square_coords(-84.478, 45.640, -84.477, 45.641)
        )),
        crs = 4326
      )
    )
  }

  namespace <- asNamespace("gloBFPr")
  old_read_gba_buildings <- get("read_gba_buildings", envir = namespace)
  unlockBinding("read_gba_buildings", namespace)
  assign("read_gba_buildings", mocked_read_gba_buildings, envir = namespace)
  lockBinding("read_gba_buildings", namespace)
  on.exit({
    unlockBinding("read_gba_buildings", namespace)
    assign("read_gba_buildings", old_read_gba_buildings, envir = namespace)
    lockBinding("read_gba_buildings", namespace)
  }, add = TRUE)

  buildings <- gloBFPr::search_3dglobdf(
    bbox = c(-84.485519, 45.636118, -84.462774, 45.650639),
    data_source = "gba",
    keep_source_id = TRUE,
    crop = FALSE
  )

  testthat::expect_s3_class(buildings, "sf")
  testthat::expect_named(buildings, c("Height", "geometry", "source_id", "group_id", "id"))
  testthat::expect_equal(as.character(sf::st_geometry_type(buildings)), "POLYGON")
  testthat::expect_equal(nrow(buildings), 1)
  testthat::expect_equal(buildings$Height, 7.5)
  testthat::expect_equal(buildings$source_id, "gba-source-1")
  testthat::expect_equal(buildings$group_id, 1)
  testthat::expect_equal(buildings$id, 1)
})

testthat::test_that("building geometry normalization preserves rows", {
  make_square_coords <- function(xmin, ymin, xmax, ymax) {
    list(matrix(
      c(xmin, ymin,
        xmin, ymax,
        xmax, ymax,
        xmax, ymin,
        xmin, ymin),
      ncol = 2,
      byrow = TRUE
    ))
  }

  buildings <- sf::st_sf(
    Height = c(4, 6),
    geometry = sf::st_sfc(
      sf::st_multipolygon(list(make_square_coords(0, 0, 1, 1))),
      sf::st_multipolygon(list(
        make_square_coords(2, 0, 3, 1),
        make_square_coords(4, 0, 5, 1)
      )),
      crs = 4326
    )
  )

  normalized <- getFromNamespace("normalize_building_geometries", "gloBFPr")(buildings)

  testthat::expect_s3_class(normalized, "sf")
  testthat::expect_equal(nrow(normalized), 2)
  testthat::expect_equal(unique(as.character(sf::st_geometry_type(normalized))), "POLYGON")
  testthat::expect_s3_class(sf::st_geometry(normalized), "sfc_POLYGON")
  testthat::expect_equal(normalized$Height, c(4, 6))
})

testthat::test_that("internal source id is optional in returned data", {
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
    Height = 4,
    .source_id = "a",
    geometry = sf::st_sfc(
      make_square(0, 0, 1, 1),
      crs = 4326
    )
  )

  cleaned <- getFromNamespace("drop_internal_building_columns", "gloBFPr")(buildings)
  kept <- getFromNamespace("drop_internal_building_columns", "gloBFPr")(
    buildings,
    keep_source_id = TRUE
  )

  testthat::expect_false(".source_id" %in% names(cleaned))
  testthat::expect_false(".source_id" %in% names(kept))
  testthat::expect_equal(kept$source_id, "a")
})

testthat::test_that("GBA reader falls back to bbox tiling", {
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

  calls <- 0L
  mocked_read_gba_bbox <- function(bbox) {
    calls <<- calls + 1L
    if (calls == 1L) {
      stop("HTTP error code : 504", call. = FALSE)
    }
    bb <- sf::st_bbox(bbox)
    sf::st_sf(
      id = paste0("tile-", calls),
      height = calls,
      geometry = sf::st_sfc(
        make_square(bb[["xmin"]], bb[["ymin"]], bb[["xmin"]] + 0.001, bb[["ymin"]] + 0.001),
        crs = 4326
      )
    )
  }

  namespace <- asNamespace("gloBFPr")
  old_read_gba_bbox <- get("read_gba_bbox", envir = namespace)
  unlockBinding("read_gba_bbox", namespace)
  assign("read_gba_bbox", mocked_read_gba_bbox, envir = namespace)
  lockBinding("read_gba_bbox", namespace)
  on.exit({
    unlockBinding("read_gba_bbox", namespace)
    assign("read_gba_bbox", old_read_gba_bbox, envir = namespace)
    lockBinding("read_gba_bbox", namespace)
  }, add = TRUE)

  result <- getFromNamespace("read_gba_buildings", "gloBFPr")(
    sf::st_as_sfc(sf::st_bbox(c(xmin = -83.1, ymin = 42.3, xmax = -83.0, ymax = 42.4), crs = 4326)),
    quiet = TRUE
  )

  testthat::expect_gt(calls, 1L)
  testthat::expect_s3_class(result, "sf")
  testthat::expect_true(all(c("Height", ".source_id", "geometry") %in% names(result)))
})

testthat::test_that("assign_building_group_id groups touching polygons", {
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
    Height = c(4, 5, 6),
    geometry = sf::st_sfc(
      make_square(0, 0, 1, 1),
      make_square(1, 0, 2, 1),
      make_square(4, 0, 5, 1),
      crs = 3857
    )
  )

  grouped <- getFromNamespace("assign_building_group_id", "gloBFPr")(buildings)

  testthat::expect_equal(grouped$group_id, c(1L, 1L, 2L))
})

testthat::test_that("search_3dglobdf validates data_source", {
  testthat::expect_error(
    gloBFPr::search_3dglobdf(
      bbox = c(-84.485519, 45.636118, -84.462774, 45.650639),
      data_source = "other"
    ),
    'Invalid data_source specified. Use "GBF" or "GBA".',
    fixed = TRUE
  )
})
