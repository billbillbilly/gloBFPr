testthat::test_that("prepare_openfoam_inputs creates a local metric case", {
  square <- sf::st_polygon(list(matrix(
    c(-84.50, 42.30,
      -84.50, 42.301,
      -84.499, 42.301,
      -84.499, 42.30,
      -84.50, 42.30),
    ncol = 2,
    byrow = TRUE
  )))
  buildings <- sf::st_sf(
    id = 1L,
    Height = 12,
    geometry = sf::st_sfc(square, crs = 4326)
  )
  case_dir <- tempfile("openfoam-case-")
  on.exit(unlink(case_dir, recursive = TRUE, force = TRUE), add = TRUE)

  result <- prepare_openfoam_inputs(
    case_dir = case_dir,
    buildings_list = list(poly = buildings, binary = NULL, graduated = NULL),
    include_fused_dsm = FALSE,
    include_tree_canopy = FALSE,
    include_morphology = FALSE,
    include_neighbors = FALSE,
    include_greenspace = FALSE,
    domain_buffer = 25,
    zmax_buffer = 20,
    quiet = TRUE
  )

  testthat::expect_equal(result$n_buildings, 1L)
  # Translation into a local engineering coordinate system intentionally drops
  # the source CRS; `result$crs` retains the projected source reference.
  testthat::expect_true(is.na(sf::st_crs(result$data$buildings)))
  testthat::expect_false(sf::st_is_longlat(result$crs))
  testthat::expect_equal(unname(result$origin["z"]), 0)
  testthat::expect_equal(result$domain$xmin, 0)
  testthat::expect_equal(result$domain$ymin, 0)
  testthat::expect_equal(result$domain$zmax, 32)
  testthat::expect_true(file.exists(result$files$building_gpkg))
  testthat::expect_true(file.exists(result$files$building_stl))
  testthat::expect_true(file.exists(result$files$metadata_rds))
  testthat::expect_match(readLines(result$files$building_stl, n = 1),
                         "^solid buildings$")
  testthat::expect_null(result$files$building_binary_raster)
  testthat::expect_null(result$files$building_height_raster)
})

testthat::test_that("prepare_openfoam_inputs validates case parameters", {
  testthat::expect_error(
    prepare_openfoam_inputs(tempfile(), cell_size = 0),
    "positive"
  )
  testthat::expect_error(
    prepare_openfoam_inputs(tempfile(), landcover_year = 2019),
    "2020 or 2021"
  )
})
