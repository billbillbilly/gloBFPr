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
  # Translation into a local engineering coordinate system keeps the projected
  # metric CRS so downstream spatial outputs remain self-describing.
  testthat::expect_false(is.na(sf::st_crs(result$data$buildings)))
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

testthat::test_that("prepare_foam_case writes the wind-only buoyant PIMPLE case", {
  case_dir <- tempfile("foam-wind-case-")
  tri_dir <- file.path(case_dir, "constant", "triSurface")
  dir.create(tri_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(case_dir, recursive = TRUE, force = TRUE), add = TRUE)

  stl_file <- file.path(tri_dir, "buildings.stl")
  writeLines(c(
    "solid buildings",
    "facet normal 0 0 1",
    "outer loop",
    "vertex 10 10 0",
    "vertex 20 10 0",
    "vertex 10 20 10",
    "endloop",
    "endfacet",
    "endsolid buildings"
  ), stl_file)

  domain <- list(xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                 zmin = 0, zmax = 50)
  case <- prepare_foam_case(
    case_dir = case_dir,
    stl_file = stl_file,
    domain = domain,
    inlet_velocity = c(3, 3, 0),
    sim_hours = 0.01,
    quiet = TRUE
  )

  control <- readLines(file.path(case_dir, "system", "controlDict"),
                       warn = FALSE)
  turbulence <- readLines(file.path(case_dir, "constant",
                                    "turbulenceProperties"), warn = FALSE)
  fields <- list.files(file.path(case_dir, "0"))

  testthat::expect_true(any(grepl("application\\s+buoyantBoussinesqPimpleFoam",
                                  control)))
  testthat::expect_true(any(grepl("RASModel\\s+kOmegaSST", turbulence)))
  testthat::expect_setequal(fields,
                            c("U", "p_rgh", "T", "k", "omega", "nut",
                              "alphat"))
  testthat::expect_equal(unname(case$params$patch_roles),
                         c("inlet", "outlet", "inlet", "outlet"))
})

testthat::test_that("prepare_foam_case fails on explicitly supplied missing geometry", {
  case_dir <- tempfile("foam-wind-case-")
  tri_dir <- file.path(case_dir, "constant", "triSurface")
  dir.create(tri_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(case_dir, recursive = TRUE, force = TRUE), add = TRUE)

  stl_file <- file.path(tri_dir, "buildings.stl")
  writeLines(c(
    "solid buildings",
    "facet normal 0 0 1",
    "outer loop",
    "vertex 10 10 0",
    "vertex 20 10 0",
    "vertex 10 20 10",
    "endloop",
    "endfacet",
    "endsolid buildings"
  ), stl_file)

  domain <- list(xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                 zmin = 0, zmax = 50)

  testthat::expect_error(
    prepare_foam_case(case_dir, stl_file, domain,
                      terrain_stl = file.path(case_dir, "missing-terrain.stl"),
                      quiet = TRUE),
    "terrain_stl"
  )
  testthat::expect_error(
    prepare_foam_case(case_dir, stl_file, domain,
                      canopy_stl = file.path(case_dir, "missing-canopy.stl"),
                      quiet = TRUE),
    "canopy_stl"
  )
  testthat::expect_error(
    prepare_foam_case(case_dir, stl_file, domain,
                      terrain_dem = file.path(case_dir, "missing-dem.tif"),
                      quiet = TRUE),
    "terrain_dem"
  )
})

testthat::test_that("legacy OpenFOAM generators are not exported", {
  exports <- getNamespaceExports("gloBFPr")
  testthat::expect_false("prepare_openfoam_case" %in% exports)
  testthat::expect_false("prepare_nocturnal_case" %in% exports)
  testthat::expect_true("prepare_foam_case" %in% exports)
})

testthat::test_that("plot_foam_map can add orientation annotations", {
  testthat::skip_if_not_installed("ggplot2")

  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 400,
                   ymin = 0, ymax = 300, crs = "")
  terra::values(r) <- seq_len(terra::ncell(r))
  names(r) <- "U_mag"

  p <- plot_foam_map(r, scalebar = TRUE, scalebar_unit = "m",
                     north_arrow = TRUE)
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_gte(length(p$layers), 8)
  testthat::expect_no_error(ggplot2::ggplot_build(p))

  p_plain <- plot_foam_map(r, scalebar = FALSE, north_arrow = FALSE)
  testthat::expect_s3_class(p_plain, "ggplot")
  testthat::expect_no_error(ggplot2::ggplot_build(p_plain))
})
