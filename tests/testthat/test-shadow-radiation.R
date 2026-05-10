testthat::test_that("shadow and radiation functions handle missing input", {
  testthat::expect_null(gloBFPr::get_shadow_footprint())
  testthat::expect_null(gloBFPr::get_shadow_height())
  testthat::expect_null(gloBFPr::get_radiation(solar_normal = 1, solar_diffuse = 1))
})

testthat::test_that("solar position validation is strict", {
  validate_solar_pos <- getFromNamespace("validate_solar_pos", "gloBFPr")
  validate_azimuth_elevation <- getFromNamespace("validate_azimuth_elevation", "gloBFPr")
  testthat::expect_equal(
    validate_solar_pos(data.frame(az = 180, elev = 45, extra = 1)),
    matrix(c(180, 45), nrow = 1, dimnames = list(NULL, c("azimuth", "elevation")))
  )
  testthat::expect_equal(
    validate_azimuth_elevation(list(90, 180), c(45, 30)),
    matrix(c(90, 180, 45, 30), ncol = 2, dimnames = list(NULL, c("azimuth", "elevation")))
  )
  testthat::expect_error(validate_solar_pos(matrix(180, ncol = 1)), "two columns")
  testthat::expect_error(validate_solar_pos(matrix(c(180, NA), ncol = 2)), "numeric")
  testthat::expect_error(validate_azimuth_elevation(c(90, 180), 45), "same length")
})

testthat::test_that("solar time uses time zone and takes precedence", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 4326)
  )
  resolve_solar_inputs <- getFromNamespace("resolve_solar_inputs", "gloBFPr")

  result <- resolve_solar_inputs(
    building,
    azimuth = 1,
    elevation = 1,
    solar_time = list("2026-06-21 12:00:00", "2026-06-21 13:00:00"),
    time_zone = "UTC",
    solar_pos = matrix(c(1, 1), ncol = 2)
  )

  testthat::expect_equal(nrow(result), 2)
  testthat::expect_false(all(result[, "azimuth"] == 1 & result[, "elevation"] == 1))
  testthat::expect_error(
    resolve_solar_inputs(building, solar_time = "2026-06-21 12:00:00"),
    "Both `solar_time` and `time_zone`"
  )
  testthat::expect_error(
    resolve_solar_inputs(building, solar_time = "2026-06-21 12:00:00", time_zone = c("UTC", "MST")),
    "`time_zone` must be a single"
  )
  testthat::expect_error(
    resolve_solar_inputs(
      building,
      solar_time = as.POSIXct("2026-06-21 12:00:00", tz = "UTC"),
      time_zone = "UTC"
    ),
    "`solar_time` must be a character"
  )
})

testthat::test_that("get_shadow_height accepts deprecated location alias", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )
  point <- sf::st_as_sf(data.frame(x = 5, y = 5), coords = c("x", "y"), crs = 3857)

  testthat::expect_warning(
    result <- gloBFPr::get_shadow_height(
      building,
      location = point,
      azimuth = 90,
      elevation = 45,
      quiet = TRUE
    ),
    "deprecated"
  )
  testthat::expect_true(is.matrix(result))
})

testthat::test_that("get_shadow_height rejects both location names", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )
  point <- sf::st_as_sf(data.frame(x = 5, y = 5), coords = c("x", "y"), crs = 3857)

  testthat::expect_error(
    gloBFPr::get_shadow_height(
      building,
      shadow_locations = point,
      location = point,
      azimuth = 90,
      elevation = 45,
      quiet = TRUE
    ),
    "only one"
  )
})

testthat::test_that("get_shadow_height does not count a building as shading itself", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )
  points <- sf::st_as_sf(
    data.frame(
      id = c("roof", "shadow"),
      x = c(5, -5),
      y = c(5, 5)
    ),
    coords = c("x", "y"),
    crs = 3857
  )

  result <- gloBFPr::get_shadow_height(
    building,
    shadow_locations = points,
    azimuth = 90,
    elevation = 45,
    quiet = TRUE
  )

  testthat::expect_true(is.na(result[points$id == "roof", 1]))
  testthat::expect_true(result[points$id == "shadow", 1] > 0)
})

testthat::test_that("shadow footprints support multiple solar positions", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )
  solar_pos <- rbind(c(90, 45), c(180, 35), c(270, 25))

  result <- gloBFPr::get_shadow_footprint(
    building,
    azimuth = solar_pos[, 1],
    elevation = solar_pos[, 2],
    quiet = TRUE
  )

  testthat::expect_s3_class(result, "sf")
  testthat::expect_equal(sort(unique(result$solar_index)), 1:3)
  testthat::expect_equal(nrow(result), 3)
})

testthat::test_that("shadow footprints use solar_time as sun_id and can plot", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(-83.001, 42.001,
        -83.001, 42.002,
        -83.000, 42.002,
        -83.000, 42.001,
        -83.001, 42.001),
      ncol = 2,
      byrow = TRUE
    ))), crs = 4326)
  )
  solar_time <- c("2026-06-21 09:00:00", "2026-06-21 12:00:00")

  result <- gloBFPr::get_shadow_footprint(
    building,
    solar_time = solar_time,
    time_zone = "America/Detroit",
    quiet = TRUE
  )

  testthat::expect_equal(sort(unique(result$sun_id)), sort(solar_time))

  shadow_plot_cols <- getFromNamespace("shadow_plot_cols", "gloBFPr")
  cols <- shadow_plot_cols(result, plot_overlap_gradient = TRUE)
  testthat::expect_equal(unname(cols), rep(grDevices::adjustcolor("grey20", alpha.f = 0.22), 2))

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- gloBFPr::get_shadow_footprint(
    building,
    solar_time = solar_time,
    time_zone = "America/Detroit",
    plot = TRUE,
    plot_overlap_gradient = TRUE,
    quiet = TRUE
  )
  testthat::expect_s3_class(plotted, "sf")
})

testthat::test_that("shadow footprints can combine overlapping solar positions", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )
  solar_pos <- rbind(c(90, 45), c(180, 35), c(270, 25))

  result <- gloBFPr::get_shadow_footprint(
    building,
    azimuth = solar_pos[, 1],
    elevation = solar_pos[, 2],
    overlap_shadow = TRUE,
    quiet = TRUE
  )

  testthat::expect_s3_class(result, "sf")
  testthat::expect_equal(result$shadow_source, "building")
  testthat::expect_equal(result$solar_count, 3)
})

testthat::test_that("shadow footprints accept deprecated combine alias", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )
  solar_pos <- rbind(c(90, 45), c(180, 35))

  testthat::expect_warning(
    result <- gloBFPr::get_shadow_footprint(
      building,
      solar_pos = solar_pos,
      combine = TRUE,
      quiet = TRUE
    ),
    "deprecated"
  )
  testthat::expect_equal(result$solar_count, 2)
})

testthat::test_that("canopy shadows reduce direct radiation", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )

  chm <- terra::rast(
    xmin = -30, xmax = 30,
    ymin = -20, ymax = 20,
    resolution = 5,
    crs = "EPSG:3857"
  )
  terra::values(chm) <- 0
  chm[terra::cellFromXY(chm, matrix(c(20, 2.5), ncol = 2))] <- 30

  solar_pos <- matrix(c(90, 45), ncol = 2)
  radiation_plain <- gloBFPr::get_radiation(
    building,
    azimuth = solar_pos[, 1],
    elevation = solar_pos[, 2],
    solar_normal = 800,
    solar_diffuse = 100,
    grid_res = 10,
    quiet = TRUE
  )
  radiation_canopy <- gloBFPr::get_radiation(
    building,
    azimuth = solar_pos[, 1],
    elevation = solar_pos[, 2],
    solar_normal = 800,
    solar_diffuse = 100,
    canopy_height = chm,
    canopy_transmissivity = 0.1,
    grid_res = 10,
    quiet = TRUE
  )

  roof_plain <- radiation_plain$total[radiation_plain$surface == "roof"]
  roof_canopy <- radiation_canopy$total[radiation_canopy$surface == "roof"]

  testthat::expect_lt(min(roof_canopy), max(roof_plain))
})

testthat::test_that("get_radiation can plot and still returns sf", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- gloBFPr::get_radiation(
    building,
    azimuth = 90,
    elevation = 45,
    solar_normal = 800,
    solar_diffuse = 100,
    grid_res = 10,
    plot = TRUE,
    quiet = TRUE
  )

  testthat::expect_s3_class(result, "sf")
  testthat::expect_true(all(c("direct", "diffuse", "total") %in% names(result)))
})

testthat::test_that("shadow footprints include canopy obstacles", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )

  chm <- terra::rast(
    xmin = -30, xmax = 30,
    ymin = -20, ymax = 20,
    resolution = 5,
    crs = "EPSG:3857"
  )
  terra::values(chm) <- 0
  chm[terra::cellFromXY(chm, matrix(c(20, 2.5), ncol = 2))] <- 30

  result <- gloBFPr::get_shadow_footprint(
    building,
    azimuth = 90,
    elevation = 45,
    canopy_height = chm,
    quiet = TRUE
  )

  testthat::expect_s3_class(result, "sf")
  testthat::expect_setequal(result$shadow_source, c("building", "canopy"))
  canopy_bbox <- sf::st_bbox(result[result$shadow_source == "canopy", ])
  testthat::expect_lte(unname(canopy_bbox["xmin"]), -7)

  canopy_cells_sf <- getFromNamespace("canopy_cells_sf", "gloBFPr")
  canopy_height_cols <- getFromNamespace("canopy_height_cols", "gloBFPr")
  order_canopy_shadows <- getFromNamespace("order_canopy_shadows", "gloBFPr")
  shadow_source_plot_style <- getFromNamespace("shadow_source_plot_style", "gloBFPr")
  canopy <- list(
    xy = matrix(c(20, 2.5), ncol = 2),
    height = 30,
    cell_size = 5
  )
  canopy_cells <- canopy_cells_sf(canopy, crs = sf::st_crs(building))
  testthat::expect_s3_class(canopy_cells, "sf")
  testthat::expect_equal(nrow(canopy_cells), 1)
  testthat::expect_equal(
    length(unique(canopy_height_cols(c(5, 30)))),
    2
  )
  ordered_canopy <- order_canopy_shadows(
    sf::st_sf(
      canopy_height = c(30, 5),
      geometry = sf::st_sfc(
        sf::st_point(c(0, 0)),
        sf::st_point(c(1, 1)),
        crs = sf::st_crs(building)
      )
    )
  )
  testthat::expect_equal(ordered_canopy$canopy_height, c(5, 30))
  source_style <- shadow_source_plot_style(result, plot_overlap_gradient = TRUE)
  testthat::expect_equal(source_style$fill[["canopy"]], source_style$fill[["building"]])

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- gloBFPr::get_shadow_footprint(
    building,
    azimuth = c(90, 120),
    elevation = c(45, 35),
    canopy_height = chm,
    plot = TRUE,
    plot_overlap_gradient = TRUE,
    quiet = TRUE
  )
  testthat::expect_s3_class(plotted, "sf")
  testthat::expect_setequal(plotted$shadow_source, c("building", "canopy"))
})

testthat::test_that("shadow functions can retrieve canopy and DEM internally", {
  building <- sf::st_sf(
    Height = 10,
    geometry = sf::st_sfc(sf::st_polygon(list(matrix(
      c(0, 0,
        0, 10,
        10, 10,
        10, 0,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ))), crs = 3857)
  )

  mocked_get_chm <- function(bbox, min_height, datasource = "metachm") {
    testthat::expect_true(tolower(datasource) %in% c("metachm", "ethchm"))
    testthat::expect_true(all(abs(bbox[c(1, 3)]) <= 180))
    testthat::expect_true(all(abs(bbox[c(2, 4)]) <= 90))
    chm <- terra::rast(
      xmin = -30, xmax = 30,
      ymin = -20, ymax = 20,
      resolution = 5,
      crs = "EPSG:3857"
    )
    terra::values(chm) <- 0
    chm[terra::cellFromXY(chm, matrix(c(20, 2.5), ncol = 2))] <- min_height + 20
    list(chm, terra::ifel(chm >= min_height, 1, 0))
  }
  mocked_get_dem <- function(bbox, key) {
    dem <- mocked_get_chm(bbox, 2)[[1]]
    terra::values(dem) <- 100
    dem
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

  result <- gloBFPr::get_shadow_footprint(
    building,
    azimuth = 90,
    elevation = 45,
    datasource_canopy_height = "ethCHM",
    key = "test-key",
    min_tree_height = 2,
    quiet = TRUE
  )

  testthat::expect_setequal(result$shadow_source, c("building", "canopy"))
})

testthat::test_that("dsmSearch bounding boxes are normalized to EPSG:4326", {
  as_wgs84_bbox_vector <- getFromNamespace("as_wgs84_bbox_vector", "gloBFPr")

  projected_bbox <- sf::st_as_sfc(
    sf::st_bbox(
      c(xmin = 276000, ymin = 4683000, xmax = 277000, ymax = 4684000),
      crs = 32617
    )
  )
  bbox <- as_wgs84_bbox_vector(projected_bbox)

  testthat::expect_named(bbox, c("xmin", "ymin", "xmax", "ymax"))
  testthat::expect_true(all(abs(bbox[c("xmin", "xmax")]) <= 180))
  testthat::expect_true(all(abs(bbox[c("ymin", "ymax")]) <= 90))
  testthat::expect_error(
    as_wgs84_bbox_vector(c(276000, 4683000, 277000, 4684000)),
    "EPSG:4326"
  )
})
