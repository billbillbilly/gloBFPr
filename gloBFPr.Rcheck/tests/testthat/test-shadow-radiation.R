testthat::test_that("shadow and radiation functions handle missing input", {
  testthat::expect_null(gloBFPr::get_shadow_footprint())
  testthat::expect_null(gloBFPr::get_shadow_height())
  testthat::expect_null(gloBFPr::get_radiation(solar_normal = 1, solar_diffuse = 1))
})

testthat::test_that("solar position validation is strict", {
  validate_azimuth_elevation <- getFromNamespace("validate_azimuth_elevation", "gloBFPr")
  testthat::expect_equal(
    validate_azimuth_elevation(list(90, 180), c(45, 30)),
    matrix(c(90, 180, 45, 30), ncol = 2, dimnames = list(NULL, c("azimuth", "elevation")))
  )
  testthat::expect_error(validate_azimuth_elevation(NULL, 45), "Provide either")
  testthat::expect_error(validate_azimuth_elevation(180, NA), "numeric")
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
    time_zone = "UTC"
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

testthat::test_that("svf returns a raster with lower values near buildings", {
  building <- sf::st_sf(
    Height = 20,
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

  result <- gloBFPr::svf(
    building,
    grid_res = 10,
    extent_buffer = 30,
    res_angle = 45,
    max_distance = 60,
    quiet = TRUE
  )

  testthat::expect_s4_class(result, "SpatRaster")
  vals <- terra::values(result, mat = FALSE)
  testthat::expect_true(all(vals[!is.na(vals)] >= 0 & vals[!is.na(vals)] <= 1))
  near <- terra::extract(result, matrix(c(15, 5), ncol = 2))[, 1]
  far <- terra::extract(result, matrix(c(35, 25), ncol = 2))[, 1]
  testthat::expect_lt(near, far)
})

testthat::test_that("get_shadow_height uses shadow_locations and rejects old location alias", {
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

  result <- gloBFPr::get_shadow_height(
    building,
    shadow_locations = point,
    azimuth = 90,
    elevation = 45,
    quiet = TRUE
  )
  testthat::expect_true(is.matrix(result))
  testthat::expect_error(
    gloBFPr::get_shadow_height(
      building,
      location = point,
      azimuth = 90,
      elevation = 45,
      quiet = TRUE
    ),
    "unused argument"
  )
})

testthat::test_that("get_shadow_height rejects unknown arguments", {
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
      unexpected = point,
      azimuth = 90,
      elevation = 45,
      quiet = TRUE
    ),
    "unused argument"
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

testthat::test_that("shadow footprints use overlap_shadow and reject old combine alias", {
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

  result <- gloBFPr::get_shadow_footprint(
    building,
    azimuth = solar_pos[, 1],
    elevation = solar_pos[, 2],
    overlap_shadow = TRUE,
    quiet = TRUE
  )
  testthat::expect_equal(result$solar_count, 2)
  testthat::expect_error(
    gloBFPr::get_shadow_footprint(
      building,
      azimuth = solar_pos[, 1],
      elevation = solar_pos[, 2],
      combine = TRUE,
      quiet = TRUE
    ),
    "unused argument"
  )
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

testthat::test_that("get_radiation plots canopy impact when canopy is supplied", {
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

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- gloBFPr::get_radiation(
    building,
    azimuth = 90,
    elevation = 45,
    solar_normal = 800,
    solar_diffuse = 100,
    canopy_height = chm,
    canopy_transmissivity = 0.1,
    grid_res = 10,
    plot = TRUE,
    quiet = TRUE
  )

  testthat::expect_s3_class(result, "sf")
  testthat::expect_true(all(c("direct", "diffuse", "total") %in% names(result)))
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
    plot_3d = TRUE,
    quiet = TRUE
  )

  testthat::expect_s3_class(result, "sf")
  testthat::expect_true(all(c("direct", "diffuse", "total") %in% names(result)))
  projection <- getFromNamespace("project_radiation_3d", "gloBFPr")(
    sf::st_coordinates(result)[, "X"],
    sf::st_coordinates(result)[, "Y"],
    result$z
  )
  testthat::expect_equal(length(projection$x), nrow(result))
})

testthat::test_that("get_radiation plots ground, facade, and roof surfaces in 2D", {
  building <- sf::st_sf(
    Height = c(10, 8),
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(
        c(0, 0,
          0, 10,
          10, 10,
          10, 0,
          0, 0),
        ncol = 2,
        byrow = TRUE
      ))),
      sf::st_polygon(list(matrix(
        c(30, 0,
          30, 10,
          40, 10,
          40, 0,
          30, 0),
        ncol = 2,
        byrow = TRUE
      ))),
      crs = 3857
    )
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
    ground = TRUE,
    ground_res = 10,
    plot = TRUE,
    quiet = TRUE
  )

  testthat::expect_s3_class(result, "sf")
  testthat::expect_true(all(c("ground", "facade", "roof") %in% result$surface))
})

testthat::test_that("get_radiation computes diffuse radiation from surface SVF", {
  buildings <- sf::st_sf(
    Height = c(20, 8),
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(
        c(0, 0,
          0, 20,
          20, 20,
          20, 0,
          0, 0),
        ncol = 2,
        byrow = TRUE
      ))),
      sf::st_polygon(list(matrix(
        c(35, 0,
          35, 15,
          50, 15,
          50, 0,
          35, 0),
        ncol = 2,
        byrow = TRUE
      )))
    ),
    crs = 3857
  )

  result <- gloBFPr::get_radiation(
    buildings,
    azimuth = 225,
    elevation = 45,
    solar_normal = 850,
    solar_diffuse = 120,
    grid_res = 10,
    svf_res_angle = 45,
    quiet = TRUE
  )

  testthat::expect_true(all(result$svf >= 0 & result$svf <= 1))
  testthat::expect_gt(length(unique(round(result$svf, 4))), 2)
  testthat::expect_gt(length(unique(round(result$diffuse, 4))), 2)
  testthat::expect_equal(result$diffuse, result$svf * 120, tolerance = 1e-8)
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

testthat::test_that("facade normals point outward for both CW and CCW polygons", {
  facade_grid_sf <- getFromNamespace("facade_grid_sf", "gloBFPr")

  make_building <- function(coords) {
    sf::st_sf(
      Height = 10,
      geometry = sf::st_sfc(
        sf::st_polygon(list(matrix(coords, ncol = 2, byrow = TRUE))),
        crs = 3857
      )
    )
  }

  # CCW exterior ring (positive signed area, GeoJSON/OGC standard)
  b_ccw <- make_building(c(
    0, 0,
    10, 0,
    10, 10,
    0, 10,
    0, 0
  ))
  # CW exterior ring (negative signed area, Shapefile/legacy convention)
  b_cw <- make_building(c(
    0, 0,
    0, 10,
    10, 10,
    10, 0,
    0, 0
  ))

  for (b in list(b_ccw, b_cw)) {
    grid <- facade_grid_sf(b, "Height", grid_res = 5, offset = 0.01)
    # All facade normals must be unit vectors in XY plane (nz = 0)
    norms <- sqrt(grid$nx^2 + grid$ny^2)
    testthat::expect_true(all(abs(norms - 1) < 1e-9))
    # Each normal must point AWAY from the building centroid
    centroid <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(b)))
    pts <- sf::st_coordinates(grid)
    outward <- (pts[, "X"] - centroid[, "X"]) * grid$nx +
               (pts[, "Y"] - centroid[, "Y"]) * grid$ny
    testthat::expect_true(all(outward > 0),
      info = paste("Some facade normals point inward for",
                   if (identical(b, b_ccw)) "CCW" else "CW", "polygon"))
  }
})

testthat::test_that("sun-facing facades receive positive direct radiation", {
  # East-facing facade should receive direct radiation when sun is in the east.
  # Uses a CCW and a CW building to catch orientation-dependent normal bugs.
  make_building <- function(coords) {
    sf::st_sf(
      Height = 10,
      geometry = sf::st_sfc(
        sf::st_polygon(list(matrix(coords, ncol = 2, byrow = TRUE))),
        crs = 3857
      )
    )
  }

  b_ccw <- make_building(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0))
  b_cw  <- make_building(c(0, 0, 0, 10, 10, 10, 10, 0, 0, 0))

  for (b in list(b_ccw, b_cw)) {
    result <- gloBFPr::get_radiation(
      b,
      azimuth = 90,
      elevation = 45,
      solar_normal = 800,
      solar_diffuse = 0,
      grid_res = 5,
      quiet = TRUE
    )
    facade <- result[result$surface == "facade", ]
    # East-facing facades (nx > 0.5) should have positive direct radiation
    east_facing <- facade[facade$nx > 0.5, ]
    testthat::expect_true(
      nrow(east_facing) > 0,
      info = "No east-facing facade points found"
    )
    testthat::expect_true(
      any(east_facing$direct > 0),
      info = paste("East-facing facades have zero direct radiation for",
                   if (identical(b, b_ccw)) "CCW" else "CW", "polygon")
    )
  }
})
