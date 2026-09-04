make_noise_square <- function(xmin, ymin, xmax, ymax) {
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

testthat::test_that("infer_osm_traffic maps highway classes and OSM speed tags", {
  roads <- sf::st_sf(
    highway = c("primary", "residential", "service"),
    maxspeed = c("35 mph", NA, "20"),
    lanes = c("4", "2", NA),
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(0, 10, 100, 10), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(0, 20, 100, 20), ncol = 2, byrow = TRUE)),
      crs = 3857
    )
  )

  out <- gloBFPr::infer_osm_traffic(roads, use_maxspeed = TRUE, use_lanes = TRUE)

  testthat::expect_equal(out$highway_class, c("primary", "residential", "service"))
  testthat::expect_equal(round(out$speed_kmh[1]), 56)
  testthat::expect_equal(out$speed_kmh[3], 20)
  testthat::expect_equal(round(out$light_veh_h[1], 1), 949.9)
  testthat::expect_true(all(c("LV_D", "HGV_D", "MV_D", "WAV_D", "WBV_D", "LV_SPD_D", "HGV_SPD_N") %in% names(out)))
  testthat::expect_equal(out$MV_D, c(0, 0, 0))
  testthat::expect_equal(out$PVMT, c("NL08", "NL08", "NL08"))
})

testthat::test_that("OSM traffic defaults match NoiseModelling Import_OSM categories", {
  defaults <- gloBFPr::osm_noise_traffic_defaults()
  primary <- defaults[defaults$highway == "primary", ]
  motorway <- defaults[defaults$highway == "motorway", ]
  service <- defaults[defaults$highway == "service", ]

  testthat::expect_equal(motorway$speed_kmh, 130)
  testthat::expect_equal(service$speed_kmh, 30)
  testthat::expect_equal(round(primary$LV_D, 3), round((1 - 0.2) * 7124 / 12, 3))
  testthat::expect_equal(round(primary$HGV_D, 3), round(0.2 * 7124 / 12, 3))
  testthat::expect_equal(round(primary$LV_E, 3), round((1 - 0.15) * 1069 / 4, 3))
  testthat::expect_equal(round(primary$HGV_N, 3), round(0.1 * 712 / 8, 3))
})

testthat::test_that("prepare_noisemodelling_inputs creates core layers", {
  buildings <- sf::st_sf(
    id = 1,
    Height = 12,
    geometry = sf::st_sfc(make_noise_square(20, -5, 40, 15), crs = 3857)
  )
  roads <- sf::st_sf(
    highway = "residential",
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      crs = 3857
    )
  )
  greenspace <- sf::st_sf(
    geometry = sf::st_sfc(make_noise_square(60, -10, 90, 20), crs = 3857)
  )

  out <- gloBFPr::prepare_noisemodelling_inputs(
    x = buildings,
    roads = roads,
    greenspace = greenspace,
    receiver = "grid",
    resolution = 20
  )

  testthat::expect_s3_class(out$buildings, "sf")
  testthat::expect_s3_class(out$roads, "sf")
  testthat::expect_s3_class(out$ground, "sf")
  testthat::expect_s3_class(out$receivers, "sf")
  testthat::expect_true("HEIGHT" %in% names(out$buildings))
  testthat::expect_true(all(c("PK", "LV_D", "MV_D", "PVMT") %in% names(out$roads)))
  testthat::expect_true("PK" %in% names(out$receivers))
  testthat::expect_true(any(out$ground$G == 1))
  testthat::expect_true(nrow(out$receivers) > 0)
  testthat::expect_true("Z" %in% colnames(sf::st_coordinates(out$receivers)))
})

testthat::test_that("prepare_noisemodelling_inputs preserves building population as POP", {
  buildings <- sf::st_sf(
    id = 1:2,
    Height = c(12, 8),
    residents = c(21, 0),
    geometry = sf::st_sfc(
      make_noise_square(20, -5, 40, 15),
      make_noise_square(45, -5, 65, 15),
      crs = 3857
    )
  )
  roads <- sf::st_sf(
    highway = "residential",
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      crs = 3857
    )
  )

  out <- gloBFPr::prepare_noisemodelling_inputs(
    x = buildings,
    roads = roads,
    population_field = "residents",
    receiver = "grid",
    resolution = 20
  )

  testthat::expect_true("POP" %in% names(out$buildings))
  testthat::expect_equal(out$buildings$POP, c(21, 0))
})

testthat::test_that("prepare_noisemodelling_inputs uses get_pop output when requested", {
  buildings <- sf::st_sf(
    id = 1,
    Height = 12,
    pop_total = 42.5,
    geometry = sf::st_sfc(make_noise_square(20, -5, 40, 15), crs = 3857)
  )
  roads <- sf::st_sf(
    highway = "residential",
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      crs = 3857
    )
  )

  out <- gloBFPr::prepare_noisemodelling_inputs(
    x = buildings,
    roads = roads,
    population = TRUE,
    receiver = "grid",
    resolution = 20
  )

  testthat::expect_equal(out$buildings$POP, 42.5)
})

testthat::test_that("prepare_noisemodelling_inputs reports missing population field", {
  buildings <- sf::st_sf(
    id = 1,
    Height = 12,
    geometry = sf::st_sfc(make_noise_square(20, -5, 40, 15), crs = 3857)
  )
  roads <- sf::st_sf(
    highway = "residential",
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      crs = 3857
    )
  )

  testthat::expect_error(
    gloBFPr::prepare_noisemodelling_inputs(
      x = buildings,
      roads = roads,
      population_field = "residents",
      receiver = "grid",
      resolution = 20
    ),
    "`population_field` was not found"
  )
})

testthat::test_that("prepare_noisemodelling_inputs can defer roads to NoiseModelling OSM import", {
  buildings <- sf::st_sf(
    id = 1,
    Height = 12,
    geometry = sf::st_sfc(make_noise_square(20, -5, 40, 15), crs = 3857)
  )

  out <- gloBFPr::prepare_noisemodelling_inputs(
    x = buildings,
    roads = NULL,
    download_roads = FALSE,
    receiver = "grid",
    resolution = 20
  )

  testthat::expect_s3_class(out$buildings, "sf")
  testthat::expect_null(out$roads)
  testthat::expect_null(out$receivers)
  testthat::expect_s3_class(out$ground, "sf")
})

testthat::test_that("NoiseModelling input writer preserves 3D receivers", {
  buildings <- sf::st_sf(
    id = 1,
    Height = 12,
    geometry = sf::st_sfc(make_noise_square(20, -5, 40, 15), crs = 3857)
  )
  roads <- sf::st_sf(
    highway = "residential",
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      crs = 3857
    )
  )
  inputs <- gloBFPr::prepare_noisemodelling_inputs(
    x = buildings,
    roads = roads,
    receiver = "grid",
    resolution = 20
  )
  write_inputs <- getFromNamespace("write_noisemodelling_shapefiles", "gloBFPr")
  paths <- write_inputs(inputs, tempdir())

  testthat::expect_true(all(grepl("\\.geojson$", paths)))
  testthat::expect_true(file.exists(paths[["RECEIVERS"]]))
  receivers <- sf::st_read(paths[["RECEIVERS"]], quiet = TRUE)
  testthat::expect_true("Z" %in% colnames(sf::st_coordinates(receivers)))
})

testthat::test_that("NoiseModelling receiver levels join back to receiver geometries", {
  receivers <- sf::st_sf(
    PK = c(1L, 2L),
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0, 4)),
      sf::st_point(c(10, 0, 4)),
      crs = 3857
    )
  )
  levels <- data.frame(
    IDRECEIVER = c(1L, 1L, 2L, 2L),
    PERIOD = c("D", "N", "D", "N"),
    LAEQ = c(55, 47, 61, 50),
    LEQ = c(56, 48, 62, 51)
  )
  make_map <- getFromNamespace("receivers_level_to_noise_map", "gloBFPr")

  out <- make_map(levels, receivers)

  testthat::expect_s3_class(out, "sf")
  testthat::expect_equal(nrow(out), 2)
  testthat::expect_true(all(c("LAEQ_D", "LAEQ_N", "LEQ_D", "LEQ_N") %in% names(out)))
  testthat::expect_equal(out$LAEQ_D[match(2L, out$PK)], 61)
  testthat::expect_true("Z" %in% colnames(sf::st_coordinates(out)))
})

testthat::test_that("NoiseModelling no-result sentinel levels are mapped to NA", {
  receivers <- sf::st_sf(
    PK = 1L,
    geometry = sf::st_sfc(sf::st_point(c(0, 0, 4)), crs = 3857)
  )
  levels <- data.frame(
    IDRECEIVER = 1L,
    PERIOD = "DEN",
    HZ63 = -99,
    HZ125 = -99,
    HZ250 = -99,
    HZ500 = -99,
    HZ1000 = -99,
    HZ2000 = -99,
    HZ4000 = -99,
    HZ8000 = -99,
    LAEQ = -92.01287,
    LEQ = -89.9691
  )
  make_map <- getFromNamespace("receivers_level_to_noise_map", "gloBFPr")

  out <- make_map(levels, receivers)

  testthat::expect_true(is.na(out$LAEQ_DEN))
  testthat::expect_true(is.na(out$LEQ_DEN))
})

testthat::test_that("receiver noise map can become a contour surface", {
  receivers <- sf::st_sf(
    PK = seq_len(9),
    geometry = sf::st_sfc(
      lapply(seq_len(9), function(i) {
        xy <- expand.grid(x = c(0, 10, 20), y = c(0, 10, 20))[i, ]
        sf::st_point(c(xy$x, xy$y, 4))
      }),
      crs = 3857
    )
  )
  receivers$LAEQ_DEN <- c(42, 45, 48, 50, 55, 60, 62, 67, 72)
  make_surface <- getFromNamespace("noise_map_to_noise_surface", "gloBFPr")

  surface <- make_surface(
    receivers,
    inputs = list(receivers = receivers),
    period = "DEN",
    field = "LAEQ",
    resolution = 10,
    breaks = c(45, 55, 65)
  )

  testthat::expect_s3_class(surface, "gloBFPr_noise_surface")
  testthat::expect_s4_class(surface$raster, "SpatRaster")
  testthat::expect_s3_class(surface$bands, "sf")
  testthat::expect_true(all(c("noise_class", "level_min", "level_max", "label") %in% names(surface$bands)))
  testthat::expect_true(nrow(surface$bands) > 0)
})

testthat::test_that("NoiseModelling isophone legend labels use dB(A) ranges", {
  isophones <- sf::st_sf(
    ISOLVL = c(1L, 2L, 7L),
    geometry = sf::st_sfc(
      make_noise_square(0, 0, 1, 1),
      make_noise_square(1, 0, 2, 1),
      make_noise_square(2, 0, 3, 1),
      crs = 3857
    )
  )
  class_values <- getFromNamespace("isophone_class_values", "gloBFPr")(isophones)
  legend_info <- getFromNamespace("isophone_legend_info", "gloBFPr")(isophones, class_values)

  testthat::expect_equal(class_values, c(1L, 2L, 3L))
  testthat::expect_true(any(grepl("dB", legend_info$label, fixed = TRUE)))
  testthat::expect_true("< 35 dB(A)" %in% legend_info$label)
})

testthat::test_that("NoiseModelling zero-based isophone classes are mapped to dB(A) ranges", {
  isophones <- sf::st_sf(
    ISOLVL = c(0L, 1L, 10L),
    geometry = sf::st_sfc(
      make_noise_square(0, 0, 1, 1),
      make_noise_square(1, 0, 2, 1),
      make_noise_square(2, 0, 3, 1),
      crs = 3857
    )
  )
  class_values <- getFromNamespace("isophone_class_values", "gloBFPr")(isophones)
  legend_info <- getFromNamespace("isophone_legend_info", "gloBFPr")(isophones, class_values)

  testthat::expect_equal(class_values, c(1L, 2L, 3L))
  testthat::expect_true("< 35 dB(A)" %in% legend_info$label)
  testthat::expect_true("35-40 dB(A)" %in% legend_info$label)
  testthat::expect_true(">= 80 dB(A)" %in% legend_info$label)
  testthat::expect_false("0-1 dB(A)" %in% legend_info$label)
})

testthat::test_that("NoiseModelling text isophone classes are parsed", {
  isophones <- sf::st_sf(
    ISOLVL = c("class 1", "class 2", "class 7"),
    geometry = sf::st_sfc(
      make_noise_square(0, 0, 1, 1),
      make_noise_square(1, 0, 2, 1),
      make_noise_square(2, 0, 3, 1),
      crs = 3857
    )
  )
  class_values <- getFromNamespace("isophone_class_values", "gloBFPr")(isophones)
  legend_info <- getFromNamespace("isophone_legend_info", "gloBFPr")(isophones, class_values)

  testthat::expect_equal(class_values, c(1L, 2L, 3L))
  testthat::expect_true("< 35 dB(A)" %in% legend_info$label)
})

testthat::test_that("Java version requests are detected", {
  is_java_version_request <- getFromNamespace("is_java_version_request", "gloBFPr")
  parse_java_major_version <- getFromNamespace("parse_java_major_version", "gloBFPr")

  testthat::expect_true(is_java_version_request(17))
  testthat::expect_true(is_java_version_request("17"))
  testthat::expect_false(is_java_version_request("/usr/bin/java"))
  testthat::expect_equal(parse_java_major_version('java version "17.0.12"'), 17L)
  testthat::expect_equal(parse_java_major_version('openjdk version "1.8.0_392"'), 8L)
})

testthat::test_that("NoiseModelling WPS argument helpers expose acoustic controls", {
  args <- getFromNamespace("noise_wps_args_vector", "gloBFPr")(
    confDiffVertical = TRUE,
    confDiffHorizontal = FALSE,
    confHumidity = 70,
    confTemperature = 15,
    confFavourableOccurrencesDefault = getFromNamespace("format_favourable_occurrences", "gloBFPr")(0.5),
    confRaysName = NULL
  )

  testthat::expect_true("-confDiffVertical" %in% args)
  testthat::expect_true("true" %in% args)
  testthat::expect_true("-confDiffHorizontal" %in% args)
  testthat::expect_true("false" %in% args)
  testthat::expect_true("-confHumidity" %in% args)
  testthat::expect_true("-confFavourableOccurrencesDefault" %in% args)
  testthat::expect_false("-confRaysName" %in% args)
  testthat::expect_equal(
    getFromNamespace("format_favourable_occurrences", "gloBFPr")(c(0, rep(0.5, 15))),
    paste(c(0, rep(0.5, 15)), collapse = ",")
  )
})

testthat::test_that("prepare_noisemodelling_inputs can fetch roads and greenspace internally", {
  buildings <- sf::st_sf(
    id = 1,
    Height = 12,
    geometry = sf::st_sfc(make_noise_square(20, -5, 40, 15), crs = 3857)
  )
  mocked_roads <- function(x, overpass_url = NULL, timeout = NULL) {
    sf::st_sf(
      highway = "secondary",
      geometry = sf::st_sfc(
        sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
        crs = sf::st_crs(x)
      )
    )
  }
  mocked_greenspace <- function(bbox = NULL, buffer = NULL, type = NULL,
                                zoom = 17, year = NULL, min_tree_height = 2) {
    template <- terra::rast(
      xmin = 0,
      xmax = 100,
      ymin = -20,
      ymax = 30,
      resolution = 10,
      crs = sf::st_crs(buildings)$wkt
    )
    terra::values(template) <- NA
    template[terra::cellFromXY(template, matrix(c(70, 0), ncol = 2))] <- 1
    template
  }

  namespace <- asNamespace("gloBFPr")
  old_roads <- get("get_osm_roads_for_buildings", envir = namespace)
  old_greenspace <- get("get_greenspace", envir = namespace)
  unlockBinding("get_osm_roads_for_buildings", namespace)
  unlockBinding("get_greenspace", namespace)
  assign("get_osm_roads_for_buildings", mocked_roads, envir = namespace)
  assign("get_greenspace", mocked_greenspace, envir = namespace)
  lockBinding("get_osm_roads_for_buildings", namespace)
  lockBinding("get_greenspace", namespace)
  on.exit({
    unlockBinding("get_osm_roads_for_buildings", namespace)
    unlockBinding("get_greenspace", namespace)
    assign("get_osm_roads_for_buildings", old_roads, envir = namespace)
    assign("get_greenspace", old_greenspace, envir = namespace)
    lockBinding("get_osm_roads_for_buildings", namespace)
    lockBinding("get_greenspace", namespace)
  }, add = TRUE)

  out <- gloBFPr::prepare_noisemodelling_inputs(
    x = buildings,
    datasource_greenspace = "esri",
    receiver = "grid",
    resolution = 20
  )

  testthat::expect_s3_class(out$roads, "sf")
  testthat::expect_equal(out$roads$highway_class, "secondary")
  testthat::expect_true(any(out$ground$source == "greenspace"))
})

testthat::test_that("noise raster input bbox is valid for projected buildings", {
  buildings <- sf::st_sf(
    id = 1,
    Height = 12,
    geometry = sf::st_sfc(make_noise_square(20, -5, 40, 15), crs = 3857)
  )
  roads <- sf::st_sf(
    highway = "residential",
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      crs = 3857
    )
  )
  captured_bbox <- NULL
  mocked_get_chm <- function(bbox, min_height, datasource = "metachm") {
    captured_bbox <<- bbox
    template <- terra::rast(
      xmin = 0,
      xmax = 100,
      ymin = -20,
      ymax = 30,
      resolution = 10,
      crs = sf::st_crs(buildings)$wkt
    )
    terra::values(template) <- 0
    list(template, template)
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

  out <- gloBFPr::prepare_noisemodelling_inputs(
    x = buildings,
    roads = roads,
    datasource_canopy_height = "metachm",
    receiver = "grid",
    resolution = 20
  )

  testthat::expect_s3_class(out$buildings, "sf")
  testthat::expect_length(captured_bbox, 4)
  testthat::expect_false(anyNA(captured_bbox))
  testthat::expect_true(captured_bbox[["xmin"]] >= -180 && captured_bbox[["xmax"]] <= 180)
  testthat::expect_true(captured_bbox[["ymin"]] >= -90 && captured_bbox[["ymax"]] <= 90)
})

testthat::test_that("greenspace download passes unnamed bbox to greenSD", {
  captured_bbox <- NULL
  mocked_get_tile_green <- function(bbox, zoom, provider, year = NULL) {
    captured_bbox <<- bbox
    template <- terra::rast(
      xmin = -83.1,
      xmax = -83.0,
      ymin = 42.3,
      ymax = 42.4,
      resolution = 0.05,
      crs = "EPSG:4326"
    )
    terra::values(template) <- 1
    list(green = template)
  }

  namespace <- asNamespace("greenSD")
  old_get_tile_green <- get("get_tile_green", envir = namespace)
  unlockBinding("get_tile_green", namespace)
  assign("get_tile_green", mocked_get_tile_green, envir = namespace)
  lockBinding("get_tile_green", namespace)
  on.exit({
    unlockBinding("get_tile_green", namespace)
    assign("get_tile_green", old_get_tile_green, envir = namespace)
    lockBinding("get_tile_green", namespace)
  }, add = TRUE)

  bbox <- c(xmin = -83.1, ymin = 42.3, xmax = -83.0, ymax = 42.4)
  out <- getFromNamespace("get_greenspace", "gloBFPr")(
    bbox = bbox,
    type = "esri",
    zoom = 14
  )

  testthat::expect_s4_class(out, "SpatRaster")
  testthat::expect_length(captured_bbox, 4)
  testthat::expect_null(names(captured_bbox))
  testthat::expect_false(anyNA(captured_bbox))
})

testthat::test_that("greenspace download reports corrupt greenSD installs clearly", {
  mocked_get_tile_green <- function(...) {
    stop("lazy-load database '/path/greenSD.rdb' is corrupt: internal error 1 in R_decompress1 with libdeflate")
  }

  namespace <- asNamespace("greenSD")
  old_get_tile_green <- get("get_tile_green", envir = namespace)
  unlockBinding("get_tile_green", namespace)
  assign("get_tile_green", mocked_get_tile_green, envir = namespace)
  lockBinding("get_tile_green", namespace)
  on.exit({
    unlockBinding("get_tile_green", namespace)
    assign("get_tile_green", old_get_tile_green, envir = namespace)
    lockBinding("get_tile_green", namespace)
  }, add = TRUE)

  testthat::expect_error(
    getFromNamespace("fetch_greenspace_tile", "gloBFPr")(
      bbox = c(-83.1, 42.3, -83.0, 42.4),
      zoom = 14,
      provider = "esri"
    ),
    "installed `greenSD` package appears to be corrupt",
    fixed = TRUE
  )
})

testthat::test_that("get_noise_map reports missing NoiseModelling runner", {
  buildings <- sf::st_sf(
    id = 1,
    Height = 12,
    geometry = sf::st_sfc(make_noise_square(20, -5, 40, 15), crs = 3857)
  )
  roads <- sf::st_sf(
    highway = "residential",
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      crs = 3857
    )
  )

  testthat::expect_error(
    gloBFPr::get_noise_map(
      x = buildings,
      roads = roads,
      run = TRUE,
      nm_path = tempfile(),
      download_nm = FALSE
    ),
    "NoiseModelling headless runner was not found"
  )
})

testthat::test_that("get_noise_map requires one road source", {
  buildings <- sf::st_sf(
    id = 1,
    Height = 12,
    geometry = sf::st_sfc(make_noise_square(20, -5, 40, 15), crs = 3857)
  )
  roads <- sf::st_sf(
    highway = "residential",
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)),
      crs = 3857
    )
  )

  testthat::expect_error(
    gloBFPr::get_noise_map(
      x = buildings,
      roads = roads,
      osm_file = tempfile(fileext = ".osm.pbf"),
      run = FALSE
    ),
    "Use either `roads` or `osm_file`"
  )
})

testthat::test_that("Java version parsing supports legacy and modern version strings", {
  parse_java_major_version <- getFromNamespace("parse_java_major_version", "gloBFPr")
  java_class_file_major <- getFromNamespace("java_class_file_major", "gloBFPr")
  testthat::expect_equal(parse_java_major_version('java version "1.8.0_402"'), 8)
  testthat::expect_equal(parse_java_major_version('openjdk version "11.0.22" 2024-01-16'), 11)
  testthat::expect_equal(parse_java_major_version('openjdk version "17.0.10" 2024-01-16'), 17)
  testthat::expect_equal(parse_java_major_version('openjdk version "26.0.1" 2026-04-21'), 26)
  testthat::expect_equal(java_class_file_major(26), 70)
})

testthat::test_that("Road emission WPS script patch replaces brittle setLength calls", {
  patch_script <- getFromNamespace("patch_road_emission_wps_script", "gloBFPr")
  script <- tempfile(fileext = ".groovy")
  writeLines(c(
    "createTableQuery.setLength(createTableQuery.length() - 2)",
    "preparedInsertQuery.setLength(preparedInsertQuery.length() - 2)"
  ), script)

  patched <- patch_script(script, tempdir())
  patched_lines <- readLines(patched, warn = FALSE)

  testthat::expect_false(identical(script, patched))
  testthat::expect_true(any(grepl("createTableQuery.delete", patched_lines, fixed = TRUE)))
  testthat::expect_true(any(grepl("preparedInsertQuery.delete", patched_lines, fixed = TRUE)))
})
