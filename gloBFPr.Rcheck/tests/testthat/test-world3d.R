# Helpers -------------------------------------------------------------------

unit_square <- function() {
  matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE)
}

l_shape <- function() {
  # Concave (L-shaped) polygon
  matrix(c(0, 0, 2, 0, 2, 1, 1, 1, 1, 2, 0, 2), ncol = 2, byrow = TRUE)
}

mesh_volume <- function(mesh) {
  # Signed volume via divergence theorem; positive when faces point outward
  v <- mesh$vertices; f <- mesh$faces
  vol <- 0
  for (i in seq_len(nrow(f))) {
    p1 <- v[f[i, 1], ]; p2 <- v[f[i, 2], ]; p3 <- v[f[i, 3], ]
    vol <- vol + (p1[1] * (p2[2] * p3[3] - p3[2] * p2[3]) -
                    p2[1] * (p1[2] * p3[3] - p3[2] * p1[3]) +
                    p3[1] * (p1[2] * p2[3] - p2[2] * p1[3])) / 6
  }
  unname(vol)
}

toy_buildings <- function() {
  sq <- function(x0, y0, s) {
    sf::st_polygon(list(matrix(
      c(x0, y0, x0 + s, y0, x0 + s, y0 + s, x0, y0 + s, x0, y0),
      ncol = 2, byrow = TRUE)))
  }
  # Coordinates near (500000, 4000000) in EPSG:32617 => lon ~ -81, zone 17
  sf::st_sf(
    Height = c(10, 25),
    id = 1:2,
    geometry = sf::st_sfc(sq(500000, 4000000, 10),
                          sq(500020, 4000000, 15), crs = 32617)
  )
}

# Triangulation and extrusion ------------------------------------------------

testthat::test_that("mesh_triangulate handles convex and concave polygons", {
  mesh_triangulate <- getFromNamespace("mesh_triangulate", "gloBFPr")

  tri_sq <- mesh_triangulate(unit_square())
  testthat::expect_equal(nrow(tri_sq), 2)

  tri_l <- mesh_triangulate(l_shape())
  testthat::expect_equal(nrow(tri_l), 4)  # n - 2 triangles

  # Triangulated area equals polygon area (concave case would fail with a fan)
  area <- sum(apply(tri_l, 1, function(idx) {
    p <- l_shape()[idx, ]
    abs((p[2, 1] - p[1, 1]) * (p[3, 2] - p[1, 2]) -
          (p[3, 1] - p[1, 1]) * (p[2, 2] - p[1, 2])) / 2
  }))
  testthat::expect_equal(area, 3)
})

testthat::test_that("mesh_extrude_polygon produces a closed prism with correct volume", {
  mesh_extrude_polygon <- getFromNamespace("mesh_extrude_polygon", "gloBFPr")

  m <- mesh_extrude_polygon(unit_square(), 0, 5)
  testthat::expect_equal(nrow(m$vertices), 8)
  testthat::expect_equal(nrow(m$faces), 12)  # 8 walls + 2 roof + 2 floor
  testthat::expect_equal(mesh_volume(m), 5, tolerance = 1e-8)

  # Concave prism volume: L-shape area 3, height 2 -> 6
  ml <- mesh_extrude_polygon(l_shape(), 0, 2)
  testthat::expect_equal(mesh_volume(ml), 6, tolerance = 1e-8)

  # Clockwise input ring must give the same (positive) volume
  mc <- mesh_extrude_polygon(unit_square()[4:1, ], 0, 5)
  testthat::expect_equal(mesh_volume(mc), 5, tolerance = 1e-8)

  testthat::expect_null(mesh_extrude_polygon(unit_square(), 0, 0))
  testthat::expect_null(mesh_extrude_polygon(unit_square(), 0, NA))
})

testthat::test_that("mesh_combine offsets face indices", {
  mesh_extrude_polygon <- getFromNamespace("mesh_extrude_polygon", "gloBFPr")
  mesh_combine <- getFromNamespace("mesh_combine", "gloBFPr")

  m1 <- mesh_extrude_polygon(unit_square(), 0, 1)
  m2 <- mesh_extrude_polygon(unit_square() + 5, 0, 2)
  cm <- mesh_combine(list(m1, m2, NULL))
  testthat::expect_equal(nrow(cm$vertices), 16)
  testthat::expect_equal(nrow(cm$faces), 24)
  testthat::expect_equal(max(cm$faces), 16)
  testthat::expect_equal(sort(unique(cm$face_groups)), c(1, 2))
  testthat::expect_equal(mesh_volume(cm), 3, tolerance = 1e-8)
})

# Terrain and trees -----------------------------------------------------------

testthat::test_that("mesh_terrain triangulates a raster grid", {
  testthat::skip_if_not_installed("terra")
  mesh_terrain <- getFromNamespace("mesh_terrain", "gloBFPr")

  r <- terra::rast(nrows = 4, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 4)
  terra::values(r) <- runif(20, 100, 110)
  m <- mesh_terrain(r)
  testthat::expect_equal(nrow(m$vertices), 20)
  testthat::expect_equal(nrow(m$faces), 2 * 3 * 4)
  testthat::expect_true(all(m$vertices[, 3] >= 100 & m$vertices[, 3] <= 110))

  # NA cells are skipped
  terra::values(r)[1] <- NA
  m2 <- mesh_terrain(r)
  testthat::expect_lt(nrow(m2$faces), nrow(m$faces))
  testthat::expect_false(anyNA(m2$vertices))
})

testthat::test_that("mesh_tree builds blocky arnis-style trunk and crown", {
  mesh_tree <- getFromNamespace("mesh_tree", "gloBFPr")
  tr <- mesh_tree(0, 0, 100, 12, variant = 0)
  testthat::expect_true(all(c("trunk", "crown") %in% names(tr)))
  testthat::expect_gt(nrow(tr$trunk$faces), 0)
  testthat::expect_gt(nrow(tr$crown$faces), 0)
  testthat::expect_true(all(tr$trunk$vertices[, 3] >= 100))
  # Apex cap top sits exactly at z_base + height
  testthat::expect_equal(max(tr$crown$vertices[, 3]), 112)
  # Voxel construction: every part is an axis-aligned box (12 triangles each)
  testthat::expect_equal(nrow(tr$trunk$faces) %% 12, 0)
  testthat::expect_equal(nrow(tr$crown$faces) %% 12, 0)
  # Canopy footprint stays within the ring radius (3 blocks + half a block)
  s <- 12 / 11  # oak_standard: canopy_top = 10 -> block size = height / 11
  testthat::expect_lte(max(abs(tr$crown$vertices[, 1])), 3.5 * s + 1e-9)

  # Deterministic for a given position, variants differ
  tr2 <- mesh_tree(0, 0, 100, 12, variant = 0)
  testthat::expect_identical(tr, tr2)
  spruce <- mesh_tree(0, 0, 100, 12, variant = 3)
  testthat::expect_false(nrow(spruce$crown$faces) == nrow(tr$crown$faces))

  testthat::expect_null(mesh_tree(0, 0, 0, -1))
})

testthat::test_that("tree variant follows CHM height", {
  mesh_tree <- getFromNamespace("mesh_tree", "gloBFPr")
  # At (0, 0) the position hash is 0
  short <- mesh_tree(0, 0, 0, 5)     # < 8 m  -> oak_compact (2)
  mid   <- mesh_tree(0, 0, 0, 12)    # < 14 m -> oak_bushy (1)
  tall  <- mesh_tree(0, 0, 0, 30)    # >= 22 m -> spruce_standard (3)
  testthat::expect_identical(short, mesh_tree(0, 0, 0, 5, variant = 2))
  testthat::expect_identical(mid, mesh_tree(0, 0, 0, 12, variant = 1))
  testthat::expect_identical(tall, mesh_tree(0, 0, 0, 30, variant = 3))
  # All heights still match the CHM value exactly
  testthat::expect_equal(max(short$crown$vertices[, 3]), 5)
  testthat::expect_equal(max(tall$crown$vertices[, 3]), 30)
  # Towering spruce reachable via position jitter
  towering <- mesh_tree(1, 0, 0, 30)
  testthat::expect_identical(towering, mesh_tree(1, 0, 0, 30, variant = 4))
})

testthat::test_that("mesh_voxel_columns merges runs and keeps volume", {
  mesh_voxel_columns <- getFromNamespace("mesh_voxel_columns", "gloBFPr")

  # 3 x 2 block slab, z 0..2 -> two row-merged boxes, volume 12
  cells <- expand.grid(x = c(0.5, 1.5, 2.5), y = c(0.5, 1.5))
  cells$z0 <- 0; cells$z1 <- 2
  m <- mesh_voxel_columns(cells, 1)
  testthat::expect_equal(nrow(m$faces), 2 * 12)  # 2 boxes
  testthat::expect_equal(mesh_volume(m), 12, tolerance = 1e-9)

  # Different column heights are not merged
  cells2 <- data.frame(x = c(0.5, 1.5), y = 0.5, z0 = 0, z1 = c(2, 3))
  m2 <- mesh_voxel_columns(cells2, 1)
  testthat::expect_equal(nrow(m2$faces), 2 * 12)
  testthat::expect_equal(mesh_volume(m2), 5, tolerance = 1e-9)

  testthat::expect_null(mesh_voxel_columns(
    data.frame(x = 0.5, y = 0.5, z0 = 0, z1 = 0), 1))
})

testthat::test_that("mesh_voxel_heightmap builds stepped terrain", {
  mesh_voxel_heightmap <- getFromNamespace("mesh_voxel_heightmap", "gloBFPr")

  # Single cell: top + 4 boundary walls = 5 quads = 10 triangles
  m1 <- mesh_voxel_heightmap(matrix(2), xs = 0.5, ys = 0.5, s = 1)
  testthat::expect_equal(nrow(m1$faces), 10)

  # 1 x 2 grid with a step: higher cell walls the shared edge
  m2 <- mesh_voxel_heightmap(matrix(c(3, 2), nrow = 1),
                             xs = c(0.5, 1.5), ys = 0.5, s = 1)
  testthat::expect_equal(nrow(m2$faces), 18)

  # Closed volume: heightmap walls + bottom cap enclose the exact block volume
  tops <- matrix(c(2, 2, 2, 3), nrow = 2, byrow = TRUE)
  m3 <- mesh_voxel_heightmap(tops, xs = c(0.5, 1.5), ys = c(1.5, 0.5),
                             s = 1, base = 0)
  bottom <- list(
    vertices = cbind(c(0, 2, 2, 0), c(0, 0, 2, 2), 0),
    faces = matrix(c(1L, 3L, 2L, 1L, 4L, 3L), ncol = 3, byrow = TRUE)
  )
  closed <- list(
    vertices = rbind(m3$vertices, bottom$vertices),
    faces = rbind(m3$faces, bottom$faces + nrow(m3$vertices))
  )
  testthat::expect_equal(mesh_volume(closed), 9, tolerance = 1e-9)

  # NA cells are skipped
  testthat::expect_null(mesh_voxel_heightmap(matrix(NA_real_), 0.5, 0.5, 1))
})

testthat::test_that("road widths and groups follow the arnis mapping", {
  road_half_width <- getFromNamespace("road_half_width", "gloBFPr")
  road_group <- getFromNamespace("road_group", "gloBFPr")

  testthat::expect_equal(
    road_half_width(c("motorway", "secondary", "tertiary", "residential",
                      "footway", "nonsense")),
    c(5, 4, 3, 2, 1, 2))
  testthat::expect_equal(
    road_group(c("primary", "footway", "path", "cycleway"),
               c(NA, NA, NA, NA)),
    c("roads", "sidewalks", "paths", "sidewalks"))
  # Overture sidewalk subclass wins over class
  testthat::expect_equal(road_group("footway", "sidewalk"), "sidewalks")
})

testthat::test_that("roads paint terrain surface classes from sf input", {
  toy_roads <- sf::st_sf(
    class = c("residential", "footway"),
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(499995, 4000005, 500040, 4000005),
                               ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(499995, 4000011, 500040, 4000011),
                               ncol = 2, byrow = TRUE)),
      crs = 32617
    )
  )
  out_dir <- tempfile("world3d")
  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    roads = toy_roads, format = "obj", out_dir = out_dir, quiet = TRUE
  )
  tm <- world$meshes$terrain
  # Terrain exists in flat mode when surface classes are present, with
  # roads/sidewalks as face classes of the terrain itself
  testthat::expect_false(is.null(tm))
  testthat::expect_true(all(c("roads", "sidewalks", "terrain") %in%
                              unique(tm$face_class)))
  # One continuous surface at z = 0 (no floating ribbons)
  testthat::expect_true(all(abs(tm$vertices[, 3]) < 1e-9))
  obj <- readLines(file.path(out_dir, "world.obj"))
  testthat::expect_true(any(obj == "o terrain"))
  testthat::expect_true(any(obj == "g roads"))
  testthat::expect_true(any(obj == "usemtl roads"))
  testthat::expect_true(any(obj == "g sidewalks"))
  # Road faces sit where the residential ribbon runs (y ~ 3..7 local)
  road_faces <- which(tm$face_class == "roads")
  vy <- tm$vertices[unique(as.vector(tm$faces[road_faces, ])), 2]
  testthat::expect_true(all(vy >= 1 & vy <= 9))

  # Voxel mode: classes on the block terrain
  world_vox <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    roads = toy_roads, all_vox = TRUE, format = "obj",
    out_dir = tempfile("world3d"), quiet = TRUE
  )
  testthat::expect_true("roads" %in% world_vox$meshes$terrain$face_class)

  testthat::expect_error(
    gloBFPr::get_3d_world(x = toy_buildings(), terrain = FALSE, canopy = NULL,
                          roads = 42, out_dir = tempfile()),
    "roads"
  )
})

testthat::test_that("bridge deck profile ramps from ground to deck", {
  bridge_deck_profile <- getFromNamespace("bridge_deck_profile", "gloBFPr")

  co <- cbind(seq(0, 100, by = 5), 0)
  g <- rep(0, nrow(co))
  z <- bridge_deck_profile(co, g, deck_z = 6, ramp = 20)
  testthat::expect_equal(z[1], 0)                 # starts at grade
  testthat::expect_equal(z[length(z)], 0)         # ends at grade
  testthat::expect_equal(z[11], 6)                # flat deck in the middle
  testthat::expect_true(all(diff(z[1:5]) > 0))    # ramps up
  testthat::expect_true(all(diff(z[17:21]) < 0))  # ramps down
  # No ramps requested: flat deck throughout
  testthat::expect_true(all(bridge_deck_profile(co, g, 6, ramp = 0) == 6))
})

testthat::test_that("bridge deck and supports build solid geometry", {
  mesh_bridge_deck <- getFromNamespace("mesh_bridge_deck", "gloBFPr")
  mesh_bridge_supports <- getFromNamespace("mesh_bridge_supports", "gloBFPr")

  co <- cbind(seq(0, 100, by = 5), 0)
  z <- rep(6, nrow(co))
  dk <- mesh_bridge_deck(co, z, half_width = 5, thickness = 0.6,
                         railing = 0)
  # 20 segments x 12 triangles, volume = 100 x 10 x 0.6
  testthat::expect_equal(nrow(dk$faces), 20 * 12)
  testthat::expect_equal(mesh_volume(dk), 600, tolerance = 1e-6)
  testthat::expect_equal(max(dk$vertices[, 3]), 6)
  testthat::expect_equal(min(dk$vertices[, 3]), 5.4)

  # Railings add geometry above the deck
  dk_r <- mesh_bridge_deck(co, z, half_width = 5, thickness = 0.6,
                           railing = 1)
  testthat::expect_gt(nrow(dk_r$faces), nrow(dk$faces))
  testthat::expect_equal(max(dk_r$vertices[, 3]), 7)

  # Pillars every 25 m, reaching from ground up to the deck underside
  sp <- mesh_bridge_supports(co, z - 0.6, rep(0, nrow(co)),
                             interval = 25, radius = 0.6)
  testthat::expect_equal(nrow(sp$faces) %% 12, 0)
  testthat::expect_equal(max(sp$vertices[, 3]), 5.4)
  testthat::expect_lte(min(sp$vertices[, 3]), 0)
  # Too little clearance: no pillars
  testthat::expect_null(
    mesh_bridge_supports(co, rep(0.5, nrow(co)), rep(0, nrow(co)),
                         interval = 25))
  # Span shorter than two intervals: no interior pillars, and no error
  short <- cbind(seq(0, 30, by = 5), 0)
  testthat::expect_null(
    mesh_bridge_supports(short, rep(6, nrow(short)), rep(0, nrow(short)),
                         interval = 25))
  testthat::expect_null(
    mesh_bridge_supports(short, rep(6, nrow(short)), rep(0, nrow(short)),
                         interval = 100))
})

testthat::test_that("bridges are elevated and excluded from the ground surface", {
  toy_roads <- sf::st_sf(
    class = c("residential", "primary"),
    is_bridge = c(FALSE, TRUE),
    is_tunnel = c(FALSE, FALSE),
    level = c(0L, 1L),
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(499995, 4000005, 500040, 4000005),
                               ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(499995, 4000030, 500060, 4000030),
                               ncol = 2, byrow = TRUE)),
      crs = 32617
    )
  )
  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    roads = toy_roads, format = "obj",
    out_dir = tempfile("world3d"), quiet = TRUE
  )
  testthat::expect_equal(world$n_bridges, 1)
  bm <- world$meshes$bridges
  testthat::expect_false(is.null(bm))
  # Deck reaches level (1) x bridge_level_height (6) above flat ground
  testthat::expect_equal(max(bm$vertices[, 3]), 6 + 1)  # incl. railing
  testthat::expect_false(is.null(world$meshes$bridge_supports))

  # The bridge is NOT painted on the terrain: road faces only along the
  # surface street at y ~ 5 local, not the bridge at y ~ 30
  tm <- world$meshes$terrain
  road_y <- tm$vertices[unique(as.vector(
    tm$faces[which(tm$face_class == "roads"), ])), 2]
  testthat::expect_lt(max(road_y), 20)

  # Tunnels are dropped from the surface entirely
  tun <- toy_roads
  tun$is_bridge <- FALSE
  tun$is_tunnel <- c(FALSE, TRUE)
  w2 <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    roads = tun, format = "obj", out_dir = tempfile("world3d"), quiet = TRUE)
  testthat::expect_equal(w2$n_bridges, 0)
  ty <- w2$meshes$terrain$vertices[unique(as.vector(
    w2$meshes$terrain$faces[
      which(w2$meshes$terrain$face_class == "roads"), ])), 2]
  testthat::expect_lt(max(ty), 20)

  # bridges = FALSE paints them flat instead
  w3 <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    roads = toy_roads, bridges = FALSE, format = "obj",
    out_dir = tempfile("world3d"), quiet = TRUE)
  testthat::expect_equal(w3$n_bridges, 0)
  testthat::expect_true(is.null(w3$meshes$bridges))
})

testthat::test_that("bridges are voxelized in all_vox mode", {
  testthat::skip_if_not_installed("terra")
  toy_roads <- sf::st_sf(
    class = "primary",
    is_bridge = TRUE, is_tunnel = FALSE, level = 1L,
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(499995, 4000030, 500060, 4000030),
                               ncol = 2, byrow = TRUE)), crs = 32617)
  )
  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    roads = toy_roads, all_vox = TRUE, vox_size = 1, format = "obj",
    out_dir = tempfile("world3d"), quiet = TRUE
  )
  bm <- world$meshes$bridges
  testthat::expect_false(is.null(bm))
  # Deck elevations snap to the block grid
  testthat::expect_true(all(abs(bm$vertices[, 3] -
                                  round(bm$vertices[, 3])) < 1e-9))
  testthat::expect_equal(nrow(bm$faces) %% 12, 0)
})

testthat::test_that("water flattens terrain and paints water plus sand fringe", {
  pond <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(matrix(
    c(500022, 4000002, 500032, 4000002, 500032, 4000012,
      500022, 4000012, 500022, 4000002),
    ncol = 2, byrow = TRUE))), crs = 32617))

  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    water = pond, sand_buffer = 3, format = "obj",
    out_dir = tempfile("world3d"), quiet = TRUE
  )
  tm <- world$meshes$terrain
  testthat::expect_true(all(c("water", "sand") %in% unique(tm$face_class)))
  # Water surface flattened 0.2 m below the flat banks (arnis method)
  water_faces <- which(tm$face_class == "water")
  wz <- tm$vertices[unique(as.vector(tm$faces[water_faces, ])), 3]
  testthat::expect_true(any(abs(wz + 0.2) < 1e-6))
  # Banks stay at 0
  ground_faces <- which(tm$face_class == "terrain")
  gz <- tm$vertices[unique(as.vector(tm$faces[ground_faces, ])), 3]
  testthat::expect_true(any(abs(gz) < 1e-9))

  testthat::expect_error(
    gloBFPr::get_3d_world(x = toy_buildings(), terrain = FALSE, canopy = NULL,
                          water = 42, out_dir = tempfile()),
    "water"
  )
})

testthat::test_that("tree detection window is resolution-independent", {
  testthat::skip_if_not_installed("terra")
  detect_tree_tops <- getFromNamespace("detect_tree_tops", "gloBFPr")

  # Two trees 8 m apart on a 1 m CHM (unequal heights, so the shorter one
  # can be suppressed by a wide window; equal heights both tie the local max)
  chm1 <- terra::rast(xmin = 0, xmax = 40, ymin = 0, ymax = 40,
                      resolution = 1, crs = "EPSG:32617")
  terra::values(chm1) <- 0
  chm1[terra::cellFromXY(chm1, cbind(10.5, 20.5))] <- 10
  chm1[terra::cellFromXY(chm1, cbind(18.5, 20.5))] <- 8
  testthat::expect_equal(nrow(detect_tree_tops(chm1, window = 5)), 2)
  # A window wider than the spacing keeps only the taller tree
  testthat::expect_equal(nrow(detect_tree_tops(chm1, window = 20)), 1)
  testthat::expect_equal(detect_tree_tops(chm1, window = 20)$height, 10)

  # Same 5 m window on a 5 m CHM must not collapse everything: cells, not
  # metres, would make this a 25 m window and drop the second tree
  chm5 <- terra::rast(xmin = 0, xmax = 40, ymin = 0, ymax = 40,
                      resolution = 5, crs = "EPSG:32617")
  terra::values(chm5) <- 0
  chm5[terra::cellFromXY(chm5, cbind(2.5, 22.5))] <- 10
  chm5[terra::cellFromXY(chm5, cbind(22.5, 22.5))] <- 8
  testthat::expect_equal(nrow(detect_tree_tops(chm5, window = 5)), 2)

  # min_height filter
  chm1[terra::cellFromXY(chm1, cbind(30.5, 30.5))] <- 1.5
  testthat::expect_equal(
    nrow(detect_tree_tops(chm1, window = 5, min_height = 2)), 2)

  # Window at or below the cell size: one tree per canopy cell, rather than
  # a forced 3-cell (here 15 m) minimum that would drop most of the canopy
  chm5[terra::cellFromXY(chm5, cbind(7.5, 22.5))] <- 9
  testthat::expect_equal(nrow(detect_tree_tops(chm5, window = 5)), 3)
})

testthat::test_that("greenspace keeps grass under tree canopy by default", {
  testthat::skip_if_not_installed("terra")
  green <- terra::rast(xmin = 499990, xmax = 500060,
                       ymin = 3999990, ymax = 4000060,
                       resolution = 1, crs = "EPSG:32617")
  terra::values(green) <- 0
  green[terra::cells(green, terra::ext(500000, 500020, 4000000, 4000010))] <- 1
  chm <- terra::rast(green)
  terra::values(chm) <- 0
  chm[terra::cells(chm, terra::ext(500000, 500010, 4000000, 4000010))] <- 10

  n_lawn <- function(w) {
    tm <- w$meshes$terrain
    length(which(tm$face_class == "greenspace"))
  }

  # Default: trees stand on grass, so the whole green patch stays lawn
  keep <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    canopy_height = chm, greenspace = green,
    format = "obj", out_dir = tempfile("world3d"), quiet = TRUE)
  # Opt-in: canopy cells cut out of the lawn
  cut <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    canopy_height = chm, greenspace = green,
    greenspace_under_canopy = FALSE,
    format = "obj", out_dir = tempfile("world3d"), quiet = TRUE)

  testthat::expect_gt(n_lawn(keep), n_lawn(cut))
  # Grass under the canopy half (local x < 10) only in the default mode
  kl <- keep$meshes$terrain
  lawn_x <- kl$vertices[unique(as.vector(
    kl$faces[which(kl$face_class == "greenspace"), ])), 1]
  testthat::expect_lt(min(lawn_x), 5)
})

testthat::test_that("max_trees cap warns and is reported", {
  testthat::skip_if_not_installed("terra")
  chm <- terra::rast(xmin = 499990, xmax = 500060,
                     ymin = 3999990, ymax = 4000060,
                     resolution = 1, crs = "EPSG:32617")
  terra::values(chm) <- 0
  for (p in list(c(500045.5, 4000045.5), c(500055.5, 4000045.5),
                 c(500045.5, 4000055.5))) {
    chm[terra::cellFromXY(chm, rbind(p))] <- c(20, 15, 10)[
      which(sapply(list(c(500045.5, 4000045.5), c(500055.5, 4000045.5),
                        c(500045.5, 4000055.5)),
                   function(q) isTRUE(all.equal(q, p))))]
  }
  testthat::expect_warning(
    world <- gloBFPr::get_3d_world(
      x = toy_buildings(), terrain = FALSE, canopy = NULL,
      canopy_height = chm, max_trees = 2,
      format = "obj", out_dir = tempfile("world3d"), quiet = TRUE),
    "trees detected"
  )
  testthat::expect_equal(world$n_trees, 2)
  testthat::expect_equal(world$n_trees_detected, 3)
})

testthat::test_that("greenspace class keeps lawns and excludes canopy", {
  testthat::skip_if_not_installed("terra")
  # Binary green raster covering x 500000..500020, y 4000000..4000010
  green <- terra::rast(xmin = 499990, xmax = 500060,
                       ymin = 3999990, ymax = 4000060,
                       resolution = 1, crs = "EPSG:32617")
  terra::values(green) <- 0
  green[terra::cells(green, terra::ext(500000, 500020, 4000000, 4000010))] <- 1
  # CHM: a 10 m tree patch inside the green area (x 500000..500010)
  chm <- terra::rast(green)
  terra::values(chm) <- 0
  chm[terra::cells(chm, terra::ext(500000, 500010, 4000000, 4000010))] <- 10

  out_dir <- tempfile("world3d")
  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    canopy_height = chm, greenspace = green,
    greenspace_under_canopy = FALSE,
    format = "obj", out_dir = out_dir, quiet = TRUE
  )
  tm <- world$meshes$terrain
  testthat::expect_true("greenspace" %in% unique(tm$face_class))
  # Canopy half excluded: lawn faces only where x >= ~10 local
  lawn_faces <- which(tm$face_class == "greenspace")
  lx <- tm$vertices[unique(as.vector(tm$faces[lawn_faces, ])), 1]
  testthat::expect_gte(min(lx), 8)
  obj <- readLines(file.path(out_dir, "world.obj"))
  testthat::expect_true(any(obj == "g greenspace"))
  testthat::expect_true(any(obj == "usemtl greenspace"))

  # Voxel mode: grass class on the block terrain
  world_vox <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    canopy_height = chm, greenspace = green, all_vox = TRUE,
    format = "obj", out_dir = tempfile("world3d"), quiet = TRUE
  )
  testthat::expect_true(
    "greenspace" %in% world_vox$meshes$terrain$face_class)

  # Greenspace uses greenSD map tiles only; the CHM source is rejected
  testthat::expect_error(
    gloBFPr::get_3d_world(x = toy_buildings(), terrain = FALSE, canopy = NULL,
                          greenspace = "metachm", out_dir = tempfile()),
    "map tiles"
  )
})

# Writers ---------------------------------------------------------------------

testthat::test_that("write_stl_binary writes a valid binary STL", {
  mesh_extrude_polygon <- getFromNamespace("mesh_extrude_polygon", "gloBFPr")
  write_stl_binary <- getFromNamespace("write_stl_binary", "gloBFPr")

  m <- mesh_extrude_polygon(unit_square(), 0, 5)
  path <- tempfile(fileext = ".stl")
  write_stl_binary(m, path, name = "test")

  testthat::expect_equal(file.size(path), 84 + 50 * nrow(m$faces))
  con <- file(path, "rb")
  on.exit(close(con))
  invisible(readBin(con, "raw", 80))
  n <- readBin(con, "integer", 1, size = 4, endian = "little")
  testthat::expect_equal(n, nrow(m$faces))
})

testthat::test_that("write_world_obj writes OBJ with objects, groups, and MTL", {
  mesh_extrude_polygon <- getFromNamespace("mesh_extrude_polygon", "gloBFPr")
  write_world_obj <- getFromNamespace("write_world_obj", "gloBFPr")

  m <- mesh_extrude_polygon(unit_square(), 0, 5)
  path <- tempfile(fileext = ".obj")
  layers <- list(buildings = list(
    mesh = m,
    groups = list(list(name = "building_1", faces = seq_len(nrow(m$faces)),
                       material = "buildings"))
  ))
  out <- write_world_obj(layers, path,
                         materials = list(buildings = c(0.7, 0.7, 0.7)))
  lines <- readLines(path)
  testthat::expect_equal(sum(startsWith(lines, "v ")), nrow(m$vertices))
  testthat::expect_equal(sum(startsWith(lines, "f ")), nrow(m$faces))
  testthat::expect_true(any(lines == "o buildings"))
  testthat::expect_true(any(lines == "g building_1"))
  testthat::expect_true(any(lines == "usemtl buildings"))
  testthat::expect_true(file.exists(sub("\\.obj$", ".mtl", path)))

  # Face indices must be within vertex count
  f_idx <- as.integer(unlist(strsplit(sub("^f ", "", lines[startsWith(lines, "f ")]), " ")))
  testthat::expect_true(all(f_idx >= 1 & f_idx <= nrow(m$vertices)))
})

testthat::test_that("write_world_obj supports UV-textured layers", {
  mesh_extrude_polygon <- getFromNamespace("mesh_extrude_polygon", "gloBFPr")
  write_world_obj <- getFromNamespace("write_world_obj", "gloBFPr")

  ground <- list(
    vertices = cbind(c(0, 10, 10, 0), c(0, 0, 10, 10), 0),
    faces = matrix(c(1L, 2L, 3L, 1L, 3L, 4L), ncol = 3, byrow = TRUE)
  )
  uv <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
  bld <- mesh_extrude_polygon(unit_square(), 0, 5)
  path <- tempfile(fileext = ".obj")
  write_world_obj(
    list(terrain = list(mesh = ground, material = "terrain_basemap", uv = uv),
         buildings = list(mesh = bld, material = "buildings")),
    path,
    materials = list(
      terrain_basemap = list(Kd = c(1, 1, 1), map_Kd = "basemap.jpg"),
      buildings = c(0.7, 0.7, 0.7)
    )
  )
  lines <- readLines(path)
  testthat::expect_equal(sum(startsWith(lines, "vt ")), 4)
  # Textured faces use v/vt indices; untextured faces stay plain
  testthat::expect_true(any(grepl("^f \\d+/\\d+ \\d+/\\d+ \\d+/\\d+$", lines)))
  testthat::expect_true(any(grepl("^f \\d+ \\d+ \\d+$", lines)))
  # Building v-indices continue after terrain, but vt-indices don't leak
  f_plain <- lines[grepl("^f \\d+ \\d+ \\d+$", lines)]
  idx <- as.integer(unlist(strsplit(sub("^f ", "", f_plain), " ")))
  testthat::expect_true(all(idx > 4))
  mtl <- readLines(sub("\\.obj$", ".mtl", path))
  testthat::expect_true(any(mtl == "map_Kd basemap.jpg"))
})

# Integration -----------------------------------------------------------------

testthat::test_that("get_3d_world exports a flat-ground scene from sf input", {
  out_dir <- tempfile("world3d")
  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    out_dir = out_dir, quiet = TRUE
  )

  testthat::expect_equal(world$n_buildings, 2)
  testthat::expect_true(file.exists(file.path(out_dir, "world.obj")))
  testthat::expect_true(file.exists(file.path(out_dir, "buildings.stl")))
  testthat::expect_true(file.exists(file.path(out_dir, "world_metadata.json")))

  # Local origin: coordinates start at 0
  vmin <- apply(world$meshes$buildings$vertices, 2, min)
  testthat::expect_equal(unname(vmin[1]), 0)
  testthat::expect_equal(unname(vmin[2]), 0)
  testthat::expect_equal(world$origin$epsg, 32617)

  # Total volume = 10*10*10 + 15*15*25
  testthat::expect_equal(mesh_volume(world$meshes$buildings),
                         1000 + 5625, tolerance = 1e-6)

  # Both buildings appear as OBJ groups
  lines <- readLines(file.path(out_dir, "world.obj"))
  testthat::expect_true(any(lines == "g building_1"))
  testthat::expect_true(any(lines == "g building_2"))

  # MTL only contains materials for layers actually in the scene
  mtl <- readLines(file.path(out_dir, "world.mtl"))
  testthat::expect_true(any(mtl == "newmtl buildings"))
  testthat::expect_false(any(grepl(
    "^newmtl (roads|sidewalks|paths|greenspace|terrain|canopy)", mtl)))
})

testthat::test_that("get_3d_world validates inputs", {
  testthat::expect_error(
    gloBFPr::get_3d_world(x = toy_buildings(), terrain = TRUE, key = NULL),
    "OpenTopography"
  )
  bad <- toy_buildings(); names(bad)[1] <- "h"
  testthat::expect_error(
    gloBFPr::get_3d_world(x = bad, terrain = FALSE, canopy = NULL),
    "Height"
  )
  testthat::expect_error(
    gloBFPr::get_3d_world(x = toy_buildings(), terrain = FALSE, canopy = NULL,
                          color_by = "nope", out_dir = tempfile()),
    "color_by"
  )
})

testthat::test_that("get_3d_world accepts dem and canopy_height rasters", {
  testthat::skip_if_not_installed("terra")
  dem <- terra::rast(xmin = 499990, xmax = 500060,
                     ymin = 3999990, ymax = 4000060,
                     resolution = 5, crs = "EPSG:32617")
  terra::values(dem) <- 100
  chm <- terra::rast(xmin = 499990, xmax = 500060,
                     ymin = 3999990, ymax = 4000060,
                     resolution = 1, crs = "EPSG:32617")
  terra::values(chm) <- 0
  chm[terra::cellFromXY(chm, cbind(500050.5, 4000050.5))] <- 10

  out_dir <- tempfile("world3d")
  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = TRUE, dem = dem,
    canopy_height = chm, canopy = NULL,
    out_dir = out_dir, quiet = TRUE
  )
  # No API key needed; buildings sit on the 100 m ground, walls 1 m below
  testthat::expect_equal(min(world$meshes$buildings$vertices[, 3]), 99)
  testthat::expect_false(is.null(world$meshes$terrain))
  # The single CHM tree is detected and stands on the terrain
  testthat::expect_equal(world$n_trees, 1)
  testthat::expect_equal(max(world$meshes$canopy_crowns$vertices[, 3]), 110)
})

testthat::test_that("get_3d_world all_vox produces a voxel world", {
  testthat::skip_if_not_installed("terra")
  out_dir <- tempfile("world3d")
  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    all_vox = TRUE, vox_size = 1, format = "obj",
    out_dir = out_dir, quiet = TRUE
  )
  testthat::expect_equal(world$n_buildings, 2)
  bm <- world$meshes$buildings
  # Every part is an axis-aligned box
  testthat::expect_equal(nrow(bm$faces) %% 12, 0)
  # Columns extend one block below flat ground; tops quantized to Height
  testthat::expect_equal(min(bm$vertices[, 3]), -1)
  testthat::expect_equal(max(bm$vertices[, 3]), 25)
  # Voxelized volume: 10x10x(10+1) + 15x15x(25+1)
  testthat::expect_equal(mesh_volume(bm), 100 * 11 + 225 * 26, tolerance = 0.05)
  # Flat slab ground plane exists
  testthat::expect_false(is.null(world$meshes$terrain))

  testthat::expect_error(
    gloBFPr::get_3d_world(x = toy_buildings(), terrain = FALSE, canopy = NULL,
                          all_vox = TRUE, vox_size = -1, out_dir = tempfile()),
    "vox_size"
  )
})

testthat::test_that("get_3d_world supports facade_palette", {
  out_dir <- tempfile("world3d")
  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    facade_palette = TRUE, format = "obj", out_dir = out_dir, quiet = TRUE
  )
  mtl <- readLines(file.path(out_dir, "world.mtl"))
  testthat::expect_true(any(grepl("^newmtl facade_", mtl)))
  obj <- readLines(file.path(out_dir, "world.obj"))
  # Deterministic assignment: id 1 -> facade_2, id 2 -> facade_3
  testthat::expect_true(any(obj == "usemtl facade_2"))
  testthat::expect_true(any(obj == "usemtl facade_3"))

  # Custom palette
  out_dir2 <- tempfile("world3d")
  gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    facade_palette = c("red", "blue"), format = "obj",
    out_dir = out_dir2, quiet = TRUE
  )
  mtl2 <- readLines(file.path(out_dir2, "world.mtl"))
  testthat::expect_equal(sum(grepl("^newmtl facade_", mtl2)), 2)
})

testthat::test_that("get_3d_world supports color_by and format subsetting", {
  out_dir <- tempfile("world3d")
  world <- gloBFPr::get_3d_world(
    x = toy_buildings(), terrain = FALSE, canopy = NULL,
    color_by = "Height", format = "obj", out_dir = out_dir, quiet = TRUE
  )
  testthat::expect_false(file.exists(file.path(out_dir, "buildings.stl")))
  mtl <- readLines(file.path(out_dir, "world.mtl"))
  testthat::expect_true(any(grepl("^newmtl bldg_ramp_", mtl)))
})
