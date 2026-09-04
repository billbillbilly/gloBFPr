# ===========================================================================
# Terrain, canopy and terrain-based building geometry for OpenFOAM cases
# ===========================================================================
# prepare_openfoam_inputs() downloads terrain, building and canopy layers but
# writes them only as rasters; the CFD case has always been meshed on a flat
# blockMesh floor with a single extruded-footprint building STL.  This file
# closes that gap.  Everything here is optional and degrades gracefully: pass
# what you have, and the parts you cannot supply are skipped.
#
# The three layers are treated differently on purpose:
#
#   TERRAIN  -> solid geometry (STL).  Slope is what makes cold air drain, so
#               it must be in the mesh, not in a roughness length.
#
#   BUILDINGS-> solid geometry (STL), re-based onto the terrain.  Kept as
#               extruded polygons rather than taken from the DSM: a 10-30 m
#               DSM turns building walls into stair-steps, whereas footprint
#               extrusion gives clean vertical faces.
#
#   CANOPY   -> POROUS DRAG, not solid.  A tree crown passes and drags air; it
#               does not block it.  Meshing crowns as solid over-blocks the
#               flow.  This matches the intent already stated in
#               prepare_openfoam_inputs() ("continuous CHM ... needed for
#               porous-zone / drag-coefficient definitions in OpenFOAM") and in
#               get_roughness_raster(), which masks tree cover "so that
#               porous-zone drag is not accumulated on top of a roughness
#               penalty".
#
# Terrain is not taken from the fused DSM directly.  The fused DSM is
# terrain + buildings + canopy, so meshing it as ground underneath a separate
# building STL would double-count every building.  Because get_fused_dsm()
# builds it as
#
#     surface    <- max(chm, building_height)
#     ground_dsm <- dem + surface
#     dsm        <- over buildings: base_elev(centroid) + height
#
# the bare earth is recoverable exactly:
#
#     dem = fused_dsm - max(building_height, canopy_height, 0)
#
# over ground (surface = 0), over canopy (surface = chm) and over buildings
# (dsm - height = base_elev).  No extra API call and no change to
# get_fused_dsm() is needed.
# ===========================================================================


#' Recover bare-earth terrain from the fused DSM
#'
#' @param fused_dsm SpatRaster or path. Fused DSM (terrain + buildings + canopy)
#'   in the case's local coordinate system.
#' @param building_height SpatRaster or path or NULL. Building heights above
#'   ground.
#' @param canopy_height SpatRaster or path or NULL. Continuous canopy height
#'   model (CHM).
#' @param smooth_iter Integer. Passes of 3x3 focal-median smoothing applied to
#'   the recovered terrain.  Subtracting rasters that were resampled onto
#'   different grids leaves small seams at building and canopy edges; a light
#'   smooth removes them without flattening real slope.  Default 1.
#' @param quiet Logical.
#'
#' @return A single-layer SpatRaster named \code{"terrain_m"}.
#' @noRd
foam_derive_terrain <- function(fused_dsm,
                                building_height = NULL,
                                canopy_height   = NULL,
                                smooth_iter     = 1L,
                                quiet           = FALSE) {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required.", call. = FALSE)

  as_rast <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "SpatRaster")) return(x)
    if (is.character(x) && file.exists(x)) return(terra::rast(x))
    NULL
  }

  dsm <- as_rast(fused_dsm)
  if (is.null(dsm))
    stop("`fused_dsm` must be a SpatRaster or an existing file path.",
         call. = FALSE)
  if (terra::nlyr(dsm) > 1L) dsm <- dsm[[1L]]

  bh  <- as_rast(building_height)
  chm <- as_rast(canopy_height)

  # Resample the subtrahends onto the DSM grid.  bilinear for the continuous
  # canopy surface; buildings use bilinear too because the DSM itself was built
  # from a bilinearly-resampled stack.
  align <- function(r) {
    if (is.null(r)) return(NULL)
    if (terra::nlyr(r) > 1L) r <- r[[1L]]
    if (!terra::same.crs(r, dsm) && !is.na(terra::crs(r)) && nzchar(terra::crs(r)))
      r <- terra::project(r, dsm, method = "bilinear")
    else if (!terra::compareGeom(r, dsm, stopOnError = FALSE))
      r <- terra::resample(r, dsm, method = "bilinear")
    terra::ifel(is.na(r), 0, r)
  }

  bh_a  <- align(bh)
  chm_a <- align(chm)

  surface <- if (!is.null(bh_a) && !is.null(chm_a)) {
    terra::ifel(chm_a > bh_a, chm_a, bh_a)
  } else if (!is.null(bh_a)) {
    bh_a
  } else if (!is.null(chm_a)) {
    chm_a
  } else {
    if (!isTRUE(quiet))
      warning("Neither `building_height` nor `canopy_height` supplied: the ",
              "fused DSM will be used as terrain, so buildings and canopy are ",
              "baked into the ground surface. Expect double-counted buildings ",
              "if a separate building STL is also meshed.", call. = FALSE)
    NULL
  }

  dem <- if (is.null(surface)) dsm else dsm - terra::ifel(surface < 0, 0, surface)

  if (smooth_iter > 0L)
    for (i in seq_len(as.integer(smooth_iter)))
      dem <- terra::focal(dem, w = 3, fun = "median", na.policy = "all",
                          na.rm = TRUE)

  # Fill any residual holes so the surface is closed.
  if (any(is.na(terra::values(dem, mat = FALSE)))) {
    med <- stats::median(terra::values(dem, mat = FALSE), na.rm = TRUE)
    dem <- terra::ifel(is.na(dem), med, dem)
  }

  names(dem) <- "terrain_m"
  if (!isTRUE(quiet)) {
    rng <- range(terra::values(dem, mat = FALSE), na.rm = TRUE)
    message(sprintf(
      "  Terrain recovered: %.1f - %.1f m (relief %.1f m) at %g m resolution",
      rng[1], rng[2], diff(rng), terra::res(dem)[1]))
  }
  dem
}


#' Watertight solid mesh from a terrain raster
#'
#' Top surface follows the raster; four side walls drop to \code{base_z} and a
#' flat bottom closes the volume.  snappyHexMesh needs a closed surface to
#' decide inside from outside reliably - an open heightmap sheet leaves cells
#' below the terrain in the mesh.
#'
#' @param dem SpatRaster (single layer, local metric coordinates).
#' @param base_z Numeric. Elevation of the flat bottom; must be below
#'   \code{min(dem)}.
#' @return mesh list with \code{vertices} and \code{faces}.
#' @noRd
foam_terrain_solid_mesh <- function(dem, base_z) {
  z <- terra::as.matrix(dem, wide = TRUE)   # row 1 = north
  nr <- nrow(z); nc <- ncol(z)
  if (nr < 2 || nc < 2)
    stop("Terrain raster needs at least 2x2 cells.", call. = FALSE)
  if (anyNA(z)) z[is.na(z)] <- stats::median(z, na.rm = TRUE)
  if (base_z >= min(z))
    stop(sprintf("`base_z` (%g) must be below the terrain minimum (%g).",
                 base_z, min(z)), call. = FALSE)

  xs <- terra::xFromCol(dem, seq_len(nc))
  ys <- terra::yFromRow(dem, seq_len(nr))

  # Top vertices, row-major (row 1 = north)
  vid_top <- matrix(seq_len(nr * nc), nrow = nr, byrow = TRUE)
  v_top   <- cbind(rep(xs, times = nr), rep(ys, each = nc), as.vector(t(z)))

  # Bottom vertices: same x/y, flat z. Only the boundary ring is used, but
  # carrying the full grid keeps indexing trivial and costs only memory.
  off      <- nr * nc
  vid_bot  <- vid_top + off
  v_bot    <- cbind(v_top[, 1], v_top[, 2], base_z)

  vertices <- rbind(v_top, v_bot)

  faces <- vector("list", 0L)

  # --- top surface (normals +z) ---
  tf <- matrix(0L, nrow = 2L * (nr - 1L) * (nc - 1L), ncol = 3L)
  k <- 0L
  for (i in seq_len(nr - 1L)) {
    for (j in seq_len(nc - 1L)) {
      v00 <- vid_top[i, j];      v01 <- vid_top[i, j + 1L]
      v10 <- vid_top[i + 1L, j]; v11 <- vid_top[i + 1L, j + 1L]
      k <- k + 1L; tf[k, ] <- c(v00, v10, v11)
      k <- k + 1L; tf[k, ] <- c(v00, v11, v01)
    }
  }
  faces[[length(faces) + 1L]] <- tf

  # --- boundary ring, one consistent circuit -------------------------------
  # north row west->east, east column north->south, south row east->west,
  # west column south->north.  Traversed in this order every wall triangle
  # comes out with an outward normal, and every ring edge is the reverse of
  # the matching top-surface edge - which is what makes the solid watertight.
  ring <- c(vid_top[1L, ],
            vid_top[2:nr, nc],
            vid_top[nr, (nc - 1L):1L],
            if (nr > 2L) vid_top[(nr - 1L):2L, 1L] else integer(0))
  nring <- length(ring)
  nxt   <- c(ring[-1L], ring[1L])

  faces[[length(faces) + 1L]] <- rbind(
    cbind(ring, nxt,       nxt + off),
    cbind(ring, nxt + off, ring + off)
  )

  # --- flat bottom: fan over the same ring ---------------------------------
  # Fanning over the full ring (not just the four corners) means every wall
  # bottom edge is shared with exactly one bottom triangle.  The ring circuit
  # above is clockwise seen from +z, so the plain fan order already yields a
  # -z (outward) normal; reversing it points the bottom into the solid and
  # inflates the divergence-theorem volume by twice the bottom's flux.
  rb <- ring + off
  faces[[length(faces) + 1L]] <- cbind(rb[1L], rb[2:(nring - 1L)], rb[3:nring])

  faces <- do.call(rbind, faces)
  storage.mode(faces) <- "integer"
  list(vertices = vertices, faces = faces)
}


#' Canopy volume mesh: one closed box per canopy cell
#'
#' Used only to select cells for the leaf-area-density field; it never becomes
#' meshed geometry.  Resample the CHM before calling - one box per cell at 1 m
#' resolution produces millions of facets.
#'
#' @param dem SpatRaster. Terrain, already aligned to \code{chm}.
#' @param chm SpatRaster. Canopy height above ground.
#' @param min_height Numeric. Cells below this canopy height are ignored.
#' @return mesh list, or NULL when no cell qualifies.
#' @noRd
foam_canopy_boxes_mesh <- function(dem, chm, min_height = 2) {
  if (!terra::compareGeom(dem, chm, stopOnError = FALSE))
    dem <- terra::resample(dem, chm, method = "bilinear")

  h    <- terra::values(chm, mat = FALSE)
  base <- terra::values(dem, mat = FALSE)
  ok   <- which(is.finite(h) & h >= min_height & is.finite(base))
  if (length(ok) == 0L) return(NULL)

  xy <- terra::xyFromCell(chm, ok)
  rs <- terra::res(chm)
  hx <- rs[1] / 2; hy <- rs[2] / 2

  n  <- length(ok)
  # 8 vertices and 12 triangles per box
  vertices <- matrix(0, nrow = 8L * n, ncol = 3L)
  faces    <- matrix(0L, nrow = 12L * n, ncol = 3L)

  # Unit box triangulation (outward normals), vertex order:
  # 1..4 bottom (nw,ne,se,sw), 5..8 top (nw,ne,se,sw)
  box_faces <- rbind(
    c(1, 3, 2), c(1, 4, 3),        # bottom (-z)
    c(5, 6, 7), c(5, 7, 8),        # top (+z)
    c(1, 2, 6), c(1, 6, 5),        # north
    c(2, 3, 7), c(2, 7, 6),        # east
    c(3, 4, 8), c(3, 8, 7),        # south
    c(4, 1, 5), c(4, 5, 8)         # west
  )

  for (i in seq_len(n)) {
    x0 <- xy[i, 1] - hx; x1 <- xy[i, 1] + hx
    y0 <- xy[i, 2] - hy; y1 <- xy[i, 2] + hy
    z0 <- base[ok[i]];   z1 <- z0 + h[ok[i]]
    vi <- (i - 1L) * 8L
    vertices[vi + 1:8, ] <- rbind(
      c(x0, y1, z0), c(x1, y1, z0), c(x1, y0, z0), c(x0, y0, z0),
      c(x0, y1, z1), c(x1, y1, z1), c(x1, y0, z1), c(x0, y0, z1)
    )
    faces[(i - 1L) * 12L + 1:12, ] <- box_faces + vi
  }
  storage.mode(faces) <- "integer"
  list(vertices = vertices, faces = faces, n_boxes = n)
}


#' Rewrite the building STL with each building based on the terrain
#'
#' prepare_openfoam_inputs() extrudes every footprint from z = 0
#' (\code{polygon_to_stl_facets}: \code{z0 <- 0}).  On sloping ground that
#' leaves buildings floating above or buried inside the terrain.  Each building
#' is re-based to the terrain elevation at its centroid - the same convention
#' get_fused_dsm() uses when it flattens roofs - and extruded from slightly
#' below that so the footprint always intersects the ground surface cleanly.
#'
#' @param buildings sf polygons in local coordinates.
#' @param height_col Character. Height column name.
#' @param dem SpatRaster. Terrain in the same local coordinates.
#' @param path Character. Output STL path (ASCII, matching the existing writer).
#' @param bury Numeric. Metres to extend each building below its base so that
#'   snappyHexMesh gets a clean intersection rather than a coincident face.
#'   Default 2.
#' @return invisible list(path, base_elev)
#' @noRd
foam_write_buildings_stl_on_terrain <- function(buildings, height_col, dem,
                                                path, bury = 2, quiet = FALSE) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required.", call. = FALSE)

  buildings <- sf::st_make_valid(buildings)
  buildings <- suppressWarnings(sf::st_collection_extract(buildings, "POLYGON"))
  buildings <- suppressWarnings(sf::st_cast(buildings, "POLYGON"))
  if (nrow(buildings) == 0)
    stop("No valid building polygons available for STL export.", call. = FALSE)

  ctr <- suppressWarnings(
    sf::st_coordinates(sf::st_centroid(sf::st_geometry(buildings))))
  base_elev <- terra::extract(dem, ctr[, 1:2, drop = FALSE],
                              method = "bilinear")[, 1]
  dem_med <- stats::median(terra::values(dem, mat = FALSE), na.rm = TRUE)
  base_elev[!is.finite(base_elev)] <- dem_med

  tri_normal <- function(p1, p2, p3) {
    u <- p2 - p1; v <- p3 - p1
    n <- c(u[2] * v[3] - u[3] * v[2],
           u[3] * v[1] - u[1] * v[3],
           u[1] * v[2] - u[2] * v[1])
    len <- sqrt(sum(n^2))
    if (!is.finite(len) || len == 0) return(c(0, 0, 0))
    n / len
  }
  facet <- function(p1, p2, p3) {
    n <- tri_normal(p1, p2, p3)
    paste0(
      "  facet normal ", paste(format(n, scientific = FALSE), collapse = " "), "\n",
      "    outer loop\n",
      "      vertex ", paste(format(p1, scientific = FALSE), collapse = " "), "\n",
      "      vertex ", paste(format(p2, scientific = FALSE), collapse = " "), "\n",
      "      vertex ", paste(format(p3, scientific = FALSE), collapse = " "), "\n",
      "    endloop\n  endfacet\n")
  }

  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines("solid buildings", con)

  for (i in seq_len(nrow(buildings))) {
    geom   <- sf::st_geometry(buildings[i, ])[[1]]
    coords <- as.matrix(geom[[1]][, 1:2, drop = FALSE])
    if (nrow(coords) > 1 && all(coords[1, ] == coords[nrow(coords), ]))
      coords <- coords[-nrow(coords), , drop = FALSE]
    if (nrow(coords) < 3) next

    h  <- as.numeric(buildings[[height_col]][i])
    if (!is.finite(h) || h <= 0) next
    z0 <- base_elev[i] - bury          # extend below ground for a clean cut
    z1 <- base_elev[i] + h

    fc <- character(0)
    n  <- nrow(coords)
    for (a in seq_len(n)) {
      b  <- if (a == n) 1L else a + 1L
      p1 <- c(coords[a, 1], coords[a, 2], z0)
      p2 <- c(coords[b, 1], coords[b, 2], z0)
      p3 <- c(coords[b, 1], coords[b, 2], z1)
      p4 <- c(coords[a, 1], coords[a, 2], z1)
      fc <- c(fc, facet(p1, p2, p3), facet(p1, p3, p4))
    }
    for (a in 2:(n - 1L))                                    # roof
      fc <- c(fc, facet(c(coords[1, 1], coords[1, 2], z1),
                        c(coords[a, 1], coords[a, 2], z1),
                        c(coords[a + 1L, 1], coords[a + 1L, 2], z1)))
    for (a in 2:(n - 1L))                                    # base cap
      fc <- c(fc, facet(c(coords[1, 1], coords[1, 2], z0),
                        c(coords[a + 1L, 1], coords[a + 1L, 2], z0),
                        c(coords[a, 1], coords[a, 2], z0)))
    writeLines(fc, con)
  }
  writeLines("endsolid buildings", con)

  if (!isTRUE(quiet))
    message(sprintf(
      "  Buildings re-based onto terrain: %d polygons, base %.1f - %.1f m",
      nrow(buildings), min(base_elev), max(base_elev)))

  invisible(list(path = normalizePath(path, mustWork = FALSE),
                 base_elev = base_elev))
}


# ---------------------------------------------------------------------------
# Dictionary fragments
# ---------------------------------------------------------------------------

#' system/setFieldsDict - paint leaf area density inside the canopy volume
#' @noRd
make_foam_set_fields_dict <- function(canopy_stl_rel, outside_point,
                                      lad = 0.4, plant_cd = 0.2) {
  paste0(
    noc_foam_header("dictionary", "setFieldsDict"),
    "// Leaf area density and the drag coefficient are read by\n",
    "// atmPlantCanopyUSource / atmPlantCanopyTurbSource.  Drag is proportional\n",
    "// to LAD, so cells left at zero feel no canopy force - which is why the\n",
    "// fvOption can use selectionMode all and needs no cellSet.\n",
    "//\n",
    "// Both field names changed between releases:\n",
    "//     leaf area density   LAD (v2506)  ->  leafAreaDensity (v2606)\n",
    "//     drag coefficient    Cd  (v2506)  ->  plantCd         (v2606)\n",
    "// All four are set so the dictionary works on either.\n",
    "defaultFieldValues\n(\n",
    "    volScalarFieldValue LAD             0\n",
    "    volScalarFieldValue leafAreaDensity 0\n",
    "    volScalarFieldValue Cd              0\n",
    "    volScalarFieldValue plantCd         0\n",
    ");\n\n",
    "regions\n(\n",
    "    surfaceToCell\n    {\n",
    sprintf("        file            \"%s\";\n", canopy_stl_rel),
    sprintf("        outsidePoints   ( ( %g %g %g ) );\n",
            outside_point[1], outside_point[2], outside_point[3]),
    # includeCut TRUE matters. The canopy box cloud is built at the CHM's
    # resolution (a few metres) while background cells are ~base_cell_size, so
    # most boxes contain no cell centre at all: with includeCut false only the
    # cells whose centre happens to fall inside a box are selected, which on a
    # real case picked up ~1/4 of the canopy volume. Canopy edges are diffuse
    # anyway, so counting every cell the crown touches is the better error.
    "        includeCut      true;\n",
    "        includeInside   true;\n",
    "        includeOutside  false;\n",
    "        nearDistance    -1;\n",
    "        curvature       -100;\n\n",
    "        fieldValues\n        (\n",
    sprintf("            volScalarFieldValue LAD             %g\n", lad),
    sprintf("            volScalarFieldValue leafAreaDensity %g\n", lad),
    sprintf("            volScalarFieldValue Cd              %g\n", plant_cd),
    sprintf("            volScalarFieldValue plantCd         %g\n", plant_cd),
    "        );\n",
    "    }\n",
    ");\n"
  )
}

# (0/leafAreaDensity, 0/plantCd and the atmPlantCanopy* fvOptions entries are
#  written by prepare_foam_case() in openfoam_case.R, so that every 0/ field
#  and every fvOptions source is emitted from one place.)


# ===========================================================================
# Public orchestrator
# ===========================================================================

#' Build terrain, canopy and terrain-based building geometry for a case
#'
#' @description
#' Turns the rasters written by \code{\link{prepare_openfoam_inputs}} into the
#' STL geometry \code{\link{prepare_foam_case}} needs.  Everything is optional:
#' supply what you have and the rest is skipped.
#'
#' Terrain is recovered from the fused DSM rather than used directly, because
#' the fused DSM is terrain + buildings + canopy and meshing it as ground under
#' a separate building STL would double-count every building.  Given how
#' \code{get_fused_dsm()} builds it, bare earth comes back exactly as
#' \code{fused_dsm - max(building_height, canopy_height, 0)}.
#'
#' @param case_dir Character. OpenFOAM case directory.
#' @param fused_dsm,building_height,canopy_height SpatRaster or path or NULL.
#'   Typically \code{foam_inputs$files$fused_dsm},
#'   \code{...$building_height_raster}, \code{...$canopy_height_raster}.
#' @param buildings sf polygons in local coordinates, or NULL to leave the
#'   existing building STL alone.  Usually
#'   \code{readRDS(foam_inputs$files$buildings_rds)}.
#' @param height_col Character. Building height column.
#' @param domain Named list with xmin/xmax/ymin/ymax/zmin/zmax; used to crop.
#' @param terrain_res Numeric. Resolution (m) to resample terrain to before
#'   triangulating.  One vertex per cell, so a 1 m DEM over 1 km2 is a million
#'   vertices; default \code{max(5, base_cell_size / 2)}.
#' @param canopy_res Numeric. Resolution (m) for the canopy box cloud.
#'   Default 5.
#' @param min_tree_height Numeric. Canopy cells below this are ignored.
#'   Default 2.
#' @param base_cell_size Numeric. Background mesh cell size, used only to pick
#'   \code{terrain_res}. Default 10.
#' @param skirt Numeric. Metres the terrain solid extends below its minimum.
#'   Default 30.
#' @param quiet Logical.
#'
#' @return Invisibly, a list with \code{terrain_stl}, \code{canopy_stl},
#'   \code{building_stl}, \code{dem} (SpatRaster) and \code{dem_file} - feed
#'   these straight to \code{\link{prepare_foam_case}}.
#' @export
prepare_foam_geometry <- function(case_dir,
                                  fused_dsm       = NULL,
                                  building_height = NULL,
                                  canopy_height   = NULL,
                                  buildings       = NULL,
                                  height_col      = NULL,
                                  domain          = NULL,
                                  terrain_res     = NULL,
                                  canopy_res      = 5,
                                  min_tree_height = 2,
                                  base_cell_size  = 10,
                                  skirt           = 30,
                                  quiet           = FALSE) {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required.", call. = FALSE)

  msg <- function(...) if (!isTRUE(quiet)) message(...)
  tri_dir <- file.path(case_dir, "constant", "triSurface")
  dir.create(tri_dir, recursive = TRUE, showWarnings = FALSE)
  out <- list(terrain_stl = NULL, canopy_stl = NULL, building_stl = NULL,
              dem = NULL, dem_file = NULL)

  if (is.null(fused_dsm)) {
    msg("No `fused_dsm` supplied - no terrain geometry will be built.")
    return(invisible(out))
  }

  msg("Recovering bare terrain from the fused DSM ...")
  dem <- foam_derive_terrain(fused_dsm, building_height, canopy_height,
                             quiet = quiet)

  if (!is.null(domain)) {
    ext <- terra::ext(domain$xmin, domain$xmax, domain$ymin, domain$ymax)
    dem <- tryCatch(terra::crop(dem, ext), error = function(e) dem)
  }

  if (is.null(terrain_res)) terrain_res <- max(5, base_cell_size / 2)
  if (terra::res(dem)[1] < terrain_res) {
    fact <- max(1L, round(terrain_res / terra::res(dem)[1]))
    if (fact > 1L) {
      dem <- terra::aggregate(dem, fact = fact, fun = "mean", na.rm = TRUE)
      msg(sprintf("  Terrain resampled to %.1f m for meshing", terra::res(dem)[1]))
    }
  }

  out$dem      <- dem
  out$dem_file <- file.path(case_dir, "constant", "gloBFPr", "rasters",
                            "terrain_bare.tif")
  dir.create(dirname(out$dem_file), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(dem, out$dem_file, overwrite = TRUE)

  msg("Meshing terrain solid ...")
  zmin <- min(terra::values(dem, mat = FALSE), na.rm = TRUE)
  tm   <- foam_terrain_solid_mesh(dem, base_z = zmin - skirt)
  out$terrain_stl <- file.path(tri_dir, "terrain.stl")
  write_stl_binary(tm, out$terrain_stl, name = "terrain")
  msg(sprintf("  terrain.stl: %d vertices, %d facets",
              nrow(tm$vertices), nrow(tm$faces)))

  # -- Canopy ---------------------------------------------------------------
  if (!is.null(canopy_height)) {
    chm <- if (inherits(canopy_height, "SpatRaster")) canopy_height
           else if (is.character(canopy_height) && file.exists(canopy_height))
             terra::rast(canopy_height) else NULL
    if (!is.null(chm)) {
      if (terra::nlyr(chm) > 1L) chm <- chm[[1L]]
      if (!is.null(domain))
        chm <- tryCatch(terra::crop(chm, terra::ext(
          domain$xmin, domain$xmax, domain$ymin, domain$ymax)),
          error = function(e) chm)
      if (terra::res(chm)[1] < canopy_res) {
        fact <- max(1L, round(canopy_res / terra::res(chm)[1]))
        if (fact > 1L) chm <- terra::aggregate(chm, fact = fact, fun = "max",
                                               na.rm = TRUE)
      }
      cb <- foam_canopy_boxes_mesh(dem, chm, min_height = min_tree_height)
      if (is.null(cb)) {
        msg("  No canopy cells above min_tree_height - skipping canopy.")
      } else {
        out$canopy_stl <- file.path(tri_dir, "canopy.stl")
        write_stl_binary(cb, out$canopy_stl, name = "canopy")
        msg(sprintf("  canopy.stl: %d boxes, %d facets at %.1f m",
                    cb$n_boxes, nrow(cb$faces), terra::res(chm)[1]))
      }
    }
  }

  # -- Buildings re-based onto terrain --------------------------------------
  if (!is.null(buildings)) {
    if (is.null(height_col)) {
      cand <- intersect(c("Height", "height", "height_m", "mean_height"),
                        names(buildings))
      height_col <- if (length(cand)) cand[1] else
        stop("`height_col` could not be guessed; supply it.", call. = FALSE)
    }
    msg("Re-basing buildings onto terrain ...")
    out$building_stl <- file.path(tri_dir, "buildings.stl")
    foam_write_buildings_stl_on_terrain(buildings, height_col, dem,
                                        out$building_stl, quiet = quiet)
  }

  invisible(out)
}
