# Internal triangle-mesh helpers shared by get_3d_world() (and reusable by the
# OpenFOAM STL pipeline). A "mesh" is a list with:
#   vertices : numeric matrix (n x 3), columns x, y, z
#   faces    : integer matrix (k x 3), 1-based vertex indices, CCW = outward

#' Ear-clipping triangulation of a simple polygon
#'
#' @param coords numeric matrix (n x 2) of an open ring (no repeated last
#' vertex), any orientation. Holes are not supported.
#' @return integer matrix (k x 3) of 1-based vertex indices, CCW orientation.
#' @noRd
mesh_triangulate <- function(coords) {
  n <- nrow(coords)
  if (n < 3) return(matrix(integer(0), ncol = 3))
  if (n == 3) return(matrix(1:3, ncol = 3))

  signed_area <- function(idx) {
    x <- coords[idx, 1]; y <- coords[idx, 2]
    sum(x * y[c(2:length(idx), 1)] - x[c(2:length(idx), 1)] * y) / 2
  }
  cross_z <- function(a, b, c) {
    (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
  }
  point_in_tri <- function(p, a, b, c) {
    d1 <- cross_z(a, b, p); d2 <- cross_z(b, c, p); d3 <- cross_z(c, a, p)
    !((d1 < 0 || d2 < 0 || d3 < 0) && (d1 > 0 || d2 > 0 || d3 > 0))
  }

  idx <- seq_len(n)
  # Work in CCW orientation
  reversed <- signed_area(idx) < 0
  if (reversed) idx <- rev(idx)

  tris <- matrix(integer(0), ncol = 3)
  guard <- 0L
  max_iter <- 2L * n * n

  while (length(idx) > 3 && guard < max_iter) {
    guard <- guard + 1L
    m <- length(idx)
    ear_found <- FALSE
    for (i in seq_len(m)) {
      i_prev <- idx[if (i == 1L) m else i - 1L]
      i_curr <- idx[i]
      i_next <- idx[if (i == m) 1L else i + 1L]
      a <- coords[i_prev, ]; b <- coords[i_curr, ]; c <- coords[i_next, ]
      # Reflex or degenerate vertex cannot be an ear
      if (cross_z(a, b, c) <= 0) next
      # No other polygon vertex may lie inside the candidate ear
      others <- setdiff(idx, c(i_prev, i_curr, i_next))
      inside <- FALSE
      for (j in others) {
        if (point_in_tri(coords[j, ], a, b, c)) { inside <- TRUE; break }
      }
      if (inside) next
      tris <- rbind(tris, c(i_prev, i_curr, i_next))
      idx <- idx[-i]
      ear_found <- TRUE
      break
    }
    if (!ear_found) {
      # Degenerate/self-intersecting ring: fall back to fan triangulation
      break
    }
  }

  if (length(idx) == 3) {
    tris <- rbind(tris, idx)
  } else if (length(idx) > 3) {
    # Fan fallback for the (possibly degenerate) remainder
    for (i in 2:(length(idx) - 1)) {
      tris <- rbind(tris, c(idx[1], idx[i], idx[i + 1]))
    }
  }
  storage.mode(tris) <- "integer"
  tris
}

#' Extrude a polygon ring into a closed prism mesh
#'
#' @param coords numeric matrix (n x 2), open ring.
#' @param z0,z1 numeric. Bottom and top elevation.
#' @return mesh list (vertices, faces) with outward-facing triangles.
#' @noRd
mesh_extrude_polygon <- function(coords, z0, z1) {
  n <- nrow(coords)
  if (n < 3 || !is.finite(z0) || !is.finite(z1) || z1 <= z0) {
    return(NULL)
  }
  # Ensure CCW so that wall winding faces outward
  sa <- sum(coords[, 1] * coords[c(2:n, 1), 2] -
              coords[c(2:n, 1), 1] * coords[, 2]) / 2
  if (sa < 0) coords <- coords[n:1, , drop = FALSE]

  bottom <- cbind(coords, z0)
  top    <- cbind(coords, z1)
  vertices <- rbind(bottom, top)  # 1..n bottom, n+1..2n top

  faces <- matrix(integer(0), ncol = 3)
  # Walls
  for (i in seq_len(n)) {
    j <- if (i == n) 1L else i + 1L
    faces <- rbind(faces,
                   c(i, j, n + j),
                   c(i, n + j, n + i))
  }
  # Caps
  cap <- mesh_triangulate(coords)
  if (nrow(cap) > 0) {
    faces <- rbind(faces, cap + n)              # roof, CCW from above = up
    faces <- rbind(faces, cap[, c(1, 3, 2), drop = FALSE])  # floor, flipped
  }
  storage.mode(faces) <- "integer"
  list(vertices = vertices, faces = faces)
}

#' Combine meshes into one, offsetting face indices
#' @param meshes list of mesh lists (NULLs are dropped).
#' @return single mesh, plus `face_groups`: integer vector mapping each face to
#' the index of its source mesh.
#' @noRd
mesh_combine <- function(meshes) {
  meshes <- Filter(Negate(is.null), meshes)
  if (length(meshes) == 0) return(NULL)
  n_vert <- vapply(meshes, function(m) nrow(m$vertices), integer(1))
  n_face <- vapply(meshes, function(m) nrow(m$faces), integer(1))
  offsets <- cumsum(c(0L, n_vert[-length(n_vert)]))
  vertices <- do.call(rbind, lapply(meshes, `[[`, "vertices"))
  faces <- do.call(rbind, Map(function(m, off) m$faces + off, meshes, offsets))
  storage.mode(faces) <- "integer"
  list(vertices = vertices, faces = faces,
       face_groups = rep(seq_along(meshes), n_face))
}

#' Triangulate a raster surface into a terrain mesh
#' @param r SpatRaster (single layer, projected CRS in metres).
#' @return mesh list with `face_cell` (linear raster cell index, row-major,
#' of each face's top-left cell), or NULL when fewer than 2x2 finite cells
#' exist.
#' @noRd
mesh_terrain <- function(r) {
  z <- terra::as.matrix(r, wide = TRUE)  # rows = raster rows (north first)
  nr <- nrow(z); nc <- ncol(z)
  if (is.null(nr) || nr < 2 || nc < 2) return(NULL)
  xs <- terra::xFromCol(r, seq_len(nc))
  ys <- terra::yFromRow(r, seq_len(nr))

  # Vertex per cell centre, row-major (row 1 = north)
  vid <- matrix(seq_len(nr * nc), nrow = nr, byrow = TRUE)
  vertices <- cbind(
    rep(xs, times = nr),
    rep(ys, each = nc),
    as.vector(t(z))
  )

  faces <- vector("list", (nr - 1) * (nc - 1))
  cells <- integer((nr - 1) * (nc - 1))
  k <- 0L
  for (i in seq_len(nr - 1)) {
    for (j in seq_len(nc - 1)) {
      v00 <- vid[i, j];     v01 <- vid[i, j + 1]
      v10 <- vid[i + 1, j]; v11 <- vid[i + 1, j + 1]
      zs <- vertices[c(v00, v01, v10, v11), 3]
      if (anyNA(zs)) next
      k <- k + 1L
      # Row 1 is north (larger y): this winding faces +z
      faces[[k]] <- rbind(c(v00, v10, v11), c(v00, v11, v01))
      cells[k] <- v00  # top-left cell, row-major linear index
    }
  }
  if (k == 0L) return(NULL)
  face_cell <- rep(cells[seq_len(k)], each = 2L)
  faces <- do.call(rbind, faces[seq_len(k)])
  storage.mode(faces) <- "integer"

  # Drop unused (NA) vertices
  used <- sort(unique(as.vector(faces)))
  remap <- integer(nrow(vertices)); remap[used] <- seq_along(used)
  faces <- matrix(remap[faces], ncol = 3)
  storage.mode(faces) <- "integer"
  list(vertices = vertices[used, , drop = FALSE], faces = faces,
       face_cell = face_cell)
}

#' Axis-aligned box mesh
#' @noRd
mesh_box <- function(x0, x1, y0, y1, z0, z1) {
  ring <- matrix(c(x0, y0, x1, y0, x1, y1, x0, y1), ncol = 2, byrow = TRUE)
  mesh_extrude_polygon(ring, z0, z1)
}

#' Blocky tree specs, ported from the arnis project
#' (<https://github.com/louis-e/arnis>, src/element_processing/tree.rs).
#'
#' Each spec describes a voxel tree in block units (1 block = 1 unit):
#' a trunk column, four inner leaf columns plus an apex cap (`fill`:
#' dx, dz, y_from, y_to), and up to three concentric canopy rings placed at
#' the listed y offsets (`rounds`, matching ROUND1/2/3 patterns in arnis).
#' @noRd
tree_spec <- function(variant) {
  round1 <- matrix(c(-2, 0, 2, 0, 0, -2, 0, 2,
                     -1, -1, 1, 1, 1, -1, -1, 1), ncol = 2, byrow = TRUE)
  round2 <- matrix(c(3, 0, 2, -1, 2, 1, 1, -2, 1, 2, -3, 0,
                     -2, -1, -2, 1, -1, 2, -1, -2, 0, -3, 0, 3),
                   ncol = 2, byrow = TRUE)
  round3 <- matrix(c(3, -1, 3, 1, 2, -2, 2, 2, 1, -3, 1, 3,
                     -3, -1, -3, 1, -2, -2, -2, 2, -1, 3, -1, -3),
                   ncol = 2, byrow = TRUE)
  specs <- list(
    oak_standard = list(
      log_height = 8,
      fill = list(c(-1, 0, 3, 9), c(1, 0, 3, 9), c(0, -1, 3, 9),
                  c(0, 1, 3, 9), c(0, 0, 9, 10)),
      rounds = list(3:8, 4:7, 5:6)
    ),
    oak_bushy = list(
      log_height = 6,
      fill = list(c(-1, 0, 3, 7), c(1, 0, 3, 7), c(0, -1, 3, 7),
                  c(0, 1, 3, 7), c(0, 0, 7, 8)),
      rounds = list(3:7, 3:6, 4:5)
    ),
    oak_compact = list(
      log_height = 5,
      fill = list(c(-1, 0, 2, 5), c(1, 0, 2, 5), c(0, -1, 2, 5),
                  c(0, 1, 2, 5), c(0, 0, 5, 6)),
      rounds = list(2:5, 3:4, integer(0))
    ),
    spruce_standard = list(
      log_height = 9,
      fill = list(c(-1, 0, 3, 10), c(0, -1, 3, 10), c(1, 0, 3, 10),
                  c(0, 1, 3, 10), c(0, 0, 11, 11)),
      rounds = list(c(9, 7, 6, 4, 3), c(6, 3), integer(0))
    ),
    spruce_towering = list(
      log_height = 12,
      fill = list(c(-1, 0, 4, 13), c(0, -1, 4, 13), c(1, 0, 4, 13),
                  c(0, 1, 4, 13), c(0, 0, 14, 14)),
      rounds = list(c(12, 10, 8, 6, 4), c(9, 6, 4), integer(0))
    )
  )
  spec <- specs[[(as.integer(variant) %% length(specs)) + 1L]]
  spec$patterns <- list(round1, round2, round3)
  spec$canopy_top <- max(vapply(spec$fill, function(f) f[4], numeric(1)))
  spec
}

#' Choose a tree variant from CHM-measured height
#'
#' Short trees get squat/compact proportions, mid-height trees the broadleaf
#' oaks, tall trees the conifer specs. A position hash adds variety between
#' the two candidates of each height class, keeping selection deterministic.
#' Spec indices: 0 oak_standard, 1 oak_bushy, 2 oak_compact,
#' 3 spruce_standard, 4 spruce_towering.
#' @noRd
tree_variant_for_height <- function(x, y, height) {
  jitter <- (abs(round(x)) + abs(round(y))) %% 2L
  if (height < 8) {
    2L                                    # oak_compact
  } else if (height < 14) {
    if (jitter == 0L) 1L else 0L          # oak_bushy / oak_standard
  } else if (height < 22) {
    if (jitter == 0L) 0L else 3L          # oak_standard / spruce_standard
  } else {
    if (jitter == 0L) 3L else 4L          # spruce_standard / spruce_towering
  }
}

#' Blocky arnis-style tree object
#'
#' Builds a voxel-look tree (trunk column, leaf columns, apex cap, canopy
#' rings) from unit blocks scaled so the total tree height matches `height`.
#' Vertical runs of blocks are merged into single boxes to keep face counts
#' low (~40 boxes per tree).
#'
#' @param x,y numeric. Tree location. @param z_base numeric ground elevation.
#' @param height numeric total tree height in metres.
#' @param variant integer or NULL. Tree variant (see `tree_spec()`); when
#' NULL, chosen deterministically from the tree position.
#' @return list(trunk = mesh, crown = mesh)
#' @noRd
mesh_tree <- function(x, y, z_base, height, variant = NULL) {
  if (!is.finite(height) || height <= 0) return(NULL)
  if (is.null(variant)) {
    variant <- tree_variant_for_height(x, y, height)
  }
  spec <- tree_spec(variant)
  s <- height / (spec$canopy_top + 1)  # block edge length in metres
  half <- s / 2

  block_col <- function(dx, dz, j1, j2) {
    mesh_box(x + dx * s - half, x + dx * s + half,
             y + dz * s - half, y + dz * s + half,
             z_base + j1 * s, z_base + (j2 + 1) * s)
  }

  trunk_h <- min(spec$log_height, spec$canopy_top - 1)
  trunk <- block_col(0, 0, 0, trunk_h)

  crown_parts <- list()
  for (f in spec$fill) {
    crown_parts[[length(crown_parts) + 1L]] <- block_col(f[1], f[2], f[3], f[4])
  }
  for (ri in seq_along(spec$rounds)) {
    ys <- sort(spec$rounds[[ri]])
    if (length(ys) == 0) next
    pat <- spec$patterns[[ri]]
    runs <- split(ys, cumsum(c(1, diff(ys) != 1)))
    for (run in runs) {
      for (ci in seq_len(nrow(pat))) {
        crown_parts[[length(crown_parts) + 1L]] <-
          block_col(pat[ci, 1], pat[ci, 2], min(run), max(run))
      }
    }
  }
  crown <- mesh_combine(crown_parts)
  list(trunk = trunk,
       crown = list(vertices = crown$vertices, faces = crown$faces))
}

#' Write a mesh (or list of layer meshes) to Wavefront OBJ (+ MTL)
#'
#' @param layers named list. Each element: list(mesh =, material =, groups =,
#' uv =) where `groups` is an optional list of list(name =, faces = integer
#' face rows, material = optional override) and `uv` an optional (n x 2)
#' matrix of texture coordinates matching the layer's vertices 1:1.
#' @param path character. Output .obj path. An .mtl with the same stem is
#' written next to it when any material is defined.
#' @param materials named list. Each element is either an RGB vector in 0-1
#' (e.g. c(0.7, 0.7, 0.7)) or a list(Kd = rgb, map_Kd = "image.jpg") for
#' textured materials.
#' @return invisible character vector of written paths.
#' @noRd
write_world_obj <- function(layers, path, materials = list()) {
  mtl_path <- sub("\\.obj$", ".mtl", path)
  has_mtl <- length(materials) > 0

  obj <- character(0)
  obj <- c(obj, "# Generated by gloBFPr::get_3d_world()")
  if (has_mtl) obj <- c(obj, paste("mtllib", basename(mtl_path)))

  v_offset <- 0L
  vt_offset <- 0L
  fmt_v <- function(m) sprintf("v %.4f %.4f %.4f", m[, 1], m[, 2], m[, 3])

  for (layer_name in names(layers)) {
    layer <- layers[[layer_name]]
    mesh <- layer$mesh
    if (is.null(mesh) || nrow(mesh$faces) == 0) next
    has_uv <- !is.null(layer$uv)
    fmt_f <- function(f) {
      if (has_uv) {
        sprintf("f %d/%d %d/%d %d/%d",
                f[, 1] + v_offset, f[, 1] + vt_offset,
                f[, 2] + v_offset, f[, 2] + vt_offset,
                f[, 3] + v_offset, f[, 3] + vt_offset)
      } else {
        sprintf("f %d %d %d",
                f[, 1] + v_offset, f[, 2] + v_offset, f[, 3] + v_offset)
      }
    }
    obj <- c(obj, paste("o", layer_name), fmt_v(mesh$vertices))
    if (has_uv) {
      obj <- c(obj, sprintf("vt %.6f %.6f", layer$uv[, 1], layer$uv[, 2]))
    }
    if (!is.null(layer$groups)) {
      for (grp in layer$groups) {
        if (length(grp$faces) == 0) next
        obj <- c(obj, paste("g", grp$name))
        mat <- if (!is.null(grp$material)) grp$material else layer$material
        if (has_mtl && !is.null(mat)) obj <- c(obj, paste("usemtl", mat))
        obj <- c(obj, fmt_f(mesh$faces[grp$faces, , drop = FALSE]))
      }
    } else {
      if (has_mtl && !is.null(layer$material)) {
        obj <- c(obj, paste("usemtl", layer$material))
      }
      obj <- c(obj, fmt_f(mesh$faces))
    }
    v_offset <- v_offset + nrow(mesh$vertices)
    if (has_uv) vt_offset <- vt_offset + nrow(mesh$vertices)
  }
  writeLines(obj, path)

  written <- path
  if (has_mtl) {
    mtl <- character(0)
    for (mat_name in names(materials)) {
      mat <- materials[[mat_name]]
      if (is.list(mat)) {
        rgb <- if (!is.null(mat$Kd)) mat$Kd else c(1, 1, 1)
        map_kd <- mat$map_Kd
      } else {
        rgb <- mat
        map_kd <- NULL
      }
      mtl <- c(mtl,
               paste("newmtl", mat_name),
               sprintf("Kd %.4f %.4f %.4f", rgb[1], rgb[2], rgb[3]),
               "Ka 0 0 0", "Ks 0 0 0", "d 1", "illum 1")
      if (!is.null(map_kd)) mtl <- c(mtl, paste("map_Kd", map_kd))
      mtl <- c(mtl, "")
    }
    writeLines(mtl, mtl_path)
    written <- c(written, mtl_path)
  }
  invisible(written)
}

#' Write a mesh as binary STL
#' @param mesh mesh list. @param path output path. @param name solid name
#' stored in the 80-byte header.
#' @return invisible path.
#' @noRd
write_stl_binary <- function(mesh, path, name = "gloBFPr") {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)

  header <- charToRaw(sprintf("%-80s", substr(paste0("gloBFPr:", name), 1, 80)))
  writeBin(header[1:80], con)
  n_faces <- nrow(mesh$faces)
  writeBin(as.integer(n_faces), con, size = 4, endian = "little")

  v <- mesh$vertices
  f <- mesh$faces
  for (i in seq_len(n_faces)) {
    p1 <- v[f[i, 1], ]; p2 <- v[f[i, 2], ]; p3 <- v[f[i, 3], ]
    nvec <- c(
      (p2[2] - p1[2]) * (p3[3] - p1[3]) - (p2[3] - p1[3]) * (p3[2] - p1[2]),
      (p2[3] - p1[3]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[3] - p1[3]),
      (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
    )
    len <- sqrt(sum(nvec^2))
    if (is.finite(len) && len > 0) nvec <- nvec / len else nvec <- c(0, 0, 0)
    writeBin(as.numeric(c(nvec, p1, p2, p3)), con, size = 4, endian = "little")
    writeBin(as.integer(0), con, size = 2, endian = "little")
  }
  invisible(path)
}

#' Orient every face of a convex solid outward
#' @noRd
mesh_orient_outward <- function(mesh) {
  v <- mesh$vertices; f <- mesh$faces
  ctr <- colMeans(v)
  for (i in seq_len(nrow(f))) {
    p1 <- v[f[i, 1], ]; p2 <- v[f[i, 2], ]; p3 <- v[f[i, 3], ]
    u <- p2 - p1; w <- p3 - p1
    n <- c(u[2] * w[3] - u[3] * w[2],
           u[3] * w[1] - u[1] * w[3],
           u[1] * w[2] - u[2] * w[1])
    if (sum(n * ((p1 + p2 + p3) / 3 - ctr)) < 0) {
      f[i, ] <- f[i, c(1, 3, 2)]
    }
  }
  storage.mode(f) <- "integer"
  list(vertices = v, faces = f)
}

#' Hexahedron from a top quad and a thickness (downward)
#' @param quad numeric 4 x 3 matrix of top-face corners, in loop order.
#' @param thickness numeric. Solid depth below the top face.
#' @noRd
mesh_hexahedron <- function(quad, thickness) {
  bottom <- quad
  bottom[, 3] <- bottom[, 3] - thickness
  v <- rbind(quad, bottom)          # 1-4 top, 5-8 bottom
  f <- rbind(
    c(1, 2, 3), c(1, 3, 4),         # top
    c(5, 7, 6), c(5, 8, 7),         # bottom
    c(1, 5, 6), c(1, 6, 2),         # sides
    c(2, 6, 7), c(2, 7, 3),
    c(3, 7, 8), c(3, 8, 4),
    c(4, 8, 5), c(4, 5, 1)
  )
  storage.mode(f) <- "integer"
  mesh_orient_outward(list(vertices = v, faces = f))
}

#' Deck elevations along a bridge centreline (arnis ramp interpolation)
#'
#' The deck is flat at `deck_z`, with linear ramps from ground level at each
#' free end, mirroring `BridgeMemberInfo::y_at()` in arnis.
#'
#' @param coords numeric n x 2 centreline.
#' @param ground numeric length-n ground elevations.
#' @param deck_z numeric. Deck elevation.
#' @param ramp numeric. Ramp length in metres (0 disables ramps).
#' @return numeric length-n deck elevations.
#' @noRd
bridge_deck_profile <- function(coords, ground, deck_z, ramp = 20) {
  n <- nrow(coords)
  if (n < 2) return(rep(deck_z, max(n, 1)))
  seg <- sqrt(diff(coords[, 1])^2 + diff(coords[, 2])^2)
  d <- c(0, cumsum(seg))
  total <- d[n]
  z <- rep(deck_z, n)
  if (!is.finite(ramp) || ramp <= 0 || total <= 0) return(z)
  rl <- min(ramp, total / 2)
  if (rl <= 0) return(z)
  g0 <- ground[1]; g1 <- ground[n]
  i <- d < rl
  if (any(i)) z[i] <- g0 + (deck_z - g0) * (d[i] / rl)
  j <- (total - d) < rl
  if (any(j)) {
    z[j] <- pmin(z[j], g1 + (deck_z - g1) * ((total - d[j]) / rl))
  }
  z
}

#' Bridge deck mesh: a ribbon swept along a centreline at given elevations
#'
#' @param coords numeric n x 2 centreline (scene coordinates).
#' @param z numeric length-n deck elevations.
#' @param half_width numeric. Half deck width in metres.
#' @param thickness numeric. Deck slab thickness.
#' @param railing numeric or 0. Parapet height; 0 omits railings.
#' @return mesh list, or NULL.
#' @noRd
mesh_bridge_deck <- function(coords, z, half_width, thickness = 0.6,
                             railing = 0) {
  n <- nrow(coords)
  if (n < 2) return(NULL)
  parts <- list()
  for (i in seq_len(n - 1)) {
    p1 <- coords[i, ]; p2 <- coords[i + 1, ]
    d <- p2 - p1
    len <- sqrt(sum(d^2))
    if (!is.finite(len) || len == 0) next
    nrm <- c(-d[2], d[1]) / len * half_width
    quad <- rbind(
      c(p1 + nrm, z[i]),
      c(p2 + nrm, z[i + 1]),
      c(p2 - nrm, z[i + 1]),
      c(p1 - nrm, z[i])
    )
    parts[[length(parts) + 1L]] <- mesh_hexahedron(quad, thickness)
    if (is.numeric(railing) && railing > 0) {
      rw <- max(half_width * 0.06, 0.15)
      for (sgn in c(1, -1)) {
        off <- nrm * sgn
        inner <- off * (1 - rw / half_width)
        rq <- rbind(
          c(p1 + off, z[i] + railing),
          c(p2 + off, z[i + 1] + railing),
          c(p2 + inner, z[i + 1] + railing),
          c(p1 + inner, z[i] + railing)
        )
        parts[[length(parts) + 1L]] <- mesh_hexahedron(rq, railing)
      }
    }
  }
  if (length(parts) == 0) return(NULL)
  m <- mesh_combine(parts)
  list(vertices = m$vertices, faces = m$faces)
}

#' Bridge support pillars dropped from the deck to the ground
#'
#' @param coords numeric n x 2 centreline. @param z deck elevations.
#' @param ground numeric length-n ground elevations.
#' @param interval numeric. Spacing between pillars in metres.
#' @param radius numeric. Half-width of the (square) pillar.
#' @param min_clear numeric. Skip pillars shorter than this.
#' @return mesh list, or NULL.
#' @noRd
mesh_bridge_supports <- function(coords, z, ground, interval = 25,
                                 radius = 0.6, min_clear = 1.5) {
  n <- nrow(coords)
  if (n < 2 || !is.finite(interval) || interval <= 0) return(NULL)
  seg <- sqrt(diff(coords[, 1])^2 + diff(coords[, 2])^2)
  d <- c(0, cumsum(seg))
  total <- d[n]
  if (total <= 0) return(NULL)
  # Pillars sit at interior positions only; a span shorter than two
  # intervals gets none (its abutments carry it).
  last <- total - interval
  if (last < interval) return(NULL)
  targets <- seq(interval, last, by = interval)
  parts <- list()
  for (t in targets) {
    k <- which.min(abs(d - t))
    top <- z[k]; base <- ground[k]
    if (!is.finite(top) || !is.finite(base) || top - base < min_clear) next
    parts[[length(parts) + 1L]] <- mesh_box(
      coords[k, 1] - radius, coords[k, 1] + radius,
      coords[k, 2] - radius, coords[k, 2] + radius,
      base - 0.5, top)
  }
  if (length(parts) == 0) return(NULL)
  m <- mesh_combine(parts)
  list(vertices = m$vertices, faces = m$faces)
}

#' Stepped voxel heightmap mesh (Minecraft-style terrain)
#'
#' Emits one top quad per cell plus vertical walls where a cell is higher
#' than its neighbor (or sits on the grid edge), producing the stepped
#' block-terrain look of arnis/Minecraft.
#'
#' @param tops numeric matrix of quantized surface elevations
#' (row 1 = north), NA cells are skipped.
#' @param xs,ys numeric. Cell-centre coordinates per column/row.
#' @param s numeric. Block edge length.
#' @param base numeric. Elevation of the boundary skirt bottom.
#' @return mesh list.
#' @noRd
mesh_voxel_heightmap <- function(tops, xs, ys, s, base = NULL) {
  nr <- nrow(tops); nc <- ncol(tops)
  if (!any(is.finite(tops))) return(NULL)
  if (is.null(base)) base <- min(tops, na.rm = TRUE) - s
  half <- s / 2
  verts <- vector("list", 0L)
  faces <- vector("list", 0L)
  cells_acc <- vector("list", 0L)
  nv <- 0L
  cur_cell <- 0L
  add_quad <- function(p1, p2, p3, p4) {
    verts[[length(verts) + 1L]] <<- rbind(p1, p2, p3, p4)
    faces[[length(faces) + 1L]] <<- rbind(c(nv + 1L, nv + 2L, nv + 3L),
                                          c(nv + 1L, nv + 3L, nv + 4L))
    cells_acc[[length(cells_acc) + 1L]] <<- c(cur_cell, cur_cell)
    nv <<- nv + 4L
  }
  neighbor <- function(i, j) {
    if (i < 1 || i > nr || j < 1 || j > nc) return(NA_real_)
    tops[i, j]
  }
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      h <- tops[i, j]
      if (!is.finite(h)) next
      cur_cell <- (i - 1L) * nc + j  # row-major linear cell index
      x0 <- xs[j] - half; x1 <- xs[j] + half
      y0 <- ys[i] - half; y1 <- ys[i] + half
      # Top face (+z)
      add_quad(c(x0, y0, h), c(x1, y0, h), c(x1, y1, h), c(x0, y1, h))
      # Walls toward lower/missing neighbors (outward normals), assigned to
      # the higher (emitting) cell's class
      hn <- neighbor(i - 1L, j)  # north (+y)
      hb <- if (is.finite(hn)) hn else base
      if (hb < h) add_quad(c(x1, y1, hb), c(x0, y1, hb),
                           c(x0, y1, h), c(x1, y1, h))
      hn <- neighbor(i + 1L, j)  # south (-y)
      hb <- if (is.finite(hn)) hn else base
      if (hb < h) add_quad(c(x0, y0, hb), c(x1, y0, hb),
                           c(x1, y0, h), c(x0, y0, h))
      hn <- neighbor(i, j + 1L)  # east (+x)
      hb <- if (is.finite(hn)) hn else base
      if (hb < h) add_quad(c(x1, y0, hb), c(x1, y1, hb),
                           c(x1, y1, h), c(x1, y0, h))
      hn <- neighbor(i, j - 1L)  # west (-x)
      hb <- if (is.finite(hn)) hn else base
      if (hb < h) add_quad(c(x0, y1, hb), c(x0, y0, hb),
                           c(x0, y0, h), c(x0, y1, h))
    }
  }
  if (nv == 0L) return(NULL)
  f <- do.call(rbind, faces)
  storage.mode(f) <- "integer"
  list(vertices = do.call(rbind, verts), faces = f,
       face_cell = unlist(cells_acc))
}

#' Voxel column mesh from grid cells, with run merging along rows
#'
#' @param cells data.frame with columns x, y (grid-aligned cell centres),
#' z0, z1 (column bottom/top).
#' @param s numeric. Block edge length.
#' @return mesh list, or NULL when no valid cells.
#' @noRd
mesh_voxel_columns <- function(cells, s) {
  cells <- cells[is.finite(cells$z1) & cells$z1 > cells$z0, , drop = FALSE]
  if (nrow(cells) == 0) return(NULL)
  ord <- order(cells$y, cells$z0, cells$z1, cells$x)
  cells <- cells[ord, , drop = FALSE]
  half <- s / 2
  boxes <- list()
  run_start <- 1L
  for (k in seq_len(nrow(cells))) {
    nxt <- k + 1L
    same_run <- nxt <= nrow(cells) &&
      cells$y[nxt] == cells$y[k] &&
      abs(cells$x[nxt] - cells$x[k] - s) < s * 1e-6 &&
      cells$z0[nxt] == cells$z0[k] &&
      cells$z1[nxt] == cells$z1[k]
    if (!same_run) {
      boxes[[length(boxes) + 1L]] <- mesh_box(
        cells$x[run_start] - half, cells$x[k] + half,
        cells$y[k] - half, cells$y[k] + half,
        cells$z0[k], cells$z1[k])
      run_start <- nxt
    }
  }
  m <- mesh_combine(boxes)
  list(vertices = m$vertices, faces = m$faces)
}

#' Arnis-style road half-width in metres by Overture/OSM class
#'
#' Ported from arnis `highway_block_range()` defaults (1 block = 1 m):
#' motorway/trunk/primary 5, secondary 4, tertiary 3, residential/service 2,
#' foot/cycle/path classes 1, unknown classes 2.
#' @noRd
road_half_width <- function(cls) {
  defaults <- c(motorway = 5, trunk = 5, primary = 5, secondary = 4,
                tertiary = 3, residential = 2, unclassified = 2,
                living_street = 2, service = 2, track = 1, path = 1,
                footway = 1, pedestrian = 1, steps = 1, cycleway = 1,
                bridleway = 1, sidewalk = 1, crosswalk = 1, alley = 2)
  out <- unname(defaults[as.character(cls)])
  out[is.na(out)] <- 2
  out
}

#' Surface group for a road segment: roads / sidewalks / paths
#' @noRd
road_group <- function(cls, subclass = NA_character_) {
  cls <- as.character(cls)
  subclass <- as.character(subclass)
  ped <- c("footway", "pedestrian", "steps", "cycleway", "crosswalk",
           "sidewalk")
  dirt <- c("path", "track", "bridleway")
  ifelse(!is.na(subclass) & subclass %in% c("sidewalk", "crosswalk"),
         "sidewalks",
         ifelse(cls %in% dirt, "paths",
                ifelse(cls %in% ped, "sidewalks", "roads")))
}

#' Detect tree tops from a canopy height model
#'
#' Local-maxima detection: one tree per `window`-metre neighbourhood. The
#' window is given in metres (not cells) so results do not change when the
#' CHM resolution does - a 5-cell window means 5 m on a 1 m metaCHM but
#' 50 m on a 10 m ethCHM, which would silently discard most trees.
#'
#' @param chm SpatRaster of canopy heights (0 = no tree), projected CRS.
#' @param window numeric. Local-maximum search window in metres.
#' @param min_height numeric. Minimum tree height.
#' @return data.frame(x, y, height)
#' @noRd
detect_tree_tops <- function(chm, window = 5, min_height = 2) {
  res_m <- terra::res(chm)[1]
  window <- max(1L, as.integer(round(window / res_m)))
  if (window %% 2L == 0L) window <- window + 1L
  if (window <= 1L) {
    # Window at or below one cell (coarse CHM): every canopy cell is a tree.
    # Forcing a 3-cell minimum here would impose e.g. 30 m spacing on a 10 m
    # raster and silently drop most of the canopy.
    tops <- chm >= min_height
  } else {
    local_max <- terra::focal(chm, w = window, fun = "max",
                              na.policy = "omit", na.rm = TRUE)
    tops <- (chm >= min_height) & (chm == local_max)
  }
  cells <- which(terra::values(tops, mat = FALSE) == 1)
  if (length(cells) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), height = numeric(0)))
  }
  xy <- terra::xyFromCell(chm, cells)
  h <- terra::values(chm, mat = FALSE)[cells]
  data.frame(x = xy[, 1], y = xy[, 2], height = h)
}
