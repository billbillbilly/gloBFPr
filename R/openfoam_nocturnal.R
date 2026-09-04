# ===========================================================================
# Pedestrian-level OpenFOAM post-processing
# ===========================================================================
# The current development branch is wind-only. What remains here is the
# pedestrian-slice reader used by the wind workflow.
# ===========================================================================


#' @noRd
prepare_nocturnal_case <- function() {
  stop("`prepare_nocturnal_case()` is not part of the current wind-only ",
       "OpenFOAM workflow. Use `prepare_foam_case()` for wind simulation.",
       call. = FALSE)
}


# ===========================================================================
# Post-processing
# ===========================================================================

#' Read the pedestrian-level slice and compute wind maps
#'
#' @description
#' Reads the pedestrian-level surface sample written by the
#' \code{pedestrianSlice} function object and returns a multi-layer
#' \code{SpatRaster}. The current OpenFOAM workflow is wind-only, so the
#' velocity layers are the primary result; temperature-derived layers should be
#' approximately zero for neutral wind cases:
#'
#' \describe{
#'   \item{T_air}{Air temperature (K)}
#'   \item{U_mag}{Wind speed (m/s)}
#'   \item{T_cool}{Cooling relative to ambient, \eqn{T_{ref} - T_{air}} (K)}
#'   \item{T_cool_flux}{\eqn{\max(T_{cool},0) \times |U|} - cool-air transport}
#'   \item{Ux, Uy}{Horizontal velocity components (m/s)}
#' }
#'
#' The case is transient, so \code{time_step} selects a physical time in
#' seconds since sunset; \code{"latest"} gives the end of the run.
#'
#' @param case_dir Character. OpenFOAM case directory.
#' @param T_ref Numeric. Reference temperature (K) for \code{T_cool}, which is
#'   \code{T_ref - T}.  It has to match the temperature the case was
#'   initialised at, or every cooling number is offset by the difference and
#'   \code{T_cool} can come out negative everywhere.  The default,
#'   \code{NULL}, reads \code{TRef} from the case's
#'   \code{constant/transportProperties} and only falls back to 295 K when
#'   that file cannot be read.
#' @param time_step \code{"latest"} (default) or a number of seconds.
#' @param resolution Numeric. Output cell size (m). Default 2.
#' @param base_cell_size Numeric. Background mesh cell size used when the case
#'   was generated (\code{foam$params$base_cell_size}); sets the gap-fill
#'   window. Default 10.
#' @param agl_tol Numeric. Samples more than this many metres above the local
#'   ground surface are discarded.  Guards against the stray near-vertical
#'   sheet a distanceSurface can generate around a closed terrain STL, which
#'   otherwise folds upper-level wind into the pedestrian map. Default 5.
#' @param trim Crop the flow-adjustment zone off the returned raster.
#'   \code{"auto"} (default) removes, per side, the building-free apron plus an
#'   adjustment fetch; \code{"buildings"} removes only the apron; a number
#'   removes that many metres from every side; \code{0} returns the full domain.
#'
#'   Trimming is on by default because the outer band is not a result. The
#'   domain is larger than the built area, and even where buildings reach the
#'   boundary the ABL inlet delivers an undisturbed profile that only slows as
#'   it works into the roughness. Measured on a real case, the outer 100 m
#'   averaged 1.86 m/s against 0.74 m/s in the interior - a bright rim that is
#'   an artifact of the domain, not a feature of the city. Excluding it is
#'   standard practice in urban CFD.
#'
#'   Each side is trimmed by \code{max(apron, fetch)} - not their sum. The
#'   fetch is measured from the domain boundary and the apron is its first
#'   part, so they overlap rather than stack.
#'
#'   The fetch is ~5x the median building height (100 m when the layer has no
#'   height column), clamped to 50-150 m and to 10% of the shorter domain
#'   span. That is deliberately less than the ~15x an internal boundary layer
#'   needs to fully equilibrate: on a measured case the edge excess was
#'   concentrated in the outer 100 m, and beyond that the band-to-band
#'   variation was genuine urban structure rather than an edge effect. The
#'   amount removed on each side is reported.
#' @param canopy Canopy height raster (path or SpatRaster) used to build a
#'   canopy overlay for \code{\link{plot_foam_map}}.  Auto-detected from
#'   \code{constant/gloBFPr/rasters/canopy_height.tif} when \code{NULL}; pass
#'   \code{FALSE} to skip it.  The result is attached as the \code{"canopy"}
#'   attribute, mirroring \code{"buildings"}.
#' @param min_canopy_height Numeric. Canopy cells below this height are ignored.
#'   Default 2 m, matching \code{min_tree_height} in
#'   \code{\link{prepare_openfoam_inputs}}.
#' @param crs CRS to assign (e.g. \code{32617}). Default \code{NA}.
#' @param buildings Optional \code{sf} footprints in local coordinates;
#'   auto-detected from the case directory when \code{NULL}.
#' @param quiet Logical. Default \code{FALSE}.
#'
#' @return A \code{terra::SpatRaster} with six layers.
#' @seealso \code{\link{prepare_foam_case}}, \code{\link{run_openfoam_docker}}
#' @export
read_foam_pedestrian_slice <- function(case_dir,
                                       T_ref          = NULL,
                                       time_step      = "latest",
                                       resolution     = 2,
                                       base_cell_size = 10,
                                       agl_tol        = 5,
                                       trim           = "auto",
                                       canopy         = NULL,
                                       min_canopy_height = 2,
                                       crs            = NA,
                                       buildings      = NULL,
                                       quiet          = FALSE) {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required. install.packages('terra')", call. = FALSE)

  case_dir <- normalizePath(case_dir, mustWork = TRUE)

  if (is.null(buildings)) {
    bpath <- file.path(case_dir, "constant", "gloBFPr", "metadata",
                       "buildings_openfoam.rds")
    if (file.exists(bpath)) {
      buildings <- tryCatch(readRDS(bpath), error = function(e) NULL)
      if (!isTRUE(quiet) && !is.null(buildings))
        message("Auto-loaded ", nrow(buildings), " buildings from case directory.")
    }
  }

  pp_dir <- file.path(case_dir, "postProcessing", "pedestrianSlice")
  if (!dir.exists(pp_dir))
    stop("No postProcessing/pedestrianSlice directory in:\n  ", case_dir,
         "\nRun the simulation first with run_openfoam_docker().", call. = FALSE)

  time_dirs <- list.dirs(pp_dir, full.names = FALSE, recursive = FALSE)
  num_times <- suppressWarnings(as.numeric(time_dirs))
  num_times <- sort(num_times[!is.na(num_times)])
  if (length(num_times) == 0)
    stop("No time directories found under ", pp_dir, call. = FALSE)

  t_chosen <- if (identical(time_step, "latest")) max(num_times) else {
    ts <- suppressWarnings(as.numeric(time_step))
    if (is.na(ts) || !(ts %in% num_times))
      stop(sprintf("time_step %s not found. Available: %s",
                   time_step, paste(num_times, collapse = ", ")), call. = FALSE)
    ts
  }
  t_name <- time_dirs[which(suppressWarnings(as.numeric(time_dirs)) == t_chosen)][1L]
  t_dir  <- file.path(pp_dir, t_name)

  # Resolve the reference temperature from the case rather than a default.
  # buoyantBoussinesqPimpleFoam's TRef is the temperature the Boussinesq
  # buoyancy is measured against and the temperature prepare_foam_case()
  # initialises the air at, so it is also the only defensible zero point for
  # "cooling below ambient".  Getting it wrong is silent: on a case with
  # TRef = 295.783 the old default of 295 made T_cool negative across the
  # whole map, and T_cool_flux - which clips at zero - collapsed to nothing.
  if (is.null(T_ref)) {
    tp <- file.path(case_dir, "constant", "transportProperties")
    T_ref <- if (file.exists(tp)) {
      ln <- grep("^\\s*TRef\\s", readLines(tp, warn = FALSE), value = TRUE)
      v  <- if (length(ln))
        suppressWarnings(as.numeric(sub(".*\\]\\s*([0-9.eE+-]+)\\s*;.*", "\\1", ln[1L])))
      else NA_real_
      if (is.finite(v)) v else NA_real_
    } else NA_real_
    if (!is.finite(T_ref)) {
      T_ref <- 295
      warning("Could not read TRef from constant/transportProperties; ",
              "using T_ref = 295 K. Pass `T_ref` explicitly if that is wrong.",
              call. = FALSE)
    } else if (!isTRUE(quiet)) {
      message(sprintf("Reference temperature: %.3f K (from constant/transportProperties)",
                      T_ref))
    }
  }

  # v2506 reversed the .raw naming convention:
  #   pre-v2506: <surface>_<field>.raw   -> z1p5m_T.raw
  #   v2506+   : <field>_<surface>.raw   -> T_z1p5m.raw
  # and v2206+ may nest under surfaces/<surface>/<field>.raw.
  find_noc_raw <- function(dir, field, surf = "z1p5m") {
    cand <- c(file.path(dir, paste0(surf, "_", field, ".raw")),
              file.path(dir, paste0(field, "_", surf, ".raw")),
              file.path(dir, "surfaces", surf, paste0(field, ".raw")))
    found <- cand[file.exists(cand)]
    if (!length(found))
      stop(sprintf("Raw file for '%s' not found in:\n  %s\nTried:\n  %s",
                   field, dir, paste(cand, collapse = "\n  ")), call. = FALSE)
    found[[1L]]
  }

  # postProcessing/ survives a re-run unless the Allrun script clears it, so a
  # time directory can belong to an earlier run in a different mode, on a
  # different mesh.  Both look like results.  The mesh is rebuilt (or at least
  # touched) by the run that produced the current slices, so a slice older than
  # the mesh cannot have come from this configuration.
  mesh_f <- file.path(case_dir, "constant", "polyMesh", "owner")
  if (file.exists(mesh_f)) {
    raw_any <- list.files(t_dir, pattern = "\\.raw$", full.names = TRUE)
    if (length(raw_any) &&
        max(file.mtime(raw_any)) < file.mtime(mesh_f)) {
      warning(sprintf(paste0("The slice at t = %s was written before the ",
                             "current mesh, so it is left over from an ",
                             "earlier run - possibly in a different mode. ",
                             "Re-run the case, or delete postProcessing/ and ",
                             "re-run, before trusting this map."), t_name),
              call. = FALSE)
    }
  }

  if (!isTRUE(quiet))
    message(sprintf("Reading slice at t = %g s (%.2f h after sunset)",
                    t_chosen, t_chosen / 3600))

  read_raw <- function(path, col_names) {
    txt <- readLines(path)
    dl  <- txt[!grepl("^\\s*#", txt) & nzchar(trimws(txt))]
    if (!length(dl)) stop("No data in ", path, call. = FALSE)
    utils::read.table(text = paste(dl, collapse = "\n"), header = FALSE,
                      col.names = col_names)
  }

  T_dat <- read_raw(find_noc_raw(t_dir, "T"), c("x", "y", "z", "T"))
  U_dat <- read_raw(find_noc_raw(t_dir, "U"), c("x", "y", "z", "Ux", "Uy", "Uz"))

  # The sampled surface is MULTI-VALUED in (x, y).  The slice is a
  # distanceSurface offset from the terrain STL, and that STL is a closed
  # solid, so its iso-distance surface does not only cover the ground: it also
  # wraps the vertical skirt at the domain edge and every near-vertical face
  # inside the domain, producing sheets that span most of the domain height.
  # Measured on a real case: z ran 157-441 m where the terrain tops out at 295,
  # and those stray samples carried upper-level wind speeds down into the
  # pedestrian map as a bright false rim.
  #
  # De-duplicating on exact (x, y) does NOT work - the stray samples sit at
  # nearby, not identical, coordinates (it removed 0.2% and changed nothing).
  # What separates a real sample from a stray one is height above the terrain:
  # by construction every genuine sample sits at ONE fixed offset above it.
  #
  # So the filter needs a ground reference, and getting that reference wrong
  # does visible damage in both directions - too low and the sheets leak in,
  # too high and legitimate ground samples are thrown away, leaving voids in
  # the map.  Prefer the DEM the case was actually meshed on, which
  # prepare_foam_geometry() leaves in the case; estimate it from the sample
  # cloud only when that file is gone.
  ground_from_dem <- function() {
    dem_file <- file.path(case_dir, "constant", "gloBFPr", "rasters",
                          "terrain_bare.tif")
    if (!file.exists(dem_file)) return(NULL)
    dem <- tryCatch(terra::rast(dem_file), error = function(e) NULL)
    if (is.null(dem)) return(NULL)
    if (terra::nlyr(dem) > 1L) dem <- dem[[1L]]
    # A hole in the DEM would reject every sample above it and punch a void
    # straight through the map, so close small ones before using it.
    if (isTRUE(terra::global(is.na(dem), "sum")[[1L]] > 0))
      for (i in 1:4)
        dem <- terra::focal(dem, w = 5, fun = mean, na.policy = "only",
                            na.rm = TRUE)
    function(x, y) terra::extract(dem, cbind(x, y), method = "bilinear")[, 1L]
  }

  # Fallback.  Takes the MINIMUM z per bin, not the minimum over the bin's 3x3
  # neighbourhood, which is what an earlier version did: over a 75 m
  # neighbourhood a 10 % slope drops 7.5 m, more than `agl_tol`, so the
  # neighbourhood minimum rejected legitimate ground samples wherever the
  # terrain rose - punching voids into the map on exactly the slopes that drive
  # cold-air drainage.  A bin whose lowest sample still stands well above the
  # surrounding ground has no ground sample in it at all (it is roofed, or it
  # holds only a vertical sheet), so it is discarded and re-interpolated.  The
  # comparison is against the neighbourhood MEDIAN, which a uniform slope
  # leaves unbiased, rather than the minimum, which it does not.
  ground_from_cloud <- function(d, bin = 15, roof_tol = 8, med_w = 11) {
    tmpl <- terra::rast(xmin = min(d$x) - bin, xmax = max(d$x) + bin,
                        ymin = min(d$y) - bin, ymax = max(d$y) + bin,
                        resolution = bin, crs = "")
    g <- terra::rasterize(terra::vect(d[, c("x", "y", "z")],
                                      geom = c("x", "y"), crs = ""),
                          tmpl, field = "z", fun = "min")
    for (i in 1:2) {
      med <- terra::focal(g, w = med_w, fun = "median", na.rm = TRUE)
      g[g > med + roof_tol] <- NA
    }
    for (i in 1:8)
      g <- terra::focal(g, w = 5, fun = mean, na.policy = "only", na.rm = TRUE)
    g <- terra::focal(g, w = 3, fun = mean, na.rm = TRUE)
    function(x, y) terra::extract(g, cbind(x, y), method = "bilinear")[, 1L]
  }

  # The sampler writes its height above ground into the file NAME, not into the
  # data, so read the offset back from the samples themselves - restricted to a
  # plausible band first, so the sheets cannot drag it.
  keep_surface <- function(d, ground_at) {
    dz   <- d$z - ground_at(d$x, d$y)
    near <- is.finite(dz) & abs(dz) <= 20
    off  <- if (any(near)) stats::median(dz[near]) else 0
    d[is.finite(dz) & abs(dz - off) <= agl_tol, , drop = FALSE]
  }

  n_before  <- nrow(U_dat)
  ground_at <- ground_from_dem()
  src       <- "constant/gloBFPr/rasters/terrain_bare.tif"
  if (!is.null(ground_at)) {
    # A DEM left over from another domain, or one still in map coordinates
    # while the case is in local metres, would silently reject nearly
    # everything.  Check before trusting it.
    keep <- nrow(keep_surface(U_dat, ground_at)) / max(1L, n_before)
    if (keep < 0.2) {
      warning(sprintf(paste0("%s keeps only %.0f%% of the slice samples and ",
                             "looks unrelated to this case; estimating the ",
                             "ground from the sample cloud instead."),
                      src, 100 * keep), call. = FALSE)
      ground_at <- NULL
    }
  }
  if (is.null(ground_at)) {
    ground_at <- ground_from_cloud(U_dat)
    src       <- "estimated from the sample cloud"
  }
  if (!isTRUE(quiet)) message("Ground reference: ", src)

  T_dat <- keep_surface(T_dat, ground_at)
  U_dat <- keep_surface(U_dat, ground_at)
  if (!isTRUE(quiet) && nrow(U_dat) < n_before)
    message(sprintf("Kept %d of %d samples within %g m of the ground surface (%.0f%% dropped)",
                    nrow(U_dat), n_before, agl_tol,
                    100 * (n_before - nrow(U_dat)) / n_before))

  dat <- merge(T_dat[, c("x", "y", "T")],
               U_dat[, c("x", "y", "Ux", "Uy", "Uz")],
               by = c("x", "y"), all = TRUE)
  dat$U_mag       <- sqrt(dat$Ux^2 + dat$Uy^2 + dat$Uz^2)
  dat$T_cool      <- T_ref - dat$T
  dat$T_cool_flux <- pmax(dat$T_cool, 0) * dat$U_mag

  r_template <- terra::rast(
    xmin = min(dat$x, na.rm = TRUE), xmax = max(dat$x, na.rm = TRUE),
    ymin = min(dat$y, na.rm = TRUE), ymax = max(dat$y, na.rm = TRUE),
    resolution = resolution,
    crs = if (is.na(crs)) "" else as.character(crs))
  pts <- terra::vect(dat, geom = c("x", "y"),
                     crs = if (is.na(crs)) "" else as.character(crs))
  rasterize_mean <- function(field)
    terra::rasterize(pts, r_template, field = field, fun = "mean")

  # Gap-fill sparse CFD points, then Gaussian-smooth, so the output is a
  # continuous field rather than isolated dots.  The window must span the
  # coarsest mesh cell, hence base_cell_size / resolution.
  smooth_noc_layer <- function(r, gap_w, sigma = 2) {
    nm <- names(r)
    r <- terra::focal(r, w = gap_w, fun = mean, na.policy = "only", na.rm = TRUE)
    r <- terra::focal(r, w = gap_w, fun = mean, na.policy = "only", na.rm = TRUE)
    ksize  <- max(3L, 2L * ceiling(3 * sigma) + 1L)
    ax     <- seq(-(ksize %/% 2L), ksize %/% 2L)
    kernel <- outer(ax, ax, function(x, y) exp(-(x^2 + y^2) / (2 * sigma^2)))
    kernel <- kernel / sum(kernel)
    # Normalise by the weight that actually landed on data.
    #
    # terra::focal(w = kernel, na.rm = TRUE) sums w*x over non-NA cells but
    # keeps the full denominator, because the kernel already sums to 1. Next
    # to an NA region - i.e. next to every building - only a fraction of the
    # weight is valid, so the result is scaled toward ZERO rather than toward
    # the local mean. Feeding it a uniform 295 K field with building-shaped
    # holes returned 0.65-295 K, with 316 cells below 100 K: the "cold spots"
    # were courtyards and alleys where most of the kernel sat on buildings.
    #
    # This silently biased every layer low near buildings - including wind
    # speed, exactly where pedestrian comfort is read.
    num <- terra::focal(r, w = kernel, fun = "sum", na.rm = TRUE)
    den <- terra::focal(!is.na(r), w = kernel, fun = "sum", na.rm = TRUE)
    r   <- num / den
    names(r) <- nm
    r
  }
  gap_w <- max(as.integer(2 * ceiling(1.5 * base_cell_size / resolution) + 1L), 7L)

  raw_layers <- list(
    T_air = rasterize_mean("T"),           U_mag = rasterize_mean("U_mag"),
    T_cool = rasterize_mean("T_cool"),     T_cool_flux = rasterize_mean("T_cool_flux"),
    Ux = rasterize_mean("Ux"),             Uy = rasterize_mean("Uy"))

  result <- terra::rast(lapply(names(raw_layers), function(nm) {
    r <- raw_layers[[nm]]; names(r) <- nm
    smooth_noc_layer(r, gap_w = gap_w, sigma = 2)
  }))
  names(result) <- c("T_air", "U_mag", "T_cool", "T_cool_flux", "Ux", "Uy")

  # Re-apply the building mask AFTER gap-filling: the fill floods over small
  # buildings, so without this only the largest ones show as voids.
  if (!is.null(buildings)) {
    if (!requireNamespace("sf", quietly = TRUE))
      warning("Package 'sf' required to mask buildings; skipping.", call. = FALSE)
    else {
      bvect      <- terra::vect(sf::st_geometry(buildings))
      building_r <- terra::rasterize(bvect, result[[1L]], field = 1L,
                                     background = NA)
      result     <- terra::mask(result, building_r, maskvalues = 1L,
                                updatevalue = NA)
    }
  }
  # ---- Crop the flow-adjustment zone -------------------------------------
  # Two distinct things make the outer band unusable, and both are per-side:
  #
  #   1. The building-free apron. The domain is larger than the built area, so
  #      part of the map is simply not city. Definitional, measurable exactly.
  #   2. Adjustment fetch. Even where buildings start at the boundary, the ABL
  #      inlet delivers an undisturbed profile that only slows as it works into
  #      the roughness. Measured on a real case: the outer 100 m averaged
  #      1.86 m/s against 0.74 m/s in the interior.
  #
  # Neither is a result, so "auto" removes both. The fetch follows the usual
  # internal-boundary-layer rule of thumb, ~15x the roughness element height,
  # taken from the building heights when the layer carries them. It is clamped
  # so it can never eat the map, and the crop is reported rather than silent.
  #
  # Done AFTER gap-filling, smoothing and the building mask, so surviving edge
  # cells were computed with a full neighbourhood - cropping earlier would just
  # move the artifact inward.
  auto_trim_sides <- function(e, buildings, fetch_cap_frac = 0.10) {
    span_x <- as.numeric(e[2] - e[1]); span_y <- as.numeric(e[4] - e[3])
    margin <- c(0, 0, 0, 0)
    if (!is.null(buildings) && requireNamespace("sf", quietly = TRUE)) {
      bb <- sf::st_bbox(buildings)
      margin <- c(max(0, bb[["xmin"]] - as.numeric(e[1])),
                  max(0, as.numeric(e[2]) - bb[["xmax"]]),
                  max(0, bb[["ymin"]] - as.numeric(e[3])),
                  max(0, as.numeric(e[4]) - bb[["ymax"]]))
    }
    h <- NULL
    if (!is.null(buildings)) {
      hc <- intersect(c("Height", "height", "height_m", "mean_height"),
                      names(buildings))
      if (length(hc)) {
        hv <- suppressWarnings(as.numeric(buildings[[hc[1]]]))
        hv <- hv[is.finite(hv) & hv > 0]
        if (length(hv)) h <- stats::median(hv)
      }
    }
    fetch <- if (is.null(h)) 100 else min(max(5 * h, 50), 150)
    fetch <- min(fetch, fetch_cap_frac * min(span_x, span_y))
    # max(), NOT margin + fetch.  The adjustment fetch is measured from the
    # DOMAIN boundary, and the building-free apron is the first part of it -
    # they overlap, they do not stack.  Adding them double-counted the apron
    # and removed 42% of a real map where 19% was warranted.
    list(trim = pmax(margin, fetch), fetch = fetch, margin = margin, h = h)
  }

  if (!identical(trim, 0) && !is.null(trim)) {
    e <- terra::ext(result)
    sides <- NULL   # c(west, east, south, north)
    if (identical(trim, "auto")) {
      at <- auto_trim_sides(e, buildings)
      sides <- at$trim
      if (!isTRUE(quiet))
        message(sprintf(
          "Auto-trim: per side max(apron, fetch); apron W/E/S/N = %.0f/%.0f/%.0f/%.0f m, fetch %.0f m%s",
          at$margin[1], at$margin[2], at$margin[3], at$margin[4], at$fetch,
          if (is.null(at$h)) " (default; no height column found)"
          else sprintf(" (5 x %.1f m median building height)", at$h)))
    } else if (identical(trim, "buildings")) {
      if (is.null(buildings)) {
        warning("`trim = \"buildings\"` needs the buildings layer; not cropping.",
                call. = FALSE)
      } else {
        bb <- sf::st_bbox(buildings)
        sides <- c(max(0, bb[["xmin"]] - as.numeric(e[1])),
                   max(0, as.numeric(e[2]) - bb[["xmax"]]),
                   max(0, bb[["ymin"]] - as.numeric(e[3])),
                   max(0, as.numeric(e[4]) - bb[["ymax"]]))
      }
    } else {
      tr <- suppressWarnings(as.numeric(trim))
      if (is.na(tr) || tr < 0)
        warning("`trim` must be a non-negative number, \"auto\" or \"buildings\"; ",
                "not cropping.", call. = FALSE)
      else sides <- rep(tr, 4)
    }

    if (!is.null(sides)) {
      if (as.numeric(e[2] - e[1]) <= sides[1] + sides[2] ||
          as.numeric(e[4] - e[3]) <= sides[3] + sides[4]) {
        warning("`trim` would remove the whole domain; not cropping.",
                call. = FALSE)
      } else {
        new_e <- terra::ext(as.numeric(e[1]) + sides[1],
                            as.numeric(e[2]) - sides[2],
                            as.numeric(e[3]) + sides[3],
                            as.numeric(e[4]) - sides[4])
        result <- terra::crop(result, new_e)
        if (!is.null(buildings) && requireNamespace("sf", quietly = TRUE)) {
          # as.numeric() is required: terra::ext()[i] carries a name ("xmin"),
          # so c(xmin = new_e[1]) becomes "xmin.xmin" and st_bbox() sees NA.
          bx <- sf::st_bbox(c(xmin = as.numeric(new_e[1]),
                              xmax = as.numeric(new_e[2]),
                              ymin = as.numeric(new_e[3]),
                              ymax = as.numeric(new_e[4])),
                            crs = sf::st_crs(buildings))
          # CROP, not select-by-intersection.  Keeping whole buildings that
          # merely touch the extent leaves them hanging outside the raster,
          # which reads as a broken plot.
          buildings <- suppressWarnings(sf::st_crop(buildings, bx))
        }
        if (!isTRUE(quiet))
          message(sprintf("Trimmed to %.0f x %.0f m (removed W/E/S/N %.0f/%.0f/%.0f/%.0f m)",
                          as.numeric(new_e[2] - new_e[1]),
                          as.numeric(new_e[4] - new_e[3]),
                          sides[1], sides[2], sides[3], sides[4]))
      }
    }
  }

  # ---- Report voids --------------------------------------------------------
  # A cell with no result inside the mapped area is not a plotting artifact -
  # it is somewhere the flow never reached at this height.  Buildings are
  # expected; anything else is a structure that is in the DSM but missing from
  # the footprint layer, so it survived into the terrain and the pedestrian
  # surface climbed over it.  Say so rather than leaving an unexplained hole,
  # because it points at the input data, not at this function.
  if (!isTRUE(quiet)) {
    na_cells <- is.na(result[[1L]])
    n_void   <- as.numeric(terra::global(na_cells, "sum", na.rm = TRUE)[[1L]])
    if (!is.null(buildings) && requireNamespace("sf", quietly = TRUE)) {
      br <- terra::rasterize(terra::vect(sf::st_geometry(buildings)),
                             result[[1L]], field = 1L, background = NA)
      n_void <- as.numeric(terra::global(na_cells & is.na(br), "sum",
                                         na.rm = TRUE)[[1L]])
    }
    if (isTRUE(n_void > 0))
      message(sprintf(paste0("%d cells (%.2f%% of the map, %.0f m2) have no ",
                             "result and no building footprint - most likely ",
                             "structures present in the fused DSM but absent ",
                             "from the building layer. They plot solid, like ",
                             "buildings."),
                      n_void, 100 * n_void / terra::ncell(result),
                      n_void * prod(terra::res(result))))
  }

  # ---- Canopy overlay ------------------------------------------------------
  # Kept as cell centres rather than polygons: the CHM is a fine raster and
  # dissolving it to polygons is slow for no visual gain at plot scale.
  canopy_df <- NULL
  if (!identical(canopy, FALSE)) {
    cr <- NULL
    if (inherits(canopy, "SpatRaster")) cr <- canopy
    else if (is.character(canopy) && file.exists(canopy)) cr <- terra::rast(canopy)
    else if (is.null(canopy)) {
      cpath <- file.path(case_dir, "constant", "gloBFPr", "rasters",
                         "canopy_height.tif")
      if (file.exists(cpath)) cr <- terra::rast(cpath)
    }
    if (!is.null(cr)) {
      cr <- tryCatch({
        if (terra::nlyr(cr) > 1L) cr <- cr[[1L]]
        cr <- terra::crop(cr, terra::ext(result))
        # Coarsen to roughly the output resolution so the overlay is a
        # manageable number of tiles; "max" keeps thin tree lines visible.
        target <- max(resolution, 4)
        if (terra::res(cr)[1] < target) {
          f <- max(1L, round(target / terra::res(cr)[1]))
          if (f > 1L) cr <- terra::aggregate(cr, fact = f, fun = "max",
                                             na.rm = TRUE)
        }
        cr
      }, error = function(e) NULL)
    }
    if (!is.null(cr)) {
      v  <- terra::values(cr, mat = FALSE)
      ok <- which(is.finite(v) & v >= min_canopy_height)
      if (length(ok)) {
        xy <- terra::xyFromCell(cr, ok)
        canopy_df <- data.frame(x = xy[, 1], y = xy[, 2])
        attr(canopy_df, "res") <- terra::res(cr)[1]
        if (!isTRUE(quiet))
          message(sprintf("Canopy overlay: %d cells >= %g m at %.1f m resolution",
                          nrow(canopy_df), min_canopy_height, terra::res(cr)[1]))
      }
    }
  }

  # re-attach: terra::crop() drops user attributes
  if (!is.null(canopy_df)) attr(result, "canopy") <- canopy_df
  if (!is.null(buildings)) attr(result, "buildings") <- buildings
  attr(result, "time_s") <- t_chosen

  if (!isTRUE(quiet))
    message(sprintf("Raster: %d x %d cells at %g m | %d sample points",
                    terra::nrow(result), terra::ncol(result), resolution,
                    nrow(dat)))
  result
}
