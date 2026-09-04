# ============================================================
# ERA5 meteorological conditions for OpenFOAM boundary inputs
# ============================================================

#' Fetch ERA5 conditions for OpenFOAM boundary conditions
#'
#' @description
#' Downloads one hour of ERA5 reanalysis data (10-m wind components, 2-m
#' temperature, skin temperature) from the Copernicus Climate Data Store for a
#' single location and time step.  Returns the values pre-formatted for direct
#' use as arguments to [prepare_foam_case()] for wind simulation.
#'
#' ERA5 is a global reanalysis with 0.25 deg (~28 km) spatial resolution and
#' hourly temporal resolution from 1940 to present.  It is **not** a local
#' measurement; treat it as a representative synoptic condition, not a
#' site-specific reading.
#'
#' @section Authentication:
#' You need a free Copernicus CDS account and a personal access token.
#' 1. Register at <https://cds.climate.copernicus.eu>
#' 2. Copy your personal access token from your user profile page.
#' 3. Set it once per session: `Sys.setenv(CDS_API_KEY = "your-token-here")`
#'    or store it in `~/.Renviron` so it loads automatically.
#'
#' The legacy CDS was retired in 2024, so `ecmwfr` >= 2.0.0 is required for
#' current tokens. This function detects the installed `ecmwfr` API version and
#' adapts automatically, but on older versions the request will be rejected by
#' the server. Upgrade with `install.packages("ecmwfr")` if you hit auth errors.
#'
#' You must also accept the dataset licence once, in the browser, at
#' <https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels>
#' (the "Terms of use" tab). Requests fail with a licence error until you do.
#'
#' @param lon Numeric. Longitude of the site in decimal degrees (WGS-84).
#' @param lat Numeric. Latitude of the site in decimal degrees (WGS-84).
#' @param datetime POSIXct or character `"YYYY-MM-DD HH:MM"` in UTC.
#'   ERA5 is available at full hours; the nearest hour is used automatically.
#' @param cds_key Character. CDS personal access token.  Defaults to the
#'   `CDS_API_KEY` environment variable.
#' @param cache_dir Character. Directory for caching downloaded NetCDF files so
#'   repeated calls for the same location and time are instant.
#'   Default `tempdir()`.
#' @param quiet Logical. Suppress progress messages.  Default `FALSE`.
#'
#' @return A named list (invisibly) containing:
#' \describe{
#'   \item{`inlet_velocity`}{`c(u, v, 0)` in m/s - direct input for
#'     `prepare_foam_case(inlet_velocity = ...)`.  x = east, y = north,
#'     matching a UTM-projected domain.}
#'   \item{`z_ref`}{Reference height for the wind measurement - always 10 m.}
#'   \item{`T_ref`}{ERA5 2-m air temperature in K. In wind-only
#'     \code{prepare_foam_case()} runs, this is used as the uniform reference
#'     temperature that switches buoyancy off.}
#'   \item{`T_skin`}{ERA5 skin (land-surface) temperature in K.  Useful as a
#'     `T_ground` estimate; `NULL` if unavailable.}
#'   \item{`wind_speed_ms`}{Scalar 10-m wind speed in m/s.}
#'   \item{`wind_dir_deg`}{Meteorological wind direction in degrees (direction
#'     FROM which the wind blows; 0 = from North, 90 = from East).}
#'   \item{`u10`, `v10`}{Raw ERA5 eastward and northward wind components (m/s).}
#'   \item{`datetime`}{Rounded POSIXct of the ERA5 time step used.}
#'   \item{`lon`, `lat`}{Site coordinates as supplied.}
#' }
#'
#' @examples
#' \dontrun{
#' Sys.setenv(CDS_API_KEY = "your-token-here")
#'
#' # Summer evening in Detroit
#' met <- get_era5_met(
#'   lon               = -83.05,
#'   lat               =  42.34,
#'   datetime          = "2023-07-15 22:00"
#' )
#'
#' ## ---- Wind simulation ------------------------------------------------
#' prepare_foam_case(
#'   case_dir       = "path/to/case",
#'   stl_file       = "path/to/buildings.stl",
#'   domain         = list(xmin = 0, xmax = 500, ymin = 0, ymax = 500,
#'                         zmin = 0, zmax = 200),
#'   inlet_velocity = met$inlet_velocity,
#'   z_ref          = met$z_ref,
#'   T_ref          = met$T_ref
#' )
#' }
#'
#' @seealso [prepare_foam_case()]
#' @export
get_era5_met <- function(
    lon,
    lat,
    datetime,
    cds_key   = Sys.getenv("CDS_API_KEY"),
    cache_dir = tempdir(),
    quiet     = FALSE
) {

  # ---- dependency checks ------------------------------------------------
  if (!requireNamespace("ecmwfr", quietly = TRUE))
    stop(
      "Package 'ecmwfr' is required. Install with:\n",
      "  install.packages('ecmwfr')\n\n",
      "Then create a free account at https://cds.climate.copernicus.eu\n",
      "and set your personal access token:\n",
      "  Sys.setenv(CDS_API_KEY = 'your-token-here')",
      call. = FALSE
    )
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required.", call. = FALSE)

  if (!nzchar(trimws(cds_key)))
    stop(
      "CDS API key not found.\n",
      "  Register at https://cds.climate.copernicus.eu, then:\n",
      "  Sys.setenv(CDS_API_KEY = 'your-token-here')\n",
      "  Or pass cds_key = 'your-token-here' directly.",
      call. = FALSE
    )

  # ---- parse & round datetime to nearest full hour ---------------------
  dt <- tryCatch(
    as.POSIXct(datetime, tz = "UTC"),
    error = function(e) NA
  )
  if (is.na(dt))
    stop("`datetime` must be a valid UTC date-time, e.g. '2023-07-15 22:00'",
         call. = FALSE)

  dt <- as.POSIXct(
    round(as.numeric(dt) / 3600) * 3600,
    tz     = "UTC",
    origin = "1970-01-01"
  )

  if (!isTRUE(quiet))
    message(sprintf("ERA5: requesting data for (%.4f N, %.4f E) at %s UTC",
                    lat, lon, format(dt, "%Y-%m-%d %H:%M")))

  # ---- cache ---------------------------------------------------------------
  nc_name <- sprintf("era5_%s_lat%.3f_lon%.3f.nc",
                     format(dt, "%Y%m%d%H"), lat, lon)
  nc_path <- file.path(cache_dir, nc_name)

  if (!file.exists(nc_path)) {

    # bounding box: 1 deg buffer so the point is never on the edge
    area <- c(
      ceiling(lat  + 1),   # N
      floor(lon    - 1),   # W
      floor(lat    - 1),   # S
      ceiling(lon  + 1)    # E
    )

    request <- list(
      dataset_short_name = "reanalysis-era5-single-levels",
      product_type       = "reanalysis",
      variable           = c(
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "2m_temperature",
        "skin_temperature"
      ),
      year   = format(dt, "%Y"),
      month  = format(dt, "%m"),
      day    = format(dt, "%d"),
      time   = format(dt, "%H:00"),
      area   = area,
      target = nc_name
    )

    # ecmwfr 2.x (new CDS backend) renamed `format` -> `data_format`, and
    # dropped the `service` / `user` arguments that 1.x required.
    new_api <- !("service" %in% names(formals(ecmwfr::wf_set_key)))

    if (new_api) {
      request$data_format     <- "netcdf"
      request$download_format <- "unarchived"
    } else {
      request$format <- "netcdf"
    }

    # -- authenticate --------------------------------------------------------
    if (new_api) {
      ecmwfr::wf_set_key(key = cds_key)
    } else {
      ecmwfr::wf_set_key(key = cds_key, service = "cds")
    }

    # -- submit --------------------------------------------------------------
    req_args <- list(
      request  = request,
      transfer = TRUE,
      path     = cache_dir,
      verbose  = !quiet
    )
    if (!new_api) req_args$user <- "ecmwfr"

    tryCatch(
      do.call(ecmwfr::wf_request, req_args),
      error = function(e) {
        m <- conditionMessage(e)
        if (grepl("licence|license", m, ignore.case = TRUE))
          stop(
            "Copernicus CDS rejected the request: the ERA5 licence has not ",
            "been accepted on your account.\n\n",
            "  Fix (one time, in a browser):\n",
            "    1. Open https://cds.climate.copernicus.eu/datasets/",
            "reanalysis-era5-single-levels?tab=download#manage-licences\n",
            "    2. Scroll to 'Terms of use' and accept the licence(s)\n",
            "    3. Re-run this function\n",
            call. = FALSE
          )
        if (grepl("401|403|not authori[sz]ed|permission denied", m,
                  ignore.case = TRUE))
          stop(
            "Copernicus CDS authentication failed.\n",
            "  Check that CDS_API_KEY holds a current personal access token ",
            "from https://cds.climate.copernicus.eu (profile page), and that ",
            "ecmwfr is >= 2.0.0 (installed: ",
            as.character(utils::packageVersion("ecmwfr")), ").\n",
            call. = FALSE
          )
        stop(m, call. = FALSE)
      }
    )

    if (!file.exists(nc_path))
      stop("Download failed - NetCDF file not found at: ", nc_path, call. = FALSE)

  } else {
    if (!isTRUE(quiet))
      message("  Using cached file: ", nc_path)
  }

  # ---- extract values at the site point ------------------------------------
  r  <- terra::rast(nc_path)
  pt <- terra::vect(
    data.frame(x = lon, y = lat),
    geom = c("x", "y"),
    crs  = "EPSG:4326"
  )
  vals <- as.data.frame(terra::extract(r, pt))
  nms  <- tolower(names(vals))

  # helper: first column whose name matches any of the patterns
  pick <- function(...) {
    pats <- c(...)
    for (p in pats) {
      idx <- which(grepl(p, nms, ignore.case = TRUE))
      if (length(idx) > 0L) return(as.numeric(vals[[idx[1L]]]))
    }
    NA_real_
  }

  u10 <- pick("^u10$", "u10", "10u", "u_component_of_wind_at_10")
  v10 <- pick("^v10$", "v10", "10v", "v_component_of_wind_at_10")
  t2m <- pick("^t2m$", "t2m", "2t",  "temperature_at_2m", "2_metre_temperature")
  skt <- pick("^skt$", "skt",        "skin_temperature")

  if (any(is.na(c(u10, v10, t2m))))
    stop(
      "Could not parse u10 / v10 / t2m from the ERA5 NetCDF.\n",
      "  Available layers: ", paste(names(vals), collapse = ", "),
      call. = FALSE
    )

  # ---- derived quantities --------------------------------------------------
  wspd <- sqrt(u10^2 + v10^2)
  # Met convention: direction FROM which wind blows (0 = from N, 90 = from E)
  wdir <- (270 - atan2(v10, u10) * 180 / pi) %% 360

  # ---- build result --------------------------------------------------------
  result <- list(
    # direct args for prepare_foam_case()
    inlet_velocity = c(u10, v10, 0),
    z_ref          = 10,

    # direct arg for prepare_foam_case()
    T_ref  = t2m,
    T_skin = if (!is.na(skt)) skt else NULL,

    # diagnostics
    wind_speed_ms = wspd,
    wind_dir_deg  = wdir,
    u10           = u10,
    v10           = v10,
    datetime      = dt,
    lon           = lon,
    lat           = lat
  )

  # ---- print summary -------------------------------------------------------
  if (!isTRUE(quiet)) {
    cat(sprintf(
      "\n  ERA5 conditions - %s UTC\n",
      format(dt, "%Y-%m-%d %H:%M")
    ))
    cat(sprintf("  Location         : %.4f N, %.4f E\n", lat, lon))
    cat(sprintf("  2-m temperature  : %.1f K  (%.1f  degC)\n",
                t2m, t2m - 273.15))
    if (!is.null(result$T_skin))
      cat(sprintf("  Skin temperature : %.1f K  (%.1f  degC)\n",
                  skt, skt - 273.15))
    cat(sprintf("  10-m wind speed  : %.1f m/s\n", wspd))
    cat(sprintf("  Wind direction   : %.0f  deg  (FROM, met convention)\n", wdir))
    cat(sprintf("  u10, v10         : %+.2f, %+.2f m/s\n", u10, v10))

    cat("\n  --- Copy-paste ready ---\n")
    cat(sprintf("  inlet_velocity = c(%+.2f, %+.2f, 0)   # m/s (x=E, y=N)\n",
                u10, v10))
    cat(sprintf("  T_ref          = %.1f               # K\n", t2m))
    cat("\n")
  }

  invisible(result)
}
