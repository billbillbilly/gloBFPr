# Fetch ERA5 conditions for OpenFOAM boundary conditions

Downloads one hour of ERA5 reanalysis data (10-m wind components, 2-m
temperature, skin temperature) from the Copernicus Climate Data Store
for a single location and time step. Returns the values pre-formatted
for direct use as arguments to
[`prepare_foam_case()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_case.md)
for wind simulation.

ERA5 is a global reanalysis with 0.25 deg (~28 km) spatial resolution
and hourly temporal resolution from 1940 to present. It is **not** a
local measurement; treat it as a representative synoptic condition, not
a site-specific reading.

## Usage

``` r
get_era5_met(
  lon,
  lat,
  datetime,
  cds_key = Sys.getenv("CDS_API_KEY"),
  cache_dir = tempdir(),
  quiet = FALSE
)
```

## Arguments

- lon:

  Numeric. Longitude of the site in decimal degrees (WGS-84).

- lat:

  Numeric. Latitude of the site in decimal degrees (WGS-84).

- datetime:

  POSIXct or character `"YYYY-MM-DD HH:MM"` in UTC. ERA5 is available at
  full hours; the nearest hour is used automatically.

- cds_key:

  Character. CDS personal access token. Defaults to the `CDS_API_KEY`
  environment variable.

- cache_dir:

  Character. Directory for caching downloaded NetCDF files so repeated
  calls for the same location and time are instant. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- quiet:

  Logical. Suppress progress messages. Default `FALSE`.

## Value

A named list (invisibly) containing:

- `inlet_velocity`:

  `c(u, v, 0)` in m/s - direct input for
  `prepare_foam_case(inlet_velocity = ...)`. x = east, y = north,
  matching a UTM-projected domain.

- `z_ref`:

  Reference height for the wind measurement - always 10 m.

- `T_ref`:

  ERA5 2-m air temperature in K. In wind-only
  [`prepare_foam_case()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_case.md)
  runs, this is used as the uniform reference temperature that switches
  buoyancy off.

- `T_skin`:

  ERA5 skin (land-surface) temperature in K. Useful as a `T_ground`
  estimate; `NULL` if unavailable.

- `wind_speed_ms`:

  Scalar 10-m wind speed in m/s.

- `wind_dir_deg`:

  Meteorological wind direction in degrees (direction FROM which the
  wind blows; 0 = from North, 90 = from East).

- `u10`, `v10`:

  Raw ERA5 eastward and northward wind components (m/s).

- `datetime`:

  Rounded POSIXct of the ERA5 time step used.

- `lon`, `lat`:

  Site coordinates as supplied.

## Authentication

You need a free Copernicus CDS account and a personal access token.

1.  Register at <https://cds.climate.copernicus.eu>

2.  Copy your personal access token from your user profile page.

3.  Set it once per session:
    `Sys.setenv(CDS_API_KEY = "your-token-here")` or store it in
    `~/.Renviron` so it loads automatically.

The legacy CDS was retired in 2024, so `ecmwfr` \>= 2.0.0 is required
for current tokens. This function detects the installed `ecmwfr` API
version and adapts automatically, but on older versions the request will
be rejected by the server. Upgrade with `install.packages("ecmwfr")` if
you hit auth errors.

You must also accept the dataset licence once, in the browser, at
<https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels>
(the "Terms of use" tab). Requests fail with a licence error until you
do.

## See also

[`prepare_foam_case()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_case.md)

## Examples

``` r
if (FALSE) { # \dontrun{
Sys.setenv(CDS_API_KEY = "your-token-here")

# Summer evening in Detroit
met <- get_era5_met(
  lon               = -83.05,
  lat               =  42.34,
  datetime          = "2023-07-15 22:00"
)

## ---- Wind simulation ------------------------------------------------
prepare_foam_case(
  case_dir       = "path/to/case",
  stl_file       = "path/to/buildings.stl",
  domain         = list(xmin = 0, xmax = 500, ymin = 0, ymax = 500,
                        zmin = 0, zmax = 200),
  inlet_velocity = met$inlet_velocity,
  z_ref          = met$z_ref,
  T_ref          = met$T_ref
)
} # }
```
