# Sample OpenFOAM results at pedestrian level and return a raster

Uses OpenFOAM's `postProcess -func surfaces` utility to cut a horizontal
plane at z = 1.5 m (or any height) through the latest time-step, then
reads the resulting `.raw` files into R and returns a multi-layer
`SpatRaster`.

For the current wind workflow,
[`prepare_foam_case()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_case.md)
already writes a pedestrian slice during the solver run via the
`pedestrianSlice` function object; use
[`read_foam_pedestrian_slice`](https://billbillbilly.github.io/gloBFPr/reference/read_foam_pedestrian_slice.md)
for that default 1.5 m output. Use this helper to sample other heights
or fields after a run.

## Usage

``` r
sample_foam_slice(
  case_dir,
  fields = c("U", "p"),
  z = 1.5,
  image = "opencfd/openfoam-run:2506",
  resolution = 5,
  time_step = "latestTime",
  interpolate = c("smooth", "focal", "none"),
  idw_power = 2,
  idw_maxdist = NULL,
  buildings = NULL,
  quiet = FALSE
)
```

## Arguments

- case_dir:

  Character. OpenFOAM case directory.

- fields:

  Character vector of field names to sample. Default `c("U", "p")`.

- z:

  Numeric. Height above ground in metres. Default 1.5.

- image:

  Character. Docker image tag. Default `"opencfd/openfoam-run:2506"`.

- resolution:

  Numeric. Output raster cell size in metres. Default 5.

- time_step:

  Character or numeric. `"latestTime"` (default) or a specific time-step
  number.

- interpolate:

  Character. Post-rasterization smoothing method. `"smooth"` (default)
  fills NA gaps then applies a Gaussian focal filter to all cells,
  producing continuous gradients similar to ParaView; `"focal"` fills NA
  gaps only (faster, but blocky in open areas); `"none"` returns the raw
  rasterization without any smoothing.

- idw_power:

  Numeric. Controls the Gaussian sigma when `interpolate = "smooth"`:
  sigma = `resolution * idw_power` metres. Default 2 (10 m sigma at 5 m
  resolution). Increase for a wider blur; decrease for sharper
  transitions near buildings.

- idw_maxdist:

  Ignored (reserved for future use).

- buildings:

  Optional `sf` polygon layer of building footprints in the domain's
  local coordinate system. When `NULL` (default) the footprints saved by
  [`prepare_openfoam_inputs`](https://billbillbilly.github.io/gloBFPr/reference/prepare_openfoam_inputs.md)
  are loaded automatically from
  `<case_dir>/constant/gloBFPr/metadata/buildings_openfoam.rds`.
  Building interiors are masked to `NA`, and the layer is attached to
  the returned raster so
  [`plot_foam_map`](https://billbillbilly.github.io/gloBFPr/reference/plot_foam_map.md)
  overlays it without any extra argument. Pass
  `buildings = sf::st_sf(...)` to override, or note that masking is
  skipped when no footprints can be found.

- quiet:

  Logical. Default FALSE.

## Value

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with layers named after the sampled fields. Velocity `U` produces three
extra layers: `U_mag` (speed), `Ux`, `Uy`. Building footprints, when
available, are attached as `attr(x, "buildings")`.

## See also

[`plot_foam_map`](https://billbillbilly.github.io/gloBFPr/reference/plot_foam_map.md),
[`read_foam_pedestrian_slice`](https://billbillbilly.github.io/gloBFPr/reference/read_foam_pedestrian_slice.md)
