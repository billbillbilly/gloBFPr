# Read the pedestrian-level slice and compute wind maps

Reads the pedestrian-level surface sample written by the
`pedestrianSlice` function object and returns a multi-layer
`SpatRaster`. The current OpenFOAM workflow is wind-only, so the
velocity layers are the primary result; temperature-derived layers
should be approximately zero for neutral wind cases:

- T_air:

  Air temperature (K)

- U_mag:

  Wind speed (m/s)

- T_cool:

  Cooling relative to ambient, \\T\_{ref} - T\_{air}\\ (K)

- T_cool_flux:

  \\\max(T\_{cool},0) \times \|U\|\\ - cool-air transport

- Ux, Uy:

  Horizontal velocity components (m/s)

The case is transient, so `time_step` selects a physical time in seconds
since sunset; `"latest"` gives the end of the run.

## Usage

``` r
read_foam_pedestrian_slice(
  case_dir,
  T_ref = NULL,
  time_step = "latest",
  resolution = 2,
  base_cell_size = 10,
  agl_tol = 5,
  trim = "auto",
  canopy = NULL,
  min_canopy_height = 2,
  crs = NA,
  buildings = NULL,
  quiet = FALSE
)
```

## Arguments

- case_dir:

  Character. OpenFOAM case directory.

- T_ref:

  Numeric. Reference temperature (K) for `T_cool`, which is `T_ref - T`.
  It has to match the temperature the case was initialised at, or every
  cooling number is offset by the difference and `T_cool` can come out
  negative everywhere. The default, `NULL`, reads `TRef` from the case's
  `constant/transportProperties` and only falls back to 295 K when that
  file cannot be read.

- time_step:

  `"latest"` (default) or a number of seconds.

- resolution:

  Numeric. Output cell size (m). Default 2.

- base_cell_size:

  Numeric. Background mesh cell size used when the case was generated
  (`foam$params$base_cell_size`); sets the gap-fill window. Default 10.

- agl_tol:

  Numeric. Samples more than this many metres above the local ground
  surface are discarded. Guards against the stray near-vertical sheet a
  distanceSurface can generate around a closed terrain STL, which
  otherwise folds upper-level wind into the pedestrian map. Default 5.

- trim:

  Crop the flow-adjustment zone off the returned raster. `"auto"`
  (default) removes, per side, the building-free apron plus an
  adjustment fetch; `"buildings"` removes only the apron; a number
  removes that many metres from every side; `0` returns the full domain.

  Trimming is on by default because the outer band is not a result. The
  domain is larger than the built area, and even where buildings reach
  the boundary the ABL inlet delivers an undisturbed profile that only
  slows as it works into the roughness. Measured on a real case, the
  outer 100 m averaged 1.86 m/s against 0.74 m/s in the interior - a
  bright rim that is an artifact of the domain, not a feature of the
  city. Excluding it is standard practice in urban CFD.

  Each side is trimmed by `max(apron, fetch)` - not their sum. The fetch
  is measured from the domain boundary and the apron is its first part,
  so they overlap rather than stack.

  The fetch is ~5x the median building height (100 m when the layer has
  no height column), clamped to 50-150 m and to 10% of the shorter
  domain span. That is deliberately less than the ~15x an internal
  boundary layer needs to fully equilibrate: on a measured case the edge
  excess was concentrated in the outer 100 m, and beyond that the
  band-to-band variation was genuine urban structure rather than an edge
  effect. The amount removed on each side is reported.

- canopy:

  Canopy height raster (path or SpatRaster) used to build a canopy
  overlay for
  [`plot_foam_map`](https://billbillbilly.github.io/gloBFPr/reference/plot_foam_map.md).
  Auto-detected from `constant/gloBFPr/rasters/canopy_height.tif` when
  `NULL`; pass `FALSE` to skip it. The result is attached as the
  `"canopy"` attribute, mirroring `"buildings"`.

- min_canopy_height:

  Numeric. Canopy cells below this height are ignored. Default 2 m,
  matching `min_tree_height` in
  [`prepare_openfoam_inputs`](https://billbillbilly.github.io/gloBFPr/reference/prepare_openfoam_inputs.md).

- crs:

  CRS to assign (e.g. `32617`). Default `NA`.

- buildings:

  Optional `sf` footprints in local coordinates; auto-detected from the
  case directory when `NULL`.

- quiet:

  Logical. Default `FALSE`.

## Value

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with six layers.

## See also

[`prepare_foam_case`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_case.md),
[`run_openfoam_docker`](https://billbillbilly.github.io/gloBFPr/reference/run_openfoam_docker.md)
