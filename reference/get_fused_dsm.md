# get_fused_dsm

Generate digital surface model using multiple datasets, including
building height map, canopy height map, and terrain model. Each building
is given a single flat roof elevation (its own base ground elevation,
sampled at its centroid, plus its height) rather than following the
terrain slope beneath it pixel by pixel, so buildings on sloped ground
do not come out tilted or warped.

## Usage

``` r
get_fused_dsm(
  x = NULL,
  datasource_canopy_height = "metachm",
  min_tree_height = 2,
  resolution = NULL,
  opentopo_key = NULL,
  key = NULL,
  quiet = TRUE
)
```

## Arguments

- x:

  sf. building footprint polygon, typically output from
  [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md)

- datasource_canopy_height:

  character or `NULL`. Canopy height source. Currently supports
  `"metachm"`, `"ethCHM"`, or `NULL`.

- min_tree_height:

  numeric. Minimum canopy height threshold in meters.

- resolution:

  numeric or `NULL`. Output raster resolution in meters. If `NULL`
  (default), the finest native resolution among the downloaded DEM and
  canopy height model is used, so the output is never silently degraded
  to the coarser of the two source rasters. Set explicitly (e.g. `1`) to
  force a finer grid than the native source data (e.g. when the DEM
  falls back to coarse SRTM data), or a coarser one to speed up large
  areas.

- opentopo_key:

  character. OpenTopography API key used to download DEM data.

- key:

  Deprecated alias for `opentopo_key`.

- quiet:

  logical. If `TRUE`, suppress cli messages and progress output. Default
  is `TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
 example <- gloBFPr::globfp_example
 dsm <- get_fused_dsm(x= example, opentopo_key = 'key')
} # }
```
