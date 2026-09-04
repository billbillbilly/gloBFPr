# Building shadow and radiation calculations

`svf()` computes a Sky View Factor raster from building and optional
canopy obstacles.

`get_shadow_footprint()` computes ground shadow footprints for extruded
building polygons.

`get_shadow_height()` computes shadow height at points or across a
`terra` surface. If `shadow_locations` is omitted, a `terra` template is
generated around the buildings.

`get_radiation()` estimates direct, diffuse, and total radiation load on
roofs and facades represented by a 3D `sf` surface grid.

## Usage

``` r
svf(
  x = NULL,
  height_field = "Height",
  min_tree_height = 2,
  datasource_canopy_height = NULL,
  key = NULL,
  canopy_height = NULL,
  dem = NULL,
  raster_buffer = NULL,
  grid_res = 2,
  extent_buffer = NULL,
  res_angle = 5,
  observer_height = 1.7,
  max_distance = NULL,
  plot = FALSE,
  scalebar = TRUE,
  scalebar_unit = c("auto", "km", "m"),
  scalebar_cex = 0.7,
  north_arrow = TRUE,
  quiet = TRUE
)

get_shadow_footprint(
  x = NULL,
  solar_time = NULL,
  time_zone = NULL,
  azimuth = NULL,
  elevation = NULL,
  height_field = "Height",
  min_tree_height = 2,
  datasource_canopy_height = NULL,
  key = NULL,
  canopy_height = NULL,
  dem = NULL,
  raster_buffer = NULL,
  b = 0.01,
  overlap_shadow = FALSE,
  plot = FALSE,
  plot_overlap_gradient = FALSE,
  scalebar = TRUE,
  scalebar_unit = c("auto", "km", "m"),
  scalebar_cex = 0.7,
  north_arrow = TRUE,
  quiet = TRUE
)

get_shadow_height(
  x = NULL,
  shadow_locations = NULL,
  solar_time = NULL,
  time_zone = NULL,
  azimuth = NULL,
  elevation = NULL,
  height_field = "Height",
  min_tree_height = 2,
  datasource_canopy_height = NULL,
  key = NULL,
  raster_buffer = NULL,
  canopy_height = NULL,
  dem = NULL,
  cell_size = 2,
  extent_buffer = NULL,
  b = 0.01,
  filter_footprint = FALSE,
  quiet = TRUE
)

get_radiation(
  x = NULL,
  grid = NULL,
  solar_time = NULL,
  time_zone = NULL,
  azimuth = NULL,
  elevation = NULL,
  solar_normal,
  solar_diffuse,
  height_field = "Height",
  min_tree_height = 2,
  datasource_canopy_height = NULL,
  key = NULL,
  raster_buffer = NULL,
  canopy_transmissivity = 0.15,
  canopy_height = NULL,
  dem = NULL,
  grid_res = 2,
  ground = FALSE,
  ground_res = NULL,
  offset = 0.01,
  radius = 500,
  svf_res_angle = 15,
  return_list = FALSE,
  plot = FALSE,
  plot_3d = FALSE,
  scalebar = TRUE,
  scalebar_unit = c("auto", "km", "m"),
  scalebar_cex = 0.7,
  north_arrow = TRUE,
  quiet = TRUE
)
```

## Arguments

- x:

  An `sf` polygon object with building footprints and a height field.

- height_field:

  Character. Name of the building height column. Defaults to `"Height"`,
  matching
  [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md)
  output.

- min_tree_height:

  Numeric. Minimum canopy height, in map units, used as a tree obstacle.

- datasource_canopy_height:

  Character or `NULL`. Canopy height source to retrieve internally when
  `canopy_height` is not supplied. Currently supports `"metachm"`,
  `"ethCHM"`, or `NULL`.

- key:

  Character or `NULL`. OpenTopography API key used to retrieve DEM
  internally when `dem` is not supplied.

- canopy_height:

  Optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  canopy height map. Values are interpreted as height above ground.

- dem:

  Optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  digital elevation model. When supplied, canopy and building shadows
  are compared in absolute elevation and shadow-height outputs are
  returned above local ground.

- raster_buffer:

  Numeric or `NULL`. Buffer distance in CRS units around buildings used
  when retrieving CHM/DEM internally. If `NULL`, a buffer is estimated
  from building height and solar elevation.

- grid_res:

  Numeric surface-grid resolution in CRS units.

- extent_buffer:

  Optional numeric buffer around `x` used when creating an automatic
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  template. If omitted, a buffer is estimated from building heights and
  solar elevation.

- res_angle:

  Numeric. Azimuth sampling interval in decimal degrees for `svf()`.
  Smaller values are slower and more detailed.

- observer_height:

  Numeric. Height above local ground for SVF query locations. Defaults
  to `1.7`, representing pedestrian eye level.

- max_distance:

  Numeric. Maximum obstacle search distance in CRS units for `svf()`.

- plot:

  Logical. For `get_shadow_footprint()`, draw a base R map of the
  building footprints and shadow polygons before returning the `sf`
  result. For `get_radiation()`, draw the default 2D base R radiation
  map colored by `total`. When ground samples are included, the 2D
  layout shows separate ground, facade, and roof maps with one shared
  legend. When canopy data are supplied, a second 2D map shows canopy
  impact as `canopy - no_canopy` total-radiation difference.

- scalebar:

  Logical. Draw a distance scale bar using the shared map layout.
  Defaults to `TRUE` for plotted maps.

- scalebar_unit:

  Scale bar unit: `"auto"` (default), `"km"`, or `"m"`.

- scalebar_cex:

  Scale bar label size. Defaults to `0.7`.

- north_arrow:

  Logical. Draw a north arrow in the map panel. Defaults to `TRUE`.

- quiet:

  Logical. If `FALSE`, emit progress messages.

- solar_time:

  Character vector or list of character strings. Local solar times such
  as `"2026-06-21 15:00:00"`. If `solar_time` and `time_zone` are
  supplied, `azimuth` and `elevation` are ignored and solar position is
  estimated from time and the building-layer centroid.

- time_zone:

  Character. A single time zone used to interpret `solar_time`, for
  example `"America/Denver"` or `"UTC"`.

- azimuth:

  Numeric vector or list. Solar azimuth in decimal degrees, measured
  clockwise from north. Must have the same length as `elevation`.

- elevation:

  Numeric vector or list. Solar elevation in decimal degrees above the
  horizon. Must have the same length as `azimuth`.

- b:

  Numeric buffer tolerance used when cleaning footprint unions.

- overlap_shadow:

  Logical. For `get_shadow_footprint()`, if `TRUE`, dissolve overlapping
  shadows across all supplied solar positions by shadow source.

- plot_overlap_gradient:

  Logical. For `get_shadow_footprint()` plots with multiple `solar_time`
  values, if `TRUE`, draw all shadows in transparent gray so overlapping
  areas appear darker.

- shadow_locations:

  Optional query locations for shadow height, as an `sf` point layer or
  a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

- cell_size:

  Numeric cell resolution in CRS units when `shadow_locations` is
  omitted.

- filter_footprint:

  Ignored. Shadow footprints are always used to limit height
  calculations.

- grid:

  Optional 3D `sf` point surface grid. If omitted, it is created from
  building roofs and facades (and optionally the ground).

- solar_normal:

  Direct Normal Irradiance vector, one value per solar position.

- solar_diffuse:

  Diffuse Horizontal Irradiance vector, one value per solar position.

- canopy_transmissivity:

  Numeric from 0 to 1. Fraction of direct irradiance transmitted through
  canopy shadows in `get_radiation()`.

- ground:

  Logical. If `TRUE`, add a regular grid of ground-level sample points
  over the study-area bounding box (excluding building footprints).
  Ground points have an upward normal and receive direct radiation
  whenever they are not in a building or canopy shadow, and diffuse
  radiation scaled by their Sky View Factor. Returned rows have
  `surface = "ground"` and `building_id = NA`.

- ground_res:

  Numeric resolution for the ground sample grid in CRS units. If `NULL`,
  defaults to `grid_res`.

- offset:

  Numeric vertical offset added to generated surface-grid points.

- radius:

  Maximum obstacle search distance in CRS units for radiation Sky View
  Factor calculations. Defaults to `500`. Obstacles beyond this distance
  contribute negligibly to Sky View Factor but dominate runtime, so a
  finite radius enables spatial culling and is typically many times
  faster. Use `Inf` to consider all obstacles regardless of distance.

- svf_res_angle:

  Numeric. Azimuth sampling interval in decimal degrees used when
  estimating Sky View Factor inside `get_radiation()`.

- return_list:

  Logical. If `TRUE`, return per-timestep radiation matrices instead of
  a summed `sf` surface grid.

- plot_3d:

  Logical. For `get_radiation()`, draw a base R 3D-style view with
  separate panels for direct, diffuse, and total radiation. This is
  opt-in; `plot = TRUE` uses the 2D map layout by default.

## Value

`get_shadow_footprint()` returns an `sf` polygon layer.

`get_shadow_height()` returns a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
for `terra` locations or a numeric matrix for point locations.

`get_radiation()` returns an `sf` point layer with `svf`, `direct`,
`diffuse`, and `total` columns, unless `return_list = TRUE`.

`svf()` returns a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with Sky View Factor values from 0 to 1.

## Details

These functions are implemented directly with `sf` and `terra` using a
projected 2.5D building model.

## References

Dorman, M. et al. `shadow`: Geometric Shadow Calculations.
<https://github.com/michaeldorman/shadow>
