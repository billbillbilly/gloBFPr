# Visualize an Individual Building BGVI Viewshed

Computes and optionally plots the BGVI viewshed for one building and one
viewpoint height. The function uses the same DSM, green feature layer,
target footprint flattening, and flat-roof handling as
[`get_bgvi()`](https://billbillbilly.github.io/gloBFPr/reference/get_metrics.md),
but returns the diagnostic layers for a single building instead of
summary columns for every building.

## Usage

``` r
plot_bgvi_viewshed(
  x,
  building = 1,
  level = c("bottom", "top"),
  floor = NULL,
  height = NULL,
  orientation = NULL,
  field_of_view = 45,
  datasource_canopy_height = "metachm",
  datasource_greenspace = NULL,
  min_tree_height = 2,
  zoom = 17,
  radius = 800,
  year = NULL,
  resolution = NULL,
  key = NULL,
  plot = TRUE,
  scalebar = TRUE,
  scalebar_unit = c("auto", "km", "m"),
  scalebar_cex = 0.7,
  north_arrow = TRUE,
  quiet = FALSE,
  ...
)
```

## Arguments

- x:

  sf. Building footprint polygons, typically output from
  [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md).
  Must include a `Height` column.

- building:

  Integer row number, or a value from the `id` column when `id` is
  present.

- level:

  Character. One of `"bottom"` or `"top"`. Ignored when `floor` or
  `height` is supplied.

- floor:

  Integer floor number to visualize. Floor 1 is 1.7 m above ground;
  higher floors add 3 m each.

- height:

  Numeric observer offset above ground in metres. Overrides `level` and
  `floor` when supplied.

- orientation:

  Optional sector orientation. Supply a bearing in degrees clockwise
  from north, or one of `"north"`, `"northeast"`, `"east"`,
  `"southeast"`, `"south"`, `"southwest"`, `"west"`, or `"northwest"`.
  If `NULL`, the full viewshed is used.

- field_of_view:

  Numeric angular width in degrees for `orientation`.

- datasource_canopy_height, datasource_greenspace, min_tree_height,
  zoom, radius, year, resolution, key, quiet:

  Passed to the BGVI raster preparation workflow; see
  [`get_bgvi()`](https://billbillbilly.github.io/gloBFPr/reference/get_metrics.md).

- plot:

  Logical. If `TRUE`, draw the viewshed map.

- scalebar:

  Logical. If `TRUE`, add a scale bar using the package's shared map
  layout.

- scalebar_unit:

  Character. Unit for the scale bar: `"auto"`, `"km"`, or `"m"`.

- scalebar_cex:

  Numeric text size for the scale bar and north arrow.

- north_arrow:

  Logical. If `TRUE`, add a north arrow.

- ...:

  Additional arguments passed to the initial viewshed
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A list containing the selected `building`, `viewpoint`, observer
`height`, `gvi`, `green_area`, `viewshed`, `viewshed_raster`,
`viewshed_area`, `radius`, `visible_green`, `sector_mask`,
`plot_raster`, `dsm`, and `binary_green`.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- plot_bgvi_viewshed(
  globfp_example,
  building = 1,
  level = "top",
  orientation = "south",
  field_of_view = 60,
  datasource_canopy_height = "metachm",
  key = "YOUR_opentopography_API_KEY"
)
} # }
```
