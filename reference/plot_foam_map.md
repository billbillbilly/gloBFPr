# Plot an OpenFOAM pedestrian-level map

Convenience wrapper around `ggplot2` for visualising a single layer from
the raster returned by
[`sample_foam_slice`](https://billbillbilly.github.io/gloBFPr/reference/sample_foam_slice.md)
or
[`read_foam_pedestrian_slice`](https://billbillbilly.github.io/gloBFPr/reference/read_foam_pedestrian_slice.md).

## Usage

``` r
plot_foam_map(
  r,
  layer = "U_mag",
  title = NULL,
  palette = "YlOrRd",
  reverse = FALSE,
  buildings = NULL,
  canopy = NULL,
  legend_title = layer,
  max_u_ref = NULL,
  na_colour = "black",
  scalebar = TRUE,
  scalebar_unit = c("auto", "km", "m"),
  scalebar_cex = 0.7,
  north_arrow = TRUE
)
```

## Arguments

- r:

  A `SpatRaster` (single layer, or one layer will be selected via
  `layer`).

- layer:

  Character. Layer name to plot. Default `"U_mag"`.

- title:

  Character. Plot title. Default auto-generated.

- palette:

  Character. `hcl.colors` palette name. Default `"YlOrRd"`.

- reverse:

  Logical. Reverse palette direction. Default FALSE.

- buildings:

  Optional `sf` object of building footprints to overlay (in the same
  local coordinate system as the raster).

- canopy:

  Optional data frame of canopy cell centres (`x`, `y`) with a `"res"`
  attribute giving the cell size. Taken from the raster's `"canopy"`
  attribute when not supplied, which
  [`read_foam_pedestrian_slice`](https://billbillbilly.github.io/gloBFPr/reference/read_foam_pedestrian_slice.md)
  attaches whenever the case has a canopy height raster.

- legend_title:

  Character. Legend label. Default `layer`.

- max_u_ref:

  Numeric. If plotting `U_mag`, annotate the colour scale as a wind
  speed ratio by dividing by this reference speed. Default `NULL` (no
  ratio).

- na_colour:

  Colour for cells with no result. Inside the mapped area a void is an
  obstruction the flow never entered, so the default is the same solid
  black used for buildings and canopy.

- scalebar:

  Logical. Draw a distance scale bar in the lower-left corner. Defaults
  to `TRUE`. Distances assume projected map units are metres; for
  geographic rasters they are approximated at the map's mid-latitude.

- scalebar_unit:

  Scale bar unit: `"km"`, `"m"`, or `"auto"` to pick whichever keeps the
  label readable. Default `"auto"`.

- scalebar_cex:

  Scale bar label size multiplier. Defaults to `0.7`.

- north_arrow:

  Logical. Draw a north arrow above the scale bar. Defaults to `TRUE`.

## Value

A `ggplot` object.
