# Plot a NoiseModelling-style road-noise map

Draws NoiseModelling isosurface polygons underneath roads and buildings
when available. If `x` does not contain official isosurfaces, it falls
back to an interpolated receiver surface.

## Usage

``` r
plot_noise_map(
  x,
  period = "DEN",
  field = "LAEQ",
  resolution = NULL,
  breaks = c(35, 40, 45, 50, 55, 60, 65, 70, 75),
  nodata = -99,
  palette = noise_map_palette(),
  road_col = "white",
  road_alpha = 0.35,
  road_lwd = 1.4,
  building_col = "black",
  legend = TRUE,
  legend_width = 0.26,
  legend_cex = 0.85,
  scalebar = TRUE,
  scalebar_unit = c("auto", "km", "m"),
  scalebar_cex = 0.7,
  north_arrow = TRUE,
  mar = c(0.2, 0.2, 0.2, 0.2),
  add = FALSE,
  ...
)
```

## Arguments

- x:

  A result from `get_noise_map(run = TRUE)`, a `gloBFPr_noise_surface`
  object, or an `sf` receiver noise map.

- period:

  Noise period to plot. Defaults to `"DEN"`.

- field:

  Noise field to plot. Defaults to `"LAEQ"`.

- resolution:

  Optional raster resolution in map units when `x` is an `sf` receiver
  map.

- breaks:

  Noise class breakpoints in dB.

- nodata:

  Values at or below this threshold are treated as no-data. Defaults to
  `-99`, NoiseModelling's no-result sentinel value.

- palette:

  Fill colors from quiet to loud.

- road_col:

  Road overlay color. Defaults to white with transparency. Use `NA` to
  omit roads.

- road_alpha:

  Road overlay alpha from `0` fully transparent to `1` opaque.

- road_lwd:

  Road overlay line width.

- building_col:

  Building fill color. Use `NA` to omit buildings.

- legend:

  Logical. Draw the dB(A) class legend in a dedicated panel to the right
  of the map. Defaults to `TRUE`.

- legend_width:

  Width of the legend panel relative to the map panel. Smaller values
  give the map more room. Defaults to `0.26`.

- legend_cex:

  Legend text size. Defaults to `0.85`.

- scalebar:

  Logical. Draw a distance scale bar just below the legend (or in the
  bottom-left of the map when `legend = FALSE`). Defaults to `TRUE`.
  Distances assume a projected CRS in metres; for geographic coordinates
  they are approximated at the map's mid-latitude.

- scalebar_unit:

  Scale bar unit: `"auto"` (default), `"km"`, or `"m"` to pick whichever
  keeps the label readable.

- scalebar_cex:

  Scale bar label size. Defaults to `0.7`.

- north_arrow:

  Logical. Draw a north arrow in the lower-left map area. Defaults to
  `TRUE`.

- mar:

  Margins (in lines) around the map panel. Defaults to a tight margin so
  the map fills the device.

- add:

  Logical. If `TRUE`, add to the current plot.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns the isosurface `sf` object or fallback
`gloBFPr_noise_surface` object used for plotting.
