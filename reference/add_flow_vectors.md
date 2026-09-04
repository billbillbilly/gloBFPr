# Add wind / flow vector arrows to a foam map plot

Takes a `ggplot` object produced by
[`plot_foam_map`](https://billbillbilly.github.io/gloBFPr/reference/plot_foam_map.md)
and overlays velocity arrows sampled on a regular sub-grid.

## Usage

``` r
add_flow_vectors(
  p,
  r,
  spacing = 20,
  scale = 1,
  colour = "black",
  alpha = 0.6,
  linewidth = 0.25,
  arrow_size = 0.07,
  u_ref = NULL,
  quiet = FALSE
)
```

## Arguments

- p:

  A `ggplot` object (output of `plot_foam_map`).

- r:

  A `SpatRaster` that contains layers named `Ux` and `Uy`.

- spacing:

  Numeric. Arrow sub-grid spacing in the same units as the raster
  coordinates (usually metres). Default 20.

- scale:

  Numeric. Arrow length multiplier. Default 1.

- colour:

  Character. Arrow colour. Default `"black"`.

- alpha:

  Numeric. Arrow opacity (0-1). Default 0.6.

- linewidth:

  Numeric. Line thickness in mm. Default 0.25 (thin).

- arrow_size:

  Numeric. Arrowhead length in cm. Default 0.07 (tiny).

- u_ref:

  Numeric. Speed (m/s) the longest arrow represents. The default,
  `NULL`, normalises to the field's own maximum, which makes a
  near-stagnant field look as vigorous as a strong one; pass a fixed
  value to put several maps on one scale.

- quiet:

  Logical. Suppress the message reporting the arrow scale.

## Value

The same `ggplot` object with arrows added.
