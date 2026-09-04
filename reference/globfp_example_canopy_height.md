# Example canopy height raster for the 3D-GloBFP sample

A `terra::PackedSpatRaster` canopy height map cropped to the bounding
box of `globfp_example` and aggregated from the source resolution for
lightweight examples. Convert it with
[`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html)
before analysis.

## Usage

``` r
globfp_example_canopy_height
```

## Format

A `terra::PackedSpatRaster` with one layer named `canopy_height`.

## Source

metaCHM canopy height data via dsmSearch, downloaded for the
`globfp_example` extent.

## Examples

``` r
data(globfp_example_canopy_height)
canopy_height <- terra::rast(globfp_example_canopy_height)
canopy_height
#> class       : SpatRaster
#> size        : 99, 106, 1  (nrow, ncol, nlyr)
#> resolution  : 17.65165, 17.65165  (x, y)
#> extent      : 329783.1, 331654.1, 4688653, 4690401  (xmin, xmax, ymin, ymax)
#> coord. ref. : WGS 84 / UTM zone 17N (EPSG:32617)
#> source(s)   : memory
#> name        : canopy_height
#> min value   :             0
#> max value   :            31
```
