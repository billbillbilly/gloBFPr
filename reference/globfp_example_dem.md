# Example DEM for the 3D-GloBFP sample

A `terra::PackedSpatRaster` digital elevation model cropped to the
bounding box of `globfp_example`. Convert it with
[`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html)
before analysis.

## Usage

``` r
globfp_example_dem
```

## Format

A `terra::PackedSpatRaster` with one layer named `dem`.

## Source

OpenTopography / dsmSearch elevation data, downloaded for the
`globfp_example` extent.

## Examples

``` r
data(globfp_example_dem)
dem <- terra::rast(globfp_example_dem)
dem
#> class       : SpatRaster
#> size        : 69, 74, 1  (nrow, ncol, nlyr)
#> resolution  : 25.69228, 25.69228  (x, y)
#> extent      : 329769.3, 331670.5, 4688640, 4690413  (xmin, xmax, ymin, ymax)
#> coord. ref. : WGS 84 / UTM zone 17N (EPSG:32617)
#> source(s)   : memory
#> name        :        dem
#> min value   : 175.163025
#> max value   : 211.040695
```
