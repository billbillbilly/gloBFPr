# Test 3D-GloBFP dataset

A sample dataset containing simplified 3D building footprint information
for demonstration and testing purposes.

## Usage

``` r
globfp_example
```

## Format

A data frame with 369 rows and 3 variables:

- id:

  Numeric. Unique identifier for each building.

- Height:

  Numeric. Estimated height of the building in meters.

- geometry:

  sfc_POLYGON. The building footprint geometry in simple feature (sf)
  format.

## Source

Che Yangzi, Li Xuecao, Liu Xiaoping, Wang Yuhao, Liao Weilin, Zheng
Xianwei, Zhang Xucai, Xu Xiaocong, Shi Qian, Zhu Jiajun, Zhang Honghui,
Yuan Hua, & Dai Yongjiu (2024). 3D-GloBFP: the first global
three-dimensional building footprint dataset. Earth Syst. Sci. Data, 16,
5357-5374

## Examples

``` r
data(globfp_example)
head(globfp_example)
#> Simple feature collection with 6 features and 2 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 329807.5 ymin: 4688864 xmax: 330110.9 ymax: 4689394
#> Projected CRS: WGS 84 / UTM zone 17N
#>          Height                       geometry id
#> 520676 19.54545 POLYGON ((330069.7 4689216,...  1
#> 520682 15.20905 POLYGON ((329807.5 4688930,...  2
#> 520683 12.94461 POLYGON ((329931.8 4689273,...  3
#> 520684 19.64550 POLYGON ((330060.6 4689241,...  4
#> 520685 15.00980 POLYGON ((329906.5 4689281,...  5
#> 520686 11.59947 POLYGON ((330035.8 4689366,...  6
```
