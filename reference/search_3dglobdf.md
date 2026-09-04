# search_3dglobdf

Search and retrieve 3D building footprint data from 3D-GloBFP or
GlobalBuildingAtlas that intersect a given bounding box or area of
interest (a city), with options to return vector or raster outputs
including building polygons, binary presence rasters, and height-coded
rasters.

## Usage

``` r
search_3dglobdf(
  bbox = NULL,
  place = NULL,
  crop = FALSE,
  data_source = "GBF",
  keep_source_id = FALSE,
  out_type = "poly",
  mask = FALSE,
  cell_size = 1,
  quiet = TRUE
)
```

## Arguments

- bbox:

  `sf`, `sfc`, or a numeric vector (xmin, ymin, xmax, ymax) defining the
  area of interest. This can be ignored if `place` is specified.

- place:

  vector (optional). A single line address, e.g. ("1600 Pennsylvania Ave
  NW, Washington") or a vector of addresses (c("Madrid", "Barcelona")).

- crop:

  logical. If `TRUE`, the resulting building footprint geometries will
  be cropped to the input `bbox`. Default is `FALSE`.

- data_source:

  character. Building data source to query. Use `"GBF"` for 3D-GloBFP
  (default) or `"GBA"`/`"gba"` for GlobalBuildingAtlas.

- keep_source_id:

  logical. If `TRUE`, keep the original source feature identifier as
  `source_id` when it is available. Default is `FALSE`.

- out_type:

  character. Default is `'poly'`. Output type(s) to return. Options
  include:

  - `"poly"`: building footprints as an `sf` polygon object.

  - `"binary_rast"`: binary `terra` raster where buildings = 1.

  - `"graduated_rast"`: `terra` raster encoding building height values.

  - `"rast"`: a named list with both binary and graduated rasters.

  - `"all"`: a named list including the polygon layer and both raster
    layers.

- mask:

  logical (optional). Default is `FALSE`. If `TRUE`, masks the graduated
  raster using the building footprint layer. Only used when `out_type`
  is `"graduated_rast"`, `"rast"`, or `"all"`.

- cell_size:

  numeric (optional). Default is 1. Only used when `out_type` is
  `"graduated_rast"`, `"rast"`, or `"all"`.

- quiet:

  logical. If `TRUE`, suppress cli messages and progress output. Default
  is `TRUE`.

## Value

Varies based on `out_type`:

- If `"poly"`: an `sf` object of building footprints. `MULTIPOLYGON`
  geometries are converted to `POLYGON` geometries while preserving one
  row per source feature. Polygons that touch or intersect share a
  `group_id`, which can be used to treat fragmented rows as one building
  group.

- If `"binary_rast"`: a binary `SpatRaster` (`terra`) indicating
  building presence.

- If `"graduated_rast"`: a quantitative `SpatRaster` of building
  heights.

- If `"rast"`: a named list with two `SpatRaster` objects: `binary` and
  `graduated`.

- If `"all"`: a named list with `poly` (sf), `binary`, and `graduated`
  rasters.

## Note

The downloading process may take some time, depending on the number and
size of building footprint tiles.

This implementation for gloBFP-3D relies on the current structure of the
dataset as hosted on Figshare. It may break if the dataset owner changes
the file organization or metadata format.

The server of GlobalBuildingAtlas may have issues sometimes, so users
may need to switch over to gloBFP-3D, which means using
`data_source="GBF"`

When using `data_source = "GBA"`, the GlobalBuildingAtlas dataset does
not provide unique identifiers for individual building parcels. As a
result, the `group_id` assigned to overlapping or fragmented polygons
may not accurately reflect true building boundaries in all cases. This
limitation may affect morphological analyses at the individual building
level (e.g., footprint area, perimeter, building-level height
statistics). However, the data remains suitable for environmental
simulation purposes such as noise mapping, solar radiation analysis, and
wind flow modelling, where parcel-level identity is less critical.

## References

Che Yangzi, Li Xuecao, Liu Xiaoping, Wang Yuhao, Liao Weilin, Zheng
Xianwei, Zhang Xucai, Xu Xiaocong, Shi Qian, Zhu Jiajun, Zhang Honghui,
Yuan Hua, & Dai Yongjiu (2024). 3D-GloBFP: the first global
three-dimensional building footprint dataset. Earth Syst. Sci. Data, 16,
5357-5374

Zhu X. X., Chen S., Zhang F., Shi Y., & Wang Y. (2025).
GlobalBuildingAtlas: an open global and complete dataset of building
polygons, heights and LoD1 3D models. Earth Syst. Sci. Data, 17,
6647-6668.

## Examples

``` r
if (FALSE) { # \dontrun{
buildings <- gloBFPr::search_3dglobdf(bbox=c(-84.485519,45.636118,-84.462774,45.650639))
} # }
```
