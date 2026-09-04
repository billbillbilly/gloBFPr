# get_metrics

`get_morphology`: Computes a set of morphological properties and
geometric descriptors for each building footprint polygon. These include
area, perimeter, surface area, volume, shape compactness, elongation,
and accessibility measures. If `x` contains a `group_id` column,
features with the same `group_id` are treated as one building for
morphology metrics while the returned rows and geometries remain
unchanged.

`get_neighbors`: Computes the number of neighboring buildings and
centroid distance summaries, based on a fixed-radius buffer and
Voronoi-based adjacency.

`get_bgvi`: Calculate the Building Green View Index (BGVI) for each
building volume. Neighboring buildings in the internal DSM used for the
viewshed are given a single flat roof elevation (base ground elevation
at their centroid, plus height) rather than following the terrain slope
pixel by pixel, matching
[`get_fused_dsm()`](https://billbillbilly.github.io/gloBFPr/reference/get_fused_dsm.md).

`get_dng`: Calculate the distance to the nearest greenspace (DNG) patch
for each building volume

## Usage

``` r
get_morphology(
  x = NULL,
  metrics = c("g_area", "pmeter", "v_surf", "t_surf", "vol", "obb_vol", "pa_ratio",
    "rec", "fra", "cbn", "hem", "cnv", "me_dist", "mp_dist", "vol_exch", "elo_x",
    "elo_y", "elo_z"),
  quiet = FALSE
)

get_neighbors(x = NULL, radius = 500, quiet = FALSE)

get_bgvi(
  x = NULL,
  datasource_canopy_height = "metachm",
  datasource_greenspace = NULL,
  min_tree_height = 2,
  zoom = 17,
  radius = 800,
  year = NULL,
  floor = FALSE,
  floor_step = 3,
  short_building_threshold = 6,
  field_of_view = 45,
  directions = NULL,
  workers = NULL,
  resolution = NULL,
  key = NULL,
  quiet = FALSE
)

get_dng(
  x = NULL,
  datasource = NULL,
  min_tree_height = 2,
  zoom = 17,
  radius = 800,
  min_area = 500,
  unit = c("m2", "ha", "km2"),
  year = NULL,
  network = NULL,
  overpass_url = "https://overpass-api.de/api/interpreter",
  timeout = 180,
  workers = NULL,
  quiet = FALSE
)
```

## Arguments

- x:

  sf. building footprint polygon, typically output from
  [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md)

- metrics:

  vector. A list of metrics to be computed. All metrics will be computed
  by default.

- quiet:

  logical. If `TRUE`, suppress cli messages and progress bars.

- radius:

  numeric. (only required for `get_neighbors`, `get_gbvi`, and
  `get_dng`) A numeric value specifying the buffer radius (in meters)
  used to define the proximity of a building footprint centroid. Default
  for `get_neighbors`, and `get_gbvi` and `get_dng` is respectfully 500
  and 800.

- datasource_canopy_height:

  character or `NULL`. Canopy height source for building the DSM and
  canopy green feature layer. Currently supports `"metachm"`,
  `"ethCHM"`, or `NULL`. If `NULL`, the DSM is built from buildings and
  DEM only.

- datasource_greenspace:

  character or `NULL`. Optional 2D greenspace map tile source for the
  visible-green feature layer. Supports `"esri"` and `"sentinel2"`. If
  both canopy height and greenspace sources are supplied, the
  visible-green layer is the union of height-filtered canopy and 2D
  greenspace.

- min_tree_height:

  numeric. (only required for `get_bgvi` and `get_dng`) When
  `datasource_canopy_height` is a canopy height source, minimum height
  threshold (in meters) to classify vegetation as trees in the CHM.
  Default is 2.

- zoom:

  numeric. (only required for `get_bgvi` and `get_dng`) Zoom level of
  map tile when `datasource_greenspace = "esri"` or
  `datasource_greenspace = "sentinel2"`. The default is `17`. The higher
  level of zoom will lead to higher resolution of greenspace data for
  computing BGVI or DNG.

- year:

  numeric. The desired year for Sentinel-2 cloudless mosaic tiles. (This
  has to be specified when `datasource_greenspace = "sentinel2"`)

- floor:

  logical. (only required for `get_bgvi`) Whether to compute Building
  Green View Index (BGVI) for each floor level based on estimated number
  of floors. Default is `FALSE`.

- floor_step:

  integer. (only required for `get_bgvi` when `floor = TRUE`) Compute
  GVI every `floor_step` floors. The top estimated floor is always
  included. Default is `1`, meaning every floor.

- short_building_threshold:

  numeric. (only required for `get_bgvi`) Height threshold in meters for
  deciding whether a building should use only the bottom viewpoint in
  non-floor mode. Default is 6.

- field_of_view:

  numeric. (only required for directional `get_bgvi`) Angular field of
  view in degrees for direction-specific GVI. Used only when
  `directions` is not `NULL`. Default is 45.

- directions:

  character vector or `NULL`. Optional direction names for
  direction-specific GVI. Valid values are `"southwest"`, `"southeast"`,
  `"northeast"`, `"northwest"`, `"north"`, `"east"`, `"west"`, and
  `"south"`. If `NULL`, `field_of_view` is ignored and only
  non-directional GVI columns are returned.

- workers:

  integer. (only required for `get_dng`) Number of parallel workers for
  per-building distance computation. Defaults to one fewer than
  available cores. Use `workers = 1` to run sequentially.

- resolution:

  numeric or `NULL`. (only used by `get_bgvi`) Output raster resolution
  in meters for the internal DSM used for viewshed computation. If
  `NULL` (default), the finest native resolution among the downloaded
  DEM and canopy height model is used, matching
  [`get_fused_dsm()`](https://billbillbilly.github.io/gloBFPr/reference/get_fused_dsm.md).
  Set explicitly (e.g. `1`) to force a finer grid.

- key:

  character. (only required for `get_bgvi`) API key of OpenTopography.

- datasource:

  character. (only required for `get_bgvi` and `get_dng`) Green/canopy
  data source. Supported values are `"metachm"`, `"esri"`, and
  `"sentinel2"` for `get_dng`.

- min_area:

  numeric. (only required for `get_dng`) The minimum area (the unit is
  defined by `unit`) of greenspace patches within the proximity of a
  building. For example, when `min_area = 500`, any greenspace patches
  with area less than 500 square meter will be excluded.

- unit:

  character. (only required for `get_dng`) The unit for `min_area`:
  'm2','ha', and 'km2'.

- network:

  `NULL`, character, or `sf`. (only required for `get_dng`) Controls
  whether distances are measured along a real street network instead of
  in a straight line. Use `NULL` (the default) for straight-line
  distance, `"osm"` to download a walkable OpenStreetMap network for the
  study extent, or supply your own `sf` line layer of road/path centre
  lines. When routing is requested, the distance is the sum of the walk
  from the building centroid to the network, the shortest path along the
  network, and the walk from the network to the green-space pixel.

- overpass_url:

  character. (only required for `get_dng` when `network = "osm"`)
  Overpass API endpoint used to download the street network.

- timeout:

  numeric. (only required for `get_dng` when `network = "osm"`) Overpass
  query timeout in seconds. Default is 180.

## Value

`get_morphology` returns an `sf` object identical to input `x`, with
additional columns for:

- `g_area`: Ground area of footprint

- `pmeter`: Perimeter length of footprint

- `v_surf`: Vertical surface area (walls)

- `t_surf`: Total surface area (walls + roof)

- `vol`: Volume (area \* height)

- `obb_vol`: Volume of oriented bounding box

- `pa_ratio`: Perimeter-area ratio

- `rec`: Rectangular compactness

- `fra`: Fractal dimension ratio

- `cbn`: Cuboidness index

- `hem`: Hemisphericality index

- `cnv`: Convexity index

- `me_dist`: Mean edge accessibility distance

- `mp_dist`: Mean pairwise distance within the footprint

- `vol_exch`: Volume exchange ratio

- `elo_x/elo_y/elo_z`: Elongation ratios (x/y/z axis)

`get_neighbors` returns an `sf` object identical to input `x`, with
additional columns for:

- `n_count`: Number of adjacent buildings within the radius.

- `m_ndist`: Mean distance from the building centroid to neighboring
  building centroids.

- `min_ndist`: Minimum distance to neighboring building centroids.

- `max_ndist`: Maximum distance to neighboring building centroids.

- `sd_ndist`: Standard deviation of those distances.

`get_bgvi` returns an `sf` object identical to input `x`, with
additional columns for:

- `mean_gvi`: Mean Building Green View Index value (0 to 1) from
  selected height(s).

- `bottom_gvi`: Building Green View Index from the bottom viewpoint, 1.7
  m above ground.

- `top_gvi`: Building Green View Index from the top viewpoint. For short
  buildings, this is equal to `bottom_gvi`.

- `bottom_green_area`: Visible green area (m\\^2\\) from the bottom
  viewpoint.

- `top_green_area`: Visible green area (m\\^2\\) from the top viewpoint.

- `mean_green_area`: Mean visible green area (m\\^2\\) across
  viewpoints.

- `min_gvi`, `max_gvi`, `sd_gvi`: Minimum, maximum, and standard
  deviation of Green View Index (GVI) (if `floor = TRUE`).

- `estimated_floors`: Estimated number of floors based on building
  height (if `floor = TRUE`).

`get_dng` returns an `sf` object identical to input `x`, with additional
columns: `dng`: Distance to nearest green space (tree canopy pixel) in
meters, measured along the street network when `network` is supplied;
and `dng_method`: either `"network"` or `"euclidean"`, recording how
each value was obtained. Buildings that cannot reach any green space
through the network fall back to straight-line distance and are flagged
`"euclidean"`.

## Details

When `floor = TRUE`, `get_bgvi()` estimates the number of floors from
building height and computes GVI from selected floor viewpoints. The
`floor_step` argument controls how densely floors are sampled. If
`floor_step = 1`, GVI is computed for every estimated floor. If
`floor_step > 1`, GVI is computed every `floor_step` floors to reduce
runtime; for example, `floor_step = 3` samples floors 1, 4, 7, and so
on. The top estimated floor is always included, even when it does not
fall on the step sequence.

If `directions` is provided, `get_bgvi()` also computes
direction-specific GVI by applying the requested `field_of_view` around
each direction within each computed viewshed. The viewshed itself is
computed once per viewpoint, then reused for each direction.
Direction-specific columns are named with the direction suffix, such as
`gvi_bottom_south`, `gvi_top_south`, `mean_gvi_south`, `min_gvi_south`,
`max_gvi_south`, and `sd_gvi_south`.

## Note

`x` must include a unique `id` field. If `x` includes `group_id`,
metrics that represent a whole building are computed on temporary
grouped features and copied back to the original rows.

To request an OpenTopography API token, please visit:
<https://portal.opentopography.org/requestService?service=api>

USGS 3DEP 1m/10m raster dataset is currently restricted to academic
users. Academic users can request access to these data via the
OpenTopography portal. Non-academic users can enquire about an
enterprise API key by emailing info@opentopography.org. See
OpenTopography Terms of Use for more information on appropriate use of
the API.

## References

Anna Labetski, Stelios Vitalis, Filip Biljecki, Ken Arroyo Ohori &
Jantien Stoter (2023): 3D building metrics for urban morphology.
International Journal of Geographical Information Science, 37(1): 36-67.
DOI: 10.1080/13658816.2022.2103818

Rachid Hamaina, Thomas Leduc, Guillaume Moreau. Towards Urban Fabrics
Characterization Based on Buildings Footprints. Jerome Gensel; Didier
Josselin; Danny Vandenbroucke. Bridging the Geographic Information
Sciences- International AGILE'2012 Conference, Avignon (France), April,
24-27, 2012, Springer Berlin Heidelberg, pp.327-346, 2012,
978-3-642-29062-6. 10.1007/978-3-642-290633_18. hal-01347299

Melchiorri, M., Freire, S., Schiavina, M. et al. The Multi-temporal and
Multi-dimensional Global Urban Centre Database to Delineate and Analyse
World Cities. Sci Data 11, 82 (2024).
https://doi.org/10.1038/s41597-023-02691-1

Essential background in Pesaresi, M. et al. (2024) "Advances on the
Global Human Settlement Layer by joint assessment of Earth Observation
and population survey data", International Journal of Digital Earth,
17(1).

Qi, L., Hu, Y., Bu, R., Xiong, Z., Li, B., Zhang, C., ... & Li, C.
(2024). Spatial-temporal patterns and influencing factors of the
Building Green View Index: A new approach for quantifying 3D urban
greenery visibility. Sustainable Cities and Society, 111, 105518.

## Examples

``` r
library(gloBFPr)
data(globfp_example)
result <- gloBFPr::get_morphology(globfp_example[c(1:3),], quiet = TRUE)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE

result <- gloBFPr::get_neighbors(globfp_example[c(1:3),], radius = 100)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#> ✔ Completed. Time taken: 0 seconds.
if (FALSE) { # \dontrun{
result <- gloBFPr::get_bgvi(globfp_example[c(1:3),],
                            datasource_canopy_height = "metachm",
                            datasource_greenspace = "esri",
                            key = "YOUR_opentopography_API_KEY")
} # }

result <- gloBFPr::get_dng(#globfp_example[c(1:3),],
                           datasource = "metachm",
                           unit = "m2")
#> ℹ Please input building footprint polygon.

# Measure along real road and path centre lines instead
result <- gloBFPr::get_dng(#globfp_example[c(1:3),],
                           datasource = "metachm",
                           unit = "m2",
                           network = "osm")
#> ℹ Please input building footprint polygon.
```
