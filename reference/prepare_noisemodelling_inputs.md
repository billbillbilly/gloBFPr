# Prepare NoiseModelling-style input layers

Builds a compact set of layers for an integrated or external
NoiseModelling workflow: `BUILDINGS`, `ROADS`, `GROUND`, `RECEIVERS`,
and optional `DEM`. The returned layers can be inspected directly in R
and optionally written to a GeoPackage.

## Usage

``` r
prepare_noisemodelling_inputs(
  x = NULL,
  height_field = "Height",
  datasource_canopy_height = NULL,
  datasource_greenspace = NULL,
  greenspace_year = NULL,
  greenspace_zoom = 17,
  opentopo_key = NULL,
  canopy_height = NULL,
  roads = NULL,
  download_roads = TRUE,
  greenspace = NULL,
  dem = NULL,
  population = FALSE,
  population_field = NULL,
  population_year = 2025,
  min_tree_height = 2,
  receiver = c("grid", "none"),
  resolution = 10,
  ground_default = 0,
  green_ground = 1,
  out_dir = NULL,
  write = FALSE,
  quiet = TRUE
)
```

## Arguments

- x:

  `sf` polygon object with building footprints and a height field.

- height_field:

  Building height column name. Defaults to `"Height"`.

- datasource_canopy_height:

  Character or `NULL`. Canopy height source to retrieve internally when
  `canopy_height` is not supplied. Currently supports `"metachm"` and
  `"ethCHM"`.

- datasource_greenspace:

  character or `NULL`. Optional 2D greenspace map tile source for the
  visible-green feature layer. Supports `"esri"` and `"sentinel2"`. If
  both canopy height and greenspace sources are supplied, the
  visible-green layer is the union of height-filtered canopy and 2D
  greenspace.

- greenspace_year:

  numeric. The desired year for Sentinel-2 cloudless mosaic tiles. (This
  has to be specified when `datasource_greenspace = "sentinel2"`)

- greenspace_zoom:

  numeric. Zoom level of map tile when `datasource_greenspace = "esri"`
  or `"sentinel2"`.

- opentopo_key:

  OpenTopography API key used to retrieve DEM data internally when `dem`
  is not supplied.

- canopy_height:

  Optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  canopy height map. Cells greater than or equal to `min_tree_height`
  are treated as green ground.

- roads:

  Optional `sf` line object, usually from OSM. If `NULL`, OSM roads are
  downloaded from the bounding box of `x`. If inferred traffic columns
  are absent,
  [`infer_osm_traffic()`](https://billbillbilly.github.io/gloBFPr/reference/infer_osm_traffic.md)
  is applied.

- download_roads:

  Logical. If `TRUE` and `roads = NULL`, download roads from OSM using
  the building extent. Set to `FALSE` only when a later workflow step
  supplies roads, for example `get_noise_map(osm_file = ...)`.

- greenspace:

  Optional `sf` polygons or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  marking green areas. Green areas are translated to ground absorption
  `G = 1`.

- dem:

  Optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  terrain layer.

- population:

  Logical. If `TRUE`, assign GHSL population to buildings with the
  package population helper (or `get_pop()` when available) before
  writing NoiseModelling `BUILDINGS`. Existing `POP` or
  `population_field` values are used when available.

- population_field:

  Optional column containing building population. The value is copied to
  NoiseModelling's `POP` field.

- population_year:

  GHSL population year passed to the package population function when
  `population = TRUE` and no usable population field is already present.

- min_tree_height:

  Minimum canopy height treated as green cover.

- receiver:

  One of `"grid"` or `"none"`.

- resolution:

  Receiver grid resolution in map units.

- ground_default:

  Default ground absorption for the analysis extent.

- green_ground:

  Ground absorption assigned to greenspace/canopy polygons.

- out_dir:

  Optional output directory for a GeoPackage.

- write:

  Logical. If `TRUE`, write layers to `noise_inputs.gpkg`.

- quiet:

  Logical. If `TRUE`, suppress informational messages.

## Value

A list with prepared `sf`/`terra` layers and optional GeoPackage path.
