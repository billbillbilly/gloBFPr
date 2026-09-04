# Create a screening-level urban road-noise workflow object

Convenience wrapper around
[`prepare_noisemodelling_inputs()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_noisemodelling_inputs.md).
With `run = FALSE`, it returns reproducible input layers. With
`run = TRUE`, it runs the official headless NoiseModelling WPS scripts
and reads the `RECEIVERS_LEVEL` output table back into R. The returned
`noise_map` element joins the calculated levels to receiver geometries
for point mapping, and `isophones` contains NoiseModelling's
`CONTOURING_NOISE_MAP` polygons.

## Usage

``` r
get_noise_map(
  x = NULL,
  height_field = "Height",
  datasource_canopy_height = NULL,
  datasource_greenspace = NULL,
  greenspace_year = NULL,
  greenspace_zoom = 17,
  opentopo_key = NULL,
  roads = NULL,
  greenspace = NULL,
  canopy_height = NULL,
  dem = NULL,
  population = FALSE,
  population_field = NULL,
  population_year = 2025,
  min_tree_height = 2,
  receiver = c("grid", "none"),
  resolution = 10,
  out_dir = NULL,
  write = FALSE,
  run = FALSE,
  nm_path = NULL,
  nm_version = "5.0.1",
  download_nm = TRUE,
  osm_file = NULL,
  osm_remove_tunnels = TRUE,
  osm_eliminate_no_traffic_roads = TRUE,
  java = NULL,
  keep_files = FALSE,
  wall_alpha = 0.1,
  reflection_order = 0,
  max_src_distance = 150,
  max_reflection_distance = 50,
  thread_count = 0,
  diffraction_vertical = FALSE,
  diffraction_horizontal = FALSE,
  export_source_id = FALSE,
  humidity = NULL,
  temperature = NULL,
  favourable_occurrences = NULL,
  rays_name = NULL,
  max_error = NULL,
  frequency_field_prepend = "HZ",
  noise_wps_args = NULL,
  delaunay_max_area = NULL,
  road_width = 2,
  building_buffer = 2,
  iso_levels = c(35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 200),
  iso_field = "LAEQ",
  iso_smooth = 0.5,
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

- roads:

  Optional `sf` line object, usually from OSM. If `NULL`, OSM roads are
  downloaded from the bounding box of `x`. If inferred traffic columns
  are absent,
  [`infer_osm_traffic()`](https://billbillbilly.github.io/gloBFPr/reference/infer_osm_traffic.md)
  is applied.

- greenspace:

  Optional `sf` polygons or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  marking green areas. Green areas are translated to ground absorption
  `G = 1`.

- canopy_height:

  Optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  canopy height map. Cells greater than or equal to `min_tree_height`
  are treated as green ground.

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

- out_dir:

  Optional output directory for a GeoPackage.

- write:

  Logical. If `TRUE`, write layers to `noise_inputs.gpkg`.

- run:

  Logical. If `TRUE`, run NoiseModelling after preparing inputs.

- nm_path:

  Optional path to a `NoiseModelling_without_gui-*` directory. If
  `NULL`, the package uses `options(gloBFPr.noisemodelling.path)`,
  `NOISEMODELLING_HOME`, or a cached download.

- nm_version:

  NoiseModelling version to download when needed.

- download_nm:

  Logical. If `TRUE`, download the headless NoiseModelling runner when
  `nm_path` is not available.

- osm_file:

  Optional path to a local `.osm`, `.osm.gz`, or `.osm.pbf` extract.
  When supplied with `roads = NULL` and `run = TRUE`, NoiseModelling's
  `Import_OSM.groovy` creates the `ROADS` table internally using its OSM
  road defaults. This is different from the default R workflow, which
  downloads roads from the building bounding box.

- osm_remove_tunnels:

  Logical passed to NoiseModelling `Import_OSM` when `osm_file` is used.
  If `TRUE`, OSM roads tagged `tunnel=yes` are removed.

- osm_eliminate_no_traffic_roads:

  Logical passed to NoiseModelling `Import_OSM` when `osm_file` is used.
  If `TRUE`, keeps only road classes that NoiseModelling treats as
  traffic-bearing roads.

- java:

  Optional Java executable, Java home directory, or installed Java major
  version such as `17`. If `NULL`, the function uses `JAVA_HOME` or
  `java` on `PATH`. NoiseModelling 5.0.1 works with Java 11-21; Java 17
  is recommended.

- keep_files:

  Logical. If `TRUE`, keep the temporary NoiseModelling work directory
  and include it in the result.

- wall_alpha:

  Wall absorption coefficient passed to NoiseModelling.

- reflection_order:

  Reflection order passed to NoiseModelling.

- max_src_distance:

  Maximum source-receiver distance in meters.

- max_reflection_distance:

  Maximum reflection distance in meters.

- thread_count:

  Number of NoiseModelling worker threads. `0` lets NoiseModelling
  choose.

- diffraction_vertical:

  Logical. Passes `confDiffVertical`; enables diffraction around
  vertical edges. NoiseModelling notes that CNOSSOS-EU uses this mainly
  for rail and industrial sources.

- diffraction_horizontal:

  Logical. Passes `confDiffHorizontal`; enables diffraction over
  horizontal building/terrain edges.

- export_source_id:

  Logical. Passes `confExportSourceId`; if `TRUE`, receiver levels are
  kept by source identifier instead of merged across all sources. This
  is useful for source-contribution diagnostics and can greatly enlarge
  the output.

- humidity:

  Relative humidity percentage for atmospheric absorption. Passes
  `confHumidity`. Defaults to NoiseModelling's script default of `70`.

- temperature:

  Air temperature in degrees Celsius for atmospheric absorption. Passes
  `confTemperature`. Defaults to NoiseModelling's script default of
  `15`.

- favourable_occurrences:

  Probability of favourable propagation conditions by 16 wind-direction
  sectors, clockwise, where the north sector is the last value. Supply
  one value to recycle to all sectors or 16 values. The NoiseModelling
  default is sixteen `0.5` values.

- rays_name:

  Optional table name or file URL passed as `confRaysName`. When
  supplied, NoiseModelling exports propagation rays/attenuation details
  for advanced diagnostics. This can create very large outputs.

- max_error:

  Maximum allowed error in dB for pruning negligible source
  contributions. Passes `confMaxError`; NoiseModelling's default is
  `0.1`.

- frequency_field_prepend:

  Prefix for source spectral columns, passed as `frequencyFieldPrepend`.
  Defaults to `"HZ"` for columns such as `HZ1000`.

- noise_wps_args:

  Optional named list of additional raw arguments passed to
  `Noise_level_from_source.groovy`. Use this only for advanced
  NoiseModelling options not yet exposed directly.

- delaunay_max_area:

  Maximum Delaunay triangle area in square map units. Defaults to
  `resolution^2`. Smaller values create denser receiver/contour meshes
  and slower runs.

- road_width:

  Receiver exclusion distance around roads in meters for the
  NoiseModelling Delaunay grid.

- building_buffer:

  Receiver exclusion distance around buildings in meters for the
  NoiseModelling Delaunay grid.

- iso_levels:

  Numeric vector of isosurface breakpoints in dB passed to
  `Create_Isosurface.groovy`.

- iso_field:

  Result field used to build isosurfaces. Defaults to `"LAEQ"`.

- iso_smooth:

  Smoothing coefficient passed to `Create_Isosurface.groovy`.

- quiet:

  Logical. If `TRUE`, suppress informational messages.

## Value

Prepared noise input layers when `run = FALSE`; otherwise a list with
prepared inputs, receiver-level results, output file paths, and logs.
