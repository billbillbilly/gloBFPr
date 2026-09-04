# get_3d_world

Export a study area as a 3D scene (inspired by the
[arnis](https://github.com/louis-e/arnis) project): terrain surface,
extruded building footprints, and canopy tree objects, written as
Wavefront OBJ and/or binary STL files that load directly in Rhino3D and
Blender.

The OBJ output contains named objects per layer (`terrain`, `buildings`,
`canopy_trunks`, `canopy_crowns`) and one group per building
(`building_<id>`), so individual buildings remain selectable after
import. Ground surface classes (greenspace, roads, sidewalks, paths,
water, sand) are part of the terrain object itself, as face groups with
their own materials - the terrain is one continuous, detailed surface
with differently colored regions, not a stack of separate ribbons. STL
has no object concept, so one STL file is written per layer.

## Usage

``` r
get_3d_world(
  x = NULL,
  bbox = NULL,
  place = NULL,
  terrain = TRUE,
  canopy = "metachm",
  key = NULL,
  dem = NULL,
  canopy_height = NULL,
  height_col = "Height",
  format = c("obj", "stl"),
  out_dir = "world3d",
  min_tree_height = 2,
  tree_window = 5,
  max_trees = 20000,
  simplify_terrain = 1,
  color_by = NULL,
  facade_palette = FALSE,
  roads = NULL,
  overture_release = NULL,
  bridges = TRUE,
  bridge_level_height = 6,
  bridge_ramp = 20,
  bridge_pillar_interval = 25,
  water = NULL,
  sand_buffer = 3,
  surface_res = 2,
  greenspace = NULL,
  greenspace_under_canopy = TRUE,
  greenspace_zoom = 17,
  greenspace_year = NULL,
  all_vox = FALSE,
  vox_size = 1,
  basemap = FALSE,
  local_origin = TRUE,
  crop = FALSE,
  data_source = "GBF",
  quiet = TRUE
)
```

## Arguments

- x:

  sf. Building footprint polygons, typically the output of
  [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md)
  (must contain the column given by `height_col`). If `NULL`, buildings
  are fetched via
  [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md)
  using `bbox`/`place`.

- bbox:

  `sf`, `sfc`, or numeric vector (xmin, ymin, xmax, ymax) in WGS84.
  Ignored when `x` is provided.

- place:

  character (optional). Address or place name, passed to
  [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md).
  Ignored when `x` is provided.

- terrain:

  logical. If `TRUE` (default), download a DEM (OpenTopography via
  `dsmSearch`; requires `key`) and mesh the ground surface. Buildings
  and trees are then placed at their sampled ground elevation. If
  `FALSE`, the scene has a flat ground at z = 0 and no API key is
  needed.

- canopy:

  character or `NULL`. Canopy height source, `"metachm"` (default) or
  `"ethchm"`; `NULL` skips trees.

- key:

  character. OpenTopography API key, required when `terrain = TRUE` and
  no `dem` raster is supplied.

- dem:

  `SpatRaster` or `PackedSpatRaster` (optional). A pre-loaded digital
  elevation model (e.g. `globfp_example_dem`). When supplied with
  `terrain = TRUE`, no download or API key is needed. Ignored when
  `terrain = FALSE`.

- canopy_height:

  `SpatRaster` or `PackedSpatRaster` (optional). A pre-loaded canopy
  height model (e.g. `globfp_example_canopy_height`). When supplied, it
  is used for tree detection instead of downloading from the `canopy`
  source.

- height_col:

  character. Building height column. Default `"Height"`.

- format:

  character. Any of `"obj"`, `"stl"`. Default writes both.

- out_dir:

  character. Output directory (created if missing).

- min_tree_height:

  numeric. Minimum canopy height in metres. Default 2.

- tree_window:

  numeric. Local-maximum search window in **metres** for tree detection:
  at most one tree is placed per window. Larger values give fewer, more
  widely spaced trees; lower it (e.g. 3) if dense canopy looks too
  sparse in the model. Default 5.

- max_trees:

  integer. Cap on tree objects; the tallest trees are kept and a warning
  reports how many were detected. Default 20000.

- simplify_terrain:

  integer. Aggregation factor for the terrain mesh (1 = full DEM
  resolution). Default 1.

- color_by:

  character or `NULL`. Name of a numeric column in `x` mapped to
  per-building OBJ materials via a viridis ramp (OBJ only; STL carries
  no color). Default `NULL`.

- facade_palette:

  logical or character vector. If `TRUE`, buildings get varied muted
  facade colors (arnis-style) instead of uniform grey, assigned
  deterministically per building. Supply a character vector of R colors
  to use a custom palette. Ignored when `color_by` is set. Default
  `FALSE`.

- roads:

  `NULL`, `"overture"`, or an sf line layer. Paints road and sidewalk
  surfaces into the terrain using the arnis method (class-based widths;
  asphalt/concrete/dirt terrain colors). `"overture"` fetches
  transportation segments from Overture Maps via a native DuckDB parquet
  query (requires the suggested `duckdb` and `DBI` packages and
  internet, matching
  [`generate_block()`](https://billbillbilly.github.io/gloBFPr/reference/generate_block.md)'s
  road source). An sf layer should contain LINESTRING geometries with
  optional `class` and `subclass` columns (Overture/OSM highway
  classes); missing classes default to `residential`. Default `NULL` (no
  roads).

- overture_release:

  character or `NULL`. Overture Maps release string used when
  `roads = "overture"` or `water = "overture"`, e.g. `"2025-03-19.0"`.
  The default `NULL` (or `"auto"`) queries the public bucket for the
  newest release and caches it for the session, so the code keeps
  working as Overture publishes new releases. See
  <https://github.com/OvertureMaps/data/releases> for the list.

- bridges:

  logical. If `TRUE` (default), road segments flagged as bridges are
  built as elevated structures instead of being painted on the ground: a
  deck raised `level * bridge_level_height` metres above the highest
  terrain along the span, ramps down to grade at each end, railings, and
  support pillars. Tunnel segments are always omitted from the ground
  surface. Requires bridge attributes in the road data (Overture
  supplies them; an sf layer may carry `is_bridge`, `is_tunnel`, and
  `level`).

- bridge_level_height:

  numeric. Metres of clearance per level, the arnis `LAYER_HEIGHT_STEP`.
  Default 6.

- bridge_ramp:

  numeric. Length in metres over which a deck ramps down to ground level
  at each end. Default 20.

- bridge_pillar_interval:

  numeric. Spacing in metres between support pillars. Default 25.

- water:

  `NULL`, `"overture"`, or an sf (MULTI)POLYGON layer of water bodies.
  Rivers, lakes, and coastal water are handled with the arnis method:
  the terrain under each water body is flattened to just below its
  lowest bank and colored as water, and a sand fringe of `sand_buffer`
  metres is painted along the shore (riparian and coastal edges).
  `"overture"` fetches water polygons from the Overture Maps base theme
  via the native DuckDB query. Default `NULL`.

- sand_buffer:

  numeric. Width in metres of the sand fringe around water bodies. 0
  disables the fringe. Default 3.

- surface_res:

  numeric. Grid resolution in metres for the classified terrain surface
  in the default (non-voxel) mode. Smaller values follow road edges more
  precisely at the cost of a denser mesh. Default 2.

- greenspace:

  `NULL`, `TRUE`, `"esri"`, `"sentinel2"`, or a `SpatRaster`. Paints
  ground-level greenery (lawns) into the terrain. Greenery is classified
  from map tiles by the suggested `greenSD` package: `TRUE` or `"esri"`
  uses Esri imagery, `"sentinel2"` uses Sentinel-2. Alternatively supply
  your own binary `SpatRaster` (green = 1). Cells where the canopy
  height model is at or above `min_tree_height` are excluded, so this
  surface class shows ground-level greenery only - tree canopy is
  already represented by the 3D tree objects. Default `NULL`.

- greenspace_under_canopy:

  logical. If `TRUE` (default), ground under tree canopy stays green -
  trees stand on grass, as in arnis. Set `FALSE` to cut canopy cells out
  of the greenspace surface, leaving bare ground beneath the crowns.

- greenspace_zoom:

  integer. Map-tile zoom level for the greenspace classification; higher
  zoom gives finer lawn edges. Default 17.

- greenspace_year:

  integer or `NULL`. Imagery year for `greenspace = "sentinel2"`.

- all_vox:

  logical. If `TRUE`, voxelize the entire scene so the model reproduces
  the blocky world arnis builds, from this package's data inputs:
  terrain becomes stepped blocks (Minecraft-style heightmap), buildings
  become block columns on flattened quantized bases, and trees are
  already voxel objects. Building footprints are rasterized at
  `vox_size`, so exact footprint edges are traded for the voxel look.
  Default `FALSE`.

- vox_size:

  numeric. Block edge length in metres for `all_vox` mode (1 =
  Minecraft-like scale). Larger values give chunkier, lighter models.
  Default 1.

- basemap:

  logical. If `TRUE`, download an Esri World Imagery snapshot of the
  scene, save it next to the OBJ, and drape it on the ground via texture
  coordinates (OBJ only; requires internet). With `terrain = FALSE` a
  flat textured ground plane is added. On download failure the export
  falls back to the flat terrain color with a warning. Default `FALSE`.

- local_origin:

  logical. If `TRUE` (default), shift the scene so its minimum x/y is at
  (0, 0). The offset and CRS are stored in `world_metadata.json` and in
  the returned object, so results can be georeferenced back. Strongly
  recommended for Rhino/Blender, which lose float precision at UTM-scale
  coordinates.

- crop, data_source:

  Passed to
  [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md)
  when `x` is `NULL`.

- quiet:

  logical. Suppress progress messages. Default `TRUE`.

## Value

(invisibly) a list:

- `paths`: character vector of written files.

- `meshes`: named list of in-memory meshes (`vertices`/`faces`).

- `origin`: list with `x`, `y` offset and `epsg` of the scene CRS.

- `n_buildings`, `n_trees`: scene composition counts.

- `n_trees_detected`: trees found in the canopy height model before the
  `max_trees` cap - compare with `n_trees` to see how many were dropped.

- `n_bridges`: elevated bridge segments built.

## Details

All coordinates are in metres (scene UTM zone). Building walls extend 1
m below their sampled ground elevation when `terrain = TRUE` to avoid
gaps on slopes. Trees are blocky voxel-style objects whose shapes are
ported from the arnis project (trunk column, leaf columns, apex cap, and
concentric canopy rings). Every tree's total height comes from the
canopy height raster; the variant is chosen from that height (compact
oaks below 8 m, standard/bushy oaks to 14 m, oak/spruce to 22 m,
towering spruce above, with position-hashed variety within each class),
then the block size is scaled so the tree matches its measured height
exactly.

Ground surface classes are painted into the terrain, not built as
separate geometry: a classification grid (`surface_res` metres in
default mode, `vox_size` in voxel mode) assigns each terrain cell one of
ground, greenspace, roads, sidewalks, paths, sand, or water, and the
terrain mesh carries these as per-face material groups. Roads follow the
arnis method: each segment gets a class-based half-width
(motorway/trunk/primary 5 m, secondary 4 m, tertiary 3 m,
residential/service 2 m, foot/cycle/path classes 1 m - the arnis
`highway_block_range` defaults), is buffered into a ribbon, and painted
as `roads` (asphalt), `sidewalks` (concrete; includes Overture
`subclass = "sidewalk"`), or `paths` (dirt).

Bridges and elevated highways follow arnis's structural approach.
Neither OSM nor Overture records deck elevations, so the height is
inferred from the segment's `level`: the deck sits
`level * bridge_level_height` metres (arnis uses 6 blocks per layer)
above the highest terrain along the span, ramps linearly down to grade
over `bridge_ramp` metres at each end, and is carried by pillars dropped
to the ground every `bridge_pillar_interval` metres. Decks are solid
slabs with railings, sized by the same class-based widths as surface
roads. Bridge and tunnel segments are excluded from the terrain surface
classification, so a bridge no longer leaves a road painted on the
ground beneath it.

Riparian and coastal areas also follow arnis: the terrain under each
water body is flattened to one level just below its lowest bank (0.2 m
in default mode, one block in voxel mode), painted as water, and fringed
with a `sand_buffer`-metre sand strip along the shoreline.

Vegetation is split across two levels so canopy and lawns are never
conflated: tree canopy comes from the canopy height model and is
rendered as 3D tree objects, while the `greenspace` terrain class
carries the ground surface. Greenery comes from `greenSD`'s map-tile
classification (Esri or Sentinel-2 imagery at `greenspace_zoom`). Ground
under canopy stays grass by default, so a park reads as a continuous
lawn with trees standing on it - exactly how arnis places tree objects
on grass blocks. Use `greenspace_under_canopy = FALSE` if you instead
want canopy cells cut out of the lawn.

With `all_vox = TRUE` the whole scene is voxelized the way arnis builds
its Minecraft worlds: the DEM becomes a stepped block heightmap, each
building footprint is rasterized at `vox_size` and raised as quantized
block columns on a flattened base (extending one block below ground),
and the voxel trees stand on the quantized ground. Heights still come
from the input data (building `Height` column, CHM tree heights, DEM
terrain) - only the geometry representation changes. Expect larger files
than the default mode; increase `vox_size` to lighten them.

When `basemap = TRUE`, imagery is retrieved from the Esri World Imagery
service; check the [Esri terms of
use](https://www.esri.com/en-us/legal/terms/master-agreement) and
provide attribution (Esri, Maxar, Earthstar Geographics, and the GIS
User Community) when publishing rendered scenes. Footprint holes
(courtyards) are not carved in this version; the exterior ring is
extruded and a warning is issued when holes are dropped.

## Why some trees are missing

Three filters reduce canopy pixels to tree objects. `min_tree_height`
(default 2 m) drops low vegetation; `tree_window` keeps only one tree
per window, so closely spaced trees merge into their tallest neighbour;
and `max_trees` caps the total, keeping the tallest. Compare
`n_trees_detected` with `n_trees` in the returned object, and lower
`tree_window` or raise `max_trees` if the model looks sparser than the
canopy data. Note that `tree_window` cannot resolve trees closer
together than the CHM resolution: on a coarse (e.g. aggregated) canopy
raster, set `tree_window` at or below the cell size to place one tree
per canopy cell.

## Examples

``` r
# \donttest{
example <- gloBFPr::globfp_example
world <- get_3d_world(
  x = example, terrain = FALSE, canopy = NULL,
  out_dir = tempfile("world3d")
)
# }
```
