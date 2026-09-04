# Prepare gloBFPr spatial outputs for an OpenFOAM Docker workflow

Collects building, terrain, canopy, and ground-cover data from gloBFPr
and writes them into a structured OpenFOAM case folder ready for
external CFD simulation (via Docker).

Output layers and their intended role in OpenFOAM:

- building STL:

  Solid geometry for snappyHexMesh

- building height raster:

  Auxiliary mesh-generation reference

- binary building raster:

  Building mask / validation

- fused DSM:

  Optional ground STL source (terrain + buildings + canopy surface)

- canopy height raster:

  Standalone CHM for defining porous-zone extents and drag coefficients
  in fvOptions / topoSetDict

- ground roughness raster (z0):

  Per-cell aerodynamic roughness length derived from ESA WorldCover land
  cover, for nutURoughWallFunction. Building footprint cells and
  tree-cover cells are set to NA because they are handled by solid
  geometry and porous zones respectively.

## Usage

``` r
prepare_openfoam_inputs(
  case_dir,
  bbox = NULL,
  place = NULL,
  buildings_list = NULL,
  data_source = "GBF",
  cell_size = 1,
  crop = FALSE,
  mask = TRUE,
  include_buildings = TRUE,
  include_fused_dsm = TRUE,
  include_tree_canopy = TRUE,
  canopy_source = "metachm",
  min_tree_height = 2,
  include_morphology = TRUE,
  include_neighbors = TRUE,
  include_greenspace = TRUE,
  mask_tree_cover = TRUE,
  landcover_source = c("esa", "esri"),
  landcover_year = 2021,
  include_bgvi = FALSE,
  include_shadow = FALSE,
  include_radiation = FALSE,
  opentopo_key = NULL,
  target_crs = NULL,
  height_col = NULL,
  default_height = 10,
  min_height = 2,
  domain_buffer = 100,
  zmax_buffer = 50,
  stl_name = "buildings.stl",
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- case_dir:

  Character. OpenFOAM case directory.

- bbox:

  Numeric vector c(xmin, ymin, xmax, ymax) in WGS-84 lon/lat.

- place:

  Optional place name passed to search_3dglobdf().

- buildings_list:

  Optional precomputed output from search_3dglobdf(..., out_type =
  "all").

- data_source:

  Character. Building source passed to search_3dglobdf(). Default "GBF".

- cell_size:

  Numeric. Raster resolution in metres for building rasters.

- crop:

  Logical. Whether to crop buildings to bbox.

- mask:

  Logical. Whether to mask height raster by building footprints.

- include_buildings:

  Logical. Prepare building vector/raster/STL.

- include_fused_dsm:

  Logical. Prepare fused DSM.

- include_tree_canopy:

  Logical. Extract standalone canopy height raster (CHM) and include
  canopy in fused DSM. The CHM is written as a separate file so OpenFOAM
  porous-zone definitions can reference it directly.

- canopy_source:

  Character or NULL. "metachm", "ethCHM", or NULL.

- min_tree_height:

  Numeric. Minimum tree canopy height in metres.

- include_morphology:

  Logical. Add morphology metrics.

- include_neighbors:

  Logical. Add neighbour metrics.

- include_greenspace:

  Logical. Produce two greenspace outputs: (1) a ground roughness raster
  (z0, metres) from land-cover data for nutURoughWallFunction, and (2)
  distance-to-greenspace as a building-level attribute for
  post-processing / context.

- mask_tree_cover:

  Logical. When building the roughness raster, set tree-cover cells to
  NA so porous-zone drag is not counted twice. Default TRUE.

- landcover_source:

  Character. Land-cover dataset for the ground roughness raster. `"esa"`
  (default) uses ESA WorldCover (years 2020-2021); `"esri"` uses the
  Sentinel-2 10 m ESRI LULC Time Series (years 2017-2025), which
  provides more recent and historically consistent annual maps.

- landcover_year:

  Integer. Year of the land-cover product. For `"esa"`: 2020 or 2021.
  For `"esri"`: 2017-2025. Default 2021.

- include_bgvi:

  Logical. Add building green visibility index.

- include_shadow:

  Logical. Add shadow outputs.

- include_radiation:

  Logical. Add radiation outputs.

- opentopo_key:

  Character. OpenTopography API key for get_fused_dsm().

- target_crs:

  Optional projected target CRS. Usually leave as `NULL`: when the
  building data are in longitude/latitude, the local UTM zone is derived
  automatically from the data centroid, whichever of `bbox`, `place`, or
  `buildings_list` was supplied. Set this only to force a specific
  projection (e.g. a national grid). OpenFOAM domains must be metric, so
  a geographic CRS is rejected.

- height_col:

  Optional height column name. If NULL, guessed automatically.

- default_height:

  Numeric. Default building height if height column is missing (metres).

- min_height:

  Numeric. Minimum building height (metres).

- domain_buffer:

  Numeric. Buffer in metres around data extent.

- zmax_buffer:

  Numeric. Buffer in metres above maximum surface height.

- stl_name:

  Character. STL file name.

- overwrite:

  Logical. Whether to overwrite prepared files.

- quiet:

  Logical. Suppress messages.

## Value

A list with:

- case_dir:

  Absolute path to the OpenFOAM case directory.

- data_dir:

  Absolute path to the gloBFPr data sub-directory.

- files:

  Named list of output file paths (NULL when not produced).

- data:

  Named list of in-memory spatial objects.

- domain:

  Named list with xmin/xmax/ymin/ymax/zmin/zmax for blockMeshDict.

- origin:

  Named numeric vector (x, y, z) of the local coordinate system origin
  in the input CRS.

- crs:

  CRS of the building data.

- height_col:

  Name of the height column used.

- n_buildings:

  Number of building polygons.

- max_building_height:

  Maximum building height in metres.

- max_surface_height:

  Maximum surface height (buildings + DSM).

- include_tree_canopy:

  Logical flag as supplied.

- canopy_source:

  Canopy data source as supplied.

- min_tree_height:

  Minimum tree height as supplied.

- landcover_source:

  Land-cover dataset used (`"esa"` or `"esri"`).

- landcover_year:

  Land-cover year used.
