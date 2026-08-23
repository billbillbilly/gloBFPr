# gloBFPr 2.0.0

## New features

* Added `get_3d_world()` to export a study area (terrain, extruded buildings,
  and canopy trees) as Wavefront OBJ and/or binary STL files for direct use
  in Rhino3D and Blender.
* Added `get_era5_met()` to fetch ERA5 reanalysis meteorological conditions
  (10 m wind, 2 m and skin temperature) from the Copernicus Climate Data
  Store, pre-formatted as inputs to `prepare_openfoam_case()` and
  `prepare_nocturnal_case()`.
* Added an integrated urban noise mapping workflow: `get_noise_map()` and
  `prepare_noisemodelling_inputs()` prepare screening-level road-noise inputs
  from building height, OSM roads, canopy height, and greenspace data, and
  can run the official headless NoiseModelling WPS scripts directly
  (`run = TRUE`). `install_noisemodelling()` installs the headless runner.
* Added `infer_osm_traffic()` and `osm_noise_traffic_defaults()` to derive
  screening-level traffic speed/volume defaults from OSM road class when
  observed counts are unavailable.
* Added pedestrian-level wind and nocturnal thermal comfort simulation
  support built on OpenFOAM, run inside Docker via `run_openfoam_docker()`:
  `prepare_openfoam_case()`, `prepare_nocturnal_case()`,
  `prepare_openfoam_inputs()`, and post-processing/plotting helpers
  `sample_foam_slice()`, `read_foam_pedestrian_slice()`, `plot_foam_map()`,
  and `add_flow_vectors()`.
* Added `generate_block()` and `aggregate_block()` to aggregate
  individual-building metrics into block-level summaries for city-scale
  analysis. `aggregate_block(population = TRUE)` fetches GHSL population
  directly at block level (`pop_total`), and
  `aggregate_block(residential = TRUE)` computes a block-level residential
  built-up surface proportion (`res_prop`) directly from GHS rasters, rather
  than aggregating per-building estimates.
* Added `get_neighbors()` for neighboring-building counts and centroid
  distance summaries within a fixed-radius buffer.
* Added `get_bgvi()` (Building Green View Index, following Qi et al. 2024)
  and `get_dng()` (distance to nearest green space), built on the
  `viewscape` package (now requires `viewscape >= 2.0.1`).
* Added `get_fused_dsm()` and `get_metadata()` for retrieving and fusing
  global building height tiles and dataset metadata.

## Bug fixes

* `get_fused_dsm()` no longer silently degrades the output DSM to the DEM's
  native resolution regardless of a finer canopy height model, and no longer
  resamples continuous elevation/canopy surfaces with nearest-neighbor
  (blocky terrain). A new `resolution` argument allows an explicit override.
  Building roofs are now flattened to a single elevation per building (base
  ground elevation at the centroid, plus height) instead of following the
  terrain slope pixel by pixel.
* `get_bgvi()`'s internal DSM (used for the viewshed computation) received
  the same resolution and flat-roof fixes as `get_fused_dsm()`, via a new
  `resolution` argument.
* `prepare_openfoam_inputs(include_fused_dsm = TRUE)` called
  `get_fused_dsm(key = opentopo_key, ...)`, but `get_fused_dsm()`'s
  parameter is `opentopo_key`; this raised an "unused argument" error at
  runtime and is now fixed.


## Documentation

* Added a pkgdown site (<https://billbillbilly.github.io/gloBFPr/>) with a
  reference index grouped by workflow (data acquisition, 3D scene export,
  metrics, shadows/radiation, block analysis, noise mapping, OpenFOAM) and
  articles for each major workflow.
* Documented external software setup (Docker/OpenFOAM, Java/NoiseModelling)
  in the README.
