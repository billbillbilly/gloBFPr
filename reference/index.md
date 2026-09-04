# Package index

## Data acquisition

Search, download, and retrieve global building footprint, elevation,
canopy height, and OpenStreetMap data.

- [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md)
  : search_3dglobdf
- [`get_fused_dsm()`](https://billbillbilly.github.io/gloBFPr/reference/get_fused_dsm.md)
  : get_fused_dsm
- [`get_metadata()`](https://billbillbilly.github.io/gloBFPr/reference/get_metadata.md)
  : get_metadata
- [`install_noisemodelling()`](https://billbillbilly.github.io/gloBFPr/reference/install_noisemodelling.md)
  : Install the headless NoiseModelling runner

## 3D scene export

Export terrain, extruded buildings, and canopy trees as OBJ/STL models
for Rhino3D and Blender.

- [`get_3d_world()`](https://billbillbilly.github.io/gloBFPr/reference/get_3d_world.md)
  : get_3d_world

## Building & environmental metrics

Compute 2D/2.5D building morphology metrics and environmental context
indicators (greenery accessibility, neighbor exposure).

- [`get_morphology()`](https://billbillbilly.github.io/gloBFPr/reference/get_metrics.md)
  [`get_neighbors()`](https://billbillbilly.github.io/gloBFPr/reference/get_metrics.md)
  [`get_bgvi()`](https://billbillbilly.github.io/gloBFPr/reference/get_metrics.md)
  [`get_dng()`](https://billbillbilly.github.io/gloBFPr/reference/get_metrics.md)
  : get_metrics
- [`plot_bgvi_viewshed()`](https://billbillbilly.github.io/gloBFPr/reference/plot_bgvi_viewshed.md)
  : Visualize an Individual Building BGVI Viewshed

## Shadows and solar radiation

Sky view factor, shadow footprints/heights, and surface radiation.

- [`svf()`](https://billbillbilly.github.io/gloBFPr/reference/get_shadows.md)
  [`get_shadow_footprint()`](https://billbillbilly.github.io/gloBFPr/reference/get_shadows.md)
  [`get_shadow_height()`](https://billbillbilly.github.io/gloBFPr/reference/get_shadows.md)
  [`get_radiation()`](https://billbillbilly.github.io/gloBFPr/reference/get_shadows.md)
  : Building shadow and radiation calculations

## Block-level analysis

Aggregate and summarize building metrics at the city-block scale.

- [`generate_block()`](https://billbillbilly.github.io/gloBFPr/reference/generate_block.md)
  : generate_block
- [`aggregate_block()`](https://billbillbilly.github.io/gloBFPr/reference/aggregate_block.md)
  : aggregate_block

## Urban noise mapping

Prepare inputs for and run screening-level road-noise modelling inspired
by NoiseModelling workflows.

- [`osm_noise_traffic_defaults()`](https://billbillbilly.github.io/gloBFPr/reference/osm_noise_traffic_defaults.md)
  : Default OSM traffic assumptions for screening-level road noise
- [`infer_osm_traffic()`](https://billbillbilly.github.io/gloBFPr/reference/infer_osm_traffic.md)
  : Infer screening-level road traffic inputs from OSM road attributes
- [`prepare_noisemodelling_inputs()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_noisemodelling_inputs.md)
  : Prepare NoiseModelling-style input layers
- [`get_noise_map()`](https://billbillbilly.github.io/gloBFPr/reference/get_noise_map.md)
  : Create a screening-level urban road-noise workflow object
- [`plot_noise_map()`](https://billbillbilly.github.io/gloBFPr/reference/plot_noise_map.md)
  : Plot a NoiseModelling-style road-noise map

## Urban wind & thermal CFD (OpenFOAM)

Prepare, run, and post-process pedestrian-level wind and nocturnal
thermal comfort simulations with OpenFOAM.

- [`get_era5_met()`](https://billbillbilly.github.io/gloBFPr/reference/get_era5_met.md)
  : Fetch ERA5 conditions for OpenFOAM boundary conditions
- [`prepare_openfoam_inputs()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_openfoam_inputs.md)
  : Prepare gloBFPr spatial outputs for an OpenFOAM Docker workflow
- [`prepare_foam_geometry()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_geometry.md)
  : Build terrain, canopy and terrain-based building geometry for a case
- [`prepare_foam_case()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_case.md)
  : Prepare a transient OpenFOAM case for urban wind
- [`run_openfoam_docker()`](https://billbillbilly.github.io/gloBFPr/reference/run_openfoam_docker.md)
  : Run an OpenFOAM case via Docker
- [`sample_foam_slice()`](https://billbillbilly.github.io/gloBFPr/reference/sample_foam_slice.md)
  : Sample OpenFOAM results at pedestrian level and return a raster
- [`read_foam_pedestrian_slice()`](https://billbillbilly.github.io/gloBFPr/reference/read_foam_pedestrian_slice.md)
  : Read the pedestrian-level slice and compute wind maps
- [`plot_foam_map()`](https://billbillbilly.github.io/gloBFPr/reference/plot_foam_map.md)
  : Plot an OpenFOAM pedestrian-level map
- [`add_flow_vectors()`](https://billbillbilly.github.io/gloBFPr/reference/add_flow_vectors.md)
  : Add wind / flow vector arrows to a foam map plot

## Example data

- [`globfp_example`](https://billbillbilly.github.io/gloBFPr/reference/globfp_example.md)
  : Test 3D-GloBFP dataset
- [`globfp_example_canopy_height`](https://billbillbilly.github.io/gloBFPr/reference/globfp_example_canopy_height.md)
  : Example canopy height raster for the 3D-GloBFP sample
- [`globfp_example_dem`](https://billbillbilly.github.io/gloBFPr/reference/globfp_example_dem.md)
  : Example DEM for the 3D-GloBFP sample
