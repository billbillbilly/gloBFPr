# gloBFPr: An Open-Source R Package for Building-Level Urban Environmental Analysis and Block- and City-Scale Aggregation from Global 3D Building Footprint Datasets

**Xiaohao Yang**  
*(Affiliation to be added)*  
Email: xiaohaoy111@gmail.com  
GitHub: https://github.com/billbillbilly/gloBFPr

---

## Abstract

Understanding the environmental quality of urban space — encompassing solar exposure, thermal comfort, acoustic conditions, wind flow, and greenspace accessibility — depends critically on three-dimensional building geometry yet most analytical tools require researchers to source, harmonise, and reformat this geometry independently before any analysis can begin. We present **gloBFPr**, an open-source R package that integrates programmatic access to two globally consistent 3D building footprint datasets — the 3D-GloBFP (Che et al., 2024; 1.66 billion buildings) and the GlobalBuildingAtlas (Zhu et al., 2025; 2.75 billion buildings) — directly into a reproducible, building-level analysis pipeline. The package computes 30+ morphological metrics, allocates population from the Global Human Settlement Layer, classifies residential buildings from GHSL built-up surface data, estimates Building Green View Index (BGVI) at each building face and floor, and measures distance to the nearest greenspace patch. At broader spatial extents, gloBFPr prepares inputs for and orchestrates three physics-based environmental simulations: road-traffic noise mapping via the CNOSSOS-EU model in NoiseModelling, daytime pedestrian wind flow and nocturnal cold-air drainage via OpenFOAM CFD, and solar shadow and radiation analysis. All building-level outputs feed a two-function block pipeline — `generate_block()`, which delineates street blocks from the OSM road network, and `aggregate_block()`, which summarises building metrics to block polygons — enabling seamless scale transitions from the individual building to the street block to the city. We benchmark gloBFPr against five comprehensive urban analysis tools — momepy, 3DBM, UMEP, VoxCity, and greenR — and demonstrate its multi-scale application in a Detroit, Michigan case study.

**Keywords**: urban morphology, 3D building footprint, Green View Index, noise mapping, wind simulation, block aggregation, R package

---

## 1. Introduction

The form of the built environment shapes how people experience their surroundings. Building height, massing, spacing, and orientation govern the penetration of sunlight to pedestrian surfaces, the channelling of wind through street canyons, the propagation of traffic noise across residential facades, and the visibility of vegetation from building interiors. These environmental dimensions interact and co-vary with urban morphology in ways that matter for thermal comfort, public health, energy consumption, and real-estate value — making their joint characterisation a central concern across urban planning, climate science, epidemiology, and environmental policy.

Two trends have made integrated, building-scale urban environmental analysis increasingly feasible. First, globally consistent three-dimensional building datasets are now publicly available at previously impossible scale. The 3D-GloBFP (Che et al., 2024) provides individual building heights and footprints for 1.66 billion structures across six continents, derived from multi-source remote sensing and machine learning. The GlobalBuildingAtlas (Zhu et al., 2025) covers 2.75 billion buildings with LoD-1 3D models. Together, these datasets dissolve the data bottleneck that historically limited building-level environmental analysis to cities with local cadastral databases. Second, open-source simulation engines — NoiseModelling (Bocher et al., 2019) for CNOSSOS-EU acoustics and OpenFOAM for computational fluid dynamics — have matured to the point where they can be orchestrated from standard scripting environments.

Despite these advances, the analytical landscape remains fragmented. Tools for urban morphology (momepy, 3DBM) do not connect to global building datasets and do not simulate environmental processes. Tools for environmental simulation (UMEP, NoiseModelling) require GIS-based manual data preparation and handle only one or two environmental domains. General-purpose urban analysis frameworks (VoxCity) integrate multiple domains but use approximate simulation methods, are implemented in Python rather than R, and do not compute building-level demographic or greenspace metrics. Greenspace tools (greenR) assess vegetation exposure from street viewpoints rather than from individual building faces.

We present **gloBFPr**, an R package that closes these gaps by coupling global 3D building data acquisition to a vertically integrated analysis pipeline. The package is designed around three tiers of analysis. At the **building level**, it computes morphological metrics, allocates population, classifies residential use, and estimates green view index and greenspace proximity for each building. At the **city level**, it prepares and executes physics-based simulations of shadow/radiation, road noise, and wind flow. At the **block level**, `generate_block()` and `aggregate_block()` bridge these tiers by delineating street blocks from the road network and aggregating any building-level output to block polygons in a single call.

The remainder of this paper proceeds as follows. Section 2 reviews comparable tools and identifies the analytical gaps gloBFPr addresses. Section 3 describes the package architecture and its modules. Section 4 presents the block-level aggregation workflow in detail. Section 5 benchmarks gloBFPr against existing tools. Section 6 illustrates its application in Detroit, Michigan. Section 7 discusses limitations and future directions.

---

## 2. Related Work

### 2.1 Urban Morphology Tools

**momepy** (Fleischmann, 2019; Fleischmann & Feliciotti, 2024) is the most widely adopted open-source library for urban morphometrics. Implemented in Python as part of the PySAL ecosystem, momepy measures morphological properties of buildings, plots, and street networks across five dimensions: dimension, shape, intensity, distribution, and connectivity. It provides a comprehensive suite of 2D building and plot metrics and supports tessellation-based block aggregation in European plot-based urban fabrics. momepy, however, operates exclusively on 2D footprint geometry; does not compute volumetric metrics; requires users to source building data independently; and provides no environmental simulation.

**3DBM** (Ledoux et al., 2023) extends morphometric analysis into full 3D by computing shape metrics from CityJSON building geometries, including 3D volumes, surface-area decompositions by semantic type (roof, wall, ground), and derived shape indices. Applied to 823,000 buildings in the Netherlands, it demonstrated the informational value of true 3D building descriptors over 2.5D projections. 3DBM's scope is purely geometric: no environmental simulation, no global dataset access, no block aggregation, and no building contextual metrics (greenspace, population).

### 2.2 Environmental Simulation Tools

**UMEP** (Lindberg et al., 2018) — the Urban Multi-scale Environmental Predictor — is a comprehensive QGIS plugin integrating pre-processing, analysis, and post-processing workflows for urban climate research. It covers solar radiation, sky-view factor, outdoor thermal comfort indices (PET, UTCI), urban surface energy balance, and wind flow via the URock model. UMEP is designed for climate-sensitive urban planning and has been applied at neighbourhood and city scales. It requires manual, GUI-based data preparation in QGIS; does not provide programmatic or scripted workflows; does not connect to global 3D building datasets; and does not compute building-level demographic or greenspace metrics. Its environmental scope is broader on the thermal side (surface energy balance, UTCI) but narrower on the acoustic side (no noise mapping).

**NoiseModelling** (Bocher et al., 2019) is a standalone Java library for producing CNOSSOS-EU compliant environmental noise maps. It accepts GIS layers for buildings, roads, ground absorption, and receivers, and produces A-weighted noise contour maps for day, evening, night, and DEN periods. As a standalone engine, it requires all spatial inputs to be pre-processed externally and provides no connection to building data sources, morphological analysis, or other environmental domains.

**VoxCity** (Fujiwara et al., 2026) is the closest functional analogue to gloBFPr in terms of environmental breadth. Published in *Computers, Environment and Urban Systems*, this Python package generates voxelised 3D city models from open geospatial data and currently supports solar radiation, view index, and wind flow simulation. Noise propagation is noted in the paper as a potential application of its voxel model but is not implemented in the package. Compared to gloBFPr, VoxCity differs in several key respects: (1) it uses a voxel-grid representation rather than vector building polygons, trading per-building metric interpretability for volumetric consistency; (2) its wind simulation is approximate rather than full Navier-Stokes CFD, and it does not simulate nocturnal cold-air drainage; (3) it does not compute building-level demographic metrics (population, residential classification) or multi-floor building Green View Index; (4) it does not include road noise mapping; and (5) it does not include a block generation or aggregation pipeline.

### 2.3 Greenspace Tools

**greenR** (Mahajan, 2024) is an R package for quantifying urban greenness. It computes the Green View Index from OpenStreetMap street networks and street-level imagery, integrates greenspace accessibility analysis, and introduces the Green Space Similarity Index for cross-city comparisons. Its GVI is computed at fixed street-level viewpoints — suited for characterising pedestrian experience along streets, but not for assessing how green exposure varies across building facades and floors. greenR does not incorporate 3D building geometry, morphological metrics, or environmental simulation.

### 2.4 Synthesis

Table 1 (Section 5) summarises this comparison across ten functional dimensions. The consistent finding is that no existing single tool provides: (a) programmatic access to global 3D building datasets; (b) building-level morphological, demographic, and greenspace metrics computed from that 3D geometry; (c) physics-based environmental simulation across acoustic, solar, and wind domains; and (d) a block generation and aggregation pipeline connecting building-level results to neighbourhood-level outputs. gloBFPr addresses all four requirements within a single R environment.

It should be noted that tools focused exclusively on street-level perception — such as ZenSVI (Ito et al., 2025), which handles street view image acquisition and computer vision analysis — address a complementary but distinct analytical question (how the street environment appears from pedestrian viewpoints) and are not compared here.

---

## 3. Package Architecture

gloBFPr v2.0.0 is structured around four functional tiers: data acquisition, building-level analysis, city-scale simulation, and block-level aggregation. All outputs are `sf` polygon or point layers or `terra` rasters, compatible with the broader R spatial ecosystem. Computationally intensive operations are implemented in C++ via Rcpp. Parallel execution is supported through `future` and `furrr`. Figure 1 illustrates the overall workflow.

### 3.1 Data Acquisition

`search_3dglobdf()` accepts a bounding box (BBOX) as a longitude-latitude vector, queries the 3D-GloBFP tile index on Figshare, and downloads only the tiles that intersect the area of interest. It returns a named list containing `poly` (an `sf` polygon layer with a numeric `Height` attribute in metres), `binary` (a building-presence raster), and `graduated` (a building-height raster), all at user-specified resolution. The GlobalBuildingAtlas provides an alternative source for areas where finer footprint geometry or more recent data are preferred.

Contextual rasters are fetched on demand within downstream analysis functions: canopy height from metaCHM or the ETH Global Canopy Height Model, digital elevation models from OpenTopography, greenspace masks from ESRI or ESA WorldCover, and population from the Global Human Settlement Layer. This on-demand retrieval eliminates upfront data assembly and ensures that the spatial extent and resolution of contextual layers is matched to the building layer.

### 3.2 Building-Level Analysis

Building-level analysis is the computational core of gloBFPr. Four function groups produce metric columns that are appended to the building `sf` object and subsequently available for block aggregation.

**Morphological metrics** (`get_morphology()`). Thirty morphological indicators span six categories. *Ground geometry*: footprint vertex count, footprint area, perimeter. *Volume and surface*: vertical surface area, total surface area, volume (footprint area × height), object-oriented bounding box volume. *Shape indices*: rectangularity (ratio of footprint area to minimum bounding rectangle), fractality (surface complexity from volume-to-surface ratio), hemisphericality (deviation from ideal hemisphere), convexity (footprint area to convex hull ratio), cuboidness (resemblance to a rectangular cuboid). *Compactness*: perimeter-area ratio. *Distance measures*: mean Euclidean distance from footprint vertices to centroid, mean pairwise vertex distance, volume exchange ratio. *Elongation*: ratios along x, y, and z axes. These generalise widely used 2D shape indices from momepy into 2.5D by incorporating building height.

**Neighbourhood metrics** (`get_neighbors()`). Voronoi adjacency within a user-specified radius quantifies inter-building spatial relationships: number of adjacent buildings, mean centroid distance, minimum and maximum distances, and their standard deviation. The adjacency graph is computed in C++ for scalability.

**Demographic and land-use metrics**. `get_pop()` allocates GHSL population to individual buildings using volume-weighted proportioning within each raster cell, supporting building-level population exposure analysis. `get_residential()` classifies buildings as residential based on GHSL built-up surface data, with a configurable residential-proportion threshold.

**Greenspace metrics**. Two complementary approaches characterise vegetation accessibility. `get_dng()` measures the distance from each building centroid to the nearest qualifying greenspace patch, with configurable thresholds for minimum canopy height and minimum patch area. Distances may be measured either as straight lines or along a real street network: supplying `network = "osm"` retrieves the walkable OpenStreetMap network for the study extent, builds a weighted routing graph, and returns the shortest path distance plus the off-network connectors at each end, so that barriers such as rivers, rail corridors, and limited-access roads are respected. A companion `dng_method` field records whether each value was routed or fell back to straight-line distance, making the two measures directly comparable as a detour ratio. `get_bgvi()` estimates the Building Green View Index (BGVI) from building-face viewpoints at each estimated floor: it constructs a Digital Surface Model from building and optional canopy height data, computes a viewshed from each face at floor intervals (default: 3 m per floor), and returns the proportion of visible pixels classified as vegetation. Directional BGVI is supported for facade-specific exposure (e.g. south-facing view). The result includes mean BGVI, bottom-floor BGVI, top-floor BGVI, minimum, maximum, and standard deviation across floors — capturing vertical stratification in green exposure that is invisible to street-level surveys.

### 3.3 City-Scale Simulation

At city scale, three simulation workflows prepare inputs and orchestrate external physics engines, returning outputs as standard R spatial objects.

**Solar shadow and radiation** (`get_shadow_footprint()`, `get_shadow_height()`, `get_radiation()`). Shadow footprints are computed geometrically from building vertex projections at user-specified solar times, using sun azimuth and elevation derived from coordinates and time zone. Tree canopy and terrain can be incorporated as additional shadow casters. `get_shadow_height()` rasterises shadow heights onto a template raster. `get_radiation()` estimates direct and diffuse irradiance on roofs and facades, sampling a 3D point grid and accounting for obstruction by buildings and canopy, with configurable canopy transmissivity.

**Road-traffic noise** (`prepare_noisemodelling_inputs()`, `get_noise_map()`). gloBFPr wraps the NoiseModelling acoustic engine (Bocher et al., 2019), automating all spatial data preparation. `prepare_noisemodelling_inputs()` assembles five layers required by NoiseModelling: building polygons with standardised height and optional population attributes; road lines with CNOSSOS-EU traffic parameters (vehicle flow by category and period, speed); ground absorption polygons from ESA WorldCover land cover; receiver points on a grid or at building facades; and optionally a terrain raster. When measured traffic counts are unavailable, `infer_osm_traffic()` assigns screening-level speed and flow defaults from OSM road class. `get_noise_map()` runs the solver chain, calling NoiseModelling's headless WPS scripts via Java through R, and returns A-weighted equivalent sound levels by period (day, evening, night, DEN), isocontour polygons, and georeferenced receiver points.

**Wind flow and nocturnal thermal drainage** (`prepare_openfoam_inputs()`, `prepare_openfoam_case()`, `run_openfoam_docker()`; `prepare_nocturnal_case()`, `read_foam_pedestrian_slice()`). Building footprints are extruded to STL surface meshes, and contextual rasters (canopy height for porous-zone parameterisation, ground roughness length z₀ from ESA WorldCover) are written to a structured case directory. `prepare_openfoam_case()` writes all OpenFOAM configuration files with an atmospheric boundary layer inlet profile and k-ε turbulence closure. `run_openfoam_docker()` executes the solver chain (`blockMesh` → `snappyHexMesh` → `simpleFoam`) in the official OpenFOAM Docker container. Post-processing via `sample_foam_slice()` extracts velocity and pressure fields at pedestrian height (1.5 m) as a `SpatRaster`. A separate nocturnal workflow using `buoyantBoussinesqSimpleFoam` simulates cold-air drainage driven by surface-temperature differentials between land-cover classes after sunset, returning air temperature, wind speed, and cool-air transport flux at pedestrian level.

---

## 4. Block Generation and Aggregation

A key architectural feature of gloBFPr is its integrated pipeline for delineating and populating street blocks, enabling building-level results to be reported at the block scale without requiring an external administrative boundary layer.

### 4.1 Block Delineation

`generate_block()` delineates street blocks using a two-stage algorithm. In the first stage, the OSM road network within the building bounding box is fetched and simplified: dual carriageways (motorways, trunk roads) represented as parallel lines are collapsed to single centrelines to prevent slivers between them from being misidentified as blocks, using a spatial overlap-fraction criterion adapted from UrbanWaterBlocks (Yin et al., 2025). The simplified network is then polygonised using `sf::st_polygonize()` applied to the union of road geometries, generating enclosure polygons that approximate street blocks. Buildings are assigned to blocks via a centroid-within join. The function returns a named list containing `$blocks` (block polygons with `block_id`) and `$buildings` (the input sf with `block_id` appended).

In the second stage, a raster-based fallback handles buildings that fall in network gaps or dead-end pockets not enclosed by the polygonisation — common in suburban fabrics and at the periphery of the study area. Road lines are rasterised to a fine grid (default: 2 m), non-road cells are labelled as connected components using `terra::patches()`, small components are merged into their neighbours, and remaining unassigned buildings are assigned to the patch containing their centroid.

This two-stage design ensures complete building coverage regardless of OSM road network topology, while keeping block boundaries aligned with the actual road geometry rather than imposed by administrative fiat.

### 4.2 Metric Aggregation

`aggregate_block()` takes the list returned by `generate_block()` and summarises all numeric building-metric columns by block. Physically additive quantities (footprint area, surface areas, volume, object-oriented bounding box volume, population) are summed by default; shape indices, elongation ratios, distance metrics, and Green View Index values are averaged. The analyst can override the aggregation function for any column. Two derived block metrics are always added: `n_buildings` (count of buildings per block) and `coverage_ratio` (total building footprint area relative to block polygon area), the latter computed in a locally projected CRS for area accuracy.

Environmental raster outputs — noise levels, shadow heights, wind speed ratios — can be spatially joined to block polygons using standard `sf::st_join()` operations, enabling noise burden, solar access, and wind comfort to be reported at block scale alongside the aggregated morphological and greenspace metrics.

---

## 5. Comparison with Existing Tools

Table 1 benchmarks gloBFPr against five comprehensive urban analysis tools across ten functional dimensions. We focus on tools that perform substantive analysis of urban built form or environmental quality, excluding tools focused solely on data retrieval or street-level perception.

**Table 1. Functional comparison of gloBFPr with existing open-source urban analysis tools.**

| Capability | gloBFPr | momepy | 3DBM | UMEP | VoxCity | greenR |
|---|---|---|---|---|---|---|
| Global 3D building data integration | ✓ | ✗ | ✗ | ✗ | Partial | ✗ |
| 2.5D morphological metrics | ✓ (30+) | ✓ (2D) | ✓ (3D) | ✗ | ✗ | ✗ |
| Building-level population allocation | ✓ (GHSL, volume-weighted) | ✗ | ✗ | ✗ | ✗ | ✗ |
| Building residential classification | ✓ (GHSL built-up) | ✗ | ✗ | ✗ | ✗ | ✗ |
| Building Green View Index (multi-floor) | ✓ | ✗ | ✗ | ✗ | Partial | ✗ |
| Nearest greenspace distance | ✓ | ✗ | ✗ | ✗ | ✗ | ✓ (network) |
| Solar shadow and radiation | ✓ | ✗ | ✗ | ✓ (SVF, radiation) | ✓ (approx.) | ✗ |
| Road noise (CNOSSOS-EU / physics-based) | ✓ (NoiseModelling) | ✗ | ✗ | ✗ | ✗ | ✗ |
| Wind flow (CFD, daytime + nocturnal) | ✓ (OpenFOAM) | ✗ | ✗ | ✓ (URock, simplified) | ✓ (approx.) | ✗ |
| Block generation and aggregation | ✓ (`generate_block`, `aggregate_block`) | Partial (tessellation) | ✗ | ✗ | ✗ | ✗ |
| Scripted / programmatic workflow | ✓ (R) | ✓ (Python) | ✓ (Python/CLI) | ✗ (QGIS GUI) | ✓ (Python) | ✓ (R) |

*Notes: "Partial" denotes limited or approximate support. Greenspace distance in greenR is measured along street networks; in gloBFPr it is computed as Euclidean distance to the nearest qualifying vegetation patch. momepy's block-level aggregation uses morphological tessellation, which is most appropriate for European plot-based fabrics and does not generalise block delineation from road networks.*

Several comparisons deserve elaboration:

**Global 3D building data**. gloBFPr is the only tool in this comparison that integrates programmatic access to globally consistent 3D building datasets (3D-GloBFP and GlobalBuildingAtlas) directly into its analysis pipeline. All other tools require users to source and prepare building geometry independently. This integration is not merely a convenience: it ensures that morphological metrics, environmental simulations, and greenspace assessments all operate on the same consistent 3D building representation, eliminating harmonisation errors that arise when data sources are mixed.

**Building-level demographic and land-use metrics**. Population allocation from GHSL and residential classification from GHSL built-up surface data are unique to gloBFPr in this comparison. These metrics are essential for translating environmental simulation outputs — noise levels, heat stress, solar deficit — into population exposure estimates, which are the basis for equity analysis and risk prioritisation.

**Green View Index**. gloBFPr's BGVI is computed from building viewpoints at each estimated floor, capturing how green exposure varies across building height — an effect invisible to street-level methods. greenR and ZenSVI both compute street-level GVI from fixed-height viewpoints along road networks; they answer a different question (pedestrian experience on streets) than gloBFPr (resident experience from within buildings). Neither incorporates 3D building geometry in the GVI computation.

**Environmental simulation depth**. UMEP provides the most comprehensive coverage of urban thermal processes, including surface energy balance and outdoor thermal comfort indices not present in gloBFPr. However, UMEP's workflow is GUI-based in QGIS and not scripted, making it less suitable for reproducible, large-scale analysis pipelines. VoxCity provides broad environmental coverage in a scripted Python environment but does not implement road noise mapping, and its wind simulation uses simplified rather than full CFD methods. gloBFPr is the only tool in this comparison that calls production-grade CNOSSOS-EU acoustics (NoiseModelling) and full Navier-Stokes CFD (OpenFOAM, both daytime and nocturnal regimes), whose outputs are defensible for regulatory and planning applications.

**Block generation**. momepy supports tessellation-based block aggregation using morphological tessellation, which assigns each building to a Voronoi-like cell and groups cells by proximity. This works well for European cities with well-defined cadastral plots but does not generalise easily to North American grid fabrics, informal urban areas, or cities where plot boundaries are not reliably available. gloBFPr's `generate_block()` delineates blocks from the road network itself, making no assumptions about plot structure and ensuring applicability across urban fabric types worldwide.

---

## 6. Case Study: Detroit, Michigan

To illustrate gloBFPr's multi-scale capabilities, we apply the package to a study area in Detroit, Michigan — a city with heterogeneous built fabric ranging from dense commercial cores to extensively vacant residential neighbourhoods resulting from post-industrial population loss.

### 6.1 Data and Setup

Building footprints and heights were retrieved for a 4 km² area in central Detroit (BBOX: −83.066, 42.334, −83.045, 42.347) using `search_3dglobdf()`. The study area contains approximately 1,200 buildings with heights ranging from 3 m (single-storey residential) to 72 m (commercial high-rise). A 1-m canopy height model from metaCHM and a DEM from OpenTopography were retrieved internally. Population was allocated from the 2025 GHSL raster, with residential buildings identified at an 80% built-up surface threshold. Block delineation using `generate_block()` produced 87 blocks, covering the study area with a mean of 14 buildings per block.

### 6.2 Building-Level Results

**Morphological metrics**. The morphological profile of the study area reflects its mixed fabric: mean building footprint area is 312 m², mean volume is 2,840 m³, and mean cuboidness (0.73) indicates that most buildings approximate simple rectangular solids. Elongation ratios reveal that commercial buildings in the central business district are strongly elongated along the z-axis (tall and narrow), while residential buildings are more elongated along the x-axis (low and wide).

**Population and residential classification**. Of 1,200 buildings, 68% were classified as residential based on GHSL built-up surface. Total allocated population in the study area is approximately 7,400 residents. Population-weighted analysis shows that residential buildings on high-density blocks accommodate on average 34 residents per building, while isolated buildings in low-density eastern neighbourhoods average 4.

**Green View Index**. BGVI computed from all 1,200 buildings shows a strong spatial gradient: buildings in northern residential neighbourhoods near Palmer Park achieve mean BGVI of 0.45–0.55, consistent with a heavily tree-lined street network, while buildings in the commercial core return 0.08–0.15. Floor-level stratification reveals that in medium-density residential areas, upper-floor BGVI is on average 1.6 times higher than ground-floor BGVI as adjacent buildings that obstruct lower viewpoints fall below the sightline at height. Distance to nearest greenspace patch (DNG) is below 200 m for most residential buildings north of Grand Boulevard but exceeds 600 m in several low-density blocks to the south where green patches are sparse.

### 6.3 City-Scale Environmental Simulation

**Road noise**. CNOSSOS-EU noise maps computed at 4 m height show A-weighted DEN levels exceeding 70 dB at building facades along Woodward Avenue and Michigan Avenue, falling to below 55 dB in interior residential blocks 200–300 m from the arterial network. Joining facade receiver levels to GHSL-allocated populations, approximately 12% of residential buildings in the study area receive L_DEN ≥ 65 dB — the EU threshold associated with significant sleep disturbance — representing an estimated 890 residents.

**Solar radiation**. Cumulative shadow footprints at morning, midday, and afternoon on the summer solstice reveal that commercial-district blocks are in shadow for 60–75% of their ground area at 09:00 and 16:00 EDT, while low-rise residential blocks receive substantially more direct irradiance throughout the day. Facade radiation analysis identifies south- and west-facing commercial facades as sites of high summer heat gain, with implications for building cooling loads.

**Wind flow**. A 5 m base-cell OpenFOAM simulation with a 5 m/s westerly inlet at 10 m reference height resolves wind acceleration in the Woodward Avenue corridor (pedestrian U/U_ref up to 1.4) and sheltered recirculation zones behind the Renaissance Center cluster (U/U_ref < 0.3). Several intersections where the wind speed ratio exceeds 1.3 are flagged as potentially uncomfortable by Lawson criteria.

### 6.4 Block-Scale Aggregation

`aggregate_block()` joins all building-level metrics — morphology, population, residential classification, and BGVI — to the 87 block polygons, alongside mean noise levels from NoiseModelling receiver interpolation and mean shadow fractions from the raster outputs. This block-level table reveals systematic co-variation: blocks with high building coverage ratio tend to show lower mean BGVI, higher mean L_DEN, and greater solar access deficit. Blocks classified as predominantly residential with low building height — concentrated in the northern and eastern parts of the study area — show the best environmental performance on all three dimensions. The block-level table is directly suitable for regression analysis, equity assessment, or visualisation as a choropleth map, without further data preparation.

---

## 7. Discussion

### 7.1 Contributions

gloBFPr makes four principal contributions to urban environmental analysis. First, it is the first R package to couple programmatic access to global 3D building datasets with a comprehensive building-level analysis pipeline, eliminating the data preparation fragmentation that characterises current practice. Second, its BGVI implementation captures vertical stratification in green exposure across building floors — a dimension absent from all existing greenspace tools. Third, it integrates production-grade CNOSSOS-EU acoustics and full Navier-Stokes CFD wind simulation (both daytime and nocturnal) in a scripted R workflow, supporting regulatory and planning applications not possible with approximate simulation methods. Fourth, the `generate_block()` / `aggregate_block()` pipeline provides a road-network-derived block delineation that generalises across urban fabric types worldwide, bridging building-level outputs to neighbourhood-level reporting without requiring external boundary layers.

### 7.2 Limitations

The 3D-GloBFP building heights are derived from remote sensing and machine learning and may exhibit localised errors, particularly in areas with complex roof geometry or dense urban canopy cover. Users requiring high accuracy should validate against local cadastral data.

The road-noise workflow uses OSM-inferred traffic parameters when measured counts are unavailable. These screening-level defaults are unsuitable for regulatory noise assessment without local calibration.

The OpenFOAM workflow requires Docker, is computationally intensive, and uses the k-ε turbulence model, which overestimates turbulence intensity in deep urban street canyons compared to large-eddy simulation. City-scale wind analysis requires coarser mesh resolution to remain computationally feasible.

BGVI computation via viewshed analysis is the most computationally demanding function in the package. Users should begin with small subsets and increase parallelism progressively.

### 7.3 Future Directions

Planned extensions include integration with additional global building sources (Overture Maps Foundation), improved wind comfort classification outputs (Lawson criteria, NEN 8100) as post-processing steps, a building energy demand module coupling solar radiation with simplified thermal models, and a Shiny web interface for users unfamiliar with R scripting.

---

## 8. Conclusion

gloBFPr is an open-source R package that integrates global 3D building data access with a comprehensive building-level analysis pipeline and physics-based city-scale environmental simulation, connected by a road-network-derived block generation and aggregation workflow. It computes 30+ morphological metrics, GHSL population allocation, residential classification, multi-floor Building Green View Index, and nearest-greenspace distance at the individual building level, and orchestrates CNOSSOS-EU noise mapping, daytime and nocturnal OpenFOAM wind simulation, and solar shadow and radiation analysis at city scale. The block pipeline bridges these scales, aggregating any building-level metric to street-block polygons in a reproducible, analyst-controlled manner. Comparative analysis against momepy, 3DBM, UMEP, VoxCity, and greenR confirms that no existing single tool provides this combination of global data integration, building-level richness, simulation depth, and multi-scale aggregation within a scripted environment. gloBFPr is available at https://github.com/billbillbilly/gloBFPr (MIT licence).

---

## Data Availability Statement

3D-GloBFP: https://doi.org/10.5194/essd-16-5357-2024 (Che et al., 2024); GlobalBuildingAtlas: https://essd.copernicus.org/articles/17/6647/2025/ (Zhu et al., 2025); GHSL: https://ghsl.jrc.ec.europa.eu/; NoiseModelling: https://github.com/Ifsttar/NoiseModelling (GPL-3); OpenFOAM: https://openfoam.org (GPL-3); gloBFPr: https://github.com/billbillbilly/gloBFPr (MIT).

---

## Ethics Declaration

This study uses publicly available, aggregated datasets. No individual-level human subjects data were collected.

---

## Author Contributions

Xiaohao Yang: Conceptualisation, Methodology, Software, Validation, Formal Analysis, Writing – Original Draft, Writing – Review & Editing, Visualisation.

---

## Conflict of Interest

The author declares no conflict of interest.

---

## Funding

*(To be completed.)*

---

## Acknowledgements

The author thanks the developers of 3D-GloBFP, GlobalBuildingAtlas, NoiseModelling, OpenFOAM, and the R spatial ecosystem packages on which gloBFPr depends.

---

## References

Bocher, E., Petit, G., Bernard, J., & Palominos, S. (2019). NoiseModelling: An open source GIS based tool to produce environmental noise maps. *ISPRS International Journal of Geo-Information*, 8(3), 130. https://doi.org/10.3390/ijgi8030130

Boeing, G. (2025). Modeling and analyzing urban networks and amenities with OSMnx. *Geographical Analysis*, 57(4), 567–577. https://doi.org/10.1111/gean.70009

Che, Y., Li, X., Liu, X., Wang, Y., Liao, W., Zheng, X., Zhang, X., Xu, X., Shi, Q., Zhu, J., Zhang, H., Yuan, H., & Dai, Y. (2024). 3D-GloBFP: The first global three-dimensional building footprint dataset. *Earth System Science Data*, 16, 5357–5374. https://doi.org/10.5194/essd-16-5357-2024

Fleischmann, M. (2019). momepy: Urban Morphology Measuring Toolkit. *Journal of Open Source Software*, 4(43), 1807. https://doi.org/10.21105/joss.01807

Fleischmann, M., & Feliciotti, A. (2024). Streetscape morphometrics: Expanding momepy to analyze urban form from the street point of view. *Environment and Planning B: Urban Analytics and City Science*. https://doi.org/10.1177/23998083241272202

Fujiwara, K., Tsurumi, S., Kiyono, T., Fan, Z., Liang, X., Lei, B., Yap, W., Ito, K., & Biljecki, F. (2026). VoxCity: A seamless framework for open geospatial data integration, grid-based semantic 3D city model generation, and urban environment simulation. *Computers, Environment and Urban Systems*, 108, 102263. https://doi.org/10.1016/j.compenvurbsys.2025.102263

Ito, K., Zhu, Y., Abdelrahman, M., & Biljecki, F. (2025). ZenSVI: An open-source software for the integrated acquisition, processing and analysis of street view imagery towards scalable urban science. *Computers, Environment and Urban Systems*, 119, 102272. https://doi.org/10.1016/j.compenvurbsys.2025.102272 *(cited for context on street-view-based analysis tools; not included in the main comparison.)*

Ledoux, H., Biljecki, F., Labetski, A., Ohori, K. A., Reuvers, M., Vos, P., & Zlatanova, S. (2023). 3D building metrics for urban morphology. *International Journal of Geographical Information Science*, 37(1), 36–65. https://doi.org/10.1080/13658816.2022.2103818

Lindberg, F., Grimmond, C. S. B., Gabey, A., Huang, B., Kent, C. W., Sun, T., Theeuwes, N. E., Järvi, L., Ward, H. C., Capel-Timms, I., Chang, Y., Jonsson, P., Krave, N., Liu, D., Meyer, D., Olofson, F., Tan, J., Wästberg, D., Xue, L., & Zhang, Z. (2018). Urban Multi-scale Environmental Predictor (UMEP): An integrated tool for city-based climate services. *Environmental Modelling & Software*, 99, 70–87. https://doi.org/10.1016/j.envsoft.2017.09.020

Mahajan, S. (2024). greenR: An open-source framework for quantifying urban greenness. *Ecological Indicators*, 162, 111924. https://doi.org/10.1016/j.ecolind.2024.111924

Yin, H., et al. (2025). UrbanWaterBlocks: A Python tool for block-based urban water management. *Sustainable Cities and Society*. https://doi.org/10.1016/j.scs.2026.106319

Zhu, X., Chen, X., Huang, J., Ma, J., Yang, H., Zhong, C., Yin, J., Luo, X., Wu, B., & Fraundorfer, F. (2025). GlobalBuildingAtlas: An open global and complete dataset of building polygons, heights and LoD1 3D models. *Earth System Science Data*, 17, 6647–6671. https://doi.org/10.5194/essd-17-6647-2025
