# Prepare a transient OpenFOAM case for urban wind

Generates a complete OpenFOAM wind case using
`buoyantBoussinesqPimpleFoam` - a transient, Boussinesq, finite-volume
Navier-Stokes solve with a kOmegaSST (URANS) closure - driven by an ABL
log-law inlet.

The solver is a buoyant one running an isothermal problem: every surface
sits at `T_ref` and the walls are adiabatic, so the Boussinesq body
force is identically zero and the result is pure mechanical flow.

Terrain and canopy are used when supplied and skipped when not. Terrain
becomes solid geometry; canopy becomes distributed drag via
`atmPlantCanopy*` rather than solid blocks, because a crown passes and
drags air instead of blocking it.

## Usage

``` r
prepare_foam_case(
  case_dir,
  stl_file,
  domain,
  inlet_velocity = c(5, 0, 0),
  z_ref = 10,
  T_ref = 295,
  sim_hours = NULL,
  n_writes = 8L,
  terrain_stl = NULL,
  terrain_dem = NULL,
  canopy_stl = NULL,
  leaf_area_density = 0.4,
  canopy_heat_source = NULL,
  plant_cd = 0.2,
  z0 = 0.1,
  base_cell_size = 10,
  building_refinement = 2L,
  terrain_refinement = 0L,
  max_cells = 3000000L,
  max_co = NULL,
  n_outer_correctors = 2L,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- case_dir:

  Character. Case directory; must already contain
  `constant/triSurface/<stl_file>`.

- stl_file:

  Character. Building STL path (host path).

- domain:

  Named list with xmin/xmax/ymin/ymax/zmin/zmax (metres, local
  coordinates), from `prepare_openfoam_inputs()$domain`.

- inlet_velocity:

  Numeric length-3 (Ux, Uy, Uz) in m/s at `z_ref`. Any horizontal
  direction is supported.

- z_ref:

  Numeric. Reference height for `inlet_velocity`. Default 10.

- T_ref:

  Numeric. Reference air temperature (K). Default 295. Every surface is
  held at this value, which is what switches buoyancy off.

- sim_hours:

  Numeric. Physical hours to simulate. `NULL` (the default) computes
  three flow-through times from the domain length and the inlet speed,
  which is what a wind case needs to flush its transient.

- n_writes:

  Integer. Output times over the run. Default 8.

- terrain_stl:

  Path to a terrain STL (see
  [`prepare_foam_geometry`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_geometry.md)),
  or `NULL` for a flat floor.

- terrain_dem:

  Optional SpatRaster / path of the bare-earth DEM. Used to set the
  domain floor. Strongly recommended whenever `terrain_stl` is supplied.

- canopy_stl:

  Path to a canopy volume STL, or `NULL`.

- leaf_area_density:

  Numeric. LAD (1/m) inside the canopy. Default 0.4.

- canopy_heat_source:

  Numeric or `NULL`. Enables `atmPlantCanopyTSource` and writes
  `0/qPlant` with this uniform value. `NULL` (default) leaves the source
  out entirely.

  Off by default deliberately: `qPlant` is a canopy energy flux that
  cannot be derived from a canopy height model, and enabling the source
  with `qPlant = 0` would be inert while still adding a way for the run
  to abort. Check the units against the `atmPlantCanopyTSource`
  documentation for your OpenFOAM version before relying on a value.

- plant_cd:

  Numeric. Canopy drag coefficient. Default 0.2.

- z0:

  Numeric. Aerodynamic roughness length (m). Default 0.1.

- base_cell_size:

  Numeric. Background cell size (m). Default 10.

- building_refinement:

  Integer. snappyHexMesh level for buildings. Default 2.

- terrain_refinement:

  Integer. snappyHexMesh level for terrain. Default 0 (the background
  cell size already resolves gentle slope).

- max_cells:

  Integer. Global cell-count cap. Default 3e6.

- max_co:

  Numeric. Maximum Courant number. `NULL` (the default) gives 20: PIMPLE
  re-converges momentum and pressure within each step, so Courant well
  above 1 is stable when marching to a steady field.

- n_outer_correctors:

  Integer. PIMPLE outer correctors. Default 2.

- overwrite, quiet:

  Logical.

## Value

Invisibly, a list with `case_dir`, `files` and `params`.

## Oblique wind

`blockMesh` emits four separately named lateral patches (`xMin`, `xMax`,
`yMin`, `yMax`) and each is assigned an inlet / outlet / lateral role
from the sign of `dot(flowDir, outward_normal)`. Any wind direction
works without rotating the domain - unlike the previous generator, which
pinned the inlet to the x-min face so that a north wind injected
velocity through a face whose normal was -x and effectively nothing
entered.

## See also

[`prepare_openfoam_inputs`](https://billbillbilly.github.io/gloBFPr/reference/prepare_openfoam_inputs.md),
[`prepare_foam_geometry`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_geometry.md),
[`read_foam_pedestrian_slice`](https://billbillbilly.github.io/gloBFPr/reference/read_foam_pedestrian_slice.md)
