# Run an OpenFOAM case via Docker

Mounts the case directory into an OpenFOAM Docker container and executes
the `Allrun` script produced by
[`prepare_foam_case()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_case.md).
The generated wind case chains `blockMesh` -\> `snappyHexMesh` -\>
`buoyantBoussinesqPimpleFoam` and writes log files in `case_dir`.

Docker Desktop must be running and the `case_dir` path must be under a
directory shared with Docker (Docker Desktop -\> Settings -\> Resources
-\> File Sharing).

## Usage

``` r
run_openfoam_docker(
  case_dir,
  image = "opencfd/openfoam-run:2506",
  ncpus = foam_default_ncpus(),
  wait = TRUE,
  quiet = FALSE
)
```

## Arguments

- case_dir:

  Character. Absolute path to the OpenFOAM case directory. Must contain
  an `Allrun` script (written by
  [`prepare_foam_case()`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_case.md)).

- image:

  Character. Docker image tag. Default `"opencfd/openfoam-run:2506"`.
  Use the tag you pulled, e.g. `"opencfd/openfoam-run:2406"` for an
  older version.

- ncpus:

  Integer. Number of CPU cores to use. Defaults to all physical cores
  detected on the host via `parallel::detectCores(logical = FALSE)`.
  When `ncpus > 1` the solver runs in MPI parallel: `decomposePar`
  splits the mesh, `mpirun -np ncpus solver -parallel` runs it, and
  `reconstructPar` reassembles. Typical speed-up is 3-6x on a 4-core
  machine. Pass `ncpus = 1` to force single-core.

- wait:

  Logical. If TRUE (default), R blocks until the simulation finishes. If
  FALSE, the container is launched in the background and the function
  returns immediately.

- quiet:

  Logical. Suppress messages. Default FALSE.

## Value

Invisibly returns a list with `case_dir`, `image`, and the exit `status`
(0 = success; only meaningful when `wait = TRUE`).

## See also

[`prepare_foam_case`](https://billbillbilly.github.io/gloBFPr/reference/prepare_foam_case.md)
