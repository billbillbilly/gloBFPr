# ===========================================================================
# OpenFOAM case generator for urban wind
# ===========================================================================
# One solver, one mode: wind. The public API deliberately does not expose a
# `mode` argument in the current development branch.
#
# WHY ONE SOLVER
#   buoyantBoussinesqPimpleFoam solves the wind momentum and turbulence
#   equations while also carrying the temperature field required by the
#   Boussinesq formulation. In this wind-only setup every surface sits at T_ref
#   and walls are adiabatic, so the buoyancy body force is inactive and the
#   result is pure mechanical flow.
#
# WHY TRANSIENT
#   Urban wakes can be unsteady even under stationary boundary conditions. The
#   solver writes both instantaneous pedestrian slices and a UMean slice so
#   wind maps can use the mean field for comfort-style interpretation.
#
# OBLIQUE WIND
#   The old generator hardcoded the inlet to the x-min face and expressed
#   direction as a flowDir vector on it. For wind from the north that injects
#   velocity in +y through a face whose normal is -x, so the inflow normal
#   component is zero and effectively nothing enters. The docstring said to
#   "rotate the domain in prepare_openfoam_inputs()", but no such rotation
#   argument exists. Here blockMesh emits four separately named lateral
#   patches and each is assigned inlet / outlet / lateral from the sign of
#   dot(flowDir, outward normal), so any direction works without rotating
#   anything.
# ===========================================================================


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

#' @noRd
noc_foam_header <- function(class, object) {
  paste0(
    "/*---------------------------------------------------------------------------*\\\n",
    "| =========                 |                                                 |\n",
    "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n",
    "|  \\\\    /   O peration     | Version:  v2506                                 |\n",
    "|   \\\\  /    A nd           | Web:      www.openfoam.com                      |\n",
    "|    \\\\/     M anipulation  |                                                 |\n",
    "\\*---------------------------------------------------------------------------*/\n",
    "FoamFile\n{\n",
    "    version     2.0;\n",
    "    format      ascii;\n",
    "    class       ", class, ";\n",
    "    object      ", object, ";\n}\n",
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"
  )
}

#' @noRd
noc_write <- function(content, path) {
  writeLines(content, con = path)
  invisible(normalizePath(path, mustWork = FALSE))
}

#' Indented boundaryField entry
#' @noRd
fc_patch <- function(name, ...) {
  paste0("    ", name, "\n    {\n",
         paste0("        ", c(...), collapse = "\n"), "\n    }\n")
}

FC_LATERAL <- c("xMin", "xMax", "yMin", "yMax")


# ---------------------------------------------------------------------------
# Patch role assignment
# ---------------------------------------------------------------------------

#' Assign inlet / outlet / lateral roles to the four vertical faces
#'
#' A face is an inlet when the flow direction points into the domain through
#' it, i.e. when dot(flowDir, outward_normal) < 0.  Faces the wind runs
#' parallel to get `lateral` (slip).
#'
#' @param flow_dir Numeric length-2 unit-ish vector (x, y).
#' @param tol Numeric. Components below this are treated as parallel.
#' @return Named character vector over xMin/xMax/yMin/yMax.
#' @noRd
fc_patch_roles <- function(flow_dir = c(1, 0), tol = 1e-6) {
  fx <- flow_dir[1]; fy <- flow_dir[2]
  # Outward normals: xMin (-1,0), xMax (+1,0), yMin (0,-1), yMax (0,+1)
  dots <- c(xMin = -fx, xMax = fx, yMin = -fy, yMax = fy)
  roles <- ifelse(dots < -tol, "inlet",
                  ifelse(dots > tol, "outlet", "lateral"))
  stats::setNames(roles, FC_LATERAL)
}


# ---------------------------------------------------------------------------
# constant/
# ---------------------------------------------------------------------------

#' @noRd
make_gravity <- function() {
  paste0(noc_foam_header("uniformDimensionedVectorField", "g"),
         "dimensions  [ 0 1 -2 0 0 0 0 ];\n",
         "value       ( 0 0 -9.81 );\n")
}

#' @noRd
make_noc_transport_properties <- function(T_ref, beta = NULL, Prt = 1.0) {
  if (is.null(beta)) beta <- 1 / T_ref   # ideal-gas approximation
  paste0(
    noc_foam_header("dictionary", "transportProperties"),
    "transportModel  Newtonian;\n\n",
    "nu              nu [ 0 2 -1 0 0 0 0 ] 1.5e-05;\n\n",
    sprintf("beta            beta [ 0 0 0 -1 0 0 0 ] %.6g;\n", beta),
    sprintf("TRef            TRef [ 0 0 0  1 0 0 0 ] %g;\n", T_ref),
    "Pr              Pr   [ 0  0  0 0 0 0 0 ] 0.71;\n",
    # Prt above the neutral 0.85-0.9: turbulent heat transport is suppressed
    # more strongly than momentum transport under stable stratification.
    sprintf("Prt             Prt  [ 0  0  0 0 0 0 0 ] %g;\n", Prt)
  )
}

#' @noRd
make_noc_turbulence_properties <- function() {
  paste0(
    noc_foam_header("dictionary", "turbulenceProperties"),
    "simulationType  RAS;\n\n",
    "RAS\n{\n",
    # kOmegaSST rather than kEpsilon: its production limiter
    # (Pk <= 10*Cmu*k*omega) prevents the k blow-up kEpsilon shows in
    # buoyancy-driven flows starting from rest, and omegaWallFunction is far
    # more robust than epsilonWallFunction at low velocity.
    # atmBoundaryLayerInletOmega exists, so the ABL inlet works with SST too.
    "    RASModel        kOmegaSST;\n",
    "    turbulence      on;\n",
    "    printCoeffs     on;\n}\n"
  )
}

#' @noRd
make_noc_fv_options <- function(canopy_thermal = TRUE) {
  parts <- noc_foam_header("dictionary", "fvOptions")

  parts <- paste0(parts,
    "// Plant canopy as distributed drag rather than solid geometry: a crown\n",
    "// passes and drags air, it does not block it. Sp -= alpha*rho*Cd*LAD*|U|*U.\n",
    "// The drag is proportional to leafAreaDensity, so cells left at zero feel\n",
    "// nothing - which is why selectionMode all needs no cellSet.\n",
    "atmPlantCanopyUSource1\n{\n",
    "    type            atmPlantCanopyUSource;\n",
    "    atmPlantCanopyUSourceCoeffs\n    {\n        selectionMode   all;\n    }\n}\n\n",
    "atmPlantCanopyTurbSource1\n{\n",
    "    type            atmPlantCanopyTurbSource;\n",
    "    atmPlantCanopyTurbSourceCoeffs\n    {\n        selectionMode   all;\n    }\n}\n",
    if (isTRUE(canopy_thermal)) paste0(
      "\n// Canopy heat source. At night the crown radiates to the sky and cools\n",
      "// fast while the trunk space below is shielded, so the temperature\n",
      "// structure inside a stand differs from open grass.\n",
      "//\n",
      "// OFF BY DEFAULT. This source reads a volScalarField `qPlant` that the\n",
      "// package cannot estimate - it is a canopy energy flux, not something\n",
      "// derivable from a canopy height model. Enabled with qPlant = 0 it is\n",
      "// inert, so it would add a failure mode for no physics. Supply\n",
      "// `canopy_heat_source` to prepare_foam_case() to turn it on, and check\n",
      "// the units against the atmPlantCanopyTSource documentation for your\n",
      "// OpenFOAM version before trusting the magnitude.\n",
      "atmPlantCanopyTSource1\n{\n",
      "    type            atmPlantCanopyTSource;\n",
      "    atmPlantCanopyTSourceCoeffs\n    {\n        selectionMode   all;\n    }\n}\n")
    else "")
  parts
}


# ---------------------------------------------------------------------------
# system/
# ---------------------------------------------------------------------------

#' Pedestrian-level sample surface
#'
#' With terrain, "1.5 m above ground" is NOT z = 1.5. The mesh floor sits below
#' the DEM, so a flat cutting plane at absolute z = 1.5 lands underground and
#' the function object silently writes nothing but a warning:
#'     z1p5m : Plane (0 0 1) (0 0 1.5) does not intersect the mesh bounds
#'             (0 0 155.829) (2061.85 1944.11 445.373)
#' A flat plane could not represent 1.5 m AGL over sloping ground anyway, so
#' with terrain we use a distanceSurface offset from the terrain STL, which
#' follows the ground. Flat cases keep the cheaper cutting plane, anchored to
#' the domain floor rather than to zero.
#' @noRd
fc_slice_block <- function(z_slice, terrain_stl = NULL, z_floor = 0) {
  if (!is.null(terrain_stl)) {
    paste0(
      "            z1p5m\n            {\n",
      "                type        distanceSurface;\n",
      "                surfaceType triSurfaceMesh;\n",
      sprintf("                surfaceName %s;\n", basename(terrain_stl)),
      sprintf("                distance    %g;   // metres above ground\n", z_slice),
      # signed distance is positive on the side the STL normals face.
      # foam_terrain_solid_mesh() emits an outward-oriented closed solid
      # (verified: every directed edge has exactly one reverse twin, and the
      # divergence-theorem volume matches the analytic value), so +1.5 is
      # above the ground, not inside the rock.
      "                signed      true;\n",
      # Documented entries for distanceSurface are topology / absProximity /
      # maxDistance. `isoMethod` belongs to sampledIsoSurface and is NOT valid
      # here - OpenFOAM would likely ignore it, but silently-ignored entries
      # are how you end up believing a setting took effect.
      "            }\n")
  } else {
    paste0(
      "            z1p5m\n            {\n",
      "                type        cuttingPlane;\n",
      "                planeType   pointAndNormal;\n",
      "                pointAndNormalDict\n                {\n",
      sprintf("                    point  (0 0 %g);\n", z_floor + z_slice),
      "                    normal (0 0 1);\n                }\n",
      "                interpolate true;\n",
      "            }\n")
  }
}

#' @noRd
make_noc_control_dict <- function(end_time_s, write_interval_s,
                                  delta_t = 0.5, max_co = 5, max_delta_t = 20,
                                  z_slice = 1.5, fields = c("T", "U"),
                                  terrain_stl = NULL, z_floor = 0) {
  # One flow-through of run-in before averaging starts. sim_hours = NULL gives
  # three flow-throughs, so this discards the first and averages the last two.
  avg_start_s <- end_time_s / 3

  paste0(
    noc_foam_header("dictionary", "controlDict"),
    "application     buoyantBoussinesqPimpleFoam;\n\n",
    "startFrom       startTime;\nstartTime       0;\nstopAt          endTime;\n",
    sprintf("endTime         %g;   // seconds of physical time\n", end_time_s),
    sprintf("deltaT          %g;\n\n", delta_t),
    # adjustableRunTime keeps write times exactly on the interval even though
    # the step size floats.
    "writeControl    adjustableRunTime;\n",
    sprintf("writeInterval   %g;\n", write_interval_s),
    "purgeWrite      0;\n\n",
    "writeFormat     ascii;\nwritePrecision  8;\nwriteCompression off;\n\n",
    "timeFormat      general;\ntimePrecision   6;\n\nrunTimeModifiable true;\n\n",
    # maxCo well above 1 is safe under PIMPLE: the outer correctors re-converge
    # momentum and pressure within the step. The limit is temporal accuracy of
    # the cooling ramp, not stability.
    "adjustTimeStep  yes;\n",
    sprintf("maxCo           %g;\nmaxDeltaT       %g;\n\n", max_co, max_delta_t),
    "functions\n{\n",
    "    pedestrianSlice\n    {\n",
    "        type            surfaces;\n",
    "        libs            (\"libsampling.so\");\n\n",
    "        executeControl  writeTime;\n        writeControl    writeTime;\n\n",
    "        surfaceFormat   raw;\n",
    sprintf("        fields          ( %s );\n\n", paste(fields, collapse = " ")),
    "        surfaces\n        {\n",
    fc_slice_block(z_slice, terrain_stl, z_floor),
    "        }\n    }\n\n",
    # Time-averaged velocity.
    #
    # A URANS run's LAST time is an instantaneous snapshot, not a mean. Urban
    # wakes shed, so that snapshot carries whichever phase of the shedding the
    # run happened to stop on, and two runs differing only in endTime give
    # visibly different maps behind the large blocks. Pedestrian comfort is
    # assessed on a MEAN field, which is what this produces.
    #
    # Averaging starts one flow-through in (endTime/3, matching the default
    # "run three flow-throughs") so the startup transient is excluded rather
    # than averaged in. `timeStart` is a base functionObject option, stable
    # across releases; because BOTH objects below carry it, UMean is never
    # sampled before fieldAverage has created it.
    "    fieldAverage1\n    {\n",
    "        type            fieldAverage;\n",
    "        libs            (\"libfieldFunctionObjects.so\");\n",
    sprintf("        timeStart       %g;\n", avg_start_s),
    "        writeControl    writeTime;\n",
    # Accumulate over the whole averaging window. restartOnOutput true would
    # reset at every write, so each written UMean would cover only the last
    # write interval - a mean over 1/8 of the run, silently.
    "        restartOnRestart false;\n        restartOnOutput  false;\n\n",
    "        fields\n        (\n            U\n            {\n",
    "                mean        on;\n",
    "                prime2Mean  off;\n",
    "                base        time;\n            }\n        );\n    }\n\n",
    # Its own surfaces object, so postProcessing/pedestrianSlice is untouched:
    # read_foam_pedestrian_slice() reads that directory and matches raw files
    # by exact name, so its behaviour does not change at all. The mean lands
    # in postProcessing/pedestrianSliceMean/.
    "    pedestrianSliceMean\n    {\n",
    "        type            surfaces;\n",
    "        libs            (\"libsampling.so\");\n",
    sprintf("        timeStart       %g;\n\n", avg_start_s),
    "        executeControl  writeTime;\n        writeControl    writeTime;\n\n",
    "        surfaceFormat   raw;\n",
    "        fields          ( UMean );\n\n",
    "        surfaces\n        {\n",
    fc_slice_block(z_slice, terrain_stl, z_floor),
    "        }\n    }\n\n",
    "    minMaxU\n    {\n",
    "        type            fieldMinMax;\n",
    "        libs            (\"libfieldFunctionObjects.so\");\n",
    "        writeControl    timeStep;\n        writeInterval   100;\n",
    sprintf("        fields          ( %s );\n", paste(fields, collapse = " ")),
    "    }\n}\n"
  )
}

#' @noRd
make_noc_fv_schemes <- function(ddt_scheme = "Euler") {
  paste0(
    noc_foam_header("dictionary", "fvSchemes"),
    # Euler is first order but unconditionally bounded and safe with a floating
    # time step. 'backward' is second order but assumes constant deltaT, so it
    # also needs adjustTimeStep no.
    sprintf("ddtSchemes\n{\n    default         %s;\n}\n\n", ddt_scheme),
    "gradSchemes\n{\n    default         Gauss linear;\n",
    "    grad(U)         cellLimited Gauss linear 1;\n",
    "    grad(T)         cellLimited Gauss linear 1;\n}\n\n",
    "divSchemes\n{\n    default                             none;\n",
    # No 'bounded' prefix anywhere. 'bounded Gauss ...' adds -div(phi)*U to
    # compensate for a non-zero continuity error while a STEADY solution is
    # still converging. In a transient run continuity is enforced every step,
    # so the term is at best a no-op and at worst a spurious source.
    "    div(phi,U)                          Gauss linearUpwind grad(U);\n",
    # limitedLinear is second order but TVD-limited: keeps the temperature
    # front sharp without the overshoots pure linear gives at the inversion top.
    "    div(phi,T)                          Gauss limitedLinear 1;\n",
    "    div(phi,k)                          Gauss limitedLinear 1;\n",
    "    div(phi,omega)                      Gauss limitedLinear 1;\n",
    "    div((nuEff*dev2(T(grad(U)))))       Gauss linear;\n}\n\n",
    "laplacianSchemes\n{\n    default         Gauss linear corrected;\n}\n\n",
    "interpolationSchemes\n{\n    default         linear;\n}\n\n",
    "snGradSchemes\n{\n    default         corrected;\n}\n\n",
    # kOmegaSST needs wall distance for its F1/F2 blending; v2506 removed the
    # default and requires an explicit method.
    "wallDist\n{\n    method          meshWave;\n}\n"
  )
}

#' @noRd
make_noc_fv_solution <- function(n_outer_correctors = 2L, n_correctors = 2L) {
  paste0(
    noc_foam_header("dictionary", "fvSolution"),
    "solvers\n{\n",
    # GAMG with GaussSeidel, not DICGaussSeidel: DIC computes 1/sqrt(diag) and
    # dies with SIGFPE when a coarsened parallel matrix has a near-zero
    # diagonal, which is common in open-boundary buoyancy cases.
    #
    # Non-final pressure solves may use a loose relTol because the outer
    # correctors repeat them; only the Final solve must be tight. That is the
    # opposite of the steady case, where relTol 0 was needed on every p_rgh
    # solve because there was no outer loop to recover.
    "    p_rgh\n    {\n        solver          GAMG;\n",
    "        agglomerator    faceAreaPair;\n        mergeLevels     1;\n",
    "        nCellsInCoarsestLevel 200;\n        smoother        GaussSeidel;\n",
    "        tolerance       1e-7;\n        relTol          0.01;\n",
    "        maxIter         500;\n    }\n\n",
    "    p_rghFinal\n    {\n        $p_rgh;\n",
    "        tolerance       1e-8;\n        relTol          0;\n    }\n\n",
    "    \"(U|T|k|omega)\"\n    {\n        solver          smoothSolver;\n",
    "        smoother        symGaussSeidel;\n",
    "        tolerance       1e-8;\n        relTol          0.1;\n    }\n\n",
    # Without these Final entries OpenFOAM aborts on the first time step with
    # "solver for UFinal not found".
    "    \"(U|T|k|omega)Final\"\n    {\n        $U;\n        relTol          0;\n    }\n}\n\n",
    "PIMPLE\n{\n",
    "    momentumPredictor           yes;\n",
    sprintf("    nOuterCorrectors            %d;\n", as.integer(n_outer_correctors)),
    sprintf("    nCorrectors                 %d;\n", as.integer(n_correctors)),
    "    nNonOrthogonalCorrectors    1;\n",
    "    pRefCell                    0;\n    pRefValue                   0;\n\n",
    # Bail out of the outer loop early once the step has converged; on a slowly
    # evolving cooling ramp most steps converge in one outer iteration.
    #
    # BOTH spellings are written on purpose. ESI (openfoam.com) reads
    # `residualControl`; the OpenFOAM Foundation branch reads
    # `outerCorrectorResidualControl`. OpenFOAM ignores dictionary keys it does
    # not recognise, so writing both is portable. Get it wrong and the run does
    # not fail - it just prints
    #     "PIMPLE: no residual control data found. Calculations will employ
    #      2 corrector loops"
    # and silently does the maximum number of outer iterations every step.
    "    residualControl\n    {\n",
    "        p_rgh   { tolerance 1e-4; relTol 0; }\n",
    "        U       { tolerance 1e-5; relTol 0; }\n",
    "        T       { tolerance 1e-5; relTol 0; }\n    }\n\n",
    "    outerCorrectorResidualControl\n    {\n",
    "        p_rgh   { tolerance 1e-4; relTol 0; }\n",
    "        U       { tolerance 1e-5; relTol 0; }\n",
    "        T       { tolerance 1e-5; relTol 0; }\n    }\n}\n\n",
    # Under-relaxation applies only to non-final outer iterations; the Final
    # iteration must be unrelaxed or the time derivative is inconsistent.
    # Far gentler than the steady case needed (U 0.2, T 0.2) because the time
    # derivative itself now stabilises the system.
    "relaxationFactors\n{\n",
    "    fields\n    {\n        p_rgh       0.7;\n        p_rghFinal  1.0;\n    }\n",
    "    equations\n    {\n",
    "        \"(U|T|k|omega)\"       0.7;\n",
    "        \"(U|T|k|omega)Final\"  1.0;\n    }\n}\n"
  )
}

#' blockMeshDict with four separately named lateral patches
#' @noRd
make_noc_block_mesh_dict <- function(domain, nx, ny, nz, z_grading = 8) {
  d <- domain
  paste0(
    noc_foam_header("dictionary", "blockMeshDict"),
    "scale 1;\n\nvertices\n(\n",
    sprintf("    (%g %g %g)\n", d$xmin, d$ymin, d$zmin),
    sprintf("    (%g %g %g)\n", d$xmax, d$ymin, d$zmin),
    sprintf("    (%g %g %g)\n", d$xmax, d$ymax, d$zmin),
    sprintf("    (%g %g %g)\n", d$xmin, d$ymax, d$zmin),
    sprintf("    (%g %g %g)\n", d$xmin, d$ymin, d$zmax),
    sprintf("    (%g %g %g)\n", d$xmax, d$ymin, d$zmax),
    sprintf("    (%g %g %g)\n", d$xmax, d$ymax, d$zmax),
    sprintf("    (%g %g %g)\n", d$xmin, d$ymax, d$zmax),
    ");\n\nblocks\n(\n",
    # Graded in z: the breeze / drainage layer is only tens of metres deep, so
    # cells cluster near the ground. Ratio is top cell / bottom cell.
    sprintf("    hex (0 1 2 3 4 5 6 7) (%d %d %d) simpleGrading (1 1 %g)\n",
            nx, ny, nz, z_grading),
    ");\n\nedges\n();\n\n",
    # Four named lateral patches, so any wind direction can assign inlet and
    # outlet roles without rotating the domain.
    "boundary\n(\n",
    "    xMin\n    {\n        type patch;\n        faces ((0 4 7 3));\n    }\n",
    "    xMax\n    {\n        type patch;\n        faces ((1 2 6 5));\n    }\n",
    "    yMin\n    {\n        type patch;\n        faces ((0 1 5 4));\n    }\n",
    "    yMax\n    {\n        type patch;\n        faces ((3 7 6 2));\n    }\n",
    "    ground\n    {\n        type wall;\n        faces ((3 2 1 0));\n    }\n",
    "    top\n    {\n        type patch;\n        faces ((4 5 6 7));\n    }\n",
    ");\n\nmergePatchPairs\n();\n"
  )
}

#' snappyHexMeshDict with optional terrain geometry
#' @noRd
make_noc_snappy_dict <- function(stl_name, loc_x, loc_y, loc_z, ref,
                                 terrain_stl      = NULL,
                                 terrain_ref      = 0L,
                                 max_local_cells  = 1000000L,
                                 max_global_cells = 4000000L) {
  patch_name <- tools::file_path_sans_ext(stl_name)
  geom <- paste0("    ", stl_name, "\n    {\n        type  triSurfaceMesh;\n",
                 "        name  ", patch_name, ";\n    }\n")
  refs <- paste0("        ", patch_name, "\n        {\n",
                 sprintf("            level (%d %d);\n", ref, ref), "        }\n")
  if (!is.null(terrain_stl)) {
    tname <- tools::file_path_sans_ext(basename(terrain_stl))
    geom <- paste0(geom, "    ", basename(terrain_stl),
                   "\n    {\n        type  triSurfaceMesh;\n",
                   "        name  ", tname, ";\n    }\n")
    refs <- paste0(refs, "        ", tname, "\n        {\n",
                   sprintf("            level (%d %d);\n", terrain_ref, terrain_ref),
                   "        }\n")
  }
  paste0(
    noc_foam_header("dictionary", "snappyHexMeshDict"),
    "castellatedMesh true;\nsnap            true;\naddLayers       false;\n\n",
    "geometry\n{\n", geom, "}\n\n",
    "castellatedMeshControls\n{\n",
    sprintf("    maxLocalCells           %d;\n", as.integer(max_local_cells)),
    sprintf("    maxGlobalCells         %d;\n", as.integer(max_global_cells)),
    "    minRefinementCells          10;\n    maxLoadUnbalance          0.10;\n",
    "    nCellsBetweenLevels           3;\n\n    features ();\n\n",
    "    refinementSurfaces\n    {\n", refs, "    }\n\n",
    "    resolveFeatureAngle 30;\n\n    refinementRegions {}\n\n",
    # locationInMesh must sit in the AIR: above the terrain and outside every
    # building, or snappy keeps the wrong side of the surfaces.
    sprintf("    locationInMesh (%g %g %g);\n", loc_x, loc_y, loc_z),
    "    allowFreeStandingZoneFaces true;\n}\n\n",
    "snapControls\n{\n    nSmoothPatch              3;\n    tolerance                 2.0;\n",
    "    nSolveIter               30;\n    nRelaxIter                5;\n",
    "    nFeatureSnapIter         10;\n    implicitFeatureSnap   false;\n",
    "    explicitFeatureSnap   false;\n    multiRegionFeatureSnap false;\n}\n\n",
    "addLayersControls\n{\n    relativeSizes           true;\n    layers {}\n",
    "    expansionRatio          1.3;\n    finalLayerThickness     0.3;\n",
    "    minThickness            0.1;\n    nGrow                     0;\n",
    "    featureAngle             60;\n    nRelaxIter                3;\n",
    "    nSmoothSurfaceNormals     1;\n    nSmoothNormals            3;\n",
    "    nSmoothThickness         10;\n    maxFaceThicknessRatio   0.5;\n",
    "    maxThicknessToMedialRatio 0.3;\n    minMedialAxisAngle       90;\n",
    "    nBufferCellsNoExtrude     0;\n    nLayerIter               50;\n}\n\n",
    "meshQualityControls\n{\n    maxNonOrtho              70;\n",
    "    maxBoundarySkewness      20;\n    maxInternalSkewness       4;\n",
    "    maxConcave               80;\n    minVol                 1e-13;\n",
    "    minTetQuality           1e-9;\n    minArea                  -1;\n",
    "    minTwist               0.02;\n    minDeterminant          0.001;\n",
    "    minFaceWeight           0.05;\n    minVolRatio             0.01;\n",
    "    minTriangleTwist         -1;\n    nSmoothScale              4;\n",
    "    errorReduction         0.75;\n}\n\ndebug 0;\nmergeTolerance 1e-6;\n"
  )
}


# ---------------------------------------------------------------------------
# Domain adequacy pre-flight
# ---------------------------------------------------------------------------
# The bright rim on pedestrian maps is a DOMAIN SIZING problem, and nothing
# downstream can undo it: by the time run_openfoam_docker() is called the mesh
# and the boundary conditions already exist. Over the building-free apron the
# ABL inlet profile is simply undisturbed, so that band shows what the wind
# would be with no city.
#
# It cannot be removed, only made small enough to discard cheaply. The usual
# guidance (COST Action 732; AIJ, Tominaga et al. 2008) is, with H the tallest
# building in the area of interest:
#     >= 5H  upwind and lateral
#     >= 15H downwind
# `prepare_openfoam_inputs(domain_buffer = )` sets this, and its default of
# 100 m is far short of that for anything but low-rise. Warning at case-write
# time costs nothing; discovering it after an hour of solving costs an hour.

#' Bounding box and height scale of an STL (ASCII or binary)
#' @noRd
fc_stl_bounds <- function(path) {
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  hdr <- readBin(con, "raw", 84L)
  is_ascii <- grepl("^solid", rawToChar(hdr[1:5]))
  if (is_ascii) {
    v <- grep("vertex", readLines(path, warn = FALSE), value = TRUE)
    if (!length(v)) return(NULL)
    m <- utils::read.table(text = sub("^\\s*vertex\\s+", "", v),
                           col.names = c("x", "y", "z"))
  } else {
    n <- readBin(hdr[81:84], "integer", size = 4L, endian = "little")
    if (!is.finite(n) || n <= 0) return(NULL)
    raw <- readBin(con, "raw", n * 50L)
    if (length(raw) < n * 50L) return(NULL)
    idx <- rep((seq_len(n) - 1L) * 50L, each = 9L) +
           rep(12L + seq_len(9L) * 4L - 4L, times = n)
    vals <- vapply(idx, function(i)
      readBin(raw[(i + 1L):(i + 4L)], "numeric", size = 4L, endian = "little"),
      numeric(1))
    m <- data.frame(x = vals[c(TRUE, FALSE, FALSE)],
                    y = vals[c(FALSE, TRUE, FALSE)],
                    z = vals[c(FALSE, FALSE, TRUE)])
  }
  i <- which.max(m$z)
  list(xmin = min(m$x), xmax = max(m$x),
       ymin = min(m$y), ymax = max(m$y),
       zmin = min(m$z), zmax = max(m$z),
       peak_x = m$x[i], peak_y = m$y[i])
}

#' Warn when the domain is too small for the tallest building
#' @noRd
fc_check_domain <- function(domain, stl_file, fd, ground_top = NULL,
                            dem = NULL, quiet = FALSE) {
  b <- tryCatch(fc_stl_bounds(stl_file), error = function(e) NULL)
  if (is.null(b)) return(invisible(NULL))

  # Building height above ITS OWN ground. Subtracting the terrain MAXIMUM is
  # wrong whenever the tallest building does not stand on the highest ground -
  # on a real site that reported a 35 m tallest building where the true figure
  # was 98.5 m. Sample the DEM directly beneath the highest STL vertex.
  base <- NA_real_
  if (!is.null(dem) && requireNamespace("terra", quietly = TRUE) &&
      !is.null(b$peak_x)) {
    base <- tryCatch(
      terra::extract(dem, cbind(b$peak_x, b$peak_y), method = "bilinear")[, 1],
      error = function(e) NA_real_)
  }
  if (!is.finite(base))
    base <- if (!is.null(ground_top) && is.finite(ground_top)) ground_top else b$zmin
  H <- b$zmax - base
  if (!is.finite(H) || H <= 0) return(invisible(NULL))

  gap <- c(W = b$xmin - domain$xmin, E = domain$xmax - b$xmax,
           S = b$ymin - domain$ymin, N = domain$ymax - b$ymax)
  # Classify each side against the flow. Downwind needs the long fetch.
  role <-
    c(W = if (fd[1] > 0.1) "upwind" else if (fd[1] < -0.1) "downwind" else "lateral",
      E = if (fd[1] < -0.1) "upwind" else if (fd[1] > 0.1) "downwind" else "lateral",
      S = if (fd[2] > 0.1) "upwind" else if (fd[2] < -0.1) "downwind" else "lateral",
      N = if (fd[2] < -0.1) "upwind" else if (fd[2] > 0.1) "downwind" else "lateral")
  need <- ifelse(role == "downwind", 15 * H, 5 * H)
  short <- gap < need

  if (any(short) && !isTRUE(quiet)) {
    warning(sprintf(
      paste0("Domain may be too small for the tallest building (%.0f m).\n",
             "  guidance (COST 732 / AIJ): >=5H upwind and lateral, >=15H downwind\n",
             "%s",
             "%s",
             "  To reduce it, raise `domain_buffer` in prepare_openfoam_inputs()\n",
             "  (currently giving %.0f m)."),
      H,
      paste(sprintf("  %s (%s): have %.0f m, want %.0f m%s\n",
                    names(gap), role, gap, need,
                    ifelse(short, "   <-- short", "")), collapse = ""),
      paste0("  Consequence: the building-free apron carries the undisturbed\n",
             "  inlet profile, so pedestrian maps show a bright rim there.\n",
             "  read_foam_pedestrian_slice(trim=) crops it but cannot undo it.\n"),
      min(gap)), call. = FALSE)
  } else if (!isTRUE(quiet)) {
    message(sprintf("Domain check: tallest building %.0f m, margins W/E/S/N = %.0f/%.0f/%.0f/%.0f m - adequate",
                    H, gap[1], gap[2], gap[3], gap[4]))
  }
  invisible(list(H = H, gap = gap, need = need, role = role))
}


# ---------------------------------------------------------------------------
# 0/ fields - role driven
# ---------------------------------------------------------------------------

#' @noRd
fc_abl_entries <- function(u_mag, fd, z_ref, z0, z_ground) {
  c(sprintf("flowDir         (%g %g 0);", fd[1], fd[2]),
    "zDir            (0 0 1);",
    sprintf("Uref            %g;", u_mag),
    sprintf("Zref            %g;", z_ref),
    sprintf("z0              uniform %g;", z0),
    "d               uniform 0;",
    sprintf("zGround         uniform %g;", z_ground))
}

#' @noRd
make_noc_u_field <- function(roles, walls, u_mag, fd, z_ref, z0, z_ground,
                             top_open) {
  lat <- vapply(names(roles), function(p) {
    switch(roles[[p]],
      inlet   = fc_patch(p, "type            atmBoundaryLayerInletVelocity;",
                         fc_abl_entries(u_mag, fd, z_ref, z0, z_ground),
                         "value           $internalField;"),
      outlet  = fc_patch(p, "type            pressureInletOutletVelocity;",
                         "value           uniform (0 0 0);"),
      lateral = fc_patch(p, "type            slip;"),
      open    = fc_patch(p, "type            pressureInletOutletVelocity;",
                         "value           uniform (0 0 0);"))
  }, character(1))

  paste0(
    noc_foam_header("volVectorField", "U"),
    "dimensions      [ 0 1 -1 0 0 0 0 ];\n\n",
    sprintf("internalField   uniform (%g %g 0);\n\n", u_mag * fd[1], u_mag * fd[2]),
    "boundaryField\n{\n", paste0(lat, collapse = ""),
    paste0(vapply(walls, function(w) fc_patch(w, "type            noSlip;"),
                  character(1)), collapse = ""),
    if (top_open) fc_patch("top", "type            pressureInletOutletVelocity;",
                           "value           uniform (0 0 0);")
    else fc_patch("top", "type            slip;"),
    "}\n")
}

#' @noRd
make_noc_p_rgh_field <- function(roles, walls, top_open) {
  # fixedValue 0 on outflow/open faces provides the pressure reference and is
  # the correct partner for pressureInletOutletVelocity on U.
  # fixedFluxPressure would be WRONG there: it adjusts the pressure gradient to
  # match U, while pressureInletOutletVelocity derives U FROM the pressure
  # gradient - circular, and continuity explodes. fixedFluxPressure belongs on
  # walls and slip faces, where the normal flux is fixed by U (zero).
  lat <- vapply(names(roles), function(p) {
    switch(roles[[p]],
      inlet   = fc_patch(p, "type            fixedFluxPressure;",
                         "value           uniform 0;"),
      outlet  = fc_patch(p, "type            fixedValue;", "value           uniform 0;"),
      lateral = fc_patch(p, "type            fixedFluxPressure;",
                         "value           uniform 0;"),
      open    = fc_patch(p, "type            fixedValue;", "value           uniform 0;"))
  }, character(1))

  paste0(
    noc_foam_header("volScalarField", "p_rgh"),
    "dimensions      [ 0 2 -2 0 0 0 0 ];\n\ninternalField   uniform 0;\n\n",
    "boundaryField\n{\n", paste0(lat, collapse = ""),
    paste0(vapply(walls, function(w)
      fc_patch(w, "type            fixedFluxPressure;", "value           uniform 0;"),
      character(1)), collapse = ""),
    if (top_open) fc_patch("top", "type            fixedValue;", "value           uniform 0;")
    else fc_patch("top", "type            fixedFluxPressure;", "value           uniform 0;"),
    "}\n")
}

#' @noRd
make_noc_t_field <- function(roles, ground_patches, buildings_patch, top_open,
                             T_ref) {
  # Every surface sits at T_ref and the walls are adiabatic, so the Boussinesq
  # body force is identically zero and the solve is pure mechanical flow. The
  # field is still written because buoyantBoussinesqPimpleFoam requires it.
  amb <- c("type            inletOutlet;",
           sprintf("inletValue      uniform %g;", T_ref),
           sprintf("value           uniform %g;", T_ref))
  lat <- vapply(names(roles), function(p) {
    switch(roles[[p]],
      inlet   = fc_patch(p, "type            fixedValue;",
                         sprintf("value           uniform %g;", T_ref)),
      outlet  = fc_patch(p, amb),
      lateral = fc_patch(p, "type            zeroGradient;"))
  }, character(1))

  ground_block <- paste0(vapply(ground_patches, function(gp)
    fc_patch(gp, "type            zeroGradient;"), character(1)), collapse = "")
  bld_block <- fc_patch(buildings_patch, "type            zeroGradient;")

  paste0(
    noc_foam_header("volScalarField", "T"),
    "dimensions      [ 0 0 0 1 0 0 0 ];\n\n",
    sprintf("internalField   uniform %g;\n\n", T_ref),
    "boundaryField\n{\n", paste0(lat, collapse = ""), ground_block, bld_block,
    if (top_open) fc_patch("top", amb) else fc_patch("top", "type            zeroGradient;"),
    "}\n")
}

#' @noRd
make_noc_turb_scalar_field <- function(object, dims, init, roles, walls,
                                       wall_fn, top_open, abl_type = NULL,
                                       u_mag = 0, fd = c(1, 0), z_ref = 10,
                                       z0 = 0.1, z_ground = 0) {
  # inletOutlet rather than zeroGradient on open faces: when buoyancy-driven
  # reverse flow re-enters the domain it must carry the safe reference value,
  # not an unbounded extrapolation (which blows up exponentially).
  io <- c("type            inletOutlet;",
          sprintf("inletValue      uniform %g;", init),
          sprintf("value           uniform %g;", init))
  lat <- vapply(names(roles), function(p) {
    switch(roles[[p]],
      inlet   = if (is.null(abl_type)) fc_patch(p, io) else
                fc_patch(p, sprintf("type            %s;", abl_type),
                         fc_abl_entries(u_mag, fd, z_ref, z0, z_ground),
                         "value           $internalField;"),
      outlet  = fc_patch(p, io),
      lateral = fc_patch(p, "type            zeroGradient;"),
      open    = fc_patch(p, io))
  }, character(1))

  paste0(
    noc_foam_header("volScalarField", object),
    sprintf("dimensions      %s;\n\n", dims),
    sprintf("internalField   uniform %g;\n\n", init),
    "boundaryField\n{\n", paste0(lat, collapse = ""),
    paste0(vapply(walls, function(w)
      fc_patch(w, sprintf("type            %s;", wall_fn),
               sprintf("value           uniform %g;", init)),
      character(1)), collapse = ""),
    if (top_open) fc_patch("top", io) else fc_patch("top", "type            zeroGradient;"),
    "}\n")
}

#' @noRd
make_noc_nut_field <- function(roles, ground_patches, buildings_patch, z0,
                               top_open) {
  ks <- z0 * 20   # Ks ~ 20*z0, standard urban-terrain approximation
  calc <- c("type            calculated;", "value           uniform 0;")
  lat <- vapply(names(roles), function(p) fc_patch(p, calc), character(1))
  paste0(
    noc_foam_header("volScalarField", "nut"),
    "dimensions      [ 0 2 -1 0 0 0 0 ];\n\ninternalField   uniform 0;\n\n",
    "boundaryField\n{\n", paste0(lat, collapse = ""),
    paste0(vapply(ground_patches, function(gp)
      fc_patch(gp, "type            nutkRoughWallFunction;",
               sprintf("Ks              uniform %g;", ks),
               "Cs              uniform 0.5;",   # v2506: Cs is Field<scalar>
               "value           uniform 0;"), character(1)), collapse = ""),
    fc_patch(buildings_patch, "type            nutkWallFunction;",
             "value           uniform 0;"),
    fc_patch("top", calc), "}\n")
}

#' @noRd
make_noc_alphat_field <- function(roles, walls, top_open, Prt = 1.0) {
  calc <- c("type            calculated;", "value           uniform 0;")
  lat <- vapply(names(roles), function(p) fc_patch(p, calc), character(1))
  paste0(
    noc_foam_header("volScalarField", "alphat"),
    "dimensions      [ 0 2 -1 0 0 0 0 ];\n\ninternalField   uniform 0;\n\n",
    "boundaryField\n{\n", paste0(lat, collapse = ""),
    paste0(vapply(walls, function(w)
      fc_patch(w, "type            alphatJayatillekeWallFunction;",
               sprintf("Prt             %g;", Prt),
               "value           uniform 0;"), character(1)), collapse = ""),
    fc_patch("top", calc), "}\n")
}

#' @noRd
make_foam_canopy_field <- function(object, dims, value = 0) {
  paste0(noc_foam_header("volScalarField", object),
         sprintf("dimensions      %s;\n\ninternalField   uniform %g;\n\n",
                 dims, value),
         "boundaryField\n{\n    \".*\"\n    {\n        type            zeroGradient;\n    }\n}\n")
}


# ---------------------------------------------------------------------------
# Allrun
# ---------------------------------------------------------------------------

#' @noRd
make_noc_allrun <- function(set_fields = FALSE) {
  paste0(
    "#!/bin/sh\n",
    "# Allrun - generated by gloBFPr::prepare_foam_case()\n",
    "# Pass NPROC (docker -e NPROC=N) to run in parallel.\n",
    "cd \"${0%/*}\" || exit 1\n\n",
    "NPROC=${NPROC:-1}\n",
    "SOLVER=buoyantBoussinesqPimpleFoam\n\n",
    "rm -rf processor* log.*\n",
    # `! -name '0.orig'` is load-bearing: 0.orig matches [0-9]* and is not
    # named '0', so without it this line deletes the pristine copy the very
    # first time it runs, and the restore below silently does nothing.
    "find . -maxdepth 1 -name '[0-9]*' ! -name '0' ! -name '0.orig' -exec rm -rf {} +\n",
    # postProcessing/ was NOT cleared here, and that is a quiet correctness
    # trap rather than untidiness.  Sampled surfaces are written under
    # postProcessing/<function>/<time>/, so a case re-run in a different mode -
    # or on a re-meshed geometry - leaves the old run's slices sitting beside
    # the new ones with nothing to tell them apart.  Measured on a real case:
    # one pedestrianSlice directory held writes from three different runs
    # spanning two days and two meshes, and reading "the" t = 3600 s slice
    # returned a wind-driven field of 4.7 m/s from the previous day while the
    # current wind run had different velocity scales. Both looked like valid
    # results.
    "rm -rf postProcessing\n\n",
    # setFields REWRITES 0/ in place, turning uniform fields into nonuniform
    # lists sized to the snappyHexMesh cell count. Re-running then rebuilds a
    # smaller blockMesh-only mesh and decomposePar dies with
    #     Size 1733307 is not equal to the expected length 1170585
    #     file: 0/Cd/internalField
    # Restoring the pristine copy is the standard OpenFOAM 0.orig idiom and is
    # what makes the case re-runnable.
    "# setFields rewrites 0/ in place; restore the pristine copy each run.\n",
    "[ -d 0.orig ] && { rm -rf 0; cp -r 0.orig 0; }\n\n",
    # Incremental meshing. Meshing dominates wall-clock (minutes on a million
    # cells) while the solver's startup failures surface in seconds, so
    # re-meshing unconditionally makes every debug cycle cost a full rebuild.
    # Re-mesh only when something the mesh actually depends on is newer than
    # the mesh itself. FOAM_FORCE_MESH=1 forces a rebuild.
    "NEED_MESH=1\n",
    "if [ -f constant/polyMesh/owner ] && [ -z \"$FOAM_FORCE_MESH\" ]; then\n",
    "    NEED_MESH=0\n",
    "    for f in system/blockMeshDict system/snappyHexMeshDict constant/triSurface/*; do\n",
    "        [ -e \"$f\" ] || continue\n",
    "        [ \"$f\" -nt constant/polyMesh/owner ] && NEED_MESH=1\n",
    "    done\n",
    "fi\n\n",
    "command -v $SOLVER >/dev/null 2>&1 || {\n",
    "    echo \"ERROR: $SOLVER not found in this OpenFOAM build.\"\n",
    "    echo \"       Recent ESI releases may ship only buoyantPimpleFoam.\"\n",
    "    echo \"       Check with: buoyantPimpleFoam -help\"\n",
    "    exit 1\n}\n\n",
    "if [ \"$NEED_MESH\" -eq 0 ]; then\n",
    "    echo '[1/4] blockMesh ... reusing existing mesh'\n",
    "    echo '[2/4] snappyHexMesh ... reusing existing mesh'\n",
    "    echo '      (FOAM_FORCE_MESH=1 to rebuild)'\n",
    "else\n",
    "echo '[1/4] blockMesh ...'\n",
    "rm -rf constant/polyMesh\n",
    "blockMesh > log.blockMesh 2>&1 && echo '      OK' \\\n",
    "  || { echo '      FAILED - see log.blockMesh'; exit 1; }\n\n",
    "echo '[2/4] snappyHexMesh ...'\n",
    "if [ \"$NPROC\" -gt 1 ]; then\n",
    "    sed -i \"s/numberOfSubdomains.*/numberOfSubdomains $NPROC;/\" system/decomposeParDict\n",
    "    decomposePar -force > log.decomposeMesh 2>&1 \\\n",
    "      || { echo '      FAILED - see log.decomposeMesh'; exit 1; }\n",
    "    mpirun -np \"$NPROC\" --allow-run-as-root snappyHexMesh -overwrite -parallel > log.snappyHexMesh 2>&1 \\\n",
    "      || { echo '      FAILED - see log.snappyHexMesh'; exit 1; }\n",
    "    reconstructParMesh -constant -mergeTol 1e-6 > log.reconstructMesh 2>&1 \\\n",
    "      || { echo '      FAILED - see log.reconstructMesh'; exit 1; }\n",
    "    rm -rf processor*\n    echo '      OK (parallel)'\n",
    "else\n",
    "    snappyHexMesh -overwrite > log.snappyHexMesh 2>&1 && echo '      OK' \\\n",
    "      || { echo '      FAILED - see log.snappyHexMesh'; exit 1; }\n",
    "fi\n",
    "fi\n\n",
    "echo '[3/4] setFields ...'\n",
    if (set_fields) paste0(
      "setFields > log.setFields 2>&1 && echo '      OK' \\\n",
      "  || { echo '      FAILED - see log.setFields'; exit 1; }\n\n")
    else "echo '      skipped (no canopy)'\n\n",
    "echo \"[4/4] $SOLVER (transient) ...\"\n",
    "if [ \"$NPROC\" -gt 1 ]; then\n",
    "    sed -i \"s/numberOfSubdomains.*/numberOfSubdomains $NPROC;/\" system/decomposeParDict\n",
    "    decomposePar -force > log.decompose 2>&1 \\\n",
    "      || { echo '      FAILED - see log.decompose'; exit 1; }\n",
    "    mpirun -np \"$NPROC\" --allow-run-as-root $SOLVER -parallel > log.solver 2>&1\n",
    "    STATUS=$?\n",
    "    reconstructPar > log.reconstruct 2>&1\n",
    "    [ $STATUS -eq 0 ] && echo '      OK' \\\n",
    "      || { echo '      FAILED - see log.solver'; exit $STATUS; }\n",
    "else\n",
    "    $SOLVER > log.solver 2>&1 && echo '      OK' \\\n",
    "      || { echo '      FAILED - see log.solver'; exit 1; }\n",
    "fi\n\n",
    "echo 'Done. Check log.solver for Courant number and time-step history.'\n")
}


# ===========================================================================
# Public function
# ===========================================================================

#' Prepare a transient OpenFOAM case for urban wind
#'
#' @description
#' Generates a complete OpenFOAM wind case using
#' \code{buoyantBoussinesqPimpleFoam} - a transient, Boussinesq, finite-volume
#' Navier-Stokes solve with a kOmegaSST (URANS) closure - driven by an ABL
#' log-law inlet.
#'
#' The solver is a buoyant one running an isothermal problem: every surface
#' sits at \code{T_ref} and the walls are adiabatic, so the Boussinesq body
#' force is identically zero and the result is pure mechanical flow.
#'
#' Terrain and canopy are used when supplied and skipped when not.  Terrain
#' becomes solid geometry; canopy becomes distributed drag
#' via \code{atmPlantCanopy*} rather than solid blocks, because a crown passes
#' and drags air instead of blocking it.
#'
#' @section Oblique wind:
#' \code{blockMesh} emits four separately named lateral patches
#' (\code{xMin}, \code{xMax}, \code{yMin}, \code{yMax}) and each is assigned an
#' inlet / outlet / lateral role from the sign of
#' \code{dot(flowDir, outward_normal)}.  Any wind direction works without
#' rotating the domain - unlike the previous generator, which pinned the inlet
#' to the x-min face so that a north wind injected velocity through a face
#' whose normal was -x and effectively nothing entered.
#'
#' @param case_dir Character. Case directory; must already contain
#'   \code{constant/triSurface/<stl_file>}.
#' @param stl_file Character. Building STL path (host path).
#' @param domain Named list with xmin/xmax/ymin/ymax/zmin/zmax (metres, local
#'   coordinates), from \code{prepare_openfoam_inputs()$domain}.
#' @param inlet_velocity Numeric length-3 (Ux, Uy, Uz) in m/s at \code{z_ref}.
#'   Any horizontal direction is supported.
#' @param z_ref Numeric. Reference height for \code{inlet_velocity}. Default 10.
#' @param T_ref Numeric. Reference air temperature (K). Default 295.  Every
#'   surface is held at this value, which is what switches buoyancy off.
#' @param sim_hours Numeric. Physical hours to simulate.  \code{NULL} (the
#'   default) computes three flow-through times from the domain length and the
#'   inlet speed, which is what a wind case needs to flush its transient.
#' @param n_writes Integer. Output times over the run. Default 8.
#' @param terrain_stl Path to a terrain STL (see
#'   \code{\link{prepare_foam_geometry}}), or \code{NULL} for a flat floor.
#' @param terrain_dem Optional SpatRaster / path of the bare-earth DEM.  Used
#'   to set the domain floor.  Strongly recommended whenever
#'   \code{terrain_stl} is supplied.
#' @param canopy_stl Path to a canopy volume STL, or \code{NULL}.
#' @param leaf_area_density Numeric. LAD (1/m) inside the canopy. Default 0.4.
#' @param canopy_heat_source Numeric or \code{NULL}. Enables
#'   \code{atmPlantCanopyTSource} and writes \code{0/qPlant} with this uniform
#'   value.  \code{NULL} (default) leaves the source out entirely.
#'
#'   Off by default deliberately: \code{qPlant} is a canopy energy flux that
#'   cannot be derived from a canopy height model, and enabling the source with
#'   \code{qPlant = 0} would be inert while still adding a way for the run to
#'   abort.  Check the units against the \code{atmPlantCanopyTSource}
#'   documentation for your OpenFOAM version before relying on a value.
#' @param plant_cd Numeric. Canopy drag coefficient. Default 0.2.
#' @param z0 Numeric. Aerodynamic roughness length (m). Default 0.1.
#' @param base_cell_size Numeric. Background cell size (m). Default 10.
#' @param building_refinement Integer. snappyHexMesh level for buildings.
#'   Default 2.
#' @param terrain_refinement Integer. snappyHexMesh level for terrain.
#'   Default 0 (the background cell size already resolves gentle slope).
#' @param max_cells Integer. Global cell-count cap. Default 3e6.
#' @param max_co Numeric. Maximum Courant number. \code{NULL} (the default)
#'   gives 20: PIMPLE re-converges momentum and pressure within each step, so
#'   Courant well above 1 is stable when marching to a steady field.
#' @param n_outer_correctors Integer. PIMPLE outer correctors. Default 2.
#' @param overwrite,quiet Logical.
#'
#' @return Invisibly, a list with \code{case_dir}, \code{files} and
#'   \code{params}.
#'
#' @seealso \code{\link{prepare_openfoam_inputs}},
#'   \code{\link{prepare_foam_geometry}},
#'   \code{\link{read_foam_pedestrian_slice}}
#' @export
prepare_foam_case <- function(
    case_dir, stl_file, domain,
    inlet_velocity       = c(5, 0, 0),
    z_ref                = 10,
    T_ref                = 295,
    sim_hours            = NULL,
    n_writes             = 8L,
    terrain_stl          = NULL,
    terrain_dem          = NULL,
    canopy_stl           = NULL,
    leaf_area_density    = 0.4,
    canopy_heat_source   = NULL,
    plant_cd             = 0.2,
    z0                   = 0.1,
    base_cell_size       = 10,
    building_refinement  = 2L,
    terrain_refinement   = 0L,
    max_cells            = 3000000L,
    max_co               = NULL,
    n_outer_correctors   = 2L,
    overwrite            = FALSE,
    quiet                = FALSE
) {
  # This generator is wind-only. Fixed here rather than exposed as an argument,
  # and kept as a local because it is reported in the return value.
  mode <- "wind"

  if (missing(case_dir) || !nzchar(case_dir))
    stop("`case_dir` must be provided.", call. = FALSE)
  if (missing(stl_file) || !file.exists(stl_file))
    stop("`stl_file` not found: ", stl_file, call. = FALSE)
  if (missing(domain) || !all(c("xmin","xmax","ymin","ymax","zmin","zmax")
                              %in% names(domain)))
    stop("`domain` must be a list with xmin/xmax/ymin/ymax/zmin/zmax.",
         call. = FALSE)

  case_dir   <- normalizePath(case_dir, mustWork = FALSE)
  system_dir <- file.path(case_dir, "system")
  ic_dir     <- file.path(case_dir, "0")
  const_dir  <- file.path(case_dir, "constant")

  if (file.exists(file.path(system_dir, "controlDict")) && !overwrite)
    stop("Case files already exist in ", case_dir,
         ".\nUse `overwrite = TRUE` to replace them.", call. = FALSE)
  for (d in c(system_dir, ic_dir, const_dir))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  msg <- function(...) if (!isTRUE(quiet)) message(...)

  # -- Wind direction and patch roles --------------------------------------
  u_horiz <- as.numeric(inlet_velocity[1:2])
  u_mag   <- sqrt(sum(as.numeric(inlet_velocity)^2))
  if (sqrt(sum(u_horiz^2)) < 1e-9)
    stop("Horizontal component of `inlet_velocity` is zero; ",
         "there is no wind direction to impose.", call. = FALSE)
  fd <- u_horiz / sqrt(sum(u_horiz^2))

  roles <- fc_patch_roles(fd)
  # The top is a slip lid, as in the classic ABL setup.
  top_open <- FALSE

  # -- Terrain -------------------------------------------------------------
  dem_r <- NULL
  if (!is.null(terrain_dem)) {
    if (inherits(terrain_dem, "SpatRaster")) dem_r <- terrain_dem
    else if (is.character(terrain_dem) && file.exists(terrain_dem) &&
             requireNamespace("terra", quietly = TRUE))
      dem_r <- terra::rast(terrain_dem)
    else
      stop("`terrain_dem` was supplied but cannot be read.", call. = FALSE)
  }
  use_terrain <- !is.null(terrain_stl)
  if (use_terrain && !(is.character(terrain_stl) && file.exists(terrain_stl)))
    stop("`terrain_stl` was supplied but cannot be read: ", terrain_stl,
         call. = FALSE)

  # With terrain the blockMesh floor must sit BELOW the lowest ground so the
  # terrain solid is fully inside the background mesh; snappy then carves the
  # rock away and the flat `ground` patch ends up with zero faces (harmless -
  # its boundaryField entry is still written).
  if (use_terrain && !is.null(dem_r)) {
    zr <- range(terra::values(dem_r, mat = FALSE), na.rm = TRUE)
    domain$zmin <- zr[1] - max(20, 0.1 * diff(zr))
    domain$zmax <- max(domain$zmax, zr[2] + 150)
    msg(sprintf("Terrain %.1f - %.1f m; domain floor lowered to %.1f m, top %.1f m",
                zr[1], zr[2], domain$zmin, domain$zmax))
  }
  ground_patches  <- if (use_terrain)
    c("ground", tools::file_path_sans_ext(basename(terrain_stl))) else "ground"
  z_ground_abl <- if (!is.null(dem_r))
    min(terra::values(dem_r, mat = FALSE), na.rm = TRUE) else domain$zmin

  stl_name        <- basename(stl_file)
  buildings_patch <- tools::file_path_sans_ext(stl_name)
  walls           <- c(ground_patches, buildings_patch)

  # A pure wind case has no time-dependent forcing, so it only has to run long
  # enough to flush the initial transient out of the domain - about three
  # flow-through times. An arbitrary "1 hour" is unrelated to that and, on a
  # 2 km domain in a 2.4 m/s wind, was ~4x longer than needed.
  if (is.null(sim_hours)) {
    streamwise <- abs(fd[1]) * (domain$xmax - domain$xmin) +
                  abs(fd[2]) * (domain$ymax - domain$ymin)
    sim_hours  <- 3 * streamwise / max(u_mag, 0.1) / 3600
    msg(sprintf("Wind case: %.0f m streamwise / %.2f m/s = %.0f s per flow-through; running 3 (%.0f s)",
                streamwise, u_mag, streamwise / max(u_mag, 0.1), sim_hours * 3600))
  }
  if (sim_hours <= 0) stop("`sim_hours` must be positive.", call. = FALSE)

  # The case is only being marched to a steady field, so the limit is
  # stability, and PIMPLE's outer correctors handle Co well above 1. This is
  # the single biggest lever on wall-clock after cell size: dt scales with
  # maxCo, and the run was pinned at dt = 0.84 s by maxCo 5.
  if (is.null(max_co)) max_co <- 20

  end_time_s       <- sim_hours * 3600
  write_interval_s <- end_time_s / max(1L, as.integer(n_writes))

  # -- Turbulence initialisation -------------------------------------------
  kappa <- 0.41; Cmu <- 0.09
  u_star <- max(u_mag * kappa / log(max(z_ref, 1.001 * z0) / z0), 0.05)
  k_init   <- max(u_star^2 / sqrt(Cmu), 1e-4)
  eps_init <- max(u_star^3 / (kappa * max(z0, 0.01)), 1e-6)
  # omega = epsilon / (Cmu * k). The Cmu factor is NOT optional: omitting it
  # under-predicts omega by 1/0.09 ~ 11x and inflates the initial eddy
  # viscosity nut = k/omega by the same factor.
  omega_init <- eps_init / (Cmu * k_init)

  # -- Mesh sizing ---------------------------------------------------------
  nx <- max(1L, ceiling((domain$xmax - domain$xmin) / base_cell_size))
  ny <- max(1L, ceiling((domain$ymax - domain$ymin) / base_cell_size))
  nz <- max(1L, ceiling((domain$zmax - domain$zmin) / base_cell_size))
  base_cells       <- as.numeric(nx) * ny * nz
  max_global_cells <- min(ceiling(base_cells * 4), as.integer(max_cells))
  max_local_cells  <- ceiling(max_global_cells / 2)

  loc_x <- (domain$xmin + domain$xmax) / 2
  loc_y <- (domain$ymin + domain$ymax) / 2
  # Must be in the air: above the highest terrain and outside every building.
  loc_z <- domain$zmin + 0.9 * (domain$zmax - domain$zmin)

  msg(sprintf("Solver buoyantBoussinesqPimpleFoam (wind)  |  %.2f h = %g s",
              sim_hours, end_time_s))
  msg(sprintf("Wind: %.2f m/s at %g m, dir (%.3f %.3f)  |  patches: %s",
              u_mag, z_ref, fd[1], fd[2],
              paste(sprintf("%s=%s", names(roles), roles), collapse = " ")))
  msg(sprintf("Domain %.0f x %.0f x %.0f m  |  mesh %d x %d x %d (base %g m)",
              domain$xmax - domain$xmin, domain$ymax - domain$ymin,
              domain$zmax - domain$zmin, nx, ny, nz, base_cell_size))

  fc_check_domain(domain, stl_file, fd,
                  ground_top = if (!is.null(dem_r))
                    max(terra::values(dem_r, mat = FALSE), na.rm = TRUE) else NULL,
                  dem = dem_r, quiet = quiet)

  files <- list()

  # Stale boundaryData from a case generated before this package became
  # wind-only: OpenFOAM would still read it if a boundary condition referred
  # to it, so it does not belong in a wind case directory.
  bd_root <- file.path(case_dir, "constant", "boundaryData")
  if (dir.exists(bd_root)) {
    unlink(bd_root, recursive = TRUE)
    msg("  Removed stale constant/boundaryData")
  }

  use_canopy <- !is.null(canopy_stl)
  if (use_canopy && !(is.character(canopy_stl) && file.exists(canopy_stl)))
    stop("`canopy_stl` was supplied but cannot be read: ", canopy_stl,
         call. = FALSE)

  # -- Write ---------------------------------------------------------------
  msg("Writing Allrun ...")
  files$allrun <- noc_write(make_noc_allrun(set_fields = use_canopy),
                            file.path(case_dir, "Allrun"))
  Sys.chmod(files$allrun, mode = "0755")

  msg("Writing system/ ...")
  files$decompose_par_dict <- noc_write(make_decompose_par_dict(1L),
    file.path(system_dir, "decomposeParDict"))
  files$block_mesh_dict <- noc_write(
    make_noc_block_mesh_dict(domain, nx, ny, nz),
    file.path(system_dir, "blockMeshDict"))
  files$snappy_hex_mesh_dict <- noc_write(
    make_noc_snappy_dict(stl_name, loc_x, loc_y, loc_z,
                         as.integer(building_refinement),
                         terrain_stl = if (use_terrain) terrain_stl else NULL,
                         terrain_ref = as.integer(terrain_refinement),
                         max_local_cells = max_local_cells,
                         max_global_cells = max_global_cells),
    file.path(system_dir, "snappyHexMeshDict"))
  files$control_dict <- noc_write(
    make_noc_control_dict(end_time_s, write_interval_s, max_co = max_co,
                          terrain_stl = if (use_terrain) terrain_stl else NULL,
                          z_floor     = domain$zmin),
    file.path(system_dir, "controlDict"))
  files$fv_schemes  <- noc_write(make_noc_fv_schemes(),
    file.path(system_dir, "fvSchemes"))
  files$fv_solution <- noc_write(
    make_noc_fv_solution(n_outer_correctors = n_outer_correctors),
    file.path(system_dir, "fvSolution"))
  if (use_canopy)
    files$set_fields_dict <- noc_write(
      make_foam_set_fields_dict(
        file.path("constant", "triSurface", basename(canopy_stl)),
        outside_point = c(loc_x, loc_y, loc_z),
        lad = leaf_area_density, plant_cd = plant_cd),
      file.path(system_dir, "setFieldsDict"))

  msg("Writing 0/ ...")
  files$U <- noc_write(
    make_noc_u_field(roles, walls, u_mag, fd, z_ref, z0, z_ground_abl, top_open),
    file.path(ic_dir, "U"))
  files$p_rgh <- noc_write(make_noc_p_rgh_field(roles, walls, top_open),
    file.path(ic_dir, "p_rgh"))
  files$T <- noc_write(
    make_noc_t_field(roles, ground_patches, buildings_patch, top_open, T_ref),
    file.path(ic_dir, "T"))
  files$k <- noc_write(
    make_noc_turb_scalar_field("k", "[ 0 2 -2 0 0 0 0 ]", k_init, roles, walls,
      "kqRWallFunction", top_open,
      abl_type = "atmBoundaryLayerInletK",
      u_mag, fd, z_ref, z0, z_ground_abl),
    file.path(ic_dir, "k"))
  files$omega <- noc_write(
    make_noc_turb_scalar_field("omega", "[ 0 0 -1 0 0 0 0 ]", omega_init,
      roles, walls, "omegaWallFunction", top_open,
      abl_type = "atmBoundaryLayerInletOmega",
      u_mag, fd, z_ref, z0, z_ground_abl),
    file.path(ic_dir, "omega"))
  files$nut <- noc_write(
    make_noc_nut_field(roles, ground_patches, buildings_patch, z0, top_open),
    file.path(ic_dir, "nut"))
  files$alphat <- noc_write(make_noc_alphat_field(roles, walls, top_open),
    file.path(ic_dir, "alphat"))
  if (use_canopy) {
    # atmPlantCanopy* renamed BOTH of its input fields between releases:
    #
    #   field              v2506     v2606 (documented)
    #   leaf area density  LAD       leafAreaDensity
    #   drag coefficient   Cd        plantCd
    #
    # All four are written, identical in pairs. The solver reads only the two
    # it wants; the others are inert volScalarFields. Getting a name wrong is
    # fatal and late - the case meshes, decomposes, constructs the turbulence
    # model and activates the fvOption, and only then every rank aborts with
    #     cannot find file "/case/processorN/0/LAD"
    # so there is no cheap way to discover it except by running.
    for (nm in c("LAD", "leafAreaDensity"))
      files[[nm]] <- noc_write(
        make_foam_canopy_field(nm, "[ 0 -1 0 0 0 0 0 ]"),
        file.path(ic_dir, nm))
    for (nm in c("Cd", "plantCd"))
      files[[nm]] <- noc_write(
        make_foam_canopy_field(nm, "[ 0 0 0 0 0 0 0 ]"),
        file.path(ic_dir, nm))
    # atmPlantCanopyTSource reads qPlant. Written only when that source is
    # actually enabled - otherwise every rank aborts with
    #     cannot find file "/case/processorN/0/qPlant"
    # after a fully successful startup.
    if (!is.null(canopy_heat_source))
      files$qPlant <- noc_write(
        make_foam_canopy_field("qPlant", "[ 1 -1 -3 0 0 0 0 ]",
                               value = as.numeric(canopy_heat_source)),
        file.path(ic_dir, "qPlant"))
  }

  # Remove stale fields left by a previous generator.  This is not tidiness:
  # OpenFOAM reads EVERY file in 0/, so a leftover 0/epsilon or 0/p from an
  # older simpleFoam case still carries the old patch names and kills the run
  # at decomposePar with
  #     "Cannot find patchField entry for xMin ... file: 0/epsilon"
  # long after blockMesh and snappyHexMesh have reported success.
  written_fields <- c("U", "p_rgh", "T", "k", "omega", "nut", "alphat",
                      if (use_canopy) c("LAD", "leafAreaDensity",
                                        "Cd", "plantCd",
                                        if (!is.null(canopy_heat_source)) "qPlant"))
  present <- list.files(ic_dir, full.names = FALSE)
  present <- present[!dir.exists(file.path(ic_dir, present))]
  stale   <- setdiff(present, written_fields)
  if (length(stale)) {
    file.remove(file.path(ic_dir, stale))
    msg("  Removed stale field(s) from 0/: ", paste(stale, collapse = ", "),
        " (left by a previous case generator)")
  }

  # Keep a pristine copy. setFields rewrites 0/ in place - uniform fields
  # become nonuniform lists sized to the snappyHexMesh cell count - so without
  # this the case runs exactly once. The second run rebuilds a smaller
  # blockMesh-only mesh and decomposePar dies with
  #     Size 1733307 is not equal to the expected length 1170585
  # Allrun restores 0/ from 0.orig/ before every run.
  orig_dir <- file.path(case_dir, "0.orig")
  unlink(orig_dir, recursive = TRUE)
  dir.create(orig_dir, showWarnings = FALSE)
  file.copy(list.files(ic_dir, full.names = TRUE), orig_dir, overwrite = TRUE)
  files$zero_orig <- orig_dir

  msg("Writing constant/ ...")
  files$g <- noc_write(make_gravity(), file.path(const_dir, "g"))
  files$transport_props <- noc_write(make_noc_transport_properties(T_ref),
    file.path(const_dir, "transportProperties"))
  files$turbulence_props <- noc_write(make_noc_turbulence_properties(),
    file.path(const_dir, "turbulenceProperties"))
  # fvOptions exists only for the canopy now. It MUST be removed when there is
  # no canopy: a stale one from an earlier canopy run still declares
  # atmPlantCanopyUSource, which reads 0/Cd and 0/LAD that are no longer
  # written, and every rank aborts with
  #     cannot find file "/case/processorN/0/Cd"
  # long after meshing and decomposition have reported success.
  fv_options_path <- file.path(const_dir, "fvOptions")
  if (use_canopy) {
    files$fv_options <- noc_write(
      make_noc_fv_options(canopy_thermal = !is.null(canopy_heat_source)),
      fv_options_path)
  } else if (file.exists(fv_options_path)) {
    file.remove(fv_options_path)
    msg("  Removed stale constant/fvOptions (no canopy in this case)")
  }

  msg("Case files written to: ", case_dir)
  msg(sprintf("Next step: run_openfoam_docker(\"%s\")", case_dir))

  invisible(list(
    case_dir = case_dir, files = files,
    params = list(
      solver = "buoyantBoussinesqPimpleFoam", mode = mode,
      patch_roles = roles, flow_dir = fd, u_mag = u_mag, z_ref = z_ref,
      T_ref = T_ref, sim_hours = sim_hours,
      end_time_s = end_time_s, write_interval_s = write_interval_s,
      domain = domain,
      terrain = use_terrain, canopy = use_canopy,
      ground_patches = ground_patches,
      z0 = z0, base_cell_size = base_cell_size,
      k_init = k_init, omega_init = omega_init, max_co = max_co)))
}
