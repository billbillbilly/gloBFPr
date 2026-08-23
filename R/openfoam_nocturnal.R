# ===========================================================================
# Nocturnal cold-air drainage OpenFOAM case generator
# ===========================================================================
# Physics: buoyantBoussinesqSimpleFoam — steady-state, Boussinesq buoyancy,
#          energy equation.  No forced inlet; flow is driven purely by
#          surface-temperature differentials (vegetation < roads < buildings).
#
# Outputs at pedestrian level (z = 1.5 m):
#   - T   — air temperature (differentials reveal heat islands / cool pools)
#   - |U| — wind speed driven by density gradients
#   - T_cool * |U| — cool-air transport flux
#   - U vectors — flow direction
# ===========================================================================

# ---------------------------------------------------------------------------
# Internal helpers — one per OpenFOAM file
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

# -- constant/g ---------------------------------------------------------------
#' @noRd
make_gravity <- function() {
  paste0(
    noc_foam_header("uniformDimensionedVectorField", "g"),
    "dimensions  [ 0 1 -2 0 0 0 0 ];\n",
    "value       ( 0 0 -9.81 );\n"
  )
}

# -- constant/transportProperties (Boussinesq) --------------------------------
#' @noRd
make_noc_transport_properties <- function(T_ref, beta = NULL) {
  if (is.null(beta)) beta <- 1 / T_ref   # ideal-gas approximation
  paste0(
    noc_foam_header("dictionary", "transportProperties"),
    "transportModel  Newtonian;\n\n",
    "// Kinematic viscosity of air at ~20 C\n",
    "nu              nu [ 0 2 -1 0 0 0 0 ] 1.5e-05;\n\n",
    "// Boussinesq parameters\n",
    sprintf("beta            beta [ 0 0 0 -1 0 0 0 ] %.6g;\n", beta),
    sprintf("TRef            TRef [ 0 0 0  1 0 0 0 ] %g;\n", T_ref),
    "Pr              Pr   [ 0  0  0 0 0 0 0 ] 0.71;\n",
    "Prt             Prt  [ 0  0  0 0 0 0 0 ] 0.9;\n"
  )
}

# -- constant/turbulenceProperties --------------------------------------------
#' @noRd
make_noc_turbulence_properties <- function() {
  paste0(
    noc_foam_header("dictionary", "turbulenceProperties"),
    "simulationType  RAS;\n\n",
    "RAS\n{\n",
    # kOmegaSST is used instead of kEpsilon because:
    #  - It has a built-in production limiter (Pk <= 10*Cmu*k*omega) that
    #    prevents the k blow-up seen with kEpsilon in buoyancy-driven cases.
    #  - omegaWallFunction provides analytical near-wall omega without needing
    #    epsilonWallFunction, which is numerically fragile at low velocities.
    "    RASModel        kOmegaSST;\n",
    "    turbulence      on;\n",
    "    printCoeffs     on;\n",
    "}\n"
  )
}

# -- system/controlDict -------------------------------------------------------
#' @noRd
make_noc_control_dict <- function(n_iterations, write_interval, z_slice = 1.5) {
  paste0(
    noc_foam_header("dictionary", "controlDict"),
    "application     buoyantBoussinesqSimpleFoam;\n\n",
    "startFrom       startTime;\n",
    "startTime       0;\n",
    "stopAt          endTime;\n",
    sprintf("endTime         %d;\n", as.integer(n_iterations)),
    "deltaT          1;\n\n",
    "writeControl    timeStep;\n",
    sprintf("writeInterval   %d;\n", as.integer(write_interval)),
    "purgeWrite      2;\n\n",
    "writeFormat     ascii;\n",
    "writePrecision  8;\n",
    "writeCompression off;\n\n",
    "timeFormat      general;\n",
    "timePrecision   6;\n\n",
    "runTimeModifiable true;\n\n",
    # Pedestrian-level slice: writes raw x-y-z-value tables under
    # postProcessing/pedestrianSlice/<time>/{T,U}_z1p5m.raw  (v2506+)
    # or z1p5m_{T,U}.raw (pre-v2506). read_foam_pedestrian_slice() tries both.
    # These are read back into R by read_foam_pedestrian_slice()
    "functions\n{\n",
    "    pedestrianSlice\n    {\n",
    "        type            surfaces;\n",
    "        libs            (\"libsampling.so\");\n\n",
    "        executeControl  writeTime;\n",
    "        writeControl    writeTime;\n\n",
    "        surfaceFormat   raw;\n",
    "        fields          ( T U );\n\n",
    "        surfaces\n        {\n",
    "            z1p5m\n            {\n",
    "                type        cuttingPlane;\n",
    "                planeType   pointAndNormal;\n",
    "                pointAndNormalDict\n                {\n",
    sprintf("                    point  (0 0 %g);\n", z_slice),
    "                    normal (0 0 1);\n",
    "                }\n",
    "                interpolate true;\n",
    "            }\n",
    "        }\n",
    "    }\n",
    "}\n"
  )
}

# -- system/fvSchemes ---------------------------------------------------------
#' @noRd
make_noc_fv_schemes <- function() {
  paste0(
    noc_foam_header("dictionary", "fvSchemes"),
    "ddtSchemes\n{\n    default         steadyState;\n}\n\n",
    "gradSchemes\n{\n",
    "    default         Gauss linear;\n",
    "    grad(U)         cellLimited Gauss linear 1;\n",
    "    grad(T)         cellLimited Gauss linear 1;\n",
    "}\n\n",
    "divSchemes\n{\n",
    "    default                             none;\n",
    # bounded: prevents U from going non-monotone during early buoyancy transient
    "    div(phi,U)                          bounded Gauss linearUpwind grad(U);\n",
    "    div(phi,T)                          bounded Gauss linearUpwind default;\n",
    # upwind for k/epsilon: eliminates oscillation-driven production blow-up
    # (linearUpwind allows local overshoots that multiply with nut → runaway)
    "    div(phi,k)                          bounded Gauss upwind;\n",
    "    div(phi,omega)                      bounded Gauss upwind;\n",
    "    div((nuEff*dev2(T(grad(U)))))       Gauss linear;\n",
    "}\n\n",
    "laplacianSchemes\n{\n    default         Gauss linear corrected;\n}\n\n",
    "interpolationSchemes\n{\n    default         linear;\n}\n\n",
    "snGradSchemes\n{\n    default         corrected;\n}\n\n",
    # kOmegaSST uses wall distance for the F1/F2 blending functions.
    # OpenFOAM v2506 removed the default and requires 'method' to be explicit.
    "wallDist\n{\n    method          meshWave;\n}\n"
  )
}

# -- system/fvSolution --------------------------------------------------------
#' @noRd
make_noc_fv_solution <- function() {
  paste0(
    noc_foam_header("dictionary", "fvSolution"),
    "solvers\n{\n",
    # GAMG with GaussSeidel (NOT DICGaussSeidel) — DIC computes 1/sqrt(diag)
    # and crashes with SIGFPE if the coarsened parallel matrix has a zero or
    # near-zero diagonal, which is common in open-boundary buoyancy cases.
    # GaussSeidel has no such reciprocal step and is the standard smoother
    # used in all OpenFOAM buoyantBoussinesqSimpleFoam tutorials.
    #
    # CRITICAL: relTol must be 0 for p_rgh.
    # relTol 0.01 stops the solver when residual drops by 99% of the initial
    # value.  If the initial residual is large (e.g. 0.88 when buoyancy kicks
    # in), the solver stops at 0.88*0.01 = 0.0088 — far above the absolute
    # tolerance, leaving a massive divergence that makes the next step worse.
    # relTol 0 forces the solver to always reach the absolute tolerance.
    "    p_rgh\n    {\n",
    "        solver          GAMG;\n",
    "        agglomerator    faceAreaPair;\n",
    "        mergeLevels     1;\n",
    "        nCellsInCoarsestLevel 200;\n",
    "        smoother        GaussSeidel;\n",
    "        tolerance       1e-6;\n",
    "        relTol          0;\n",
    "        maxIter         500;\n",
    "    }\n\n",
    "    p_rghFinal\n    {\n",
    "        solver          GAMG;\n",
    "        agglomerator    faceAreaPair;\n",
    "        mergeLevels     1;\n",
    "        nCellsInCoarsestLevel 200;\n",
    "        smoother        GaussSeidel;\n",
    "        tolerance       1e-7;\n",
    "        relTol          0;\n",
    "        maxIter         500;\n",
    "    }\n\n",
    "    \"(U|T|k|omega)\"\n    {\n",
    "        solver          smoothSolver;\n",
    "        smoother        symGaussSeidel;\n",
    "        tolerance       1e-7;\n",
    "        relTol          0.1;\n",
    "    }\n}\n\n",
    "SIMPLE\n{\n",
    "    nNonOrthogonalCorrectors 1;\n",
    "    pRefCell    0;\n",
    "    pRefValue   0;\n",
    "    residualControl\n    {\n",
    "        p_rgh           1e-4;\n",
    "        U               1e-4;\n",
    "        T               1e-4;\n",
    "        \"(k|omega)\"     1e-4;\n",
    "    }\n}\n\n",
    # Relaxation factors tuned for buoyancy-driven SIMPLE stability.
    # Key insight: T is the primary instability driver.  Once T converges,
    # the resulting buoyancy force creates a large U correction, which
    # requires a large p_rgh correction. If any of these corrections are
    # too large per iteration, the chain explodes.
    # Reference: OpenFOAM hotRoom tutorial (buoyantBoussinesqSimpleFoam):
    #   p_rgh=0.7, rho=1.0, U=0.2, T=0.2
    "relaxationFactors\n{\n",
    "    fields\n    {\n",
    "        p_rgh   0.3;\n",  # increase from 0.2: too-low p relaxation starves divergence correction
    "        rho     1.0;\n",  # no rho under-relaxation: rho should track T immediately
    "    }\n",
    "    equations\n    {\n",
    "        U       0.2;\n",  # reduce from 0.4: slow down velocity updates
    "        T       0.2;\n",  # reduce from 0.5: T is the buoyancy driver — must be conservative
    "        k       0.3;\n",
    "        omega   0.3;\n",
    "    }\n}\n"
  )
}

# -- 0/U (no inlet velocity — purely buoyancy driven) -------------------------
#' @noRd
make_noc_u_field <- function(patch_name) {
  paste0(
    noc_foam_header("volVectorField", "U"),
    "dimensions      [ 0 1 -1 0 0 0 0 ];\n\n",
    "internalField   uniform (0 0 0);\n\n",
    "boundaryField\n{\n",
    # All lateral faces: allow flow driven by buoyancy to exit/enter freely
    "    inlet\n    {\n",
    "        type            pressureInletOutletVelocity;\n",
    "        value           uniform (0 0 0);\n",
    "    }\n",
    "    outlet\n    {\n",
    "        type            pressureInletOutletVelocity;\n",
    "        value           uniform (0 0 0);\n",
    "    }\n",
    "    ground\n    {\n        type            noSlip;\n    }\n",
    "    ", patch_name, "\n    {\n        type            noSlip;\n    }\n",
    "    top\n    {\n",
    "        type            pressureInletOutletVelocity;\n",
    "        value           uniform (0 0 0);\n",
    "    }\n",
    "    sides\n    {\n",
    "        type            slip;\n",
    "    }\n",
    "}\n"
  )
}

# -- 0/p_rgh (pressure minus hydrostatic: p - rho*g*h) ------------------------
#' @noRd
make_noc_p_rgh_field <- function(patch_name) {
  paste0(
    noc_foam_header("volScalarField", "p_rgh"),
    "dimensions      [ 0 2 -2 0 0 0 0 ];\n\n",
    "internalField   uniform 0;\n\n",
    # Open faces (inlet/outlet/top): fixedValue 0 provides the pressure
    # reference and is the correct pair for pressureInletOutletVelocity on U.
    # fixedFluxPressure on open faces is WRONG — it adjusts the pressure
    # gradient to match U, but pressureInletOutletVelocity derives U FROM the
    # pressure gradient: circular → continuity explosion.
    # fixedFluxPressure is only correct on no-slip walls and slip faces where
    # the normal flux is fixed by U (zero), not derived from pressure.
    "boundaryField\n{\n",
    "    inlet\n    {\n        type            fixedValue;\n        value           uniform 0;\n    }\n",
    "    outlet\n    {\n        type            fixedValue;\n        value           uniform 0;\n    }\n",
    "    ground\n    {\n        type            fixedFluxPressure;\n        value           uniform 0;\n    }\n",
    "    ", patch_name, "\n    {\n        type            fixedFluxPressure;\n        value           uniform 0;\n    }\n",
    "    top\n    {\n        type            fixedValue;\n        value           uniform 0;\n    }\n",
    "    sides\n    {\n        type            fixedFluxPressure;\n        value           uniform 0;\n    }\n",
    "}\n"
  )
}

# -- 0/T (temperature field) ---------------------------------------------------
#' @noRd
make_noc_t_field <- function(T_ref, T_ground, T_buildings, patch_name) {
  paste0(
    noc_foam_header("volScalarField", "T"),
    "dimensions      [ 0 0 0 1 0 0 0 ];\n\n",
    sprintf("internalField   uniform %g;\n\n", T_ref),
    "boundaryField\n{\n",
    # Lateral faces: ambient temperature (inletOutlet allows outflow)
    "    inlet\n    {\n",
    "        type            inletOutlet;\n",
    sprintf("        inletValue      uniform %g;\n", T_ref),
    sprintf("        value           uniform %g;\n", T_ref),
    "    }\n",
    "    outlet\n    {\n",
    "        type            inletOutlet;\n",
    sprintf("        inletValue      uniform %g;\n", T_ref),
    sprintf("        value           uniform %g;\n", T_ref),
    "    }\n",
    # Ground: cooled surface (land-use weighted mean)
    "    ground\n    {\n",
    "        type            fixedValue;\n",
    sprintf("        value           uniform %g;\n", T_ground),
    "    }\n",
    # Buildings: slower cooling than vegetation, faster than water
    "    ", patch_name, "\n    {\n",
    "        type            fixedValue;\n",
    sprintf("        value           uniform %g;\n", T_buildings),
    "    }\n",
    # Top: ambient air temperature
    "    top\n    {\n",
    "        type            inletOutlet;\n",
    sprintf("        inletValue      uniform %g;\n", T_ref),
    sprintf("        value           uniform %g;\n", T_ref),
    "    }\n",
    "    sides\n    {\n        type            zeroGradient;\n    }\n",
    "}\n"
  )
}

# -- 0/k ----------------------------------------------------------------------
#' @noRd
make_noc_k_field <- function(k_init, patch_name) {
  # Open boundaries (inlet/outlet/top) use inletOutlet so that when
  # buoyancy-driven reverse flow re-enters the domain it carries the safe
  # reference value k_init rather than the unbounded extrapolated value that
  # zeroGradient would allow (which causes exponential blow-up).
  paste0(
    noc_foam_header("volScalarField", "k"),
    "dimensions      [ 0 2 -2 0 0 0 0 ];\n\n",
    sprintf("internalField   uniform %g;\n\n", k_init),
    "boundaryField\n{\n",
    "    inlet\n    {\n",
    "        type            inletOutlet;\n",
    sprintf("        inletValue      uniform %g;\n", k_init),
    sprintf("        value           uniform %g;\n", k_init),
    "    }\n",
    "    outlet\n    {\n",
    "        type            inletOutlet;\n",
    sprintf("        inletValue      uniform %g;\n", k_init),
    sprintf("        value           uniform %g;\n", k_init),
    "    }\n",
    "    ground\n    {\n",
    "        type            kqRWallFunction;\n",
    sprintf("        value           uniform %g;\n", k_init),
    "    }\n",
    "    ", patch_name, "\n    {\n",
    "        type            kqRWallFunction;\n",
    sprintf("        value           uniform %g;\n", k_init),
    "    }\n",
    "    top\n    {\n",
    "        type            inletOutlet;\n",
    sprintf("        inletValue      uniform %g;\n", k_init),
    sprintf("        value           uniform %g;\n", k_init),
    "    }\n",
    "    sides\n    {\n        type            zeroGradient;\n    }\n",
    "}\n"
  )
}

# -- 0/epsilon ----------------------------------------------------------------
#' @noRd
make_noc_omega_field <- function(omega_init, patch_name) {
  paste0(
    noc_foam_header("volScalarField", "omega"),
    "dimensions      [ 0 0 -1 0 0 0 0 ];\n\n",
    sprintf("internalField   uniform %g;\n\n", omega_init),
    "boundaryField\n{\n",
    "    inlet\n    {\n",
    "        type            inletOutlet;\n",
    sprintf("        inletValue      uniform %g;\n", omega_init),
    sprintf("        value           uniform %g;\n", omega_init),
    "    }\n",
    "    outlet\n    {\n",
    "        type            inletOutlet;\n",
    sprintf("        inletValue      uniform %g;\n", omega_init),
    sprintf("        value           uniform %g;\n", omega_init),
    "    }\n",
    # omegaWallFunction provides the analytical near-wall omega = 6*nu/(beta1*y^2)
    # without needing epsilon — numerically much more robust at low velocities.
    "    ground\n    {\n",
    "        type            omegaWallFunction;\n",
    sprintf("        value           uniform %g;\n", omega_init),
    "    }\n",
    "    ", patch_name, "\n    {\n",
    "        type            omegaWallFunction;\n",
    sprintf("        value           uniform %g;\n", omega_init),
    "    }\n",
    "    top\n    {\n",
    "        type            inletOutlet;\n",
    sprintf("        inletValue      uniform %g;\n", omega_init),
    sprintf("        value           uniform %g;\n", omega_init),
    "    }\n",
    "    sides\n    {\n        type            zeroGradient;\n    }\n",
    "}\n"
  )
}

# -- 0/nut --------------------------------------------------------------------
#' @noRd
make_noc_nut_field <- function(z0, patch_name) {
  ks <- z0 * 20
  paste0(
    noc_foam_header("volScalarField", "nut"),
    "dimensions      [ 0 2 -1 0 0 0 0 ];\n\n",
    "internalField   uniform 0;\n\n",
    "boundaryField\n{\n",
    "    inlet\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "    outlet\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "    ground\n    {\n",
    "        type            nutkRoughWallFunction;\n",
    sprintf("        Ks              uniform %g;\n", ks),
    "        Cs              uniform 0.5;\n",  # v2506: Cs is Field<scalar>
    "        value           uniform 0;\n",
    "    }\n",
    "    ", patch_name, "\n    {\n",
    "        type            nutkWallFunction;\n",
    "        value           uniform 0;\n",
    "    }\n",
    "    top\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "    sides\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "}\n"
  )
}

# -- 0/alphat (turbulent thermal diffusivity) ----------------------------------
#' @noRd
make_noc_alphat_field <- function(patch_name) {
  paste0(
    noc_foam_header("volScalarField", "alphat"),
    "dimensions      [ 0 2 -1 0 0 0 0 ];\n\n",
    "internalField   uniform 0;\n\n",
    "boundaryField\n{\n",
    "    inlet\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "    outlet\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "    ground\n    {\n",
    "        type            alphatJayatillekeWallFunction;\n",
    "        Prt             0.9;\n",
    "        value           uniform 0;\n",
    "    }\n",
    "    ", patch_name, "\n    {\n",
    "        type            alphatJayatillekeWallFunction;\n",
    "        Prt             0.9;\n",
    "        value           uniform 0;\n",
    "    }\n",
    "    top\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "    sides\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "}\n"
  )
}

# -- system/blockMeshDict (same geometry helpers as wind case) ----------------
#' @noRd
make_noc_block_mesh_dict <- function(domain, nx, ny, nz) {
  d <- domain
  paste0(
    noc_foam_header("dictionary", "blockMeshDict"),
    "scale 1;\n\n",
    "vertices\n(\n",
    sprintf("    (%g %g %g)\n", d$xmin, d$ymin, d$zmin),  # 0
    sprintf("    (%g %g %g)\n", d$xmax, d$ymin, d$zmin),  # 1
    sprintf("    (%g %g %g)\n", d$xmax, d$ymax, d$zmin),  # 2
    sprintf("    (%g %g %g)\n", d$xmin, d$ymax, d$zmin),  # 3
    sprintf("    (%g %g %g)\n", d$xmin, d$ymin, d$zmax),  # 4
    sprintf("    (%g %g %g)\n", d$xmax, d$ymin, d$zmax),  # 5
    sprintf("    (%g %g %g)\n", d$xmax, d$ymax, d$zmax),  # 6
    sprintf("    (%g %g %g)\n", d$xmin, d$ymax, d$zmax),  # 7
    ");\n\n",
    "blocks\n(\n",
    sprintf("    hex (0 1 2 3 4 5 6 7) (%d %d %d) simpleGrading (1 1 1)\n", nx, ny, nz),
    ");\n\n",
    "edges\n();\n\n",
    "boundary\n(\n",
    "    inlet\n    {\n        type patch;\n        faces ((0 4 7 3));\n    }\n",
    "    outlet\n    {\n        type patch;\n        faces ((1 2 6 5));\n    }\n",
    "    ground\n    {\n        type wall;\n        faces ((3 2 1 0));\n    }\n",
    "    top\n    {\n        type patch;\n        faces ((4 5 6 7));\n    }\n",
    "    sides\n    {\n        type patch;\n        faces\n        (\n",
    "            (0 1 5 4)\n",
    "            (3 7 6 2)\n",
    "        );\n    }\n",
    ");\n\n",
    "mergePatchPairs\n();\n"
  )
}

# -- system/snappyHexMeshDict (same as wind case) ----------------------------
#' @noRd
make_noc_snappy_dict <- function(stl_name, loc_x, loc_y, loc_z, ref,
                                  max_local_cells  = 1000000L,
                                  max_global_cells = 4000000L) {
  patch_name <- tools::file_path_sans_ext(stl_name)
  paste0(
    noc_foam_header("dictionary", "snappyHexMeshDict"),
    "castellatedMesh true;\n",
    "snap            true;\n",
    "addLayers       false;\n\n",
    "geometry\n{\n",
    "    ", stl_name, "\n    {\n",
    "        type  triSurfaceMesh;\n",
    "        name  ", patch_name, ";\n",
    "    }\n}\n\n",
    "castellatedMeshControls\n{\n",
    sprintf("    maxLocalCells           %d;\n", as.integer(max_local_cells)),
    sprintf("    maxGlobalCells         %d;\n", as.integer(max_global_cells)),
    "    minRefinementCells          10;\n",
    "    maxLoadUnbalance          0.10;\n",
    "    nCellsBetweenLevels           3;\n\n",
    "    features ();\n\n",
    "    refinementSurfaces\n    {\n",
    "        ", patch_name, "\n        {\n",
    sprintf("            level (%d %d);\n", ref, ref),
    "        }\n    }\n\n",
    "    resolveFeatureAngle 30;\n\n",
    "    refinementRegions {}\n\n",
    sprintf("    locationInMesh (%g %g %g);\n", loc_x, loc_y, loc_z),
    "    allowFreeStandingZoneFaces true;\n}\n\n",
    "snapControls\n{\n",
    "    nSmoothPatch              3;\n",
    "    tolerance                 2.0;\n",
    "    nSolveIter               30;\n",
    "    nRelaxIter                5;\n",
    "    nFeatureSnapIter         10;\n",
    "    implicitFeatureSnap   false;\n",
    "    explicitFeatureSnap   false;\n",
    "    multiRegionFeatureSnap false;\n}\n\n",
    "addLayersControls\n{\n",
    "    relativeSizes           true;\n",
    "    layers {}\n",
    "    expansionRatio          1.3;\n",
    "    finalLayerThickness     0.3;\n",
    "    minThickness            0.1;\n",
    "    nGrow                     0;\n",
    "    featureAngle             60;\n",
    "    nRelaxIter                3;\n",
    "    nSmoothSurfaceNormals     1;\n",
    "    nSmoothNormals            3;\n",
    "    nSmoothThickness         10;\n",
    "    maxFaceThicknessRatio   0.5;\n",
    "    maxThicknessToMedialRatio 0.3;\n",
    "    minMedialAxisAngle       90;\n",
    "    nBufferCellsNoExtrude     0;\n",
    "    nLayerIter               50;\n}\n\n",
    "meshQualityControls\n{\n",
    "    maxNonOrtho              70;\n",
    "    maxBoundarySkewness      20;\n",
    "    maxInternalSkewness       4;\n",
    "    maxConcave               80;\n",
    "    minVol                 1e-13;\n",
    "    minTetQuality           1e-9;\n",
    "    minArea                  -1;\n",
    "    minTwist               0.02;\n",
    "    minDeterminant          0.001;\n",
    "    minFaceWeight           0.05;\n",
    "    minVolRatio             0.01;\n",
    "    minTriangleTwist         -1;\n",
    "    nSmoothScale              4;\n",
    "    errorReduction         0.75;\n}\n\n",
    "debug 0;\n",
    "mergeTolerance 1e-6;\n"
  )
}

# -- Allrun -------------------------------------------------------------------
#' @noRd
make_noc_allrun <- function() {
  paste0(
    "#!/bin/sh\n",
    "# Allrun — nocturnal cold-air drainage case\n",
    "# Generated by gloBFPr::prepare_nocturnal_case()\n",
    "# Pass NPROC env variable (via docker -e NPROC=N) to enable MPI parallel.\n",
    "cd \"${0%/*}\" || exit 1\n\n",
    "NPROC=${NPROC:-1}\n\n",
    "rm -rf constant/polyMesh processor* log.*\n",
    "find . -maxdepth 1 -name '[0-9]*' ! -name '0' -exec rm -rf {} +\n\n",
    "echo '[1/3] blockMesh ...'\n",
    "blockMesh > log.blockMesh 2>&1 \\\n",
    "  && echo '      OK' \\\n",
    "  || { echo '      FAILED — see log.blockMesh'; exit 1; }\n\n",
    "echo '[2/3] snappyHexMesh ...'\n",
    "if [ \"$NPROC\" -gt 1 ]; then\n",
    "    sed -i \"s/numberOfSubdomains.*/numberOfSubdomains $NPROC;/\" system/decomposeParDict\n",
    "    decomposePar -force > log.decomposeMesh 2>&1 \\\n",
    "      || { echo '      FAILED — see log.decomposeMesh'; exit 1; }\n",
    "    mpirun -np \"$NPROC\" --allow-run-as-root snappyHexMesh -overwrite -parallel > log.snappyHexMesh 2>&1 \\\n",
    "      || { echo '      FAILED — see log.snappyHexMesh'; exit 1; }\n",
    "    reconstructParMesh -constant -mergeTol 1e-6 > log.reconstructMesh 2>&1 \\\n",
    "      || { echo '      FAILED — see log.reconstructMesh'; exit 1; }\n",
    "    rm -rf processor*\n",
    "    echo '      OK (parallel)'\n",
    "else\n",
    "    snappyHexMesh -overwrite > log.snappyHexMesh 2>&1 \\\n",
    "      && echo '      OK' \\\n",
    "      || { echo '      FAILED — see log.snappyHexMesh'; exit 1; }\n",
    "fi\n\n",
    "echo '[3/3] buoyantBoussinesqSimpleFoam ...'\n",
    "if [ \"$NPROC\" -gt 1 ]; then\n",
    "    sed -i \"s/numberOfSubdomains.*/numberOfSubdomains $NPROC;/\" system/decomposeParDict\n",
    "    decomposePar -force > log.decompose 2>&1 \\\n",
    "      || { echo '      FAILED — see log.decompose'; exit 1; }\n",
    "    mpirun -np \"$NPROC\" --allow-run-as-root buoyantBoussinesqSimpleFoam -parallel > log.simpleFoam 2>&1\n",
    "    STATUS=$?\n",
    "    reconstructPar -latestTime > log.reconstruct 2>&1\n",
    "    [ $STATUS -eq 0 ] && echo '      OK' \\\n",
    "      || { echo '      FAILED — see log.simpleFoam'; exit $STATUS; }\n",
    "else\n",
    "    buoyantBoussinesqSimpleFoam > log.simpleFoam 2>&1 \\\n",
    "      && echo '      OK' \\\n",
    "      || { echo '      FAILED — see log.simpleFoam'; exit 1; }\n",
    "fi\n\n",
    "echo 'Done. Check log.simpleFoam for convergence.'\n"
  )
}


# ===========================================================================
# Public function
# ===========================================================================

#' Prepare OpenFOAM case for nocturnal cold-air drainage simulation
#'
#' @description
#' Generates all OpenFOAM configuration files for a buoyancy-driven nighttime
#' urban cooling simulation using \code{buoyantBoussinesqSimpleFoam}.
#'
#' No wind inlet is applied.  Flow is driven entirely by surface-temperature
#' differentials: vegetation cools fastest after sunset, buildings retain heat
#' longest.  The resulting density gradients create cold-air pooling and
#' drainage patterns at pedestrian level.
#'
#' Outputs to post-process at z = 1.5 m:
#' \describe{
#'   \item{T}{Air temperature — reveals heat islands and cool-air pools}
#'   \item{|U|}{Density-gradient-driven air speed}
#'   \item{T_cool × |U|}{Cool-air transport flux}
#'   \item{U vectors}{Flow direction at pedestrian level}
#' }
#'
#' @param case_dir Character. OpenFOAM case directory (must already contain
#'   \code{constant/triSurface/<stl_file>} written by
#'   \code{prepare_openfoam_inputs()}).
#' @param stl_file Character. Path to building STL file (host path).
#' @param domain Named list with \code{xmin}, \code{xmax}, \code{ymin},
#'   \code{ymax}, \code{zmin}, \code{zmax} (metres, local coordinates).
#'   Returned by \code{prepare_openfoam_inputs()$domain}.
#' @param T_ref Numeric. Ambient (reference) air temperature in Kelvin.
#'   Default 295 K (~22°C).
#' @param T_ground Numeric. Mean ground surface temperature in Kelvin at the
#'   simulation time (e.g. two hours after sunset).  Defaults to
#'   \code{T_ref - 3}, representing net radiative cooling of grass/road mix.
#'   Use the mean of \code{surface_temps} if you supply that argument.
#' @param T_buildings Numeric. Mean building wall/roof temperature in Kelvin.
#'   Buildings retain heat longer than open ground.  Default \code{T_ref - 1}.
#' @param surface_temps Optional named numeric vector of per-land-cover
#'   surface temperatures (K), used only to compute a weighted \code{T_ground}
#'   when \code{T_ground} is not supplied explicitly.  Names should match ESA
#'   WorldCover class labels: \code{"tree"}, \code{"shrub"}, \code{"grass"},
#'   \code{"crop"}, \code{"built"}, \code{"bare"}, \code{"water"}.
#' @param hours_after_sunset Numeric. Hours since sunset (used to document the
#'   case; does not alter BCs automatically).  Default 2.
#' @param z0 Numeric. Aerodynamic roughness length in metres (for
#'   \code{nutkRoughWallFunction} on the ground patch).  Default 0.1 m.
#' @param roughness_raster Character path or \code{NULL}. Path to the ground
#'   roughness raster produced by \code{\link{prepare_openfoam_inputs}}.
#'   Not used during case generation, but stored in \code{files$roughness_raster}
#'   of the return value so downstream code (e.g. \code{plot_foam_map}) can
#'   access it uniformly via \code{foam_inputs$files$roughness_raster}.
#' @param base_cell_size Numeric. Background mesh cell size in metres.
#'   Default 10 m. Use 5 m for finer detail (increases solve time ~8×).
#'   Default 5 m.
#' @param building_refinement Integer. snappyHexMesh refinement level for
#'   buildings.  Default 2.
#' @param n_iterations Integer. Solver iterations.  Buoyancy-driven cases
#'   typically need more than wind cases; default 2000.
#' @param write_interval Integer. Write every N iterations.  Default 200.
#' @param overwrite Logical.  Default FALSE.
#' @param quiet Logical.  Default FALSE.
#'
#' @return Invisibly, a list with \code{case_dir}, \code{files} (paths to
#'   written files), and \code{params} (derived simulation parameters).
#'
#' @seealso \code{\link{prepare_openfoam_inputs}},
#'   \code{\link{run_openfoam_docker}}
#'
#' @export
prepare_nocturnal_case <- function(
    case_dir,
    stl_file,
    domain,
    T_ref                = 295,
    T_ground             = NULL,
    T_buildings          = NULL,
    surface_temps        = NULL,
    hours_after_sunset   = 2,
    z0                   = 0.1,
    roughness_raster     = NULL,
    base_cell_size       = 10,
    building_refinement  = 2L,
    n_iterations         = 500L,
    write_interval       = 100L,
    max_cells            = 3000000L,
    overwrite            = FALSE,
    quiet                = FALSE
) {
  # -- Validation ----------------------------------------------------------
  if (missing(case_dir) || !nzchar(case_dir))
    stop("`case_dir` must be provided.", call. = FALSE)
  if (missing(stl_file) || !file.exists(stl_file))
    stop("`stl_file` not found: ", stl_file, call. = FALSE)
  if (missing(domain) ||
      !all(c("xmin","xmax","ymin","ymax","zmin","zmax") %in% names(domain)))
    stop("`domain` must be a list with xmin/xmax/ymin/ymax/zmin/zmax.",
         call. = FALSE)

  case_dir <- normalizePath(case_dir, mustWork = FALSE)

  system_dir <- file.path(case_dir, "system")
  ic_dir     <- file.path(case_dir, "0")
  const_dir  <- file.path(case_dir, "constant")

  sentinel <- file.path(system_dir, "controlDict")
  if (file.exists(sentinel) && !overwrite)
    stop("Case files already exist in ", case_dir,
         ".\nUse `overwrite = TRUE` to replace them.", call. = FALSE)

  for (d in c(system_dir, ic_dir, const_dir))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  msg <- function(...) if (!isTRUE(quiet)) message(...)

  # -- Derived temperature parameters --------------------------------------
  # Default surface temperatures two hours after sunset (K)
  default_surf <- c(
    tree  = T_ref - 4.0,  # vegetation cools fastest via latent heat
    shrub = T_ref - 3.0,
    grass = T_ref - 3.5,
    crop  = T_ref - 3.0,
    built = T_ref - 1.0,  # sealed surfaces cool slowly
    bare  = T_ref - 2.5,
    water = T_ref - 0.5   # high thermal mass
  )

  surf <- if (!is.null(surface_temps)) {
    nm <- names(surface_temps)
    replace(default_surf, nm[nm %in% names(default_surf)],
            surface_temps[nm[nm %in% names(default_surf)]])
  } else {
    default_surf
  }

  if (is.null(T_ground))    T_ground    <- mean(surf[c("tree","grass","built","bare")])
  if (is.null(T_buildings)) T_buildings <- T_ref - 1.0

  # Turbulence initialisation: low values for calm night
  kappa  <- 0.41; Cmu <- 0.09
  u_star <- 0.1   # representative gentle drainage flow ~ 0.2 m/s
  k_init  <- max(u_star^2 / sqrt(Cmu), 1e-4)
  eps_init <- max(u_star^3 / (kappa * max(z0, 0.01)), 1e-6)

  # Mesh cell counts
  nx <- max(1L, ceiling((domain$xmax - domain$xmin) / base_cell_size))
  ny <- max(1L, ceiling((domain$ymax - domain$ymin) / base_cell_size))
  nz <- max(1L, ceiling((domain$zmax - domain$zmin) / base_cell_size))

  # snappyHexMesh cell limits (same formula as prepare_openfoam_case)
  base_cells       <- as.numeric(nx) * ny * nz
  max_global_cells <- min(ceiling(base_cells * 4), as.integer(max_cells))
  max_local_cells  <- ceiling(max_global_cells / 2)

  loc_x <- (domain$xmin + domain$xmax) / 2
  loc_y <- (domain$ymin + domain$ymax) / 2
  loc_z <- domain$zmax * 0.7

  stl_name   <- basename(stl_file)
  patch_name <- tools::file_path_sans_ext(stl_name)
  ref        <- as.integer(building_refinement)

  msg(sprintf(
    "Solver: buoyantBoussinesqSimpleFoam  |  %d iterations",
    as.integer(n_iterations)
  ))
  msg(sprintf(
    "T_ref = %g K  |  T_ground = %g K (-%g K)  |  T_buildings = %g K (-%g K)",
    T_ref, T_ground, T_ref - T_ground, T_buildings, T_ref - T_buildings
  ))
  msg(sprintf(
    "Domain: %.0f x %.0f x %.0f m  |  Mesh: %d x %d x %d cells (base %g m)",
    domain$xmax - domain$xmin, domain$ymax - domain$ymin,
    domain$zmax - domain$zmin, nx, ny, nz, base_cell_size
  ))

  # -- Write files ---------------------------------------------------------
  files <- list()

  msg("Writing Allrun ...")
  files$allrun <- noc_write(
    make_noc_allrun(), file.path(case_dir, "Allrun"))
  Sys.chmod(files$allrun, mode = "0755")

  msg("Writing system/ ...")
  files$decompose_par_dict <- noc_write(
    make_decompose_par_dict(1L),
    file.path(system_dir, "decomposeParDict"))
  files$block_mesh_dict <- noc_write(
    make_noc_block_mesh_dict(domain, nx, ny, nz),
    file.path(system_dir, "blockMeshDict"))
  files$snappy_hex_mesh_dict <- noc_write(
    make_noc_snappy_dict(stl_name, loc_x, loc_y, loc_z, ref,
                         max_local_cells  = max_local_cells,
                         max_global_cells = max_global_cells),
    file.path(system_dir, "snappyHexMeshDict"))
  files$control_dict <- noc_write(
    make_noc_control_dict(n_iterations, write_interval, z_slice = 1.5),
    file.path(system_dir, "controlDict"))
  files$fv_schemes <- noc_write(
    make_noc_fv_schemes(), file.path(system_dir, "fvSchemes"))
  files$fv_solution <- noc_write(
    make_noc_fv_solution(), file.path(system_dir, "fvSolution"))

  msg("Writing 0/ ...")
  files$U <- noc_write(
    make_noc_u_field(patch_name), file.path(ic_dir, "U"))
  files$p_rgh <- noc_write(
    make_noc_p_rgh_field(patch_name), file.path(ic_dir, "p_rgh"))
  files$T <- noc_write(
    make_noc_t_field(T_ref, T_ground, T_buildings, patch_name),
    file.path(ic_dir, "T"))
  omega_init <- eps_init / k_init   # kOmegaSST: omega = eps / k
  files$k <- noc_write(
    make_noc_k_field(k_init, patch_name), file.path(ic_dir, "k"))
  files$omega <- noc_write(
    make_noc_omega_field(omega_init, patch_name), file.path(ic_dir, "omega"))
  files$nut <- noc_write(
    make_noc_nut_field(z0, patch_name), file.path(ic_dir, "nut"))
  files$alphat <- noc_write(
    make_noc_alphat_field(patch_name), file.path(ic_dir, "alphat"))

  # Pass-through: roughness_raster is not written here (it comes from
  # prepare_openfoam_inputs), but store the path so callers can access it
  # consistently as foam_inputs$files$roughness_raster.
  if (!is.null(roughness_raster))
    files$roughness_raster <- roughness_raster

  msg("Writing constant/ ...")
  files$g <- noc_write(
    make_gravity(), file.path(const_dir, "g"))
  files$transport_props <- noc_write(
    make_noc_transport_properties(T_ref), file.path(const_dir, "transportProperties"))
  files$turbulence_props <- noc_write(
    make_noc_turbulence_properties(), file.path(const_dir, "turbulenceProperties"))

  msg("Case files written to: ", case_dir)
  msg(sprintf(
    "Next step: run_openfoam_docker(\"%s\")", case_dir))

  invisible(list(
    case_dir = case_dir,
    files    = files,
    params   = list(
      T_ref              = T_ref,
      T_ground           = T_ground,
      T_buildings        = T_buildings,
      surface_temps      = surf,
      hours_after_sunset = hours_after_sunset,
      z0                 = z0,
      k_init             = k_init,
      omega_init         = omega_init,
      n_iterations       = as.integer(n_iterations)
    )
  ))
}


# ===========================================================================
# Post-processing
# ===========================================================================

#' Read OpenFOAM pedestrian-level slice and compute thermal maps
#'
#' @description
#' Reads the horizontal surface sample at z = 1.5 m written by
#' \code{buoyantBoussinesqSimpleFoam} (via the \code{pedestrianSlice}
#' function object in \code{controlDict}) and returns a multi-layer
#' \code{SpatRaster} with the four pedestrian-level maps from
#' \code{prepare_nocturnal_case()}:
#'
#' \describe{
#'   \item{T_air}{Air temperature in K}
#'   \item{U_mag}{Wind speed magnitude in m/s (density-gradient driven)}
#'   \item{T_cool}{Cooling relative to ambient: \eqn{T_{ref} - T_{air}} (K)}
#'   \item{T_cool_flux}{\eqn{\max(T_{cool},0) \times |U|} — cool-air transport flux}
#'   \item{Ux}{East–west velocity component (m/s)}
#'   \item{Uy}{North–south velocity component (m/s)}
#' }
#'
#' @param case_dir Character. OpenFOAM case directory (same as passed to
#'   \code{prepare_nocturnal_case()}).
#' @param T_ref Numeric. Reference (ambient) temperature in K used to compute
#'   \code{T_cool = T_ref - T_air}.  Should match the value used in
#'   \code{prepare_nocturnal_case()}.  Default 295 K.
#' @param time_step Character or numeric. \code{"latest"} (default) uses the
#'   highest-numbered time directory; pass an integer to pick a specific step.
#' @param resolution Numeric. Raster cell size in metres for output grids.
#'   Default 2 m.
#' @param crs Character or numeric. CRS to assign to the returned raster
#'   (e.g. \code{32617} for UTM 17N).  The OpenFOAM output is in local
#'   coordinates (origin at domain SW corner); use \code{foam_inputs$origin}
#'   and \code{foam_inputs$crs} to georeferenc the result after rasterization.
#'   Default \code{NA} (no CRS).
#' @param buildings Optional \code{sf} object of building footprints in the
#'   same local coordinate system as the raster (pass
#'   \code{readRDS(foam_inputs\$files\$buildings_rds)}).  When supplied,
#'   building interiors are masked to \code{NA} \emph{after} gap-filling, so
#'   that small buildings which would otherwise be flooded over by the
#'   interpolation step are correctly restored as voids.  Default \code{NULL}.
#' @param quiet Logical. Default \code{FALSE}.
#'
#' @return A \code{terra::SpatRaster} with six layers:
#'   \code{T_air}, \code{U_mag}, \code{T_cool}, \code{T_cool_flux},
#'   \code{Ux}, \code{Uy}.
#'
#' @seealso \code{\link{prepare_nocturnal_case}}, \code{\link{run_openfoam_docker}}
#'
#' @export
read_foam_pedestrian_slice <- function(
    case_dir,
    T_ref      = 295,
    time_step  = "latest",
    resolution = 2,
    crs        = NA,
    buildings  = NULL,   # auto-detected from case_dir if NULL
    quiet      = FALSE
) {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required. Install it with install.packages('terra').",
         call. = FALSE)

  case_dir <- normalizePath(case_dir, mustWork = TRUE)

  # Auto-detect buildings from the case directory if not supplied.
  # prepare_openfoam_inputs() always writes buildings to this path.
  if (is.null(buildings)) {
    bpath <- file.path(case_dir, "constant", "gloBFPr", "metadata",
                       "buildings_openfoam.rds")
    if (file.exists(bpath)) {
      buildings <- tryCatch(readRDS(bpath), error = function(e) NULL)
      if (!isTRUE(quiet) && !is.null(buildings))
        message("Auto-loaded ", nrow(buildings), " buildings from case directory.")
    }
  }

  pp_dir   <- file.path(case_dir, "postProcessing", "pedestrianSlice")

  if (!dir.exists(pp_dir))
    stop("No postProcessing/pedestrianSlice directory found in:\n  ", case_dir,
         "\nRun the simulation first with run_openfoam_docker().", call. = FALSE)

  # Enumerate available time steps (numeric subdirectories)
  time_dirs <- list.dirs(pp_dir, full.names = FALSE, recursive = FALSE)
  num_times <- suppressWarnings(as.numeric(time_dirs))
  num_times <- sort(num_times[!is.na(num_times)])

  if (length(num_times) == 0)
    stop("No time-step directories found under ", pp_dir, call. = FALSE)

  t_chosen <- if (identical(time_step, "latest")) {
    max(num_times)
  } else {
    ts <- suppressWarnings(as.numeric(time_step))
    if (is.na(ts) || !(ts %in% num_times))
      stop(sprintf("time_step %s not found. Available: %s",
                   time_step, paste(num_times, collapse = ", ")), call. = FALSE)
    ts
  }

  t_dir <- file.path(pp_dir, as.character(t_chosen))

  # OpenFOAM v2506 reversed the .raw naming convention:
  #   old (pre-v2506): <surfaceName>_<field>.raw  → z1p5m_T.raw
  #   new (v2506+):    <field>_<surfaceName>.raw  → T_z1p5m.raw
  # Try both patterns so the package works across versions.
  find_noc_raw <- function(dir, field, surf = "z1p5m") {
    candidates <- c(
      file.path(dir, paste0(surf, "_", field, ".raw")),  # pre-v2506
      file.path(dir, paste0(field, "_", surf, ".raw"))   # v2506+
    )
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0L)
      stop(sprintf(
        "Raw file for field '%s' not found in:\n  %s\nTried:\n  %s",
        field, dir, paste(candidates, collapse = "\n  ")
      ), call. = FALSE)
    found[[1L]]
  }

  t_file <- find_noc_raw(t_dir, "T")
  u_file <- find_noc_raw(t_dir, "U")

  if (!isTRUE(quiet)) {
    message(sprintf("Reading slice at time step %g from %s", t_chosen, t_dir))
  }

  # raw format: 3 header comment lines, then x y z value (scalar)
  #             or x y z vx vy vz (vector)
  read_raw <- function(path, col_names) {
    # Skip lines beginning with #
    raw_text <- readLines(path)
    data_lines <- raw_text[!grepl("^\\s*#", raw_text) & nzchar(trimws(raw_text))]
    if (length(data_lines) == 0)
      stop("No data found in ", path, call. = FALSE)
    utils::read.table(
      text = paste(data_lines, collapse = "\n"),
      header = FALSE, col.names = col_names
    )
  }

  T_dat <- read_raw(t_file, c("x", "y", "z", "T"))
  U_dat <- read_raw(u_file, c("x", "y", "z", "Ux", "Uy", "Uz"))

  dat <- merge(T_dat[, c("x", "y", "T")],
               U_dat[, c("x", "y", "Ux", "Uy", "Uz")],
               by = c("x", "y"), all = TRUE)

  dat$U_mag        <- sqrt(dat$Ux^2 + dat$Uy^2 + dat$Uz^2)
  dat$T_cool       <- T_ref - dat$T
  dat$T_cool_flux  <- pmax(dat$T_cool, 0) * dat$U_mag

  # Build raster template from point extent
  r_template <- terra::rast(
    xmin       = min(dat$x, na.rm = TRUE),
    xmax       = max(dat$x, na.rm = TRUE),
    ymin       = min(dat$y, na.rm = TRUE),
    ymax       = max(dat$y, na.rm = TRUE),
    resolution = resolution,
    crs        = if (is.na(crs)) "" else as.character(crs)
  )

  pts <- terra::vect(dat, geom = c("x", "y"),
                     crs = if (is.na(crs)) "" else as.character(crs))

  rasterize_mean <- function(field)
    terra::rasterize(pts, r_template, field = field, fun = "mean")

  # Gap-fill sparse CFD points then apply Gaussian smoothing so the output
  # is a continuous field rather than isolated dots.
  # The background mesh is 'base_cell_size' m, but the raster resolution
  # may be much finer (e.g. 2 m vs 10 m mesh), leaving large empty areas.
  #
  # Strategy:
  #   1. Focal mean with a window large enough to span the coarsest mesh
  #      cell (~base_cell_size / resolution), applied twice to propagate
  #      values into corners.
  #   2. Gaussian smooth with sigma = smooth_sigma cells.
  smooth_noc_layer <- function(r, gap_w = 15L, sigma = 2) {
    nm <- names(r)
    # Pass 1 — fill NAs from neighbouring cells
    r <- terra::focal(r, w = gap_w, fun = mean, na.policy = "only",
                      na.rm = TRUE)
    # Pass 2 — fill any remaining NAs at domain edges
    r <- terra::focal(r, w = gap_w, fun = mean, na.policy = "only",
                      na.rm = TRUE)
    # Gaussian smoothing kernel
    ksize  <- max(3L, 2L * ceiling(3 * sigma) + 1L)
    ax     <- seq(-(ksize %/% 2L), ksize %/% 2L)
    kernel <- outer(ax, ax, function(x, y) exp(-(x^2 + y^2) / (2 * sigma^2)))
    kernel <- kernel / sum(kernel)
    r <- terra::focal(r, w = kernel, na.rm = TRUE)
    names(r) <- nm
    r
  }

  # Determine gap-fill window: span at least 1.5 × (base mesh / resolution)
  # Default base_cell_size for nocturnal is 10 m; resolution default 2 m.
  # ceil(1.5 * 10 / 2) * 2 + 1 = 16 → round up to odd number.
  gap_w <- as.integer(2 * ceiling(1.5 * 10 / resolution) + 1L)
  gap_w <- max(gap_w, 7L)

  raw_layers <- list(
    T_air       = rasterize_mean("T"),
    U_mag       = rasterize_mean("U_mag"),
    T_cool      = rasterize_mean("T_cool"),
    T_cool_flux = rasterize_mean("T_cool_flux"),
    Ux          = rasterize_mean("Ux"),
    Uy          = rasterize_mean("Uy")
  )

  smoothed <- lapply(names(raw_layers), function(nm) {
    r <- raw_layers[[nm]]
    names(r) <- nm
    smooth_noc_layer(r, gap_w = gap_w, sigma = 2)
  })

  result <- terra::rast(smoothed)
  names(result) <- c("T_air", "U_mag", "T_cool", "T_cool_flux", "Ux", "Uy")

  # Attach buildings as an R attribute so plot_foam_map can overlay them
  # automatically without the user needing to pass them explicitly.
  if (!is.null(buildings))
    attr(result, "buildings") <- buildings

  # Re-apply building mask AFTER gap-filling.
  # Gap-fill floods over small buildings (those narrower than the fill window),
  # so without this step only the largest buildings show as blank patches.
  # Rasterising the footprints and masking to NA restores all building voids
  # regardless of building size.
  if (!is.null(buildings)) {
    if (!requireNamespace("sf", quietly = TRUE))
      warning("Package 'sf' required to mask buildings; skipping.", call. = FALSE)
    else {
      bvect        <- terra::vect(sf::st_geometry(buildings))
      building_r   <- terra::rasterize(bvect, result[[1L]], field = 1L,
                                       background = NA)
      result       <- terra::mask(result, building_r,
                                  maskvalues = 1L, updatevalue = NA)
    }
  }

  if (!isTRUE(quiet)) {
    message(sprintf(
      "Raster: %d x %d cells at %g m resolution | %d sample points",
      terra::nrow(result), terra::ncol(result), resolution, nrow(dat)
    ))
  }

  result
}
