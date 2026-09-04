# ===========================================================================
# Internal string-building helpers
# ===========================================================================

#' @noRd
foam_header <- function(class, object) {
  paste0(
    "/*--------------------------------*- C++ -*----------------------------------*\\\n",
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
write_foam_file <- function(content, path) {
  writeLines(content, con = path)
  normalizePath(path, mustWork = FALSE)
}

# ---------------------------------------------------------------------------
# blockMeshDict
# ---------------------------------------------------------------------------
#' @noRd
make_block_mesh_dict <- function(domain, nx, ny, nz) {
  d <- domain
  paste0(
    foam_header("dictionary", "blockMeshDict"),
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
    # No curved edges - required section even when empty
    "edges\n(\n);\n\n",
    # Face vertex ordering: right-hand rule, outward normals.
    # Vertices: 0=(xmin,ymin,zmin) 1=(xmax,ymin,zmin) 2=(xmax,ymax,zmin)
    #           3=(xmin,ymax,zmin) 4=(xmin,ymin,zmax) 5=(xmax,ymin,zmax)
    #           6=(xmax,ymax,zmax) 7=(xmin,ymax,zmax)
    # Boundary face vertex ordering: right-hand rule with outward normals.
    # Vertices: 0=(xmin,ymin,zmin) 1=(xmax,ymin,zmin) 2=(xmax,ymax,zmin)
    #           3=(xmin,ymax,zmin) 4=(xmin,ymin,zmax) 5=(xmax,ymin,zmax)
    #           6=(xmax,ymax,zmax) 7=(xmin,ymax,zmax)
    # Normal direction verified via (v1-v0) x (v3-v0):
    #   inlet  (0 4 7 3): (0,0,+z)x(0,+y,0) = (-yz,0,0) -> -x OK
    #   outlet (1 2 6 5): (0,+y,0)x(0,0,+z) = (+yz,0,0) -> +x OK
    #   ground (3 2 1 0): (+x,0,0)x(0,-y,0) = (0,0,-xy)  -> -z OK
    #   top    (4 5 6 7): (+x,0,0)x(0,+y,0) = (0,0,+xy)  -> +z OK
    #   y-min  (0 1 5 4): (+x,0,0)x(0,0,+z) = (0,-xz,0)  -> -y OK
    #   y-max  (3 7 6 2): (0,0,+z)x(+x,0,0) = (0,+xz,0)  -> +y OK
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
    # Required closing section
    "mergePatchPairs\n(\n);\n"
  )
}

# ---------------------------------------------------------------------------
# snappyHexMeshDict
# ---------------------------------------------------------------------------
#' @noRd
make_snappy_hex_mesh_dict <- function(stl_name, loc_x, loc_y, loc_z,
                                       ref_min, ref_max,
                                       max_local_cells  = 1000000L,
                                       max_global_cells = 4000000L) {
  patch_name <- tools::file_path_sans_ext(stl_name)   # e.g. "buildings"
  paste0(
    foam_header("dictionary", "snappyHexMeshDict"),
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
    sprintf("            level (%d %d);\n", ref_min, ref_max),
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

# ---------------------------------------------------------------------------
# controlDict
# ---------------------------------------------------------------------------
#' @noRd
make_control_dict <- function(n_iterations) {
  paste0(
    foam_header("dictionary", "controlDict"),
    "application     simpleFoam;\n\n",
    # simpleFoam links the atmospheric models library itself, but utilities
    # such as potentialFoam and postProcess do not. Without this they fail
    # with "Unknown patchField type atmBoundaryLayerInletVelocity" when
    # reading the ABL inlet in 0/U.
    "libs            (\"libatmosphericModels.so\");\n\n",
    "startFrom       startTime;\n",
    "startTime       0;\n",
    "stopAt          endTime;\n",
    sprintf("endTime         %d;\n", as.integer(n_iterations)),
    "deltaT          1;\n\n",
    "writeControl    timeStep;\n",
    "writeInterval   100;\n",
    "purgeWrite      3;\n\n",
    "writeFormat     ascii;\n",
    "writePrecision  8;\n",
    "writeCompression off;\n\n",
    "timeFormat      general;\n",
    "timePrecision   6;\n\n",
    "runTimeModifiable true;\n\n",
    # Empty functions block - residuals and solver output are captured in
    # log.simpleFoam by the Allrun redirect (> log.simpleFoam 2>&1).
    # Omitting #includeFunc / function objects avoids path-lookup failures
    # that occur in some Docker image configurations of OpenFOAM 2506.
    "functions\n{\n}\n"
  )
}

# ---------------------------------------------------------------------------
# fvSchemes
# ---------------------------------------------------------------------------
#' @noRd
make_fv_schemes <- function() {
  paste0(
    foam_header("dictionary", "fvSchemes"),
    "ddtSchemes\n{\n    default         steadyState;\n}\n\n",
    "gradSchemes\n{\n",
    "    default         Gauss linear;\n",
    "    grad(U)         cellLimited Gauss linear 1;\n}\n\n",
    # `bounded` is required for steadyState: phi is not exactly
    # divergence-free while iterating, and the bounded variant adds the
    # -fvm::Sp(fvc::div(phi), psi) correction for that. Without it the
    # convection term drives k and epsilon negative ("bounding epsilon,
    # min: -..."), nut = Cmu*k^2/epsilon blows up, and the run dies with
    # SIGFPE a few hundred iterations later.
    # limitedLinear (not linearUpwind) on k/epsilon: TVD-limited, so it
    # cannot overshoot into negative values for these strictly-positive
    # quantities.
    "divSchemes\n{\n",
    "    default                         none;\n",
    "    div(phi,U)                      bounded Gauss linearUpwind grad(U);\n",
    # First-order upwind on k/epsilon. limitedLinear is second-order and can
    # still overshoot in the strong shear layers between buildings, which lets
    # k and epsilon drift negative and eventually run away (k ~ 1e50 around
    # iteration 1200). Upwind is unconditionally bounded; the accuracy cost
    # falls on the turbulence quantities, not on U, and this matches the
    # configuration used by the stable nocturnal case.
    "    div(phi,k)                      bounded Gauss upwind;\n",
    "    div(phi,epsilon)                bounded Gauss upwind;\n",
    "    div((nuEff*dev2(T(grad(U)))))   Gauss linear;\n}\n\n",  # dev2 required in v2206+
    "laplacianSchemes\n{\n    default         Gauss linear corrected;\n}\n\n",
    "interpolationSchemes\n{\n    default         linear;\n}\n\n",
    "snGradSchemes\n{\n    default         corrected;\n}\n"
  )
}

# ---------------------------------------------------------------------------
# fvSolution
# ---------------------------------------------------------------------------
#' @noRd
make_fv_solution <- function(consistent = FALSE, residual_control = 1e-4) {
  paste0(
    foam_header("dictionary", "fvSolution"),
    "solvers\n{\n",
    "    p\n    {\n",
    "        solver          GAMG;\n",
    # faceAreaPair + mergeLevels 1 give a well-conditioned agglomeration
    # hierarchy on snappyHexMesh's irregular cells.
    "        agglomerator    faceAreaPair;\n",
    "        mergeLevels     1;\n",
    "        cacheAgglomeration true;\n",
    # The v2506 default (10) is far too small in MPI parallel: split across
    # N ranks the coarsest level becomes degenerate, PCG's inner product
    # (sumProd) evaluates to zero, and the solver dies with SIGFPE inside
    # GAMGSolver::solveCoarsestLevel. 200 keeps that level solvable.
    "        nCellsInCoarsestLevel 200;\n",
    "        smoother        GaussSeidel;\n",
    "        tolerance       1e-7;\n",
    "        relTol          0.01;\n",
    "    }\n\n",
    "    \"(U|k|epsilon)\"\n    {\n",
    "        solver          smoothSolver;\n",
    "        smoother        symGaussSeidel;\n",
    "        tolerance       1e-7;\n",
    "        relTol          0.1;\n",
    "    }\n}\n\n",
    "SIMPLE\n{\n",
    "    nNonOrthogonalCorrectors 2;\n",
    # SIMPLEC ("consistent yes") retains the neighbour-coefficient term the
    # standard SIMPLE pressure correction drops. The correction is then
    # consistent with the momentum equation, so pressure needs no
    # under-relaxation at all and U can run near 0.9 instead of 0.5 -
    # typically 2-3x fewer iterations to reach the same residuals.
    if (isTRUE(consistent)) "    consistent      yes;\n" else "",
    # pRefCell/pRefValue anchor the pressure level - prevents drift that can
    # cause SIGFPE when p reaches machine-epsilon near walls
    "    pRefCell        0;\n",
    "    pRefValue       0;\n",
    "    residualControl\n    {\n",
    sprintf("        p               %g;\n", residual_control),
    sprintf("        U               %g;\n", residual_control),
    sprintf("        \"(k|epsilon)\"   %g;\n", residual_control),
    "    }\n}\n\n",
    "relaxationFactors\n{\n",
    if (isTRUE(consistent)) {
      # SIMPLEC: no pressure under-relaxation, aggressive equation relaxation.
      paste0(
        "    fields\n    {\n        p    1.0;\n    }\n",
        "    equations\n    {\n",
        "        U        0.9;\n",
        "        k        0.7;\n",
        "        epsilon  0.7;\n",
        "    }\n"
      )
    } else {
      # Standard SIMPLE: conservative values that tolerate the steep gradients
      # in building wakes.
      paste0(
        "    fields\n    {\n        p    0.3;\n    }\n",
        "    equations\n    {\n",
        "        U        0.5;\n",
        # k/epsilon relaxed harder than U: they are the fields that ran away,
        # and damping them costs little because they are already first-order
        # upwinded.
        "        k        0.3;\n",
        "        epsilon  0.3;\n",
        "    }\n"
      )
    },
    "}\n"
  )
}

# ---------------------------------------------------------------------------
# 0/U
# ---------------------------------------------------------------------------
#' @noRd
make_u_field <- function(u_mag, fd_x, fd_y, z_ref, z0, patch_name) {
  # Reasonable internal field: uniform inlet speed in flow direction
  paste0(
    foam_header("volVectorField", "U"),
    sprintf("dimensions      [ 0 1 -1 0 0 0 0 ];\n\n"),
    sprintf("internalField   uniform (%g %g 0);\n\n", u_mag * fd_x, u_mag * fd_y),
    "boundaryField\n{\n",
    # inlet - ABL log-law profile
    # Note: OpenFOAM v2206+ requires d (displacement height) as PatchFunction1
    "    inlet\n    {\n",
    "        type            atmBoundaryLayerInletVelocity;\n",
    sprintf("        flowDir         (%g %g 0);\n", fd_x, fd_y),
    "        zDir            (0 0 1);\n",
    sprintf("        Uref            %g;\n", u_mag),
    sprintf("        Zref            %g;\n", z_ref),
    sprintf("        z0              uniform %g;\n", z0),
    "        d               uniform 0;\n",
    "        zGround         uniform 0;\n",
    "        value           $internalField;\n",
    "    }\n",
    # outlet
    "    outlet\n    {\n        type            zeroGradient;\n    }\n",
    # ground - no-slip wall
    "    ground\n    {\n        type            noSlip;\n    }\n",
    # buildings - no-slip wall
    "    ", patch_name, "\n    {\n        type            noSlip;\n    }\n",
    # top - slip (free-slip symmetry)
    "    top\n    {\n        type            slip;\n    }\n",
    # sides - slip
    "    sides\n    {\n        type            slip;\n    }\n",
    "}\n"
  )
}

# ---------------------------------------------------------------------------
# 0/p
# ---------------------------------------------------------------------------
#' @noRd
make_p_field <- function(patch_name) {
  paste0(
    foam_header("volScalarField", "p"),
    "dimensions      [ 0 2 -2 0 0 0 0 ];\n\n",
    "internalField   uniform 0;\n\n",
    "boundaryField\n{\n",
    "    inlet\n    {\n        type            zeroGradient;\n    }\n",
    "    outlet\n    {\n        type            fixedValue;\n        value           uniform 0;\n    }\n",
    "    ground\n    {\n        type            zeroGradient;\n    }\n",
    "    ", patch_name, "\n    {\n        type            zeroGradient;\n    }\n",
    "    top\n    {\n        type            zeroGradient;\n    }\n",
    "    sides\n    {\n        type            zeroGradient;\n    }\n",
    "}\n"
  )
}

# ---------------------------------------------------------------------------
# 0/k
# ---------------------------------------------------------------------------
#' @noRd
make_k_field <- function(k_init, u_mag, fd_x, fd_y, z_ref, z0, patch_name) {
  paste0(
    foam_header("volScalarField", "k"),
    "dimensions      [ 0 2 -2 0 0 0 0 ];\n\n",
    sprintf("internalField   uniform %g;\n\n", k_init),
    "boundaryField\n{\n",
    "    inlet\n    {\n",
    "        type            atmBoundaryLayerInletK;\n",
    sprintf("        flowDir         (%g %g 0);\n", fd_x, fd_y),
    "        zDir            (0 0 1);\n",
    sprintf("        Uref            %g;\n", u_mag),
    sprintf("        Zref            %g;\n", z_ref),
    sprintf("        z0              uniform %g;\n", z0),
    "        d               uniform 0;\n",
    "        zGround         uniform 0;\n",
    "        value           $internalField;\n",
    "    }\n",
    "    outlet\n    {\n        type            zeroGradient;\n    }\n",
    "    ground\n    {\n",
    "        type            kqRWallFunction;\n",
    sprintf("        value           uniform %g;\n", k_init),
    "    }\n",
    "    ", patch_name, "\n    {\n",
    "        type            kqRWallFunction;\n",
    sprintf("        value           uniform %g;\n", k_init),
    "    }\n",
    "    top\n    {\n        type            zeroGradient;\n    }\n",
    "    sides\n    {\n        type            zeroGradient;\n    }\n",
    "}\n"
  )
}

# ---------------------------------------------------------------------------
# 0/epsilon
# ---------------------------------------------------------------------------
#' @noRd
make_epsilon_field <- function(eps_init, u_mag, fd_x, fd_y, z_ref, z0,
                               patch_name) {
  paste0(
    foam_header("volScalarField", "epsilon"),
    "dimensions      [ 0 2 -3 0 0 0 0 ];\n\n",
    sprintf("internalField   uniform %g;\n\n", eps_init),
    "boundaryField\n{\n",
    "    inlet\n    {\n",
    "        type            atmBoundaryLayerInletEpsilon;\n",
    sprintf("        flowDir         (%g %g 0);\n", fd_x, fd_y),
    "        zDir            (0 0 1);\n",
    sprintf("        Uref            %g;\n", u_mag),
    sprintf("        Zref            %g;\n", z_ref),
    sprintf("        z0              uniform %g;\n", z0),
    "        d               uniform 0;\n",
    "        zGround         uniform 0;\n",
    "        value           $internalField;\n",
    "    }\n",
    "    outlet\n    {\n        type            zeroGradient;\n    }\n",
    "    ground\n    {\n",
    "        type            epsilonWallFunction;\n",
    sprintf("        value           uniform %g;\n", eps_init),
    "    }\n",
    "    ", patch_name, "\n    {\n",
    "        type            epsilonWallFunction;\n",
    sprintf("        value           uniform %g;\n", eps_init),
    "    }\n",
    "    top\n    {\n        type            zeroGradient;\n    }\n",
    "    sides\n    {\n        type            zeroGradient;\n    }\n",
    "}\n"
  )
}

# ---------------------------------------------------------------------------
# 0/nut
# ---------------------------------------------------------------------------
#' @noRd
make_nut_field <- function(z0, patch_name) {
  # Ks (roughness height) ~ 20 * z0 is a widely used urban-terrain approximation
  # for nutkRoughWallFunction with Cs = 0.5
  ks <- z0 * 20
  paste0(
    foam_header("volScalarField", "nut"),
    "dimensions      [ 0 2 -1 0 0 0 0 ];\n\n",
    "internalField   uniform 0;\n\n",
    "boundaryField\n{\n",
    "    inlet\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "    outlet\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    # ground - rough wall function
    "    ground\n    {\n",
    "        type            nutkRoughWallFunction;\n",
    sprintf("        Ks              uniform %g;  // ~ 20 * z0 (z0 = %g m)\n", ks, z0),
    "        Cs              uniform 0.5;\n",  # v2506: Cs is Field<scalar>
    "        value           uniform 0;\n",
    "    }\n",
    # buildings - smooth wall function
    "    ", patch_name, "\n    {\n",
    "        type            nutkWallFunction;\n",
    "        value           uniform 0;\n",
    "    }\n",
    "    top\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "    sides\n    {\n        type            calculated;\n        value           uniform 0;\n    }\n",
    "}\n"
  )
}

# ---------------------------------------------------------------------------
# constant/turbulenceProperties
# ---------------------------------------------------------------------------
#' @noRd
make_turbulence_properties <- function() {
  paste0(
    foam_header("dictionary", "turbulenceProperties"),
    "simulationType  RAS;\n\n",
    "RAS\n{\n",
    "    RASModel        kEpsilon;\n",
    "    turbulence      on;\n",
    "    printCoeffs     on;\n\n",
    "    // Lower bounds prevent k/epsilon going negative -> avoids SIGFPE\n",
    "    kMin            kMin [ 0 2 -2 0 0 0 0 ] 1e-10;\n",
    "    epsilonMin      epsilonMin [ 0 2 -3 0 0 0 0 ] 1e-10;\n",
    "}\n"
  )
}

# ---------------------------------------------------------------------------
# constant/transportProperties
# ---------------------------------------------------------------------------
#' @noRd
make_transport_properties <- function() {
  paste0(
    foam_header("dictionary", "transportProperties"),
    "transportModel  Newtonian;\n\n",
    "// Kinematic viscosity of air at ~20  degC\n",
    "nu              [ 0 2 -1 0 0 0 0 ] 1.5e-05;\n"
  )
}

# ---------------------------------------------------------------------------
# decomposeParDict  (scotch - no manual layout, works for any ncpus)
# ---------------------------------------------------------------------------
#' @noRd
make_decompose_par_dict <- function(ncpus = 1L) {
  paste0(
    foam_header("dictionary", "decomposeParDict"),
    # numberOfSubdomains is patched at runtime by the Allrun sed command
    # when NPROC > 1, so the value written here only matters for ncpus = 1.
    sprintf("numberOfSubdomains  %d;\n\n", as.integer(ncpus)),
    "method          scotch;\n"
  )
}

# ---------------------------------------------------------------------------
# Allrun
# ---------------------------------------------------------------------------
#' @noRd
make_allrun <- function(potential_init = FALSE) {
  paste0(
    "#!/bin/sh\n",
    "# Allrun - generated by gloBFPr::prepare_openfoam_case()\n",
    "# Pass NPROC env variable (via docker -e NPROC=N) to enable MPI parallel.\n",
    "cd \"${0%/*}\" || exit 1\n\n",
    "NPROC=${NPROC:-1}\n\n",
    "# Clean previous run artefacts\n",
    "rm -rf constant/polyMesh processor* log.*\n",
    "find . -maxdepth 1 -name '[0-9]*' ! -name '0' -exec rm -rf {} +\n\n",
    "echo '[1/3] blockMesh ...'\n",
    "blockMesh > log.blockMesh 2>&1 \\\n",
    "  && echo '      OK' \\\n",
    "  || { echo '      FAILED - see log.blockMesh'; exit 1; }\n\n",
    "echo '[2/3] snappyHexMesh ...'\n",
    "if [ \"$NPROC\" -gt 1 ]; then\n",
    "    sed -i \"s/numberOfSubdomains.*/numberOfSubdomains $NPROC;/\" system/decomposeParDict\n",
    "    decomposePar -force > log.decomposeMesh 2>&1 \\\n",
    "      || { echo '      FAILED - see log.decomposeMesh'; exit 1; }\n",
    "    mpirun -np \"$NPROC\" --allow-run-as-root snappyHexMesh -overwrite -parallel > log.snappyHexMesh 2>&1 \\\n",
    "      || { echo '      FAILED - see log.snappyHexMesh'; exit 1; }\n",
    "    reconstructParMesh -constant -mergeTol 1e-6 > log.reconstructMesh 2>&1 \\\n",
    "      || { echo '      FAILED - see log.reconstructMesh'; exit 1; }\n",
    "    rm -rf processor*\n",
    "    echo '      OK (parallel)'\n",
    "else\n",
    "    snappyHexMesh -overwrite > log.snappyHexMesh 2>&1 \\\n",
    "      && echo '      OK' \\\n",
    "      || { echo '      FAILED - see log.snappyHexMesh'; exit 1; }\n",
    "fi\n\n",
    "echo '[3/3] simpleFoam ...'\n",
    "if [ \"$NPROC\" -gt 1 ]; then\n",
    "    sed -i \"s/numberOfSubdomains.*/numberOfSubdomains $NPROC;/\" system/decomposeParDict\n",
    "    decomposePar -force > log.decompose 2>&1 \\\n",
    "      || { echo '      FAILED - see log.decompose'; exit 1; }\n",
    # potentialFoam solves inviscid potential flow (seconds) to give SIMPLE a
    # divergence-free starting velocity field instead of uniform/zero flow.
    # This removes the violent initial transient and typically saves 30-50%
    # of SIMPLE iterations. Non-fatal: if it fails we just start from 0/.
    if (isTRUE(potential_init)) paste0(
      "    echo '      potentialFoam init ...'\n",
      "    mpirun -np \"$NPROC\" --allow-run-as-root potentialFoam ",
      "-parallel -writephi > log.potentialFoam 2>&1 \\\n",
      "      && echo '      init OK' \\\n",
      "      || echo '      init skipped (see log.potentialFoam)'\n"
    ) else "",
    "    mpirun -np \"$NPROC\" --allow-run-as-root simpleFoam -parallel > log.simpleFoam 2>&1\n",
    "    STATUS=$?\n",
    "    reconstructPar -latestTime > log.reconstruct 2>&1\n",
    "    [ $STATUS -eq 0 ] && echo '      OK' \\\n",
    "      || { echo '      FAILED - see log.simpleFoam'; exit $STATUS; }\n",
    "else\n",
    if (isTRUE(potential_init)) paste0(
      "    echo '      potentialFoam init ...'\n",
      "    potentialFoam -writephi > log.potentialFoam 2>&1 \\\n",
      "      && echo '      init OK' \\\n",
      "      || echo '      init skipped (see log.potentialFoam)'\n"
    ) else "",
    "    simpleFoam > log.simpleFoam 2>&1 \\\n",
    "      && echo '      OK' \\\n",
    "      || { echo '      FAILED - see log.simpleFoam'; exit 1; }\n",
    "fi\n\n",
    "echo 'Done. Check log.simpleFoam for convergence.'\n"
  )
}


#' Default MPI rank count for OpenFOAM runs
#'
#' @description
#' On Apple Silicon the physical cores are heterogeneous: a 10-core part is
#' typically 8 performance + 2 efficiency cores. MPI is synchronous, so every
#' iteration waits on the slowest rank - placing ranks on efficiency cores
#' throttles the whole solve to E-core speed. Using performance cores only is
#' usually faster despite the lower rank count.
#'
#' Falls back to all physical cores on other platforms.
#'
#' @return Integer number of MPI ranks.
#' @noRd
foam_default_ncpus <- function() {
  fallback <- tryCatch(
    parallel::detectCores(logical = FALSE),
    error = function(e) 1L
  )
  if (!is.finite(fallback) || fallback < 1L) fallback <- 1L

  # Apple Silicon: query the performance-core count directly.
  if (Sys.info()[["sysname"]] == "Darwin") {
    p <- suppressWarnings(tryCatch(
      as.integer(system("sysctl -n hw.perflevel0.physicalcpu",
                        intern = TRUE, ignore.stderr = TRUE)[1]),
      error = function(e) NA_integer_
    ))
    if (!is.na(p) && p >= 1L) return(p)
  }

  as.integer(fallback)
}


# ===========================================================================
# Public functions
# ===========================================================================

#' Prepare OpenFOAM case files for urban wind CFD
#'
#' @description
#' Writes a complete set of OpenFOAM case files for a steady-state urban wind
#' simulation using \code{simpleFoam} with standard k-epsilon turbulence.
#' The generated case is ready to run with \code{run_openfoam_docker()} (or
#' manually via Docker) without further editing.
#'
#' The solver chain is: \strong{blockMesh} (background Cartesian mesh) ->
#' \strong{snappyHexMesh} (building STL cut-in) ->
#' \strong{simpleFoam} (steady-state incompressible RANS).
#'
#' Boundary conditions:
#' \itemize{
#'   \item \strong{inlet} (x-min face): atmospheric boundary layer log-law
#'     profiles via \code{atmBoundaryLayerInletVelocity / K / Epsilon}.
#'   \item \strong{outlet} (x-max face): \code{zeroGradient} U,
#'     \code{fixedValue} p = 0.
#'   \item \strong{ground}: no-slip with \code{nutkRoughWallFunction}
#'     (Ks approx. 20 x z0).
#'   \item \strong{buildings}: no-slip with smooth \code{nutkWallFunction}.
#'   \item \strong{top / sides}: free-slip.
#' }
#'
#' @param case_dir Character. Path to the OpenFOAM case directory. Must already
#'   contain \code{constant/triSurface/<stl_name>} (written by
#'   \code{prepare_openfoam_inputs()}).
#' @param stl_file Character. Full path to the building STL file. The file name
#'   is used as the \code{snappyHexMesh} geometry key; the patch name in the
#'   mesh is the file name without extension (default \code{"buildings"}).
#' @param domain Named list with \code{xmin}, \code{xmax}, \code{ymin},
#'   \code{ymax}, \code{zmin}, \code{zmax} in metres (local coordinates).
#'   Pass \code{foam_inputs$domain} directly.
#' @param inlet_velocity Numeric vector of length 3 (Ux, Uy, Uz) in m/s.
#'   Default \code{c(5, 0, 0)} - 5 m/s in the +x direction. The horizontal
#'   magnitude and direction set the ABL profile; Uz is forced to 0 in the
#'   profile. The inlet is always the x-min face; rotate the domain in
#'   \code{prepare_openfoam_inputs()} for other primary wind directions.
#' @param z_ref Numeric. Reference height for \code{inlet_velocity} (metres).
#'   Default 10.
#' @param z0 Numeric. Aerodynamic roughness length (metres) for the ground
#'   patch \code{nutkRoughWallFunction}. Use the domain-average value from the
#'   roughness raster (\code{mean(values(rast(foam_inputs$files$roughness_raster)),
#'   na.rm = TRUE)}), or a category value: open country 0.03, suburban 0.1,
#'   dense urban 0.5. Default 0.1.
#' @param base_cell_size Numeric. Background cell edge length in metres for
#'   \code{blockMesh}. Default 5. Building surfaces receive
#'   \code{building_refinement} additional levels of 2x refinement.
#' @param building_refinement Integer. Number of \code{snappyHexMesh}
#'   refinement levels at building surfaces. Default 2 -> effective building
#'   resolution = \code{base_cell_size / 4}.
#' @param n_iterations Integer. \code{simpleFoam} end time (= iteration count).
#'   Default 500.
#' @param consistent Logical. Use the SIMPLEC algorithm (`consistent yes` in
#'   the SIMPLE dict) instead of standard SIMPLE. SIMPLEC keeps the
#'   neighbour-coefficient term that SIMPLE drops from the pressure
#'   correction, so no pressure under-relaxation is needed and momentum can
#'   run near 0.9 instead of 0.5. Typically converges in 2-3x fewer
#'   iterations. Default `FALSE` (conservative standard SIMPLE).
#' @param potential_init Logical. Run `potentialFoam` before the main solver
#'   to seed a divergence-free initial velocity field. Costs seconds and
#'   commonly saves 30-50% of SIMPLE iterations by removing the violent
#'   startup transient. Non-fatal: if it fails the run proceeds from `0/`.
#'   Default `FALSE`.
#' @param residual_control Numeric. Convergence threshold applied to p, U, k
#'   and epsilon. SIMPLE stops as soon as all initial residuals fall below it.
#'   Default `1e-4`. Note that steady RANS residuals hit a floor when the real
#'   flow is unsteady (vortex shedding behind buildings), so many urban cases
#'   plateau between 1e-4 and 1e-3 and can never reach `1e-4`. When that
#'   happens the early exit never fires and the run always uses the full
#'   `n_iterations`; setting `5e-4` lets such a case stop cleanly.
#' @param max_cells Integer. Upper limit on total cells in
#'   \code{snappyHexMesh} (\code{maxGlobalCells}). Prevents OOM kills when
#'   Docker has limited RAM. Default 3 000 000 (safe at 4 GB Docker RAM;
#'   raise to 8 000 000 or more if Docker is set to 8+ GB).
#' @param overwrite Logical. If TRUE, overwrite existing case files. Default
#'   FALSE.
#' @param quiet Logical. Suppress progress messages. Default FALSE.
#'
#' @return Invisibly returns a named list of paths to the written files.
#'
#' @seealso \code{\link{prepare_openfoam_inputs}}, \code{\link{run_openfoam_docker}}
#'
#' @noRd
prepare_openfoam_case <- function(
    case_dir,
    stl_file,
    domain,
    inlet_velocity      = c(5, 0, 0),
    z_ref               = 10,
    z0                  = 0.1,
    base_cell_size      = 5,
    building_refinement = 2L,
    n_iterations        = 500L,
    max_cells           = 3000000L,
    consistent          = FALSE,
    potential_init      = FALSE,
    residual_control    = 1e-4,
    overwrite           = FALSE,
    quiet               = FALSE
) {
  # SUPERSEDED by prepare_foam_case(mode = "wind"), which uses the same
  # buoyant transient solver as the thermal case. simpleFoam is isothermal, so
  # this generator cannot produce air temperature or cool-air transport, and
  # its inlet is pinned to the x-min face: a non-x-aligned inlet_velocity
  # injects momentum through a face whose normal is perpendicular to it, so
  # effectively nothing enters. prepare_foam_case() assigns inlet/outlet roles
  # to four separately named lateral patches from the wind direction instead.
  if (!isTRUE(getOption("gloBFPr.quiet_deprecation", FALSE)))
    warning("`prepare_openfoam_case()` is superseded by ",
            "`prepare_foam_case(mode = \"wind\")`, which also supports oblique ",
            "wind directions. See ?prepare_foam_case.", call. = FALSE)

  # -- Validation ----------------------------------------------------------
  if (missing(case_dir) || !nzchar(case_dir))
    stop("`case_dir` must be provided.", call. = FALSE)
  if (missing(stl_file) || !file.exists(stl_file))
    stop("`stl_file` not found: ", stl_file, call. = FALSE)
  if (missing(domain) || !all(c("xmin","xmax","ymin","ymax","zmin","zmax")
                               %in% names(domain)))
    stop("`domain` must be a list with xmin/xmax/ymin/ymax/zmin/zmax.",
         call. = FALSE)
  if (length(inlet_velocity) != 3 || !is.numeric(inlet_velocity))
    stop("`inlet_velocity` must be a numeric vector of length 3.", call. = FALSE)

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

  # -- Derived quantities --------------------------------------------------
  u_horiz <- inlet_velocity[1:2]
  u_hmag  <- sqrt(sum(u_horiz^2))
  if (u_hmag < 1e-6)
    stop("Horizontal component of `inlet_velocity` is zero. ",
         "Provide a non-zero Ux or Uy.", call. = FALSE)
  u_mag  <- sqrt(sum(inlet_velocity^2))
  fd_x   <- u_horiz[1] / u_hmag         # normalized flow direction (x)
  fd_y   <- u_horiz[2] / u_hmag         # normalized flow direction (y)

  kappa  <- 0.41
  Cmu    <- 0.09
  u_star <- u_mag * kappa / log((z_ref + z0) / z0)
  k_init <- u_star^2 / sqrt(Cmu)
  e_init <- u_star^3 / (kappa * (z_ref + z0))

  # blockMesh cell counts
  nx <- max(1L, ceiling((domain$xmax - domain$xmin) / base_cell_size))
  ny <- max(1L, ceiling((domain$ymax - domain$ymin) / base_cell_size))
  nz <- max(1L, ceiling((domain$zmax - domain$zmin) / base_cell_size))

  # snappyHexMesh cell limits
  # Empirically, snappyHexMesh produces ~3-5x base cells for typical urban
  # domains at refinement 1-2. We cap at `max_cells` (default 3M) which is
  # safe for Docker at 4 GB RAM (peak usage during flood-fill approx. 2x cell count).
  # Users with 8 GB+ Docker can raise max_cells for finer resolution.
  base_cells       <- as.numeric(nx) * ny * nz
  max_global_cells <- min(ceiling(base_cells * 4), as.integer(max_cells))
  max_local_cells  <- ceiling(max_global_cells / 2)

  # locationInMesh: domain centre at 70 % height (above buildings)
  loc_x <- (domain$xmin + domain$xmax) / 2
  loc_y <- (domain$ymin + domain$ymax) / 2
  loc_z <- domain$zmax * 0.7

  stl_name   <- basename(stl_file)
  patch_name <- tools::file_path_sans_ext(stl_name)   # e.g. "buildings"
  ref        <- as.integer(building_refinement)

  msg(sprintf(
    "Domain: %.0f x %.0f x %.0f m  |  Mesh: %d x %d x %d cells (base %g m)",
    domain$xmax - domain$xmin, domain$ymax - domain$ymin,
    domain$zmax - domain$zmin, nx, ny, nz, base_cell_size
  ))
  # -- Roughness / mesh compatibility check --------------------------------
  # nutkRoughWallFunction is only valid when the first cell centre sits above
  # the equivalent sand-grain roughness height Ks ~ 20*z0. Near buildings the
  # cell size is base_cell_size / 2^refinement, so the first cell centre is at
  # roughly half of that. Violating this produces NaN in nut, which propagates
  # into the pressure equation and kills the solver with SIGFPE inside GAMG.
  near_cell   <- base_cell_size / (2^max(ref, 0L))
  first_centre <- near_cell / 2
  ks_est      <- 20 * z0
  z0_max_safe <- first_centre / 20

  if (ks_est > first_centre) {
    warning(sprintf(
      paste0(
        "z0 = %.3f m is too large for this mesh.\n",
        "  Near-building cells are %.2f m, so the first cell centre is ~%.2f m,\n",
        "  but the rough-wall function needs it above Ks = 20*z0 = %.2f m.\n",
        "  Capping z0 at %.4f m. To use the full roughness instead, coarsen the\n",
        "  mesh (larger base_cell_size) or lower building_refinement.\n",
        "  Note: when buildings are resolved as STL geometry, ground z0 should\n",
        "  represent only the paved/vegetated surface between them (~0.01-0.05 m),\n",
        "  not an unresolved urban canopy (~0.5 m)."
      ),
      z0, near_cell, first_centre, ks_est, z0_max_safe
    ), call. = FALSE)
    z0 <- z0_max_safe
  }

  msg(sprintf(
    "Wind: %.1f m/s @ (%.2f, %.2f, 0)  |  u* = %.3f m/s  |  z0 = %g m",
    u_mag, fd_x, fd_y, u_star, z0
  ))
  msg(sprintf("Initial k = %.4f m2/s2  |  epsilon = %.4e m2/s3", k_init, e_init))
  msg(sprintf(
    "Algorithm: %s  |  potentialFoam init: %s  |  residualControl = %g",
    if (isTRUE(consistent)) "SIMPLEC (consistent)" else "SIMPLE",
    if (isTRUE(potential_init)) "yes" else "no",
    residual_control
  ))

  # -- Write files ---------------------------------------------------------
  files <- list()

  msg("Writing Allrun ...")
  files$allrun <- write_foam_file(
    make_allrun(potential_init = potential_init),
    file.path(case_dir, "Allrun")
  )
  Sys.chmod(files$allrun, mode = "0755")

  msg("Writing system/ ...")
  files$decompose_par_dict <- write_foam_file(
    make_decompose_par_dict(1L),
    file.path(system_dir, "decomposeParDict")
  )
  files$block_mesh_dict <- write_foam_file(
    make_block_mesh_dict(domain, nx, ny, nz),
    file.path(system_dir, "blockMeshDict")
  )
  files$snappy_hex_mesh_dict <- write_foam_file(
    make_snappy_hex_mesh_dict(stl_name, loc_x, loc_y, loc_z, ref, ref,
                              max_local_cells  = max_local_cells,
                              max_global_cells = max_global_cells),
    file.path(system_dir, "snappyHexMeshDict")
  )
  files$control_dict <- write_foam_file(
    make_control_dict(n_iterations),
    file.path(system_dir, "controlDict")
  )
  files$fv_schemes <- write_foam_file(
    make_fv_schemes(),
    file.path(system_dir, "fvSchemes")
  )
  files$fv_solution <- write_foam_file(
    make_fv_solution(consistent = consistent,
                     residual_control = residual_control),
    file.path(system_dir, "fvSolution")
  )

  msg("Writing 0/ ...")
  files$u <- write_foam_file(
    make_u_field(u_mag, fd_x, fd_y, z_ref, z0, patch_name),
    file.path(ic_dir, "U")
  )
  files$p <- write_foam_file(
    make_p_field(patch_name),
    file.path(ic_dir, "p")
  )
  files$k <- write_foam_file(
    make_k_field(k_init, u_mag, fd_x, fd_y, z_ref, z0, patch_name),
    file.path(ic_dir, "k")
  )
  files$epsilon <- write_foam_file(
    make_epsilon_field(e_init, u_mag, fd_x, fd_y, z_ref, z0, patch_name),
    file.path(ic_dir, "epsilon")
  )
  files$nut <- write_foam_file(
    make_nut_field(z0, patch_name),
    file.path(ic_dir, "nut")
  )

  msg("Writing constant/ ...")
  files$turbulence_properties <- write_foam_file(
    make_turbulence_properties(),
    file.path(const_dir, "turbulenceProperties")
  )
  files$transport_properties <- write_foam_file(
    make_transport_properties(),
    file.path(const_dir, "transportProperties")
  )

  msg("Case files written to: ", normalizePath(case_dir, mustWork = FALSE))
  msg("Next step: run_openfoam_docker(\"", case_dir, "\")", sep = "")

  invisible(lapply(files, normalizePath, mustWork = FALSE))
}


#' Run an OpenFOAM case via Docker
#'
#' @description
#' Mounts the case directory into an OpenFOAM Docker container and executes
#' the \code{Allrun} script produced by \code{prepare_foam_case()}.
#' The generated wind case chains \code{blockMesh} -> \code{snappyHexMesh} ->
#' \code{buoyantBoussinesqPimpleFoam} and writes log files in
#' \code{case_dir}.
#'
#' Docker Desktop must be running and the \code{case_dir} path must be under
#' a directory shared with Docker (Docker Desktop -> Settings -> Resources ->
#' File Sharing).
#'
#' @param case_dir Character. Absolute path to the OpenFOAM case directory.
#'   Must contain an \code{Allrun} script (written by
#'   \code{prepare_foam_case()}).
#' @param image Character. Docker image tag. Default
#'   \code{"opencfd/openfoam-run:2506"}. Use the tag you pulled, e.g.
#'   \code{"opencfd/openfoam-run:2406"} for an older version.
#' @param ncpus Integer. Number of CPU cores to use.  Defaults to all
#'   physical cores detected on the host via
#'   \code{parallel::detectCores(logical = FALSE)}.  When \code{ncpus > 1}
#'   the solver runs in MPI parallel: \code{decomposePar} splits the mesh,
#'   \code{mpirun -np ncpus solver -parallel} runs it, and
#'   \code{reconstructPar} reassembles.  Typical speed-up is 3-6x on a
#'   4-core machine.  Pass \code{ncpus = 1} to force single-core.
#' @param wait Logical. If TRUE (default), R blocks until the simulation
#'   finishes. If FALSE, the container is launched in the background and the
#'   function returns immediately.
#' @param quiet Logical. Suppress messages. Default FALSE.
#'
#' @return Invisibly returns a list with \code{case_dir}, \code{image}, and
#'   the exit \code{status} (0 = success; only meaningful when
#'   \code{wait = TRUE}).
#'
#' @seealso \code{\link{prepare_foam_case}}
#'
#' @export
run_openfoam_docker <- function(
    case_dir,
    image  = "opencfd/openfoam-run:2506",
    ncpus  = foam_default_ncpus(),
    wait   = TRUE,
    quiet  = FALSE
) {
  case_dir <- normalizePath(case_dir, mustWork = TRUE)
  allrun   <- file.path(case_dir, "Allrun")

  if (!file.exists(allrun))
    stop("No Allrun script found in:\n  ", case_dir,
         "\nRun prepare_foam_case() first.", call. = FALSE)

  # On macOS, Docker Desktop only shares /Users, /Volumes, and /tmp by
  # default.  Paths under /var/folders (which is where tempdir() points)
  # are silently invisible inside the container, causing "Allrun: No such
  # file or directory" with exit code 127.
  if (.Platform$OS.type == "unix" && grepl("^/var/folders", case_dir))
    stop(
      "case_dir is inside /var/folders:\n  ", case_dir, "\n",
      "Docker Desktop on macOS cannot access this path by default.\n",
      "Move the case to a directory under your home folder, e.g.:\n",
      "  prepare_foam_case(case_dir = file.path(path.expand('~'), 'openfoam_demo'), ...)",
      call. = FALSE
    )

  # Verify Docker is reachable
  if (system("docker info > /dev/null 2>&1", ignore.stdout = TRUE,
             ignore.stderr = TRUE) != 0)
    stop("Docker does not appear to be running. ",
         "Start Docker Desktop and try again.", call. = FALSE)

  # Use absolute in-container path (/case/Allrun) rather than relying on the
  # working-directory lookup.  On macOS with Docker Desktop + VirtioFS, `bash
  # Allrun` (relative) sometimes fails with exit 127 even when the volume is
  # mounted correctly, because the CWD is set after the shell is started.
  # Using the absolute path is always unambiguous.
  ncpus <- max(1L, as.integer(ncpus))
  vol_flag <- paste0("-v ", shQuote(case_dir), ":/case")
  cmd <- paste(
    "docker run --rm",
    vol_flag,
    "-w /case",
    sprintf("--cpus %d", ncpus),
    sprintf("-e NPROC=%d", ncpus),
    shQuote(image),
    "bash /case/Allrun"
  )

  # Pre-flight: verify that the container can see at least one file in /case.
  # A silent empty mount (VirtioFS sharing not enabled) would give a confusing
  # "No such file" error from bash rather than a clear mount failure.
  probe_cmd <- paste(
    "docker run --rm",
    vol_flag,
    shQuote(image),
    "ls /case/Allrun"
  )
  probe_status <- system(probe_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (probe_status != 0)
    stop(
      "Docker cannot see the case files at:\n  ", case_dir, "\n",
      "The volume mount appears to be empty inside the container.\n",
      "On macOS, open Docker Desktop -> Settings -> Resources -> File Sharing\n",
      "and add the directory (or a parent of it) to the shared paths list.\n",
      "Then click 'Apply & Restart' and retry.",
      call. = FALSE
    )

  if (!isTRUE(quiet)) {
    message("OpenFOAM Docker run")
    message("  Image     : ", image)
    message("  Case      : ", case_dir)
    message("  Log files : ", case_dir, "/log.*")
    if (!wait) message("  (running in background - R will not block)")
  }

  status <- system(cmd, wait = wait)

  if (wait && status != 0)
    warning(
      "Docker exited with status ", status, ".\n",
      "Check log files for errors:\n",
      "  log.blockMesh, log.snappyHexMesh, log.decomposeMesh,\n",
      "  log.setFields, log.decompose, log.solver\n",
      "The FOAM FATAL error is usually the LAST few lines of the newest log.",
      call. = FALSE
    )

  invisible(list(case_dir = case_dir, image = image, status = status))
}


# ===========================================================================
# Post-processing helpers
# ===========================================================================

#' Sample OpenFOAM results at pedestrian level and return a raster
#'
#' @description
#' Uses OpenFOAM's \code{postProcess -func surfaces} utility to cut a
#' horizontal plane at z = 1.5 m (or any height) through the latest
#' time-step, then reads the resulting \code{.raw} files into R and
#' returns a multi-layer \code{SpatRaster}.
#'
#' For the current wind workflow, \code{prepare_foam_case()} already writes a
#' pedestrian slice during the solver run via the \code{pedestrianSlice}
#' function object; use \code{\link{read_foam_pedestrian_slice}} for that
#' default 1.5 m output. Use this helper to sample other heights or fields
#' after a run.
#'
#' @param case_dir Character. OpenFOAM case directory.
#' @param fields Character vector of field names to sample.
#'   Default \code{c("U", "p")}.
#' @param z Numeric. Height above ground in metres.  Default 1.5.
#' @param image Character. Docker image tag.
#'   Default \code{"opencfd/openfoam-run:2506"}.
#' @param resolution Numeric. Output raster cell size in metres.  Default 5.
#' @param time_step Character or numeric. \code{"latestTime"} (default) or a
#'   specific time-step number.
#' @param interpolate Character. Post-rasterization smoothing method.
#'   \code{"smooth"} (default) fills NA gaps then applies a Gaussian focal
#'   filter to all cells, producing continuous gradients similar to ParaView;
#'   \code{"focal"} fills NA gaps only (faster, but blocky in open areas);
#'   \code{"none"} returns the raw rasterization without any smoothing.
#' @param idw_power Numeric. Controls the Gaussian sigma when
#'   \code{interpolate = "smooth"}: sigma = \code{resolution * idw_power}
#'   metres. Default 2 (10 m sigma at 5 m resolution). Increase for a wider
#'   blur; decrease for sharper transitions near buildings.
#' @param idw_maxdist Ignored (reserved for future use).
#' @param buildings Optional \code{sf} polygon layer of building footprints in
#'   the domain's local coordinate system. When \code{NULL} (default) the
#'   footprints saved by \code{\link{prepare_openfoam_inputs}} are loaded
#'   automatically from
#'   \code{<case_dir>/constant/gloBFPr/metadata/buildings_openfoam.rds}.
#'   Building interiors are masked to \code{NA}, and the layer is attached to
#'   the returned raster so \code{\link{plot_foam_map}} overlays it without
#'   any extra argument. Pass \code{buildings = sf::st_sf(...)} to override,
#'   or note that masking is skipped when no footprints can be found.
#' @param quiet Logical.  Default FALSE.
#'
#' @return A \code{terra::SpatRaster} with layers named after the sampled
#'   fields.  Velocity \code{U} produces three extra layers:
#'   \code{U_mag} (speed), \code{Ux}, \code{Uy}.  Building footprints, when
#'   available, are attached as \code{attr(x, "buildings")}.
#'
#' @seealso \code{\link{plot_foam_map}}, \code{\link{read_foam_pedestrian_slice}}
#' @export
sample_foam_slice <- function(
    case_dir,
    fields     = c("U", "p"),
    z          = 1.5,
    image      = "opencfd/openfoam-run:2506",
    resolution = 5,
    time_step  = "latestTime",
    interpolate = c("smooth", "focal", "none"),
    idw_power   = 2,
    idw_maxdist = NULL,
    buildings   = NULL,
    quiet       = FALSE
) {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required.", call. = FALSE)

  interpolate <- match.arg(interpolate)
  case_dir <- normalizePath(case_dir, mustWork = TRUE)
  msg <- function(...) if (!isTRUE(quiet)) message(...)

  # Build surface sampling dictionary inline
  surf_name <- paste0("z", gsub("\\.", "p", as.character(z)), "m")
  fields_str <- paste(fields, collapse = " ")

  dict <- paste0(
    foam_header("dictionary", "pedestrianSliceDict"),
    "type            surfaces;\n",
    "libs            (\"libsampling.so\");\n\n",
    "executeControl  onEnd;\n",
    "writeControl    onEnd;\n\n",
    "surfaceFormat   raw;\n",
    "fields          ( ", fields_str, " );\n\n",
    "surfaces\n{\n",
    "    ", surf_name, "\n    {\n",
    "        type        cuttingPlane;\n",
    "        planeType   pointAndNormal;\n",
    "        pointAndNormalDict\n        {\n",
    sprintf("            point  (0 0 %g);\n", z),
    "            normal (0 0 1);\n",
    "        }\n",
    "        interpolate true;\n",
    "    }\n",
    "}\n"
  )

  dict_path <- file.path(case_dir, "system", "pedestrianSliceDict")
  writeLines(dict, dict_path)
  msg("Wrote ", dict_path)

  # Run postProcess inside Docker
  ts_flag <- if (identical(time_step, "latestTime")) "-latestTime" else
    paste("-time", time_step)

  cmd <- paste(
    "docker run --rm",
    paste0("-v ", shQuote(case_dir), ":/case"),
    "-w /case",
    shQuote(image),
    paste0("bash -c 'postProcess -case /case -func pedestrianSliceDict ",
           ts_flag,
           " > /case/log.postProcess 2>&1'")
  )
  msg("Running postProcess ...")
  exit <- system(cmd, wait = TRUE)
  if (exit != 0) {
    log <- tryCatch(
      utils::tail(readLines(file.path(case_dir, "log.postProcess")), 20),
      error = function(e) character(0)
    )
    stop("postProcess failed.\n", paste(log, collapse = "\n"), call. = FALSE)
  }

  # Find the output directory
  pp_base <- file.path(case_dir, "postProcessing", "pedestrianSliceDict")
  if (!dir.exists(pp_base))
    stop("postProcessing output not found at:\n  ", pp_base, call. = FALSE)

  time_dirs <- sort(as.numeric(
    list.dirs(pp_base, full.names = FALSE, recursive = FALSE)
  ), na.last = NA)
  if (length(time_dirs) == 0)
    stop("No time directories in ", pp_base, call. = FALSE)

  out_dir <- file.path(pp_base, max(time_dirs))
  msg("Reading from: ", out_dir)

  # Discover all .raw files written under the time directory (recursive,
  # handles both flat layout and v2206+ subdirectory layout):
  #   flat:   <time>/<surfName>_<field>.raw
  #   nested: <time>/surfaces/<surfName>/<field>.raw   (OpenFOAM v2206+)
  all_raw <- list.files(out_dir, pattern = "\\.raw$",
                        recursive = TRUE, full.names = TRUE)
  msg("Found raw files:\n  ", paste(basename(all_raw), collapse = "\n  "))

  # Locate a .raw file for a given field, trying all known naming conventions:
  #   v2212- (old): <surfName>_<field>.raw
  #   v2506  (new): <field>_<surfName>.raw
  #   v2206+ nested: surfaces/<surfName>/<field>.raw
  find_raw <- function(fld) {
    # old flat: <surfName>_<field>.raw
    flat_old <- file.path(out_dir, paste0(surf_name, "_", fld, ".raw"))
    if (file.exists(flat_old)) return(flat_old)
    # new flat: <field>_<surfName>.raw  (OpenFOAM v2506+)
    flat_new <- file.path(out_dir, paste0(fld, "_", surf_name, ".raw"))
    if (file.exists(flat_new)) return(flat_new)
    # nested: .../surfaces/<surfName>/<field>.raw
    nested <- file.path(out_dir, "surfaces", surf_name,
                        paste0(fld, ".raw"))
    if (file.exists(nested)) return(nested)
    # nested alt: .../surfaces/<surfName>_<field>.raw
    nested2 <- file.path(out_dir, "surfaces",
                         paste0(surf_name, "_", fld, ".raw"))
    if (file.exists(nested2)) return(nested2)
    # nested alt: .../surfaces/<field>_<surfName>.raw
    nested3 <- file.path(out_dir, "surfaces",
                         paste0(fld, "_", surf_name, ".raw"))
    if (file.exists(nested3)) return(nested3)
    # Last resort: search by field name in all discovered files
    hits <- all_raw[grepl(paste0("(^|[/_])", fld, "([._]|$)"),
                          basename(all_raw))]
    if (length(hits) == 1L) return(hits)
    if (length(hits) > 1L) {
      msg("Multiple candidates for field '", fld, "': ",
          paste(hits, collapse = ", "), ". Using first.")
      return(hits[1])
    }
    stop(sprintf(
      paste0("Cannot locate .raw file for field '%s'.\n",
             "Searched under: %s\n",
             "Files found:\n  %s\n",
             "If the surface name differs, check the postProcessing folder ",
             "and set the 'z' argument to match the exact plane used."),
      fld, out_dir,
      if (length(all_raw)) paste(all_raw, collapse = "\n  ") else "(none)"
    ), call. = FALSE)
  }

  # Helper: parse a .raw file (comment lines start with #)
  read_raw <- function(path, col_names) {
    lines <- readLines(path)
    data_lines <- lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]
    utils::read.table(text = paste(data_lines, collapse = "\n"),
                      header = FALSE, col.names = col_names)
  }

  # Helper: build raster template from a data frame with x/y columns
  make_tmpl <- function(dat) {
    terra::rast(
      xmin = min(dat$x), xmax = max(dat$x),
      ymin = min(dat$y), ymax = max(dat$y),
      resolution = resolution, crs = ""
    )
  }

  # Helper: rasterize one scalar column; apply chosen interpolation.
  #
  #  "smooth" (default):
  #    1. Rasterize at the requested resolution.
  #    2. Fill NA cells caused by the coarse blockMesh background using an
  #       expanding focal mean (runs until no NAs remain).
  #    3. Apply a Gaussian focal smooth to ALL cells. This is what produces the
  #       continuous gradient ParaView shows - it is equivalent to bilinear
  #       interpolation across the underlying mesh faces.
  #  "focal":  steps 1-2 only (faster, but blocky in open areas).
  #  "none":   step 1 only (raw rasterization; useful for debugging).
  rasterize_col <- function(dat, col, tmpl) {
    pts <- terra::vect(dat, geom = c("x", "y"), crs = "")
    r   <- terra::rasterize(pts, tmpl, field = col, fun = "mean")

    if (interpolate %in% c("smooth", "focal")) {
      # --- Step 2: fill NAs from coarse background cells ---
      for (i in seq_len(40)) {
        n_before <- terra::global(is.na(r), "sum")[[1]]
        if (n_before == 0) break
        r <- terra::focal(r, w = 7, fun = "mean", na.rm = TRUE,
                          na.policy = "only")
        if (terra::global(is.na(r), "sum")[[1]] == n_before) break
      }
    }

    if (interpolate == "smooth") {
      # --- Step 3: Gaussian smooth across all cells ---
      # sigma = idw_power grid cells (repurposed as sigma multiplier)
      sigma_m <- resolution * max(idw_power, 1)
      gw <- terra::focalMat(r, d = sigma_m, type = "Gauss")
      r  <- terra::focal(r, w = gw, na.rm = TRUE)
    }

    names(r) <- col
    r
  }

  layers <- list()

  for (fld in fields) {
    if (fld == "U") {
      ufile <- find_raw("U")
      udat  <- read_raw(ufile, c("x", "y", "z", "Ux", "Uy", "Uz"))
      udat$U_mag <- sqrt(udat$Ux^2 + udat$Uy^2 + udat$Uz^2)
      tmpl_u <- make_tmpl(udat)
      msg("Interpolating U fields (", interpolate, ") ...")
      layers[["U_mag"]] <- rasterize_col(udat, "U_mag", tmpl_u)
      layers[["Ux"]]    <- rasterize_col(udat, "Ux",    tmpl_u)
      layers[["Uy"]]    <- rasterize_col(udat, "Uy",    tmpl_u)
    } else {
      sfile <- find_raw(fld)
      sdat  <- read_raw(sfile, c("x", "y", "z", fld))
      tmpl_s <- make_tmpl(sdat)
      msg("Interpolating ", fld, " (", interpolate, ") ...")
      layers[[fld]] <- rasterize_col(sdat, fld, tmpl_s)
    }
  }

  result <- terra::rast(layers)

  # -- Building footprints -------------------------------------------------
  # Auto-load from the case directory unless the caller supplied them, so
  # plot_foam_map() can overlay buildings without any extra arguments.
  if (is.null(buildings)) {
    bpath <- file.path(case_dir, "constant", "gloBFPr", "metadata",
                       "buildings_openfoam.rds")
    if (file.exists(bpath)) {
      buildings <- tryCatch(readRDS(bpath), error = function(e) NULL)
      if (!is.null(buildings))
        msg("Auto-loaded ", nrow(buildings), " buildings from case directory.")
    }
  }

  if (!is.null(buildings)) {
    # Attach as an attribute so plot_foam_map() finds them automatically.
    attr(result, "buildings") <- buildings

    # Mask building interiors to NA. The NA-fill and Gaussian smoothing above
    # flood over footprints narrower than the fill window, so without this
    # only the largest buildings appear as voids.
    if (!requireNamespace("sf", quietly = TRUE)) {
      warning("Package 'sf' is required to mask buildings; skipping mask.",
              call. = FALSE)
    } else {
      masked <- tryCatch({
        bvect      <- terra::vect(sf::st_geometry(buildings))
        building_r <- terra::rasterize(bvect, result[[1L]], field = 1L,
                                       background = NA)
        terra::mask(result, building_r, maskvalues = 1L, updatevalue = NA)
      }, error = function(e) {
        warning("Could not mask building footprints: ", conditionMessage(e),
                call. = FALSE)
        NULL
      })
      if (!is.null(masked)) {
        attr(masked, "buildings") <- buildings   # mask() drops attributes
        result <- masked
      }
    }
  }

  msg(sprintf("Raster: %d x %d cells  |  layers: %s",
              terra::nrow(result), terra::ncol(result),
              paste(names(result), collapse = ", ")))
  result
}


#' Plot an OpenFOAM pedestrian-level map
#'
#' @description
#' Convenience wrapper around \code{ggplot2} for visualising a single layer
#' from the raster returned by \code{\link{sample_foam_slice}} or
#' \code{\link{read_foam_pedestrian_slice}}.
#'
#' @param r A \code{SpatRaster} (single layer, or one layer will be selected
#'   via \code{layer}).
#' @param layer Character. Layer name to plot.  Default \code{"U_mag"}.
#' @param canopy Optional data frame of canopy cell centres (\code{x}, \code{y})
#'   with a \code{"res"} attribute giving the cell size.  Taken from the
#'   raster's \code{"canopy"} attribute when not supplied, which
#'   \code{\link{read_foam_pedestrian_slice}} attaches whenever the case has a
#'   canopy height raster.
#' @param title Character.  Plot title.  Default auto-generated.
#' @param palette Character.  \code{hcl.colors} palette name.
#'   Default \code{"YlOrRd"}.
#' @param reverse Logical.  Reverse palette direction.  Default FALSE.
#' @param buildings Optional \code{sf} object of building footprints to
#'   overlay (in the same local coordinate system as the raster).
#' @param legend_title Character.  Legend label.  Default \code{layer}.
#' @param max_u_ref Numeric.  If plotting \code{U_mag}, annotate the colour
#'   scale as a wind speed ratio by dividing by this reference speed.
#'   Default \code{NULL} (no ratio).
#' @param na_colour Colour for cells with no result.  Inside the mapped area a
#'   void is an obstruction the flow never entered, so the default is the same
#'   solid black used for buildings and canopy.
#' @param scalebar Logical. Draw a distance scale bar in the lower-left corner.
#'   Defaults to \code{TRUE}. Distances assume projected map units are metres;
#'   for geographic rasters they are approximated at the map's mid-latitude.
#' @param scalebar_unit Scale bar unit: \code{"km"}, \code{"m"}, or
#'   \code{"auto"} to pick whichever keeps the label readable. Default
#'   \code{"auto"}.
#' @param scalebar_cex Scale bar label size multiplier. Defaults to \code{0.7}.
#' @param north_arrow Logical. Draw a north arrow above the scale bar. Defaults
#'   to \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_foam_map <- function(
    r,
    layer          = "U_mag",
    title          = NULL,
    palette        = "YlOrRd",
    reverse        = FALSE,
    buildings      = NULL,
    canopy         = NULL,
    legend_title   = layer,
    max_u_ref      = NULL,
    na_colour      = "black",
    scalebar       = TRUE,
    scalebar_unit  = c("auto", "km", "m"),
    scalebar_cex   = 0.7,
    north_arrow    = TRUE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.", call. = FALSE)
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required.", call. = FALSE)
  scalebar_unit <- match.arg(scalebar_unit)

  # Keep a reference to the original raster before subsetting so we can
  # retrieve the buildings attribute (which terra drops on [[]] subsetting).
  r_orig <- r

  # Select layer
  if (terra::nlyr(r) > 1) {
    if (!layer %in% names(r))
      stop("Layer '", layer, "' not found. Available: ",
           paste(names(r), collapse = ", "), call. = FALSE)
    r <- r[[layer]]
  }

  # Convert raster to data frame.
  #
  # na.rm = FALSE is deliberate.  terra drops NA cells by default, which does
  # NOT leave them unpainted in `na_colour` - it leaves no tile at all, so the
  # panel background shows through and the map appears to have holes punched in
  # it.  Every void inside the mapped area is an obstruction (a building, or a
  # structure that is in the DSM but missing from the footprint layer), so it
  # should read as solid, not as a gap in the figure.  Keeping the NA rows is
  # what makes `na.value = na_colour` below do anything at all.
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- "value"

  # Wind speed ratio annotation
  if (!is.null(max_u_ref) && layer == "U_mag") {
    df$value <- df$value / max_u_ref
    if (is.null(legend_title)) legend_title <- "U / U_ref"
  }

  if (is.null(title)) {
    title <- switch(layer,
      U_mag       = "Wind speed at pedestrian level (1.5 m AGL)",
      p           = "Pressure at pedestrian level (1.5 m AGL)",
      T_air       = "Air temperature at pedestrian level (1.5 m AGL)",
      T_cool      = "Cooling below ambient at pedestrian level (1.5 m AGL)",
      T_cool_flux = "Cool-air transport flux at pedestrian level (1.5 m AGL)",
      layer
    )
  }

  cols <- if (reverse) rev(hcl.colors(100, palette)) else hcl.colors(100, palette)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colours = cols, name = legend_title,
                                  na.value = na_colour) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = NULL,
                  caption = title,
                  x = "Easting (m, local)",
                  y = "Northing (m, local)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = ggplot2::rel(1.15)
      )
    )

  # Fall back to buildings stored as an attribute by read_foam_pedestrian_slice
  # (or sample_foam_slice) when the caller did not pass them explicitly.
  if (is.null(buildings))
    buildings <- attr(r_orig, "buildings")

  # Canopy first, so buildings draw on top of it where they overlap.
  if (is.null(canopy)) canopy <- attr(r_orig, "canopy")
  if (!is.null(canopy) && nrow(canopy) > 0) {
    cres <- attr(canopy, "res")
    if (is.null(cres) || !is.finite(cres)) cres <- 5
    p <- p + ggplot2::geom_tile(
      data        = canopy,
      ggplot2::aes(x = .data$x, y = .data$y),
      width       = cres, height = cres,
      fill        = "black",
      inherit.aes = FALSE
    )
  }

  if (!is.null(buildings)) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      warning("Package 'sf' needed to overlay buildings.", call. = FALSE)
    } else {
      # Convert to plain data frame of polygon coordinates so that ggplot2
      # never invokes CRS transformation (the local-metre system has no EPSG).
      bpoly <- sf::st_cast(
        sf::st_make_valid(buildings), "MULTIPOLYGON",
        warn = FALSE
      )
      bpoly <- suppressWarnings(sf::st_cast(bpoly, "POLYGON"))
      bdf   <- do.call(rbind, lapply(seq_len(nrow(bpoly)), function(i) {
        coords <- as.data.frame(sf::st_coordinates(bpoly[i, ])[, 1:2])
        names(coords) <- c("x", "y")
        coords$group  <- i
        coords
      }))
      p <- p + ggplot2::geom_polygon(
        data        = bdf,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
        fill        = "black",
        colour      = "black",
        linewidth   = 0.1,
        inherit.aes = FALSE
      )
    }
  }

  foam_map_add_orientation(
    p,
    r,
    scalebar = scalebar,
    scalebar_unit = scalebar_unit,
    scalebar_cex = scalebar_cex,
    north_arrow = north_arrow
  )
}

#' @noRd
foam_scalebar_nice <- function(x) {
  if (!is.finite(x) || x <= 0) return(NA_real_)
  pow <- 10^floor(log10(x))
  frac <- x / pow
  nice <- if (frac >= 5) 5 else if (frac >= 2.5) 2.5 else if (frac >= 2) 2 else 1
  nice * pow
}

#' @noRd
foam_map_scale_context <- function(r) {
  ext <- terra::ext(r)
  xmin <- ext[1]
  xmax <- ext[2]
  ymin <- ext[3]
  ymax <- ext[4]
  w <- xmax - xmin
  h <- ymax - ymin
  if (!all(is.finite(c(xmin, xmax, ymin, ymax, w, h))) || w <= 0 || h <= 0) {
    return(NULL)
  }

  longlat <- suppressWarnings(try(terra::is.lonlat(r), silent = TRUE))
  longlat <- isTRUE(!inherits(longlat, "try-error") && isTRUE(longlat))
  metres_per_x <- 1
  if (longlat) {
    mid_lat <- (ymin + ymax) / 2
    metres_per_x <- 111320 * cos(mid_lat * pi / 180)
  }
  if (!is.finite(metres_per_x) || metres_per_x <= 0) return(NULL)

  list(
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    w = w, h = h, metres_per_x = metres_per_x
  )
}

#' @noRd
foam_map_add_orientation <- function(p,
                                     r,
                                     scalebar = TRUE,
                                     scalebar_unit = "auto",
                                     scalebar_cex = 0.7,
                                     north_arrow = TRUE) {
  if (!isTRUE(scalebar) && !isTRUE(north_arrow)) return(p)
  ctx <- foam_map_scale_context(r)
  if (is.null(ctx)) return(p)

  x0 <- ctx$xmin + 0.045 * ctx$w
  y0 <- ctx$ymin + 0.085 * ctx$h
  bar_h <- 0.012 * ctx$h
  text_y <- y0 - 0.026 * ctx$h
  label_size <- 3.1 * scalebar_cex
  col <- "#222222"
  scalebar_drawn <- FALSE
  arrow_h <- 0.075 * ctx$h
  arrow_x <- x0 + 0.025 * ctx$w
  bar_x0 <- if (isTRUE(north_arrow)) x0 + 0.075 * ctx$w else x0

  if (isTRUE(scalebar)) {
    metres <- foam_scalebar_nice(0.24 * ctx$w * ctx$metres_per_x)
    if (is.finite(metres) && metres > 0) {
      unit <- scalebar_unit
      if (identical(unit, "auto")) unit <- if (metres >= 1000) "km" else "m"
      divisor <- if (identical(unit, "km")) 1000 else 1
      value <- metres / divisor
      label <- paste(format(value, trim = TRUE, scientific = FALSE), unit)
      bar_w <- metres / ctx$metres_per_x
      bg_pad_x <- 0.018 * ctx$w
      bg_top <- y0 + max(bar_h, arrow_h) + 0.065 * ctx$h
      scalebar_drawn <- TRUE

      p <- p +
        ggplot2::annotate(
          "rect",
          xmin = x0 - bg_pad_x,
          xmax = bar_x0 + bar_w + bg_pad_x,
          ymin = text_y - 0.035 * ctx$h,
          ymax = bg_top,
          fill = grDevices::adjustcolor("white", alpha.f = 0.82),
          colour = NA
        ) +
        ggplot2::annotate(
          "rect",
          xmin = bar_x0,
          xmax = bar_x0 + bar_w / 2,
          ymin = y0,
          ymax = y0 + bar_h,
          fill = col,
          colour = col,
          linewidth = 0.2
        ) +
        ggplot2::annotate(
          "rect",
          xmin = bar_x0 + bar_w / 2,
          xmax = bar_x0 + bar_w,
          ymin = y0,
          ymax = y0 + bar_h,
          fill = "white",
          colour = col,
          linewidth = 0.2
        ) +
        ggplot2::annotate(
          "text",
          x = bar_x0,
          y = text_y,
          label = "0",
          hjust = 0.5,
          vjust = 1,
          size = label_size,
          colour = col
        ) +
        ggplot2::annotate(
          "text",
          x = bar_x0 + bar_w,
          y = text_y,
          label = label,
          hjust = 0.5,
          vjust = 1,
          size = label_size,
          colour = col
        )
    }
  }

  if (isTRUE(north_arrow)) {
    arrow_y0 <- y0
    arrow_y1 <- y0 + arrow_h
    if (!scalebar_drawn) {
      p <- p +
        ggplot2::annotate(
          "rect",
          xmin = arrow_x - 0.035 * ctx$w,
          xmax = arrow_x + 0.035 * ctx$w,
          ymin = arrow_y0 - 0.025 * ctx$h,
          ymax = arrow_y1 + 0.055 * ctx$h,
          fill = grDevices::adjustcolor("white", alpha.f = 0.82),
          colour = NA
        )
    }
    p <- p +
      ggplot2::annotate(
        "segment",
        x = arrow_x,
        xend = arrow_x,
        y = arrow_y0,
        yend = arrow_y1,
        linewidth = 0.55,
        colour = col,
        arrow = grid::arrow(
          length = grid::unit(0.10, "inches"),
          type = "closed"
        )
      ) +
      ggplot2::annotate(
        "text",
        x = arrow_x,
        y = arrow_y1 + 0.010 * ctx$h,
        label = "N",
        fontface = "bold",
        size = 3.7 * scalebar_cex,
        colour = col
      )
  }

  p
}


#' Add wind / flow vector arrows to a foam map plot
#'
#' @description
#' Takes a \code{ggplot} object produced by \code{\link{plot_foam_map}} and
#' overlays velocity arrows sampled on a regular sub-grid.
#'
#' @param p A \code{ggplot} object (output of \code{plot_foam_map}).
#' @param r A \code{SpatRaster} that contains layers named \code{Ux} and
#'   \code{Uy}.
#' @param spacing Numeric. Arrow sub-grid spacing in the same units as the
#'   raster coordinates (usually metres).  Default 20.
#' @param scale Numeric. Arrow length multiplier.  Default 1.
#' @param colour Character. Arrow colour.  Default \code{"black"}.
#' @param alpha Numeric.  Arrow opacity (0-1).  Default 0.6.
#' @param linewidth Numeric.  Line thickness in mm.  Default 0.25 (thin).
#' @param arrow_size Numeric.  Arrowhead length in cm.  Default 0.07 (tiny).
#' @param u_ref Numeric.  Speed (m/s) the longest arrow represents.  The
#'   default, \code{NULL}, normalises to the field's own maximum, which makes
#'   a near-stagnant field look as vigorous as a strong one; pass a fixed value
#'   to put several maps on one scale.
#' @param quiet Logical.  Suppress the message reporting the arrow scale.
#'
#' @return The same \code{ggplot} object with arrows added.
#' @export
add_flow_vectors <- function(
    p,
    r,
    spacing    = 20,
    scale      = 1,
    colour     = "black",
    alpha      = 0.6,
    linewidth  = 0.25,
    arrow_size = 0.07,
    u_ref      = NULL,
    quiet      = FALSE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.", call. = FALSE)
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required.", call. = FALSE)
  if (!all(c("Ux", "Uy") %in% names(r)))
    stop("Raster must have layers 'Ux' and 'Uy'.", call. = FALSE)

  # Sub-sample grid
  ext   <- terra::ext(r)
  xs    <- seq(ext[1] + spacing / 2, ext[2] - spacing / 2, by = spacing)
  ys    <- seq(ext[3] + spacing / 2, ext[4] - spacing / 2, by = spacing)
  grid  <- expand.grid(x = xs, y = ys)
  pts   <- terra::vect(grid, geom = c("x", "y"), crs = terra::crs(r))
  vals  <- terra::extract(r[[c("Ux", "Uy")]], pts, ID = FALSE)
  df    <- cbind(grid, vals)
  df    <- df[!is.na(df$Ux) & !is.na(df$Uy), ]

  # Normalise arrow length so the longest arrow = spacing * scale * 0.8.
  #
  # Normalising to the field's own maximum means a stagnant field is drawn at
  # exactly the same visual weight as a strong one: a 0.1 m/s drift with no
  # preferred direction renders as a dense mat of full-length arrows and reads
  # as vigorous, chaotic wind. It cost a real debugging session. Report the
  # scale, and let the caller fix it with `u_ref` so several maps share one.
  mag     <- sqrt(df$Ux^2 + df$Uy^2)
  mag_max <- max(mag, na.rm = TRUE)
  if (mag_max == 0) return(p)
  ref <- if (is.null(u_ref)) mag_max else u_ref
  if (!isTRUE(quiet))
    message(sprintf("Vectors: longest arrow = %.3g m/s%s (field max %.3g, mean %.3g)",
                    ref, if (is.null(u_ref)) "" else " (u_ref)",
                    mag_max, mean(mag, na.rm = TRUE)))
  fac <- spacing * scale * 0.8 / ref
  df$xend <- df$x + df$Ux * fac
  df$yend <- df$y + df$Uy * fac

  p + ggplot2::geom_segment(
    data        = df,
    ggplot2::aes(
      x = .data$x, y = .data$y,
      xend = .data$xend, yend = .data$yend
    ),
    arrow       = ggplot2::arrow(
                    length = ggplot2::unit(arrow_size, "cm"),
                    type   = "open",
                    angle  = 20
                  ),
    colour      = colour,
    alpha       = alpha,
    linewidth   = linewidth,
    inherit.aes = FALSE
  )
}
