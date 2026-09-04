#' get_3d_world
#' @description
#' Export a study area as a 3D scene (inspired by the
#' \href{https://github.com/louis-e/arnis}{arnis} project): terrain surface,
#' extruded building footprints, and canopy tree objects, written as Wavefront
#' OBJ and/or binary STL files that load directly in Rhino3D and Blender.
#'
#' The OBJ output contains named objects per layer (`terrain`, `buildings`,
#' `canopy_trunks`, `canopy_crowns`) and one group per building
#' (`building_<id>`), so individual buildings remain selectable after import.
#' Ground surface classes (greenspace, roads, sidewalks, paths, water, sand)
#' are part of the terrain object itself, as face groups with their own
#' materials - the terrain is one continuous, detailed surface with
#' differently colored regions, not a stack of separate ribbons.
#' STL has no object concept, so one STL file is written per layer.
#'
#' @param x sf. Building footprint polygons, typically the output of
#' [search_3dglobdf()] (must contain the column given by `height_col`). If
#' `NULL`, buildings are fetched via [search_3dglobdf()] using `bbox`/`place`.
#' @param bbox `sf`, `sfc`, or numeric vector (xmin, ymin, xmax, ymax) in
#' WGS84. Ignored when `x` is provided.
#' @param place character (optional). Address or place name, passed to
#' [search_3dglobdf()]. Ignored when `x` is provided.
#' @param terrain logical. If `TRUE` (default), download a DEM (OpenTopography
#' via `dsmSearch`; requires `key`) and mesh the ground surface. Buildings and
#' trees are then placed at their sampled ground elevation. If `FALSE`, the
#' scene has a flat ground at z = 0 and no API key is needed.
#' @param canopy character or `NULL`. Canopy height source, `"metachm"`
#' (default) or `"ethchm"`; `NULL` skips trees.
#' @param key character. OpenTopography API key, required when
#' `terrain = TRUE` and no `dem` raster is supplied.
#' @param dem `SpatRaster` or `PackedSpatRaster` (optional). A pre-loaded
#' digital elevation model (e.g. `globfp_example_dem`). When supplied with
#' `terrain = TRUE`, no download or API key is needed. Ignored when
#' `terrain = FALSE`.
#' @param canopy_height `SpatRaster` or `PackedSpatRaster` (optional). A
#' pre-loaded canopy height model (e.g. `globfp_example_canopy_height`).
#' When supplied, it is used for tree detection instead of downloading from
#' the `canopy` source.
#' @param height_col character. Building height column. Default `"Height"`.
#' @param format character. Any of `"obj"`, `"stl"`. Default writes both.
#' @param out_dir character. Output directory (created if missing).
#' @param min_tree_height numeric. Minimum canopy height in metres. Default 2.
#' @param tree_window numeric. Local-maximum search window in **metres** for
#' tree detection: at most one tree is placed per window. Larger values give
#' fewer, more widely spaced trees; lower it (e.g. 3) if dense canopy looks
#' too sparse in the model. Default 5.
#' @param max_trees integer. Cap on tree objects; the tallest trees are kept
#' and a warning reports how many were detected. Default 20000.
#' @param simplify_terrain integer. Aggregation factor for the terrain mesh
#' (1 = full DEM resolution). Default 1.
#' @param color_by character or `NULL`. Name of a numeric column in `x` mapped
#' to per-building OBJ materials via a viridis ramp (OBJ only; STL carries no
#' color). Default `NULL`.
#' @param facade_palette logical or character vector. If `TRUE`, buildings get
#' varied muted facade colors (arnis-style) instead of uniform grey, assigned
#' deterministically per building. Supply a character vector of R colors to
#' use a custom palette. Ignored when `color_by` is set. Default `FALSE`.
#' @param roads `NULL`, `"overture"`, or an sf line layer. Paints road and
#' sidewalk surfaces into the terrain using the arnis method (class-based
#' widths; asphalt/concrete/dirt terrain colors). `"overture"` fetches
#' transportation segments from Overture Maps via a native DuckDB parquet
#' query (requires the suggested `duckdb` and `DBI` packages and internet,
#' matching [generate_block()]'s road source). An sf layer should contain
#' LINESTRING geometries with optional `class` and `subclass` columns
#' (Overture/OSM highway classes); missing classes default to `residential`.
#' Default `NULL` (no roads).
#' @param overture_release character or `NULL`. Overture Maps release string
#' used when `roads = "overture"` or `water = "overture"`, e.g.
#' `"2025-03-19.0"`. The default `NULL` (or `"auto"`) queries the public
#' bucket for the newest release and caches it for the session, so the code
#' keeps working as Overture publishes new releases. See
#' <https://github.com/OvertureMaps/data/releases> for the list.
#' @param bridges logical. If `TRUE` (default), road segments flagged as
#' bridges are built as elevated structures instead of being painted on the
#' ground: a deck raised `level * bridge_level_height` metres above the
#' highest terrain along the span, ramps down to grade at each end, railings,
#' and support pillars. Tunnel segments are always omitted from the ground
#' surface. Requires bridge attributes in the road data (Overture supplies
#' them; an sf layer may carry `is_bridge`, `is_tunnel`, and `level`).
#' @param bridge_level_height numeric. Metres of clearance per level, the
#' arnis `LAYER_HEIGHT_STEP`. Default 6.
#' @param bridge_ramp numeric. Length in metres over which a deck ramps down
#' to ground level at each end. Default 20.
#' @param bridge_pillar_interval numeric. Spacing in metres between support
#' pillars. Default 25.
#' @param water `NULL`, `"overture"`, or an sf (MULTI)POLYGON layer of water
#' bodies. Rivers, lakes, and coastal water are handled with the arnis
#' method: the terrain under each water body is flattened to just below its
#' lowest bank and colored as water, and a sand fringe of `sand_buffer`
#' metres is painted along the shore (riparian and coastal edges).
#' `"overture"` fetches water polygons from the Overture Maps base theme via
#' the native DuckDB query. Default `NULL`.
#' @param sand_buffer numeric. Width in metres of the sand fringe around
#' water bodies. 0 disables the fringe. Default 3.
#' @param surface_res numeric. Grid resolution in metres for the classified
#' terrain surface in the default (non-voxel) mode. Smaller values follow
#' road edges more precisely at the cost of a denser mesh. Default 2.
#' @param greenspace `NULL`, `TRUE`, `"esri"`, `"sentinel2"`, or a
#' `SpatRaster`. Paints ground-level greenery (lawns) into the terrain.
#' Greenery is classified from map tiles by the suggested `greenSD` package:
#' `TRUE` or `"esri"` uses Esri imagery, `"sentinel2"` uses Sentinel-2.
#' Alternatively supply your own binary `SpatRaster` (green = 1). Cells where
#' the canopy height model is at or above `min_tree_height` are excluded, so
#' this surface class shows ground-level greenery only - tree canopy is
#' already represented by the 3D tree objects. Default `NULL`.
#' @param greenspace_under_canopy logical. If `TRUE` (default), ground under
#' tree canopy stays green - trees stand on grass, as in arnis. Set `FALSE`
#' to cut canopy cells out of the greenspace surface, leaving bare ground
#' beneath the crowns.
#' @param greenspace_zoom integer. Map-tile zoom level for the greenspace
#' classification; higher zoom gives finer lawn edges. Default 17.
#' @param greenspace_year integer or `NULL`. Imagery year for
#' `greenspace = "sentinel2"`.
#' @param all_vox logical. If `TRUE`, voxelize the entire scene so the model
#' reproduces the blocky world arnis builds, from this package's data inputs:
#' terrain becomes stepped blocks (Minecraft-style heightmap), buildings
#' become block columns on flattened quantized bases, and trees are already
#' voxel objects. Building footprints are rasterized at `vox_size`, so exact
#' footprint edges are traded for the voxel look. Default `FALSE`.
#' @param vox_size numeric. Block edge length in metres for `all_vox` mode
#' (1 = Minecraft-like scale). Larger values give chunkier, lighter models.
#' Default 1.
#' @param basemap logical. If `TRUE`, download an Esri World Imagery snapshot
#' of the scene, save it next to the OBJ, and drape it on the ground via
#' texture coordinates (OBJ only; requires internet). With `terrain = FALSE`
#' a flat textured ground plane is added. On download failure the export
#' falls back to the flat terrain color with a warning. Default `FALSE`.
#' @param local_origin logical. If `TRUE` (default), shift the scene so its
#' minimum x/y is at (0, 0). The offset and CRS are stored in
#' `world_metadata.json` and in the returned object, so results can be
#' georeferenced back. Strongly recommended for Rhino/Blender, which lose
#' float precision at UTM-scale coordinates.
#' @param crop,data_source Passed to [search_3dglobdf()] when `x` is `NULL`.
#' @param quiet logical. Suppress progress messages. Default `TRUE`.
#'
#' @return (invisibly) a list:
#' \itemize{
#'   \item `paths`: character vector of written files.
#'   \item `meshes`: named list of in-memory meshes (`vertices`/`faces`).
#'   \item `origin`: list with `x`, `y` offset and `epsg` of the scene CRS.
#'   \item `n_buildings`, `n_trees`: scene composition counts.
#'   \item `n_trees_detected`: trees found in the canopy height model before
#'   the `max_trees` cap - compare with `n_trees` to see how many were
#'   dropped.
#'   \item `n_bridges`: elevated bridge segments built.
#' }
#'
#' @section Why some trees are missing:
#' Three filters reduce canopy pixels to tree objects. `min_tree_height`
#' (default 2 m) drops low vegetation; `tree_window` keeps only one tree per
#' window, so closely spaced trees merge into their tallest neighbour; and
#' `max_trees` caps the total, keeping the tallest. Compare
#' `n_trees_detected` with `n_trees` in the returned object, and lower
#' `tree_window` or raise `max_trees` if the model looks sparser than the
#' canopy data. Note that `tree_window` cannot resolve trees closer together
#' than the CHM resolution: on a coarse (e.g. aggregated) canopy raster, set
#' `tree_window` at or below the cell size to place one tree per canopy cell.
#'
#' @details
#' All coordinates are in metres (scene UTM zone). Building walls extend 1 m
#' below their sampled ground elevation when `terrain = TRUE` to avoid gaps on
#' slopes. Trees are blocky voxel-style objects whose shapes are ported from
#' the arnis project (trunk column, leaf columns, apex cap, and concentric
#' canopy rings). Every tree's total height comes from the canopy height
#' raster; the variant is chosen from that height (compact oaks below 8 m,
#' standard/bushy oaks to 14 m, oak/spruce to 22 m, towering spruce above,
#' with position-hashed variety within each class), then the block size is
#' scaled so the tree matches its measured height exactly.
#'
#' Ground surface classes are painted into the terrain, not built as
#' separate geometry: a classification grid (`surface_res` metres in default
#' mode, `vox_size` in voxel mode) assigns each terrain cell one of ground,
#' greenspace, roads, sidewalks, paths, sand, or water, and the terrain mesh
#' carries these as per-face material groups. Roads follow the arnis method:
#' each segment gets a class-based half-width (motorway/trunk/primary 5 m,
#' secondary 4 m, tertiary 3 m, residential/service 2 m, foot/cycle/path
#' classes 1 m - the arnis `highway_block_range` defaults), is buffered into
#' a ribbon, and painted as `roads` (asphalt), `sidewalks` (concrete;
#' includes Overture `subclass = "sidewalk"`), or `paths` (dirt).
#'
#' Bridges and elevated highways follow arnis's structural approach. Neither
#' OSM nor Overture records deck elevations, so the height is inferred from
#' the segment's `level`: the deck sits `level * bridge_level_height` metres
#' (arnis uses 6 blocks per layer) above the highest terrain along the span,
#' ramps linearly down to grade over `bridge_ramp` metres at each end, and is
#' carried by pillars dropped to the ground every `bridge_pillar_interval`
#' metres. Decks are solid slabs with railings, sized by the same class-based
#' widths as surface roads. Bridge and tunnel segments are excluded from the
#' terrain surface classification, so a bridge no longer leaves a road
#' painted on the ground beneath it.
#'
#' Riparian and coastal areas also follow arnis: the terrain under each
#' water body is flattened to one level just below its lowest bank (0.2 m in
#' default mode, one block in voxel mode), painted as water, and fringed
#' with a `sand_buffer`-metre sand strip along the shoreline.
#'
#' Vegetation is split across two levels so canopy and lawns are never
#' conflated: tree canopy comes from the canopy height model and is rendered
#' as 3D tree objects, while the `greenspace` terrain class carries the
#' ground surface. Greenery comes from `greenSD`'s map-tile classification
#' (Esri or Sentinel-2 imagery at `greenspace_zoom`). Ground under canopy
#' stays grass by default, so a park reads as a continuous lawn with trees
#' standing on it - exactly how arnis places tree objects on grass blocks.
#' Use `greenspace_under_canopy = FALSE` if you instead want canopy cells cut
#' out of the lawn.
#'
#' With `all_vox = TRUE` the whole scene is voxelized the way arnis builds
#' its Minecraft worlds: the DEM becomes a stepped block heightmap, each
#' building footprint is rasterized at `vox_size` and raised as quantized
#' block columns on a flattened base (extending one block below ground), and
#' the voxel trees stand on the quantized ground. Heights still come from the
#' input data (building `Height` column, CHM tree heights, DEM terrain) -
#' only the geometry representation changes. Expect larger files than the
#' default mode; increase `vox_size` to lighten them.
#'
#' When `basemap = TRUE`, imagery is retrieved from the Esri World Imagery
#' service; check the
#' \href{https://www.esri.com/en-us/legal/terms/master-agreement}{Esri
#' terms of use} and provide attribution (Esri, Maxar, Earthstar Geographics,
#' and the GIS User Community) when publishing rendered scenes. Footprint holes (courtyards) are not carved in this version; the
#' exterior ring is extruded and a warning is issued when holes are dropped.
#'
#' @examples
#' \donttest{
#' example <- gloBFPr::globfp_example
#' world <- get_3d_world(
#'   x = example, terrain = FALSE, canopy = NULL,
#'   out_dir = tempfile("world3d")
#' )
#' }
#'
#' @importFrom grDevices hcl.colors col2rgb
#' @importFrom stats quantile median
#' @export
get_3d_world <- function(x = NULL,
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
                         quiet = TRUE) {
  format <- match.arg(tolower(format), c("obj", "stl"), several.ok = TRUE)
  if (isTRUE(terrain) && is.null(key) && is.null(dem)) {
    stop("`key` (OpenTopography API key) is required when `terrain = TRUE` ",
         "and no `dem` raster is supplied. Use `terrain = FALSE` for a ",
         "flat-ground scene.", call. = FALSE)
  }
  if (!is.null(canopy)) {
    canopy <- match.arg(tolower(as.character(canopy[1])),
                        c("metachm", "ethchm"))
  }

  # ---- Buildings -----------------------------------------------------------
  if (is.null(x)) {
    if (is.null(bbox) && is.null(place)) {
      stop("Provide `x`, `bbox`, or `place`.", call. = FALSE)
    }
    if (!quiet) cli::cli_alert_info("Fetching buildings with search_3dglobdf() ...")
    x <- search_3dglobdf(bbox = bbox, place = place, crop = crop,
                         data_source = data_source, out_type = "poly",
                         quiet = quiet)
    if (is.null(x)) stop("No buildings found for the area of interest.", call. = FALSE)
  }
  if (!inherits(x, "sf")) stop("`x` must be an sf polygon object.", call. = FALSE)
  if (!height_col %in% names(x)) {
    stop("Column `", height_col, "` not found in `x`.", call. = FALSE)
  }

  if (is.na(sf::st_crs(x))) {
    stop("`x` must have a coordinate reference system.", call. = FALSE)
  }
  bbox_wgs <- get_bbox(x)
  utm_crs <- get_utm_crs(bbox_wgs)
  if (isTRUE(sf::st_is_longlat(x))) {
    x <- sf::st_transform(x, crs = utm_crs)
  }
  x <- sf::st_make_valid(x)
  x <- suppressWarnings(sf::st_collection_extract(x, "POLYGON"))
  x <- suppressWarnings(sf::st_cast(x, "POLYGON"))
  if (nrow(x) == 0) stop("No valid building polygons in `x`.", call. = FALSE)

  bbox_vector <- bbox_poly_to_list(bbox_wgs)

  # ---- Terrain -------------------------------------------------------------
  if (isTRUE(terrain)) {
    if (is.null(dem)) {
      if (!quiet) cli::cli_alert_info("Downloading DEM ...")
      dem <- get_dem(bbox_vector, key)
    } else if (inherits(dem, "PackedSpatRaster")) {
      dem <- terra::rast(dem)
    }
    if (!terra::same.crs(dem, paste0("EPSG:", utm_crs))) {
      dem <- terra::project(dem, paste0("EPSG:", utm_crs), method = "bilinear")
    }
    scene_ext <- terra::ext(terra::vect(sf::st_transform(bbox_wgs, utm_crs)))
    dem <- terra::crop(dem, scene_ext, snap = "out")
    if (is.numeric(simplify_terrain) && simplify_terrain > 1) {
      dem <- terra::aggregate(dem, fact = round(simplify_terrain),
                              fun = "mean", na.rm = TRUE)
    }
  } else {
    dem <- NULL
  }

  # ---- Canopy --------------------------------------------------------------
  tree_pts <- NULL
  chm <- NULL
  if (!is.null(canopy_height)) {
    chm <- canopy_height
    if (inherits(chm, "PackedSpatRaster")) chm <- terra::rast(chm)
    if (!terra::same.crs(chm, paste0("EPSG:", utm_crs))) {
      chm <- terra::project(chm, paste0("EPSG:", utm_crs), method = "near")
    }
    chm <- terra::ifel(chm < min_tree_height, 0, chm)
  } else if (!is.null(canopy)) {
    if (!quiet) cli::cli_alert_info("Downloading canopy height model ...")
    chm <- tryCatch(
      get_chm(bbox_vector, min_tree_height, datasource = canopy)[[1]],
      error = function(e) {
        warning("Canopy download failed (", conditionMessage(e),
                "); continuing without trees.", call. = FALSE)
        NULL
      }
    )
  }
  n_trees_detected <- 0L
  if (!is.null(chm)) {
    tree_pts <- detect_tree_tops(chm, window = tree_window,
                                 min_height = min_tree_height)
    n_trees_detected <- nrow(tree_pts)
    if (nrow(tree_pts) > max_trees) {
      # Always warn: silently dropping trees is confusing when the rendered
      # scene looks sparser than the canopy data.
      warning(n_trees_detected, " trees detected; only the ", max_trees,
              " tallest are exported. Raise `max_trees` to keep more.",
              call. = FALSE)
      tree_pts <- tree_pts[order(-tree_pts$height)[seq_len(max_trees)], ]
    }
  }

  # ---- Surface data: roads, greenspace, water ------------------------------
  # These classes are part of the terrain surface itself: they are rasterized
  # into a classification grid and rendered as differently colored regions of
  # one detailed terrain mesh (not as separate floating geometry).
  ribbons_by_group <- list()
  bridge_lines <- NULL
  if (!is.null(roads)) {
    roads_sf <- if (inherits(roads, "sf")) {
      roads
    } else if (is.character(roads) && identical(tolower(roads[1]), "overture")) {
      fetch_overture_roads(bbox_wgs, release = overture_release, quiet = quiet)
    } else {
      stop("`roads` must be NULL, \"overture\", or an sf line layer.",
           call. = FALSE)
    }
    if (!is.null(roads_sf) && nrow(roads_sf) > 0) {
      roads_sf <- sf::st_zm(roads_sf, drop = TRUE)
      roads_sf <- sf::st_transform(roads_sf, crs = utm_crs)
      geom_type <- as.character(sf::st_geometry_type(roads_sf))
      roads_sf <- roads_sf[grepl("LINESTRING", geom_type), , drop = FALSE]
    }
    if (!is.null(roads_sf) && nrow(roads_sf) > 0) {
      cls <- if ("class" %in% names(roads_sf)) {
        as.character(roads_sf$class)
      } else {
        rep("residential", nrow(roads_sf))
      }
      cls[is.na(cls) | cls == ""] <- "residential"
      subclass <- if ("subclass" %in% names(roads_sf)) {
        as.character(roads_sf$subclass)
      } else {
        rep(NA_character_, nrow(roads_sf))
      }
      grp <- road_group(cls, subclass)
      hw <- road_half_width(cls)
      is_bridge <- if ("is_bridge" %in% names(roads_sf)) {
        isTRUE(bridges) & !is.na(roads_sf$is_bridge) & roads_sf$is_bridge
      } else {
        rep(FALSE, nrow(roads_sf))
      }
      is_tunnel <- if ("is_tunnel" %in% names(roads_sf)) {
        !is.na(roads_sf$is_tunnel) & roads_sf$is_tunnel
      } else {
        rep(FALSE, nrow(roads_sf))
      }
      # Bridges are built as elevated structures and tunnels run below
      # ground, so neither is painted onto the terrain surface.
      on_ground <- !is_bridge & !is_tunnel
      ribbons <- sf::st_buffer(sf::st_geometry(roads_sf), dist = hw,
                               endCapStyle = "FLAT")
      for (g in c("roads", "sidewalks", "paths")) {
        sel <- ribbons[grp == g & on_ground]
        if (length(sel) > 0) ribbons_by_group[[g]] <- sf::st_union(sel)
      }
      if (any(is_bridge)) {
        bridge_lines <- roads_sf[is_bridge, , drop = FALSE]
        bridge_lines$half_width <- hw[is_bridge]
        bridge_lines$vehicular <- grp[is_bridge] == "roads"
        if (!"level" %in% names(bridge_lines)) bridge_lines$level <- 1L
        bridge_lines$level[!is.finite(bridge_lines$level) |
                             bridge_lines$level < 1] <- 1L
      }
    } else if (!quiet) {
      cli::cli_alert_warning("No road segments found for the scene.")
    }
  }

  green_r <- NULL
  if (!is.null(greenspace) && !isFALSE(greenspace)) {
    green_src <- NULL
    if (!inherits(greenspace, c("SpatRaster", "PackedSpatRaster"))) {
      # Validate up front, outside tryCatch, so a bad source is an error
      # rather than a swallowed warning. Map-tile sources only: greenSD
      # classifies imagery into green / non-green. The CHM option belongs to
      # `canopy`, which drives the 3D trees.
      green_src <- if (isTRUE(greenspace)) "esri" else tolower(greenspace[1])
      if (identical(green_src, "metachm")) {
        stop("`greenspace` uses greenSD map tiles: choose \"esri\" or ",
             "\"sentinel2\" (or supply your own binary SpatRaster). ",
             "\"metachm\" is a `canopy` source - canopy is represented by ",
             "the 3D trees, not the greenspace surface.", call. = FALSE)
      }
      green_src <- match.arg(green_src, c("esri", "sentinel2"))
    }
    green_r <- tryCatch({
      if (inherits(greenspace, "SpatRaster")) {
        greenspace
      } else if (inherits(greenspace, "PackedSpatRaster")) {
        terra::rast(greenspace)
      } else {
        if (!quiet) cli::cli_alert_info("Fetching greenspace map tiles ...")
        get_greenspace(bbox = bbox_wgs, type = green_src,
                       zoom = greenspace_zoom, year = greenspace_year)
      }
    }, error = function(e) {
      warning("Greenspace retrieval failed (", conditionMessage(e),
              "); continuing without a greenspace layer.", call. = FALSE)
      NULL
    })
    if (!is.null(green_r)) {
      if (!terra::same.crs(green_r, paste0("EPSG:", utm_crs))) {
        green_r <- terra::project(green_r, paste0("EPSG:", utm_crs),
                                  method = "near")
      }
      green_r <- terra::ifel(green_r > 0, 1, NA)
      # Canopy vs ground: the 3D trees carry the canopy, but the ground
      # beneath them is still grass, so by default it stays green (as in
      # arnis, where trees stand on grass blocks). Set
      # `greenspace_under_canopy = FALSE` to cut canopy cells out instead.
      if (!is.null(chm) && !isTRUE(greenspace_under_canopy)) {
        chm_on_green <- terra::resample(chm, green_r, method = "near")
        green_r <- terra::ifel(
          !is.na(chm_on_green) & chm_on_green >= min_tree_height,
          NA, green_r)
      }
    }
  }

  water_polys <- NULL
  if (!is.null(water)) {
    water_sf <- if (inherits(water, c("sf", "sfc"))) {
      sf::st_as_sf(water)
    } else if (is.character(water) && identical(tolower(water[1]), "overture")) {
      fetch_overture_water(bbox_wgs, release = overture_release, quiet = quiet)
    } else {
      stop("`water` must be NULL, \"overture\", or an sf polygon layer.",
           call. = FALSE)
    }
    if (!is.null(water_sf) && nrow(water_sf) > 0) {
      water_sf <- sf::st_zm(water_sf, drop = TRUE)
      water_sf <- sf::st_transform(water_sf, crs = utm_crs)
      water_sf <- sf::st_make_valid(water_sf)
      geom_type <- as.character(sf::st_geometry_type(water_sf))
      water_sf <- water_sf[grepl("POLYGON", geom_type), , drop = FALSE]
      if (nrow(water_sf) > 0) {
        water_polys <- suppressWarnings(
          sf::st_cast(sf::st_geometry(water_sf), "POLYGON"))
      }
    }
    if (is.null(water_polys) && !quiet) {
      cli::cli_alert_warning("No water polygons found for the scene.")
    }
  }

  has_surface <- length(ribbons_by_group) > 0 || !is.null(green_r) ||
    !is.null(water_polys)

  # ---- Ground elevations ---------------------------------------------------
  centroids <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(x))))
  if (!is.null(dem)) {
    base_elev <- terra::extract(dem, centroids[, 1:2, drop = FALSE],
                                method = "bilinear")[, 1]
    dem_med <- stats::median(terra::values(dem, mat = FALSE), na.rm = TRUE)
    base_elev[!is.finite(base_elev)] <- dem_med
    wall_below <- 1
    if (!is.null(tree_pts) && nrow(tree_pts) > 0) {
      tz <- terra::extract(dem, as.matrix(tree_pts[, c("x", "y")]),
                           method = "bilinear")[, 1]
      tz[!is.finite(tz)] <- dem_med
      tree_pts$z <- tz
    }
  } else {
    base_elev <- rep(0, nrow(x))
    wall_below <- 0
    if (!is.null(tree_pts) && nrow(tree_pts) > 0) tree_pts$z <- 0
  }

  # ---- Local origin --------------------------------------------------------
  scene_bb <- sf::st_bbox(x)
  if (!is.null(dem)) {
    e <- terra::ext(dem)
    scene_bb["xmin"] <- min(scene_bb["xmin"], e[1])
    scene_bb["ymin"] <- min(scene_bb["ymin"], e[3])
  }
  origin <- if (isTRUE(local_origin)) {
    c(x = unname(scene_bb["xmin"]), y = unname(scene_bb["ymin"]))
  } else {
    c(x = 0, y = 0)
  }

  # ---- Build meshes --------------------------------------------------------
  if (!quiet) cli::cli_alert_info("Building meshes ...")

  heights <- as.numeric(x[[height_col]])
  terrain_mesh <- NULL

  if (isTRUE(all_vox)) {
    # -- Voxel world (arnis-style) ------------------------------------------
    if (!is.numeric(vox_size) || length(vox_size) != 1L ||
        !is.finite(vox_size) || vox_size <= 0) {
      stop("`vox_size` must be one positive number.", call. = FALSE)
    }
    s <- vox_size
    gxmin <- floor(as.numeric(scene_bb["xmin"]) / s) * s
    gymin <- floor(as.numeric(scene_bb["ymin"]) / s) * s
    gxmax <- ceiling(as.numeric(scene_bb["xmax"]) / s) * s
    gymax <- ceiling(as.numeric(scene_bb["ymax"]) / s) * s
    if (!is.null(dem)) {
      e <- terra::ext(dem)
      gxmin <- min(gxmin, floor(as.numeric(e[1]) / s) * s)
      gymin <- min(gymin, floor(as.numeric(e[3]) / s) * s)
      gxmax <- max(gxmax, ceiling(as.numeric(e[2]) / s) * s)
      gymax <- max(gymax, ceiling(as.numeric(e[4]) / s) * s)
    }
    template <- terra::rast(xmin = gxmin, xmax = gxmax,
                            ymin = gymin, ymax = gymax,
                            resolution = s,
                            crs = paste0("EPSG:", utm_crs))
    if (!is.null(dem)) {
      ground_r <- terra::resample(dem, template, method = "bilinear")
      gv <- terra::values(ground_r, mat = FALSE)
      gv[!is.finite(gv)] <- stats::median(gv, na.rm = TRUE)
      terra::values(ground_r) <- gv
    } else {
      ground_r <- template
      terra::values(ground_r) <- 0
    }
    # Arnis-style water: flatten each water body to just below its banks
    # (one block) before quantization
    if (!is.null(water_polys)) {
      ground_r <- flatten_water_elev(ground_r, water_polys, drop = s)
    }
    gv <- terra::values(ground_r, mat = FALSE)
    gv <- round(gv / s) * s
    terra::values(ground_r) <- gv

    # Buildings: rasterized footprints -> quantized block columns on a
    # flattened base (like arnis grades under each building)
    x$vox_id <- seq_len(nrow(x))
    bld_r <- terra::rasterize(terra::vect(x), template, field = "vox_id")
    bv <- terra::values(bld_r, mat = FALSE)
    building_meshes <- vector("list", nrow(x))
    for (i in seq_len(nrow(x))) {
      h <- heights[i]
      if (!is.finite(h) || h <= 0) next
      cell_idx <- which(bv == i)
      if (length(cell_idx) == 0) next
      xy <- terra::xyFromCell(template, cell_idx)
      base_flat <- min(gv[cell_idx])
      cells <- data.frame(
        x = xy[, 1] - origin["x"],
        y = xy[, 2] - origin["y"],
        z0 = base_flat - s,
        z1 = base_flat + ceiling(h / s) * s
      )
      building_meshes[[i]] <- mesh_voxel_columns(cells, s)
    }
    keep <- !vapply(building_meshes, is.null, logical(1))
    buildings_mesh <- mesh_combine(building_meshes[keep])
    if (is.null(buildings_mesh)) {
      stop("No buildings could be voxelized (check `", height_col,
           "` values, or decrease `vox_size`).", call. = FALSE)
    }
    building_ids <- if ("id" %in% names(x)) x$id[keep] else which(keep)

    # Terrain: stepped voxel heightmap with per-cell surface classes
    # (flat slab when terrain = FALSE and no surface classes are present)
    if (!is.null(dem) || has_surface) {
      tops <- matrix(gv, nrow = terra::nrow(template), byrow = TRUE)
      xs <- terra::xFromCol(template, seq_len(terra::ncol(template))) - origin["x"]
      ys <- terra::yFromRow(template, seq_len(terra::nrow(template))) - origin["y"]
      terrain_mesh <- mesh_voxel_heightmap(tops, xs, ys, s)
      if (has_surface && !is.null(terrain_mesh)) {
        cls_r <- classify_surface_raster(template, ribbons_by_group,
                                         green_r, water_polys, sand_buffer)
        cls_v <- terra::values(cls_r, mat = FALSE)
        terrain_mesh$face_class <-
          surface_class_names()[cls_v[terrain_mesh$face_cell]]
      }
    } else {
      terrain_mesh <- mesh_box(gxmin - origin["x"], gxmax - origin["x"],
                               gymin - origin["y"], gymax - origin["y"],
                               -s, 0)
    }

    # Trees stand on the quantized voxel ground
    if (!is.null(tree_pts) && nrow(tree_pts) > 0) {
      tz <- terra::extract(ground_r, as.matrix(tree_pts[, c("x", "y")]))[, 1]
      tz[!is.finite(tz)] <- stats::median(gv, na.rm = TRUE)
      tree_pts$z <- tz
    }
  } else {
    # -- True-footprint extrusion (default) ---------------------------------
    building_meshes <- vector("list", nrow(x))
    dropped_holes <- 0L
    for (i in seq_len(nrow(x))) {
      geom <- sf::st_geometry(x[i, ])[[1]]
      if (length(geom) > 1) dropped_holes <- dropped_holes + 1L
      ring <- as.matrix(geom[[1]])[, 1:2, drop = FALSE]
      if (nrow(ring) > 1 && all(ring[1, ] == ring[nrow(ring), ])) {
        ring <- ring[-nrow(ring), , drop = FALSE]
      }
      ring[, 1] <- ring[, 1] - origin["x"]
      ring[, 2] <- ring[, 2] - origin["y"]
      h <- heights[i]
      if (!is.finite(h) || h <= 0) next
      building_meshes[[i]] <- mesh_extrude_polygon(
        ring, z0 = base_elev[i] - wall_below, z1 = base_elev[i] + h)
    }
    if (dropped_holes > 0 && !quiet) {
      cli::cli_alert_warning(paste0(
        dropped_holes, " footprint(s) had interior rings (courtyards); ",
        "holes were not carved."))
    }
    keep <- !vapply(building_meshes, is.null, logical(1))
    buildings_mesh <- mesh_combine(building_meshes[keep])
    if (is.null(buildings_mesh)) {
      stop("No buildings could be meshed (check `", height_col, "` values).",
           call. = FALSE)
    }
    building_ids <- if ("id" %in% names(x)) x$id[keep] else which(keep)

    if (has_surface) {
      # Detailed classified terrain: one surface grid at `surface_res`, with
      # roads/sidewalks/paths/greenspace/water/sand as per-face materials
      if (!is.numeric(surface_res) || surface_res <= 0) {
        stop("`surface_res` must be one positive number.", call. = FALSE)
      }
      sr <- surface_res
      sxmin <- floor(as.numeric(scene_bb["xmin"]) / sr) * sr
      symin <- floor(as.numeric(scene_bb["ymin"]) / sr) * sr
      sxmax <- ceiling(as.numeric(scene_bb["xmax"]) / sr) * sr
      symax <- ceiling(as.numeric(scene_bb["ymax"]) / sr) * sr
      if (!is.null(dem)) {
        e <- terra::ext(dem)
        sxmax <- max(sxmax, ceiling(as.numeric(e[2]) / sr) * sr)
        symax <- max(symax, ceiling(as.numeric(e[4]) / sr) * sr)
      }
      surf_t <- terra::rast(xmin = sxmin, xmax = sxmax,
                            ymin = symin, ymax = symax,
                            resolution = sr,
                            crs = paste0("EPSG:", utm_crs))
      if (!is.null(dem)) {
        elev_r <- terra::resample(dem, surf_t, method = "bilinear")
        ev <- terra::values(elev_r, mat = FALSE)
        ev[!is.finite(ev)] <- stats::median(ev, na.rm = TRUE)
        terra::values(elev_r) <- ev
      } else {
        elev_r <- surf_t
        terra::values(elev_r) <- 0
      }
      # Arnis-style water: flatten each water body to just below its banks
      if (!is.null(water_polys)) {
        elev_r <- flatten_water_elev(elev_r, water_polys, drop = 0.2)
      }
      cls_r <- classify_surface_raster(surf_t, ribbons_by_group,
                                       green_r, water_polys, sand_buffer)
      elev_local <- elev_r
      e <- terra::ext(elev_local)
      terra::ext(elev_local) <- terra::ext(
        e[1] - origin["x"], e[2] - origin["x"],
        e[3] - origin["y"], e[4] - origin["y"])
      terrain_mesh <- mesh_terrain(elev_local)
      if (!is.null(terrain_mesh)) {
        cls_v <- terra::values(cls_r, mat = FALSE)
        terrain_mesh$face_class <-
          surface_class_names()[cls_v[terrain_mesh$face_cell]]
      }
    } else if (!is.null(dem)) {
      dem_local <- dem
      e <- terra::ext(dem_local)
      terra::ext(dem_local) <- terra::ext(e[1] - origin["x"], e[2] - origin["x"],
                                          e[3] - origin["y"], e[4] - origin["y"])
      terrain_mesh <- mesh_terrain(dem_local)
    }
  }

  trunks_mesh <- crowns_mesh <- NULL
  n_trees <- 0L
  if (!is.null(tree_pts) && nrow(tree_pts) > 0) {
    trunk_list <- vector("list", nrow(tree_pts))
    crown_list <- vector("list", nrow(tree_pts))
    for (i in seq_len(nrow(tree_pts))) {
      tr <- mesh_tree(tree_pts$x[i] - origin["x"],
                      tree_pts$y[i] - origin["y"],
                      tree_pts$z[i], tree_pts$height[i])
      if (is.null(tr)) next
      trunk_list[[i]] <- tr$trunk
      crown_list[[i]] <- tr$crown
    }
    trunks_mesh <- mesh_combine(trunk_list)
    crowns_mesh <- mesh_combine(crown_list)
    n_trees <- sum(!vapply(crown_list, is.null, logical(1)))
  }

  # ---- Bridges (arnis method: layer-derived deck, ramps, pillars) ----------
  bridge_mesh <- support_mesh <- NULL
  n_bridges <- 0L
  if (!is.null(bridge_lines) && nrow(bridge_lines) > 0) {
    s_vox <- if (isTRUE(all_vox)) vox_size else NA_real_
    ground_at <- function(xy) {
      if (isTRUE(all_vox)) {
        z <- terra::extract(ground_r, xy)[, 1]
        z[!is.finite(z)] <- stats::median(gv, na.rm = TRUE)
        z
      } else if (!is.null(dem)) {
        z <- terra::extract(dem, xy, method = "bilinear")[, 1]
        z[!is.finite(z)] <- stats::median(
          terra::values(dem, mat = FALSE), na.rm = TRUE)
        z
      } else {
        rep(0, nrow(xy))
      }
    }
    deck_parts <- list()
    support_parts <- list()
    for (i in seq_len(nrow(bridge_lines))) {
      co <- sf::st_coordinates(
        sf::st_segmentize(sf::st_geometry(bridge_lines[i, ]), 5))
      co <- co[, 1:2, drop = FALSE]
      if (nrow(co) < 2) next
      g <- ground_at(co)
      # Deck height: arnis raises the deck `level` steps above the highest
      # terrain along the span
      clearance <- bridge_lines$level[i] * bridge_level_height
      deck_z <- max(g, na.rm = TRUE) + clearance
      hw_i <- bridge_lines$half_width[i]
      thick <- if (isTRUE(all_vox)) s_vox else 0.6
      if (isTRUE(all_vox)) {
        deck_z <- round(deck_z / s_vox) * s_vox
        hw_i <- max(round(hw_i / s_vox), 1) * s_vox
      }
      z <- bridge_deck_profile(co, g, deck_z, ramp = bridge_ramp)
      if (isTRUE(all_vox)) z <- round(z / s_vox) * s_vox
      rail <- if (isTRUE(all_vox)) {
        s_vox
      } else if (isTRUE(bridge_lines$vehicular[i])) 1 else 0.9
      co_local <- cbind(co[, 1] - origin["x"], co[, 2] - origin["y"])
      dk <- mesh_bridge_deck(co_local, z, hw_i, thickness = thick,
                             railing = rail)
      if (is.null(dk)) next
      deck_parts[[length(deck_parts) + 1L]] <- dk
      n_bridges <- n_bridges + 1L
      sp <- mesh_bridge_supports(
        co_local, z - thick, g,
        interval = bridge_pillar_interval,
        radius = if (isTRUE(all_vox)) s_vox / 2 else max(hw_i * 0.15, 0.4))
      if (!is.null(sp)) support_parts[[length(support_parts) + 1L]] <- sp
    }
    if (length(deck_parts) > 0) {
      m <- mesh_combine(deck_parts)
      bridge_mesh <- list(vertices = m$vertices, faces = m$faces)
    }
    if (length(support_parts) > 0) {
      m <- mesh_combine(support_parts)
      support_mesh <- list(vertices = m$vertices, faces = m$faces)
    }
    if (!quiet && n_bridges > 0) {
      cli::cli_alert_info(paste(n_bridges, "bridge segment(s) built."))
    }
  }

  # ---- Materials / colors --------------------------------------------------
  materials <- list(
    terrain       = c(0.82, 0.77, 0.68),
    buildings     = c(0.75, 0.75, 0.78),
    canopy_trunks = c(0.45, 0.32, 0.20),
    canopy_crowns = c(0.30, 0.55, 0.25),
    roads         = c(0.16, 0.16, 0.18),
    sidewalks     = c(0.62, 0.62, 0.62),
    paths         = c(0.55, 0.45, 0.32),
    greenspace    = c(0.45, 0.62, 0.30),
    water         = c(0.25, 0.42, 0.65),
    sand          = c(0.86, 0.80, 0.62),
    bridges       = c(0.55, 0.55, 0.57),
    bridge_supports = c(0.42, 0.42, 0.44)
  )
  building_groups <- NULL
  if (!is.null(color_by)) {
    if (!color_by %in% names(x)) {
      stop("`color_by` column `", color_by, "` not found in `x`.", call. = FALSE)
    }
    vals <- as.numeric(sf::st_drop_geometry(x)[[color_by]])[keep]
    pal <- grDevices::hcl.colors(11, "viridis")
    brks <- unique(stats::quantile(vals, probs = seq(0, 1, length.out = 12),
                                   na.rm = TRUE))
    if (length(brks) < 2) {
      bins <- rep(1L, length(vals))
    } else {
      bins <- cut(vals, breaks = brks, include.lowest = TRUE, labels = FALSE)
    }
    bins[is.na(bins)] <- 1L
    bins <- pmin(bins, 11L)
    for (b in sort(unique(bins))) {
      rgb01 <- as.numeric(grDevices::col2rgb(pal[b])) / 255
      materials[[paste0("bldg_ramp_", b)]] <- rgb01
    }
    building_groups <- lapply(seq_along(building_ids), function(k) {
      list(name = paste0("building_", building_ids[k]),
           faces = which(buildings_mesh$face_groups == k),
           material = paste0("bldg_ramp_", bins[k]))
    })
  } else if (!isFALSE(facade_palette)) {
    pal_cols <- if (is.character(facade_palette)) {
      facade_palette
    } else {
      # Muted facade tones in the spirit of arnis' wall-color palette
      c("#BFB9AC", "#C8BEB0", "#A39E93", "#B08968",
        "#8D99AE", "#CFC8B8", "#9E948A", "#ADA79E")
    }
    for (i in seq_along(pal_cols)) {
      materials[[paste0("facade_", i)]] <-
        as.numeric(grDevices::col2rgb(pal_cols[i])) / 255
    }
    building_groups <- lapply(seq_along(building_ids), function(k) {
      list(name = paste0("building_", building_ids[k]),
           faces = which(buildings_mesh$face_groups == k),
           material = paste0(
             "facade_", (as.integer(building_ids[k]) %% length(pal_cols)) + 1L))
    })
  } else {
    building_groups <- lapply(seq_along(building_ids), function(k) {
      list(name = paste0("building_", building_ids[k]),
           faces = which(buildings_mesh$face_groups == k),
           material = "buildings")
    })
  }

  # ---- Write ---------------------------------------------------------------
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  paths <- character(0)

  # Optional basemap texture (Esri World Imagery) draped on the ground
  basemap_file <- NULL
  gx <- gy <- NULL
  if (isTRUE(basemap) && "obj" %in% format) {
    if (!is.null(dem)) {
      e <- terra::ext(dem)
      gx <- as.numeric(c(e[1], e[2])); gy <- as.numeric(c(e[3], e[4]))
    } else {
      gx <- as.numeric(c(scene_bb["xmin"], scene_bb["xmax"]))
      gy <- as.numeric(c(scene_bb["ymin"], scene_bb["ymax"]))
    }
    if (!quiet) cli::cli_alert_info("Downloading basemap imagery ...")
    basemap_file <- tryCatch(
      fetch_basemap_image(gx, gy, utm_crs, out_dir),
      error = function(e) {
        warning("Basemap download failed (", conditionMessage(e),
                "); using flat terrain color.", call. = FALSE)
        NULL
      }
    )
    if (!is.null(basemap_file) && is.null(terrain_mesh)) {
      # Flat textured ground plane for terrain = FALSE scenes
      x0 <- gx[1] - origin["x"]; x1 <- gx[2] - origin["x"]
      y0 <- gy[1] - origin["y"]; y1 <- gy[2] - origin["y"]
      terrain_mesh <- list(
        vertices = cbind(c(x0, x1, x1, x0), c(y0, y0, y1, y1), 0),
        faces = matrix(c(1L, 2L, 3L, 1L, 3L, 4L), ncol = 3, byrow = TRUE)
      )
    }
  }

  layers <- list()
  if (!is.null(terrain_mesh)) {
    if (!is.null(basemap_file)) {
      # Satellite texture takes precedence over surface-class colors
      uv <- cbind(
        (terrain_mesh$vertices[, 1] - (gx[1] - origin["x"])) / (gx[2] - gx[1]),
        (terrain_mesh$vertices[, 2] - (gy[1] - origin["y"])) / (gy[2] - gy[1])
      )
      materials$terrain_basemap <- list(Kd = c(1, 1, 1),
                                        map_Kd = basename(basemap_file))
      layers$terrain <- list(mesh = terrain_mesh,
                             material = "terrain_basemap", uv = uv)
      paths <- c(paths, basemap_file)
    } else if (!is.null(terrain_mesh$face_class)) {
      # Classified terrain: surface classes as face groups with materials
      terrain_groups <- lapply(
        intersect(surface_class_names(), unique(terrain_mesh$face_class)),
        function(nm) {
          list(name = nm,
               faces = which(terrain_mesh$face_class == nm),
               material = nm)
        })
      layers$terrain <- list(mesh = terrain_mesh, groups = terrain_groups)
    } else {
      layers$terrain <- list(mesh = terrain_mesh, material = "terrain")
    }
  }
  if (!is.null(bridge_mesh)) {
    layers$bridges <- list(mesh = bridge_mesh, material = "bridges")
  }
  if (!is.null(support_mesh)) {
    layers$bridge_supports <- list(mesh = support_mesh,
                                   material = "bridge_supports")
  }
  layers$buildings <- list(mesh = buildings_mesh, groups = building_groups)
  if (!is.null(trunks_mesh)) {
    layers$canopy_trunks <- list(mesh = trunks_mesh, material = "canopy_trunks")
    layers$canopy_crowns <- list(mesh = crowns_mesh, material = "canopy_crowns")
  }

  # Only keep materials actually referenced by a layer or building group
  used_mats <- character(0)
  for (layer in layers) {
    if (!is.null(layer$material)) used_mats <- c(used_mats, layer$material)
    if (!is.null(layer$groups)) {
      used_mats <- c(used_mats,
                     unlist(lapply(layer$groups, `[[`, "material")))
    }
  }
  materials <- materials[intersect(names(materials), unique(used_mats))]

  if ("obj" %in% format) {
    obj_path <- file.path(out_dir, "world.obj")
    paths <- c(paths, write_world_obj(layers, obj_path, materials = materials))
  }
  if ("stl" %in% format) {
    stl_layers <- list(buildings = buildings_mesh, terrain = terrain_mesh,
                       canopy_trunks = trunks_mesh, canopy_crowns = crowns_mesh,
                       bridges = bridge_mesh, bridge_supports = support_mesh)
    for (nm in names(stl_layers)) {
      if (is.null(stl_layers[[nm]])) next
      p <- file.path(out_dir, paste0(nm, ".stl"))
      write_stl_binary(stl_layers[[nm]], p, name = nm)
      paths <- c(paths, p)
    }
  }

  meta_path <- file.path(out_dir, "world_metadata.json")
  meta <- paste0(
    "{\n",
    '  "generator": "gloBFPr::get_3d_world",\n',
    '  "epsg": ', utm_crs, ",\n",
    '  "origin_x": ', format(origin["x"], scientific = FALSE), ",\n",
    '  "origin_y": ', format(origin["y"], scientific = FALSE), ",\n",
    '  "units": "m",\n',
    '  "n_buildings": ', sum(keep), ",\n",
    '  "n_trees": ', n_trees, "\n",
    "}"
  )
  writeLines(meta, meta_path)
  paths <- c(paths, meta_path)

  if (!quiet) {
    cli::cli_alert_success(paste0(
      "3D world written to ", normalizePath(out_dir, mustWork = FALSE),
      " (", sum(keep), " buildings, ", n_trees, " trees)."))
  }

  invisible(list(
    paths = normalizePath(paths, mustWork = FALSE),
    meshes = list(buildings = buildings_mesh, terrain = terrain_mesh,
                  canopy_trunks = trunks_mesh, canopy_crowns = crowns_mesh,
                  bridges = bridge_mesh, bridge_supports = support_mesh),
    origin = list(x = unname(origin["x"]), y = unname(origin["y"]),
                  epsg = utm_crs),
    n_buildings = sum(keep),
    n_trees = n_trees,
    n_trees_detected = n_trees_detected,
    n_bridges = n_bridges
  ))
}

#' Cache for the resolved Overture Maps release string
#' @noRd
.overture_cache <- new.env(parent = emptyenv())

#' Resolve a usable Overture Maps release
#'
#' Overture publishes roughly monthly; a hard-coded release eventually 404s.
#' When `release` is `NULL` or `"auto"`, list the release prefixes in the
#' public S3 bucket and take the newest. The result is cached per session.
#'
#' @param con a live DuckDB connection with httpfs loaded.
#' @param release character or `NULL`.
#' @return character release string.
#' @noRd
resolve_overture_release <- function(con, release = NULL) {
  if (!is.null(release) && !identical(tolower(release[1]), "auto")) {
    return(as.character(release[1]))
  }
  if (!is.null(.overture_cache$release)) return(.overture_cache$release)

  rel <- tryCatch({
    # Listing with a delimiter returns one row per release prefix, which is
    # far cheaper than globbing the parquet files themselves.
    res <- DBI::dbGetQuery(con, paste0(
      "SELECT DISTINCT regexp_extract(file, 'release/([^/]+)/', 1) AS rel ",
      "FROM glob('s3://overturemaps-us-west-2/release/*/theme=base/",
      "type=infrastructure/*') ORDER BY rel DESC LIMIT 1"))
    if (nrow(res) > 0 && nzchar(res$rel[1])) res$rel[1] else NULL
  }, error = function(e) NULL)

  if (is.null(rel)) {
    stop("Could not determine the latest Overture Maps release. Pass one ",
         "explicitly via `overture_release`, e.g. \"2025-03-19.0\"; see ",
         "<https://github.com/OvertureMaps/data/releases> for the current ",
         "list.", call. = FALSE)
  }
  .overture_cache$release <- rel
  rel
}

#' Surface class code -> material name mapping for classified terrain
#' @noRd
surface_class_names <- function() {
  c("terrain", "greenspace", "roads", "sidewalks", "paths", "water", "sand")
}

#' Rasterize surface classes onto a template grid
#'
#' Codes: 1 ground, 2 greenspace, 3 roads, 4 sidewalks, 5 paths, 6 water,
#' 7 sand. Painting order (later overwrites earlier): greenspace, paths,
#' sidewalks, roads, sand fringe, water.
#' @noRd
classify_surface_raster <- function(template, ribbons_by_group, green_r,
                                    water_polys, sand_buffer) {
  cls <- template
  terra::values(cls) <- 1
  if (!is.null(green_r)) {
    g <- terra::resample(green_r, template, method = "near")
    cls <- terra::ifel(!is.na(g) & g == 1, 2, cls)
  }
  codes <- c(paths = 5, sidewalks = 4, roads = 3)
  for (nm in names(codes)) {
    if (is.null(ribbons_by_group[[nm]])) next
    r <- terra::rasterize(terra::vect(ribbons_by_group[[nm]]), template,
                          field = 1)
    cls <- terra::ifel(!is.na(r), codes[[nm]], cls)
  }
  if (!is.null(water_polys)) {
    if (is.numeric(sand_buffer) && sand_buffer > 0) {
      sand <- sf::st_buffer(water_polys, sand_buffer)
      rs <- terra::rasterize(terra::vect(sand), template, field = 1)
      cls <- terra::ifel(!is.na(rs), 7, cls)
    }
    rw <- terra::rasterize(terra::vect(water_polys), template, field = 1)
    cls <- terra::ifel(!is.na(rw), 6, cls)
  }
  cls
}

#' Flatten terrain under water bodies (arnis-style riparian/coastal water)
#'
#' Each water polygon gets one flat level: the 5th percentile of the terrain
#' beneath it, minus `drop`, so the water surface sits just below its banks.
#' @noRd
flatten_water_elev <- function(elev_r, water_polys, drop = 0.2) {
  wp <- sf::st_sf(wid = seq_along(water_polys), geometry = water_polys)
  wr <- terra::rasterize(terra::vect(wp), elev_r, field = "wid")
  wv <- terra::values(wr, mat = FALSE)
  ev <- terra::values(elev_r, mat = FALSE)
  for (i in unique(wv[!is.na(wv)])) {
    idx <- which(wv == i)
    lvl <- stats::quantile(ev[idx], 0.05, na.rm = TRUE, names = FALSE)
    if (!is.finite(lvl)) lvl <- 0
    ev[idx] <- lvl - drop
  }
  terra::values(elev_r) <- ev
  elev_r
}

#' Fetch Overture Maps water polygons for a WGS84 extent
#'
#' Native DuckDB parquet query against the Overture Maps base theme
#' (type "water"), same approach as `fetch_overture_roads()`.
#' @noRd
fetch_overture_water <- function(bbox_wgs,
                                 release = NULL,
                                 quiet = TRUE) {
  if (!requireNamespace("duckdb", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'duckdb' and 'DBI' are required for ",
         "`water = \"overture\"`. Install them with: ",
         "install.packages(c('duckdb', 'DBI')), or supply `water` as an ",
         "sf polygon layer.", call. = FALSE)
  }
  if (!quiet) cli::cli_alert_info("Fetching Overture water polygons ...")

  bb <- sf::st_bbox(bbox_wgs)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  res <- tryCatch({
    DBI::dbExecute(con, "INSTALL httpfs;  LOAD httpfs;")
    DBI::dbExecute(con, "INSTALL spatial; LOAD spatial;")
    DBI::dbExecute(con, "SET s3_region = 'us-west-2';")
    release <- resolve_overture_release(con, release)
    if (!quiet) cli::cli_alert_info(paste("Overture release:", release))
    s3_path <- sprintf(
      "s3://overturemaps-us-west-2/release/%s/theme=base/type=water/*",
      release
    )
    query <- sprintf("
      SELECT
        ST_AsWKB(geometry) AS wkb,
        subtype,
        class
      FROM read_parquet('%s', hive_partitioning = 1)
      WHERE bbox.xmax >= %f AND bbox.xmin <= %f
        AND bbox.ymax >= %f AND bbox.ymin <= %f
    ", s3_path,
      bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]])
    DBI::dbGetQuery(con, query)
  }, error = function(e) {
    warning("Overture Maps water query failed (", conditionMessage(e),
            "); continuing without water.", call. = FALSE)
    NULL
  })

  if (is.null(res) || nrow(res) == 0) return(NULL)
  geom <- sf::st_as_sfc(
    structure(as.list(res$wkb), class = "WKB"),
    crs = 4326
  )
  sf::st_sf(subtype = res$subtype, class = res$class, geometry = geom)
}

#' Fetch Overture Maps transportation segments for a WGS84 extent
#'
#' Native DuckDB parquet query against the Overture Maps S3 release,
#' following the same approach as `fetch_overture_network()` (block.R) but
#' keeping all road subtypes (including footways, cycleways, and paths) and
#' the `subclass` column so sidewalks can be identified.
#'
#' @param bbox_wgs sfc polygon in EPSG:4326.
#' @param release character. Overture Maps release string.
#' @param quiet logical.
#' @return sf of road segments with `class` and `subclass` columns, or NULL
#' when the query fails or returns nothing.
#' @noRd
fetch_overture_roads <- function(bbox_wgs,
                                 release = NULL,
                                 quiet = TRUE) {
  if (!requireNamespace("duckdb", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'duckdb' and 'DBI' are required for ",
         "`roads = \"overture\"`. Install them with: ",
         "install.packages(c('duckdb', 'DBI')), or supply `roads` as an ",
         "sf line layer.", call. = FALSE)
  }
  if (!quiet) cli::cli_alert_info("Fetching Overture transportation segments ...")

  bb <- sf::st_bbox(bbox_wgs)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  res <- tryCatch({
    DBI::dbExecute(con, "INSTALL httpfs;  LOAD httpfs;")
    DBI::dbExecute(con, "INSTALL spatial; LOAD spatial;")
    DBI::dbExecute(con, "SET s3_region = 'us-west-2';")
    release <- resolve_overture_release(con, release)
    if (!quiet) cli::cli_alert_info(paste("Overture release:", release))
    s3_path <- sprintf(
      "s3://overturemaps-us-west-2/release/%s/theme=transportation/type=segment/*",
      release
    )
    # Bridge/tunnel/level attributes moved between releases, so probe the
    # schema and select only the columns that exist. DESCRIBE returns plain
    # text, avoiding the nested-type conversion that `SELECT *` would force.
    have <- DBI::dbGetQuery(con, sprintf(
      "DESCRIBE SELECT * FROM read_parquet('%s', hive_partitioning = 1)",
      s3_path))$column_name
    # Rendering STRUCT/LIST columns as text needs the json extension, which
    # the R duckdb build cannot autoload; install it once, and fall back to
    # geometry-only roads (no bridges) when that is not possible.
    has_json <- isTRUE(tryCatch({
      DBI::dbExecute(con, "INSTALL json; LOAD json;")
      TRUE
    }, error = function(e) FALSE))

    extra <- c("ST_AsWKB(geometry) AS wkb", "class")
    if ("subclass" %in% have) extra <- c(extra, "subclass")
    if (has_json) {
      for (nm in c("subclass_rules", "road_flags", "level_rules")) {
        if (nm %in% have) {
          extra <- c(extra, sprintf("CAST(%s AS VARCHAR) AS %s", nm, nm))
        }
      }
    } else if (!quiet) {
      cli::cli_alert_warning(paste(
        "DuckDB 'json' extension unavailable; bridge and tunnel attributes",
        "will be skipped."))
    }
    if ("level" %in% have) extra <- c(extra, "level")
    query <- sprintf("
      SELECT %s
      FROM read_parquet('%s', hive_partitioning = 1)
      WHERE subtype = 'road'
        AND bbox.xmax >= %f AND bbox.xmin <= %f
        AND bbox.ymax >= %f AND bbox.ymin <= %f
    ", paste(extra, collapse = ", "), s3_path,
      bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]])
    DBI::dbGetQuery(con, query)
  }, error = function(e) {
    warning("Overture Maps query failed (", conditionMessage(e),
            "); continuing without roads.", call. = FALSE)
    NULL
  })

  if (is.null(res) || nrow(res) == 0) return(NULL)
  geom <- sf::st_as_sfc(
    structure(as.list(res$wkb), class = "WKB"),
    crs = 4326
  )
  txt <- function(nm) if (nm %in% names(res)) as.character(res[[nm]]) else ""
  blob <- paste(txt("subclass"), txt("subclass_rules"), txt("road_flags"))
  lvl <- rep(0L, nrow(res))
  if ("level" %in% names(res)) {
    lvl <- suppressWarnings(as.integer(res$level))
  } else if ("level_rules" %in% names(res)) {
    # VARCHAR cast renders structs as {'value': 1, ...}; accept quoted or
    # unquoted keys so either rendering parses.
    pat <- "['\"]?value['\"]?[[:space:]]*:[[:space:]]*-?[0-9]+"
    hit <- grepl(pat, res$level_rules)
    if (any(hit)) {
      m <- regmatches(res$level_rules[hit], regexpr(pat, res$level_rules[hit]))
      lvl[hit] <- suppressWarnings(
        as.integer(sub(".*:[[:space:]]*", "", m)))
    }
  }
  lvl[!is.finite(lvl)] <- 0L

  sf::st_sf(
    class = res$class,
    subclass = if ("subclass" %in% names(res)) res$subclass else NA_character_,
    is_bridge = grepl("bridge", blob, ignore.case = TRUE),
    is_tunnel = grepl("tunnel", blob, ignore.case = TRUE),
    level = lvl,
    geometry = geom
  )
}

#' Download an Esri World Imagery snapshot for a projected extent
#' @param xr,yr numeric length-2. x/y range in the scene CRS.
#' @param epsg integer. EPSG code of the scene CRS.
#' @param out_dir character. Directory for basemap.jpg.
#' @param max_px integer. Longest image side in pixels.
#' @return path to the written JPEG.
#' @noRd
fetch_basemap_image <- function(xr, yr, epsg, out_dir, max_px = 2048) {
  w <- diff(xr); h <- diff(yr)
  if (w >= h) {
    px_w <- max_px; px_h <- max(1L, round(max_px * h / w))
  } else {
    px_h <- max_px; px_w <- max(1L, round(max_px * w / h))
  }
  url <- paste0(
    "https://services.arcgisonline.com/ArcGIS/rest/services/",
    "World_Imagery/MapServer/export",
    "?f=image&format=jpg",
    "&bbox=", paste(c(xr[1], yr[1], xr[2], yr[2]), collapse = ","),
    "&bboxSR=", epsg, "&imageSR=", epsg,
    "&size=", px_w, ",", px_h
  )
  resp <- httr2::req_perform(httr2::request(url))
  body <- httr2::resp_body_raw(resp)
  if (length(body) < 1000) {
    stop("Imagery server returned an unexpectedly small response.")
  }
  path <- file.path(out_dir, "basemap.jpg")
  writeBin(body, path)
  path
}
