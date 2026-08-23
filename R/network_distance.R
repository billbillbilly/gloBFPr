#### Street-network routing helpers ####
#
# Shared machinery for computing shortest-path distances along real road/path
# centre lines instead of straight-line ("as the crow flies") distances.
#
# The pipeline is:
#   1. `get_osm_walk_network()`   - download walkable OSM ways for a bbox.
#   2. `build_network_graph()`    - node the line work and build a weighted
#                                   `igraph` where edge weights are segment
#                                   lengths in projected metres.
#   3. `snap_to_network()`        - attach off-network points (building
#                                   centroids, green-space pixels) to their
#                                   nearest graph node and record the length of
#                                   the straight-line connector ("leg").
#   4. `network_distance_to_targets()` - one Dijkstra pass from the origin node,
#                                   then origin leg + path + target leg.

#' Download a walkable OSM street network for a bounding box
#'
#' Returns an `sf` line layer of ways that a pedestrian can plausibly use.
#' Motorways, trunk roads, and their link roads are dropped, as are ways
#' explicitly tagged `foot=no` or `foot=private`. Returns `NULL` (with a
#' warning) rather than erroring when the download fails or nothing is found,
#' so callers can fall back to Euclidean distance.
#'
#' @noRd
get_osm_walk_network <- function(bbox,
                                 target_crs,
                                 overpass_url = "https://overpass-api.de/api/interpreter",
                                 timeout = 180,
                                 quiet = FALSE) {
  bbox <- as_wgs84_bbox_vector(bbox)
  excluded <- paste(
    c("motorway", "motorway_link", "trunk", "trunk_link",
      "construction", "proposed", "raceway", "bus_guideway", "escape"),
    collapse = "|"
  )
  query <- sprintf(
    paste0(
      "[out:json][timeout:%d];",
      "(",
      "way[\"highway\"][\"highway\"!~\"^(%s)$\"][\"foot\"!~\"^(no|private)$\"](%f,%f,%f,%f);",
      ");",
      "out tags geom;"
    ),
    timeout,
    excluded,
    bbox[2], bbox[1], bbox[4], bbox[3]
  )

  roads <- tryCatch({
    response <- httr2::request(overpass_url) |>
      httr2::req_body_form(data = query) |>
      httr2::req_perform()
    body <- httr2::resp_body_json(response, simplifyVector = FALSE)
    osm_roads_json_to_sf(body, target_crs)
  }, error = function(e) {
    if (!quiet) {
      cli::cli_alert_warning(
        paste0("Could not download the OSM street network: ", conditionMessage(e))
      )
    }
    NULL
  })

  roads
}

#' Build a weighted, undirected graph from road centre lines
#'
#' Vertices of the input line work are snapped onto a grid of `10^-digits`
#' metres so that ways sharing an intersection collapse onto a single graph
#' node. Edge weights are planar segment lengths, so the input must already be
#' in a projected (metric) CRS.
#'
#' Because off-network points are attached to graph *vertices*, long unbroken
#' OSM ways would otherwise force a mid-block building to snap to a distant
#' street corner. The line work is therefore densified first, which bounds the
#' snapping error by `max_segment` metres.
#'
#' @param max_segment numeric. Maximum spacing between graph vertices in
#'   metres. Smaller values are more accurate but produce a larger graph.
#' @return A list with `graph` (igraph), `coords` (node coordinate matrix), and
#'   `nodes` (node `sfc` used for nearest-node lookups), or `NULL` if no usable
#'   line work was supplied.
#' @noRd
build_network_graph <- function(roads, digits = 1, max_segment = 20) {
  if (is.null(roads)) {
    return(NULL)
  }

  geom <- sf::st_geometry(roads)
  if (length(geom) == 0L) {
    return(NULL)
  }
  crs <- sf::st_crs(geom)
  if (is.na(crs)) {
    stop("`network` must have a coordinate reference system.", call. = FALSE)
  }

  if (is.finite(max_segment) && max_segment > 0) {
    geom <- sf::st_segmentize(geom, dfMaxLength = max_segment)
  }

  coords <- sf::st_coordinates(geom)
  if (is.null(coords) || nrow(coords) < 2) {
    return(NULL)
  }

  xy <- coords[, 1:2, drop = FALSE]
  # Trailing columns of st_coordinates() identify the parent geometry: one
  # column for LINESTRING (L1), two for MULTILINESTRING (L1, L2).
  id_cols <- coords[, -(1:2), drop = FALSE]
  line_id <- if (ncol(id_cols) == 1L) {
    as.integer(id_cols[, 1])
  } else {
    as.integer(factor(do.call(paste, as.data.frame(id_cols))))
  }

  # Node the line work: identical (rounded) coordinates become one vertex.
  key <- paste(round(xy[, 1], digits), round(xy[, 2], digits), sep = "_")
  node_id <- match(key, unique(key))
  # !duplicated() yields rows in first-appearance order, which is node id order.
  node_xy <- xy[!duplicated(node_id), , drop = FALSE]

  n_vertex <- length(node_id)
  from <- node_id[-n_vertex]
  to <- node_id[-1L]
  # Only consecutive vertices of the same line form an edge, and zero-length
  # edges (duplicate vertices) are dropped.
  keep <- (line_id[-n_vertex] == line_id[-1L]) & (from != to)
  from <- from[keep]
  to <- to[keep]
  if (length(from) == 0L) {
    return(NULL)
  }

  weight <- sqrt(
    (node_xy[from, 1] - node_xy[to, 1])^2 +
      (node_xy[from, 2] - node_xy[to, 2])^2
  )

  g <- igraph::make_empty_graph(n = nrow(node_xy), directed = FALSE)
  g <- igraph::add_edges(g, rbind(from, to), attr = list(weight = weight))

  nodes <- sf::st_as_sf(
    data.frame(node = seq_len(nrow(node_xy)), x = node_xy[, 1], y = node_xy[, 2]),
    coords = c("x", "y"),
    crs = crs
  )

  list(graph = g, coords = node_xy, nodes = nodes, crs = crs)
}

#' Snap off-network points to their nearest graph node
#'
#' @return A list with `index` (node ids) and `leg` (straight-line distance in
#'   metres from each point to its snapped node).
#' @noRd
snap_to_network <- function(graph_obj, xy) {
  xy <- as.matrix(xy)
  pts <- sf::st_as_sf(
    data.frame(x = xy[, 1], y = xy[, 2]),
    coords = c("x", "y"),
    crs = graph_obj$crs
  )
  idx <- sf::st_nearest_feature(pts, graph_obj$nodes)
  leg <- sqrt(
    (xy[, 1] - graph_obj$coords[idx, 1])^2 +
      (xy[, 2] - graph_obj$coords[idx, 2])^2
  )
  list(index = idx, leg = leg)
}

#' Reproject a coordinate matrix into the routing graph's CRS
#'
#' `get_buffer()` picks a UTM zone from each unit's own bounding box, which can
#' differ from the zone used to build the graph for buildings near a zone
#' boundary. This is a no-op in the common case where the two agree.
#' @noRd
to_graph_crs <- function(xy, from_crs, graph_obj) {
  if (is.na(from_crs) || from_crs == graph_obj$crs) {
    return(xy)
  }
  pts <- sf::st_as_sf(
    data.frame(x = xy[, 1], y = xy[, 2]),
    coords = c("x", "y"),
    crs = from_crs
  )
  sf::st_coordinates(sf::st_transform(pts, graph_obj$crs))[, 1:2, drop = FALSE]
}

#' Shortest network distance from one origin to the nearest of many targets
#'
#' Runs a single Dijkstra pass from the origin's snapped node to every node in
#' the graph, then adds the two off-network legs. Because a network distance is
#' never shorter than the corresponding straight-line distance, targets whose
#' Euclidean distance already exceeds the best network distance found for the
#' Euclidean-nearest target cannot win and are pruned before snapping.
#'
#' @param origin_xy 1x2 matrix of the origin coordinate (projected CRS).
#' @param target_xy nx2 matrix of candidate target coordinates.
#' @param max_candidates integer. Upper bound on how many targets are routed
#'   after pruning, evaluated nearest-first.
#' @return A list with `distance` (metres, `NA_real_` if unreachable) and
#'   `index` (row of `target_xy` that was nearest, `NA_integer_` if none).
#' @noRd
network_distance_to_targets <- function(graph_obj,
                                        origin_xy,
                                        target_xy,
                                        max_candidates = 5000L) {
  empty <- list(distance = NA_real_, index = NA_integer_)
  if (is.null(graph_obj) || is.null(target_xy) || nrow(target_xy) == 0L) {
    return(empty)
  }

  origin_xy <- matrix(as.numeric(origin_xy)[1:2], ncol = 2)
  origin <- snap_to_network(graph_obj, origin_xy)

  node_dist <- igraph::distances(
    graph_obj$graph,
    v = origin$index,
    weights = igraph::E(graph_obj$graph)$weight
  )[1, ]

  euclid <- sqrt(
    (target_xy[, 1] - origin_xy[1, 1])^2 +
      (target_xy[, 2] - origin_xy[1, 2])^2
  )
  ord <- order(euclid)

  # Route the Euclidean-nearest target first to obtain an upper bound.
  seed <- ord[1L]
  seed_snap <- snap_to_network(graph_obj, target_xy[seed, , drop = FALSE])
  best <- origin$leg + node_dist[seed_snap$index] + seed_snap$leg
  best_index <- if (is.finite(best)) seed else NA_integer_
  if (!is.finite(best)) best <- Inf

  candidates <- ord[euclid[ord] <= best]
  if (length(candidates) == 0L) {
    return(list(
      distance = if (is.finite(best)) as.numeric(best) else NA_real_,
      index = best_index
    ))
  }
  if (length(candidates) > max_candidates) {
    candidates <- candidates[seq_len(max_candidates)]
  }

  snapped <- snap_to_network(graph_obj, target_xy[candidates, , drop = FALSE])
  totals <- origin$leg + node_dist[snapped$index] + snapped$leg

  if (any(is.finite(totals))) {
    winner <- which.min(totals)
    if (totals[winner] < best) {
      best <- totals[winner]
      best_index <- candidates[winner]
    }
  }

  if (!is.finite(best)) {
    return(empty)
  }
  list(distance = as.numeric(best), index = as.integer(best_index))
}

#' Resolve the `network` argument of `get_dng()` into a routing graph
#'
#' Accepts `NULL` (no routing), `"osm"` (download), or a user-supplied `sf`
#' line layer. Returns `NULL` when routing is unavailable so the caller can
#' fall back to Euclidean distance.
#' @noRd
resolve_dng_network <- function(network,
                                extent,
                                utm_crs,
                                overpass_url = "https://overpass-api.de/api/interpreter",
                                timeout = 180,
                                quiet = FALSE) {
  if (is.null(network)) {
    return(NULL)
  }

  roads <- NULL
  if (inherits(network, c("sf", "sfc"))) {
    roads <- sf::st_transform(sf::st_as_sf(network), utm_crs)
  } else if (is.character(network) && identical(tolower(network[1]), "osm")) {
    if (!quiet) cli::cli_alert_info("Downloading OSM street network ...")
    roads <- get_osm_walk_network(
      bbox = get_bbox(extent),
      target_crs = utm_crs,
      overpass_url = overpass_url,
      timeout = timeout,
      quiet = quiet
    )
  } else {
    stop(
      "`network` must be NULL, \"osm\", or an `sf` object of road/path lines.",
      call. = FALSE
    )
  }

  if (is.null(roads) || length(sf::st_geometry(roads)) == 0L) {
    if (!quiet) {
      cli::cli_alert_warning(
        "No street network available; falling back to straight-line distance."
      )
    }
    return(NULL)
  }

  graph_obj <- build_network_graph(roads)
  if (is.null(graph_obj) && !quiet) {
    cli::cli_alert_warning(
      "Street network had no routable segments; falling back to straight-line distance."
    )
  }
  graph_obj
}
