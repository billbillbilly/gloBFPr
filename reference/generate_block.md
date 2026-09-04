# generate_block

Cluster given buildings into blocks based on street network. Uses a
two-stage approach: (1) vector polygonization of the road network for
well-formed areas, then (2) a raster fallback for buildings that fall in
network gaps or dead-end pockets. Blocks with no buildings are dropped.

Before polygonization, dual carriageways (motorways and trunk roads
represented as parallel lines) are simplified to single centrelines to
prevent artificially narrow slivers between them from being
misidentified as blocks. The simplification approach is conceptually
adapted from UrbanWaterBlocks (Yin et al., 2025).

## Usage

``` r
generate_block(
  x,
  network = NULL,
  network_source = c("overture", "osm"),
  overture_release = NULL,
  res = 2,
  min_block_area = 500,
  dc_highway_types = c("motorway", "trunk"),
  dc_overlap_threshold = 0.7,
  quiet = FALSE
)
```

## Arguments

- x:

  sf. Building footprint polygons, typically output from
  [`search_3dglobdf()`](https://billbillbilly.github.io/gloBFPr/reference/search_3dglobdf.md).

- network:

  sf or character. Optional road network or path to one. When supplied,
  `network_source` is ignored.

- network_source:

  character. Source for automatic network fetching when
  `network = NULL`. Either `"overture"` (Overture Maps via DuckDB
  parquet query, default, requires the `duckdb` and `DBI` packages) or
  `"osm"` (OpenStreetMap via the Overpass API, slower for large areas).

- overture_release:

  character or `NULL`. Overture Maps release string used when
  `network_source = "overture"`, e.g. `"2025-03-19.0"`. The default
  `NULL` (or `"auto"`) queries the public bucket for the newest release
  and caches it for the session. See
  <https://github.com/OvertureMaps/data/releases> for available
  releases.

- res:

  numeric. Raster resolution in metres for the fallback stage. Default
  `2`.

- min_block_area:

  numeric. Minimum block area in m^2 below which polygons are treated as
  slivers and merged into neighbours (raster stage) or dropped
  (polygonize stage). Default `500`.

- dc_highway_types:

  character vector. Highway class values treated as dual carriageway
  candidates (OSM `highway` tag or Overture `class` column). Default
  `c("motorway", "trunk")`.

- dc_overlap_threshold:

  numeric. Minimum overlap fraction (0-1) for a line to be considered a
  duplicate carriageway and removed. Default `0.7`.

- quiet:

  logical. If `TRUE`, suppress cli messages. Default `FALSE`.

## Value

A named list with two elements:

- `blocks`:

  An `sf` polygon object, one row per block, with a `block_id` column.
  CRS matches the input `x`.

- `buildings`:

  The input `x` with an added `block_id` integer column linking each
  building to its block. Buildings that cannot be assigned to any block
  receive `NA`.

## References

Yin, H., et al. (2025). UrbanWaterBlocks: A python tool for block-based
urban water management. Sustainable Cities and Society.
