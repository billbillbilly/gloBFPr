# Infer screening-level road traffic inputs from OSM road attributes

Adds traffic volume and speed columns to an OSM road `sf` object based
on the `highway` class and optional `maxspeed`, `lanes`, or `oneway`
values. These values are intended for relative/scenario noise mapping
when measured traffic counts are unavailable, not calibrated regulatory
noise maps.

## Usage

``` r
infer_osm_traffic(
  roads,
  defaults = osm_noise_traffic_defaults(),
  use_maxspeed = FALSE,
  use_lanes = FALSE,
  use_oneway = FALSE,
  quiet = TRUE
)
```

## Arguments

- roads:

  `sf` line object containing at least a `highway` column.

- defaults:

  Data frame of road-class assumptions. Defaults to
  [`osm_noise_traffic_defaults()`](https://billbillbilly.github.io/gloBFPr/reference/osm_noise_traffic_defaults.md).

- use_maxspeed:

  Logical. If `TRUE`, parse numeric speeds from `maxspeed` and use them
  when available. Defaults to `FALSE` to match NoiseModelling's OSM
  import defaults.

- use_lanes:

  Logical. If `TRUE`, scale vehicle counts by lane count relative to a
  two-lane road. Defaults to `FALSE` to match NoiseModelling's OSM
  import defaults.

- use_oneway:

  Logical. If `TRUE`, halve traffic counts on OSM one-way roads.
  Defaults to `FALSE`.

- quiet:

  Logical. If `TRUE`, suppress informational messages.

## Value

An `sf` object with inferred traffic columns, including `speed_kmh`,
`light_veh_h`, `heavy_veh_h`, and NoiseModelling-style `LV_*`, `MV_*`,
`HGV_*`, `WAV_*`, `WBV_*`, and speed columns for day/evening/night.
