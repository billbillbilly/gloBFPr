# Default OSM traffic assumptions for screening-level road noise

Returns the built-in lookup table used by
[`infer_osm_traffic()`](https://billbillbilly.github.io/gloBFPr/reference/infer_osm_traffic.md)
when observed traffic counts are unavailable. The defaults mirror the
road-category assumptions embedded in NoiseModelling's
`Import_OSM.groovy`, which cites the Good Practice Guide for Strategic
Noise Mapping and the Production of Associated Data on Noise Exposure,
version 2. They should still be replaced with local traffic counts
whenever available.

## Usage

``` r
osm_noise_traffic_defaults()
```

## Value

A data frame with OSM highway class, speed, and hourly light/heavy
vehicle assumptions.
