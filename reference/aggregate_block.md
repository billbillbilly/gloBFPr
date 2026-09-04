# aggregate_block

Aggregate building-level metrics to block level using the output of
[`generate_block()`](https://billbillbilly.github.io/gloBFPr/reference/generate_block.md).
Additive quantities (areas, volumes, population) are summed by default;
per-building indices (shape metrics, elongation ratios, etc.) are
averaged. Both defaults can be overridden per column via `.fns`. Two
derived metrics are always added: `n_buildings` (count of buildings per
block) and `coverage_ratio` (total building footprint area / block
area).

## Usage

``` r
aggregate_block(
  block_output,
  .fns = NULL,
  population = FALSE,
  population_year = 2025,
  residential = FALSE,
  residential_year = 2020,
  quiet = FALSE
)
```

## Arguments

- block_output:

  list. The named list returned by
  [`generate_block()`](https://billbillbilly.github.io/gloBFPr/reference/generate_block.md),
  containing `$blocks` (sf polygons) and `$buildings` (sf with
  `block_id`).

- .fns:

  named list. Optional overrides mapping column names to aggregation
  functions, e.g. `list(vol = max, Height = median)`. Any column not
  named here uses the built-in default (sum or mean as described above).

- population:

  logical. If `TRUE`, fetch GHSL population at block level and add a
  `pop_total` column. Default `FALSE`.

- population_year:

  integer. GHSL population year to use when `population = TRUE`. One of
  1975, 1980, ..., 2025, 2030. Default `2025`.

- residential:

  logical. If `TRUE`, fetch GHS built-up surface rasters at block level
  and compute `res_prop` (residential built-up surface fraction per
  block). Default `FALSE`. Note: per-building `res` flags are excluded
  from block-level aggregation - use this parameter for block-level
  residential proportion instead.

- residential_year:

  integer. GHS built-up surface year to use when `residential = TRUE`.
  One of 1975, 1980, ..., 2020, 2025, 2030. Default `2020`.

- quiet:

  logical. If `TRUE`, suppress cli messages. Default `FALSE`.

## Value

The `blocks` sf object from `block_output` with one column per
aggregated metric, plus `n_buildings` and `coverage_ratio`. When
`population = TRUE`, also includes `pop_total`. When
`residential = TRUE`, also includes `res_prop`.
