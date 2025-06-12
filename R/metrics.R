#' get_metrics
#' @description
#' calculate
#' @param x sf. building footprint polygon, typically output from [get_3dglobdf()]
#' @param metric_type character. Default is `'2d'`.
#' Metric type(s) to return. Options include:
#'   \itemize{
#'     \item `"2d"`: building footprints as an `sf` polygon object.
#'     \item `"3d"`: binary `terra` raster where buildings = 1.
#'     \item `"pop"`: `terra` raster encoding building height values.
#'     \item `"acc"`: a named list with both binary and graduated rasters.
#'     \item `"all"`: a named list including the polygon layer and both raster layers.
#'   }
#' @param out_type character. Default is `'poly'`
#' @param resolution numeric. Default is 100. Only required when `out_type`
#' is `"rast"`.
#'
#' @return `sf` or stacked `rasterslayer`
#'
#' @details
#'
#'
#' @reference
#' Anna Labetski, Stelios Vitalis, Filip Biljecki, Ken Arroyo Ohori &
#' Jantien Stoter (2023): 3D building metrics for urban morphology.
#' International Journal of Geographical Information Science, 37(1):
#' 36-67. DOI: 10.1080/13658816.2022.2103818
#'
#' Melchiorri, M., Freire, S., Schiavina, M. et al.
#' The Multi-temporal and Multi-dimensional Global Urban Centre Database to
#' Delineate and Analyse World Cities. Sci Data 11, 82 (2024).
#' https://doi.org/10.1038/s41597-023-02691-1

get_metrics <- function(x, metric_type = '2d', out_type = 'poly') {
  if (missing(x)) {
    stop('Please input building footprint polygon')
  }

  updated <- add_2d_3d(x, metric_type)
  if (metric_type == 'all' || metric_type == 'pop') {
    updated <- add_pop(updated)
  }
}
