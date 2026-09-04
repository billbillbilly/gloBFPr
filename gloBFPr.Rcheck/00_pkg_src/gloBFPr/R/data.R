#' Test 3D-GloBFP dataset
#'
#' A sample dataset containing simplified 3D building footprint information for
#' demonstration and testing purposes.
#'
#' @format A data frame with 369 rows and 3 variables:
#' \describe{
#'   \item{id}{Numeric. Unique identifier for each building.}
#'   \item{Height}{Numeric. Estimated height of the building in meters.}
#'   \item{geometry}{sfc_POLYGON. The building footprint geometry in simple feature (sf) format.}
#' }
#' @source
#' Che Yangzi, Li Xuecao, Liu Xiaoping, Wang Yuhao, Liao Weilin, Zheng Xianwei,
#' Zhang Xucai, Xu Xiaocong, Shi Qian, Zhu Jiajun, Zhang Honghui, Yuan Hua, &
#' Dai Yongjiu (2024). 3D-GloBFP: the first global three-dimensional building
#' footprint dataset. Earth Syst. Sci. Data, 16, 5357-5374
#'
#' @examples
#' data(globfp_example)
#' head(globfp_example)
"globfp_example"

#' Example DEM for the 3D-GloBFP sample
#'
#' A `terra::PackedSpatRaster` digital elevation model cropped to the bounding
#' box of `globfp_example`. Convert it with `terra::rast()` before analysis.
#'
#' @format A `terra::PackedSpatRaster` with one layer named `dem`.
#'
#' @source OpenTopography / dsmSearch elevation data, downloaded for the
#' `globfp_example` extent.
#'
#' @examples
#' data(globfp_example_dem)
#' dem <- terra::rast(globfp_example_dem)
#' dem
"globfp_example_dem"

#' Example canopy height raster for the 3D-GloBFP sample
#'
#' A `terra::PackedSpatRaster` canopy height map cropped to the bounding box of
#' `globfp_example` and aggregated from the source resolution for lightweight
#' examples. Convert it with `terra::rast()` before analysis.
#'
#' @format A `terra::PackedSpatRaster` with one layer named `canopy_height`.
#'
#' @source metaCHM canopy height data via dsmSearch, downloaded for the
#' `globfp_example` extent.
#'
#' @examples
#' data(globfp_example_canopy_height)
#' canopy_height <- terra::rast(globfp_example_canopy_height)
#' canopy_height
"globfp_example_canopy_height"

#' GHSL Tile Index
#'
#' A simple-feature tile index used internally to identify GHSL raster tiles
#' that intersect a requested area of interest.
#'
#' @format An `sf` object with GHSL tile identifiers and polygon geometry.
#'
#' @keywords internal
"ghsl_tiles"
