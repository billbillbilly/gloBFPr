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
