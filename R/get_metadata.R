#' get_metadata
#' @description
#' Returns a spatial grid (as an sf object) containing metadata and download URLs
#' for global 3D building footprint tiles (3D-GloBFP).
#'
#' @param test logic, Ignored during normal use; included for internal testing
#' purposes. Defaults to \code{FALSE}.
#'
#' @return sf a spatial polygon grid with attributes:
#' `id`, `gridID`, bounding box coordinates, and `download_url`.
#'
#' @details
#' The metadata of 3D Global Building Footprints (3D-GloBFP) dataset is uploaded on zenodo.
#' More detials about this dataset can to found [here](https://zenodo.org/records/15487037).
#'
#' The data is detailed in the following article
#'
#' @references
#' Che, Y., Li, X., Liu, X., Wang, Y., Liao, W., Zheng, X., Zhang, X., Xu, X.,
#' Shi, Q., Zhu, J., Zhang, H., Yuan, H., & Dai, Y. (2025).
#' 3D-GloBFP: the first global three-dimensional building footprint dataset.
#' Zenodo. https://doi.org/10.5281/zenodo.15487037
#'
#' Che Yangzi, Li Xuecao, Liu Xiaoping, Wang Yuhao, Liao Weilin, Zheng Xianwei,
#' Zhang Xucai, Xu Xiaocong, Shi Qian, Zhu Jiajun, Zhang Honghui, Yuan Hua, &
#' Dai Yongjiu (2024). 3D-GloBFP: the first global three-dimensional building
#' footprint dataset. Earth Syst. Sci. Data, 16, 5357-5374
#'
#' @examples
#' meta <- gloBFPr::get_metadata(test=TRUE)
#'
#' @importFrom dplyr "%>%"
#' @importFrom dplyr bind_rows
#' @importFrom dplyr rowwise
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom sf st_as_sf
#' @importFrom sf st_sfc
#' @importFrom sf st_polygon
#' @importFrom httr2 request
#' @importFrom httr2 request
#' @importFrom httr2 req_perform
#' @importFrom httr2 resp_body_json
#' @importFrom rlang .data
#'
#' @export

get_metadata <- function(test=FALSE) {
  if (test) {return(NULL)}
  info_extract <- function(name) {
    name <- sub("\\.zip$", "", name)
    parts <- unlist(strsplit(name, "_"))
    return(
      list(
        gridID = as.numeric(parts[1]),
        xmin = as.numeric(parts[2]),
        ymin = as.numeric(parts[3]),
        xmax = as.numeric(parts[4]),
        ymax = as.numeric(parts[5])
      )
    )
  }

  subdatasets <- c(
    28879733,
    28881749,
    28882700,
    28889813,
    28890593,
    28891631,
    28903454,
    28903853,
    28904453,
    28906499
  )
  result <- list()
  for (i in 1:length(subdatasets)) {
    article_id <- subdatasets[i]
    res <- httr2::request(paste0("https://api.figshare.com/v2/articles/", article_id)) %>%
      httr2::req_perform()
    data <- res %>% httr2::resp_body_json()
    file_info <- do.call(rbind, lapply(data$files, function(f) {
      info <- info_extract(f$name)
      data.frame(
        id = f$id,
        gridID = info['gridID'],
        bbox = c(info['xmin'], info['ymin'], info['xmax'], info['ymax']),
        download_url = f$download_url
      )
    }))
    result[[i]] <- file_info[!duplicated(file_info$gridID),]
  }
  result <- dplyr::bind_rows(result)

  grids_sf <- result %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      geometry = sf::st_sfc(
        sf::st_polygon(
          list(
            matrix(
              c(
                .data$bbox.xmin, .data$bbox.ymin,
                .data$bbox.xmin, .data$bbox.ymax,
                .data$bbox.xmax, .data$bbox.ymax,
                .data$bbox.xmax, .data$bbox.ymin,
                .data$bbox.xmin, .data$bbox.ymin
              ),
              ncol = 2,
              byrow = TRUE)
            )
          ), crs = 4326)) %>%
    dplyr::ungroup() %>%
    sf::st_as_sf()
  return(grids_sf)
}
