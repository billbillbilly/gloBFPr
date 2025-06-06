#' search_3dglobdf
#' @description
#' Search and retrieve 3D-GloBFP tiles that intersect a given bounding box or
#' area of interest, with options to return vector or raster outputs including
#' building polygons, binary presence rasters, and height-coded rasters.
#'
#' @param bbox `sf`, `sfc`, or a numeric vector (xmin, ymin, xmax, ymax)
#' defining the area of interest.
#' @param metadata sf. Typically output from [get_metadata()], containing tile
#' extents and download URLs.
#' @param crop logical. If `TRUE`, the resulting building footprint geometries
#' will be cropped to the input `bbox`. Default is `TRUE`.
#' @param out_type character. Default is `'poly'`.
#' Output type(s) to return. Options include:
#'   \itemize{
#'     \item `"poly"`: building footprints as an `sf` polygon object.
#'     \item `"binary_rast"`: binary `terra` raster where buildings = 1.
#'     \item `"graduated_rast"`: `terra` raster encoding building height values.
#'     \item `"rast"`: a named list with both binary and graduated rasters.
#'     \item `"all"`: a named list including the polygon layer and both raster layers.
#'   }
#' @param mask logical (optional). Default is `FALSE`. If `TRUE`, masks the
#' graduated raster using the building footprint layer. Only used when `out_type`
#' is `"graduated_rast"`, `"rast"`, or `"all"`.
#' @param cell_size numeric (optional). Default is 1. Only used when `out_type`
#' is `"graduated_rast"`, `"rast"`, or `"all"`.
#'
#' @return Varies based on `out_type`:
#' \itemize{
#'   \item If `"poly"`: an `sf` object of building footprints.
#'   \item If `"binary_rast"`: a binary `SpatRaster` (`terra`) indicating building presence.
#'   \item If `"graduated_rast"`: a quantitative `SpatRaster` of building heights.
#'   \item If `"rast"`: a named list with two `SpatRaster` objects: `binary` and `graduated`.
#'   \item If `"all"`: a named list with `poly` (sf), `binary`, and `graduated` rasters.
#' }
#'
#' @note
#' The downloading process may take some time, depending on the number and size
#' of building footprint tiles.
#'
#' This implementation relies on the current structure of the dataset as hosted on Figshare.
#' It may break if the dataset owner changes the file organization or metadata format.
#'
#' @references
#' Che Yangzi, Li Xuecao, Liu Xiaoping, Wang Yuhao, Liao Weilin, Zheng Xianwei,
#' Zhang Xucai, Xu Xiaocong, Shi Qian, Zhu Jiajun, Zhang Honghui, Yuan Hua, &
#' Dai Yongjiu (2024). 3D-GloBFP: the first global three-dimensional building
#' footprint dataset. Earth Syst. Sci. Data, 16, 5357-5374
#'
#' @examples
#' metadata <- gloBFPr::get_metadata(test=TRUE)
#' buildings <- gloBFPr::search_3dglobdf(bbox=c(-84.485519,45.636118,-84.462774,45.650639),
#'                      metadata=metadata)
#'
#' @importFrom sf st_bbox
#' @importFrom sf st_read
#' @importFrom sf st_crop
#' @importFrom sf st_cast
#' @importFrom sf st_is_valid
#' @importFrom sf st_intersects
#' @importFrom utils download.file
#' @importFrom utils unzip
#'
#' @export

search_3dglobdf <- function(bbox,
                            metadata,
                            crop=TRUE,
                            out_type='poly',
                            mask=FALSE,
                            cell_size=1) {
  if (missing(bbox) || missing(metadata)) {
    stop('bbox or metadata is missing')
  }

  if (inherits(metadata, "NULL")) {
    return(NULL)
  }

  # check type of bbox
  if (is.numeric(bbox) && length(bbox) == 4) {
    # convert bbox to sf when it is not a sf polygon
    bbox <- sf::st_as_sfc(
      sf::st_bbox(
        c(xmin = bbox[1],
          ymin = bbox[2],
          xmax = bbox[3],
          ymax = bbox[4]),
        crs = 4326
      )
    )
  }

  # Ensure input is in WGS84
  bbox <- sf::st_transform(bbox, 4326)

  # find all areas of spatial grid that intersect with bbox
  intersecting <- metadata[sf::st_intersects(metadata, bbox, sparse = FALSE), ]

  if (nrow(intersecting) == 0) {
    base::warning("No tiles intersect with the provided bbox.")
    return(NULL)
  }

  # download and load shapefiles
  result_list <- list()
  for (i in seq_len(nrow(intersecting))) {
    temp_zip <- tempfile(fileext = ".zip")
    utils::download.file(intersecting$download_url[i],
                  destfile = temp_zip,
                  mode = "auto",
                  quiet = TRUE)

    unzip_dir <- tempfile()
    utils::unzip(temp_zip, exdir = unzip_dir)

    # Find .shp file
    shp_files <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp_files) == 0) next

    sf_data <- tryCatch({
      sf::st_read(shp_files[1], quiet = TRUE)
    }, error = function(e) {
      base::message("Failed to read shapefile: ", shp_files[1])
      return(NULL)
    })
    sf_data <- sf::st_cast(sf_data, "POLYGON")

    if (!is.null(sf_data)) {
      result_list[[length(result_list) + 1]] <- sf_data
    }
    unlink(c(temp_zip, unzip_dir), recursive = TRUE)
  }

  # Combine all into one sf object
  all_data <- do.call(rbind, result_list)

  # crop the data if 'crop' == true
  if (crop) {
    # all_data <- sf::st_make_valid(all_data)
    valid <- sf::st_is_valid(all_data)
    all_data <- all_data[valid, ]
    all_data <- suppressWarnings(all_data[sf::st_intersects(all_data, bbox, sparse = FALSE), ])
    all_data <- suppressWarnings(sf::st_crop(all_data, bbox))
  }

  # export data as sf polygons
  if(out_type == 'poly') {
    return(all_data)
  }

  # auto-generate raster outputs
  if (out_type %in% c("binary_rast", "graduated_rast", "rast", "all")) {
    binary <- rasterize_binary(all_data, bbox, res=cell_size)
    if (isTRUE(mask)) {
      graduated <- rasterize_height(all_data, bbox, res=cell_size, mask=binary)
    } else {
      graduated <- rasterize_height(all_data, bbox, res=cell_size)
    }

    if (out_type == "binary_rast") {return(binary)}
    if (out_type == "graduated_rast") {return(graduated)}
    if (out_type == "rast") {return(list(binary = binary, graduated = graduated))}
    if (out_type == "all") {return(list(poly = all_data, binary = binary, graduated = graduated))}
  }

  stop("Invalid out_type specified.")
}
