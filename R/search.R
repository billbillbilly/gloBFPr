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
#' will be cropped to the input `bbox`. Default is `FALSE`.
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
#' @importFrom sf st_intersects
#' @importFrom utils download.file
#' @importFrom utils unzip
#' @importFrom dplyr bind_rows
#' @importFrom cli cli_alert_info
#' @importFrom cli cli_alert_success
#'
#' @export

search_3dglobdf <- function(bbox,
                            metadata,
                            crop=FALSE,
                            out_type='poly',
                            mask=FALSE,
                            cell_size=1) {
  if (missing(bbox) || missing(metadata)) {
    stop('bbox or metadata is missing')
  }

  if (inherits(metadata, "NULL")) {
    return(NULL)
  }

  start_time <- Sys.time()

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
  d_mode <- 'auto'
  # check os
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    d_mode <- 'wb'
  }

  cli::cli_alert_info('Start downloading data ...')
  # Store the original 'timeout' option and ensure it's reset upon function exit
  original_timeout <- getOption('timeout')
  on.exit(options(timeout = original_timeout), add = TRUE)
  options(timeout=9999)
  for (i in seq_len(nrow(intersecting))) {
    temp_zip <- tempfile(fileext = ".zip")
    utils::download.file(intersecting$download_url[i],
                  destfile = temp_zip,
                  # method = method,
                  mode = d_mode,
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

    if (!is.null(sf_data)) {
      result_list[[length(result_list) + 1]] <- sf_data
    }
    unlink(c(temp_zip, unzip_dir), recursive = TRUE)
  }
  cli::cli_alert_success('Finished downloading and loading')
  cli::cli_alert_info('Start processing data ...')
  result_list <- lapply(result_list, function(x) {
    #x <- sf::st_cast(x, "POLYGON")  # ensure same geometry type
    x <- x[, intersect(names(x), names(result_list[[1]]))]  # keep common columns only
    return(x)
  })
  # Combine all into one sf object
  all_data <- dplyr::bind_rows(result_list)
  all_data <- all_data[,c('Height','geometry')]

  utm_crs <- get_utm_crs(bbox)
  bbox_proj <- sf::st_transform(bbox, crs = utm_crs)
  all_data <- sf::st_transform(all_data, crs = utm_crs)
  # sf_data <- sf_data[!sf::st_is_empty(sf_data), ]
  all_data <- suppressWarnings(all_data[sf::st_intersects(all_data, bbox_proj, sparse = FALSE), ])
  cli::cli_alert_success('Found building footprints within bbox')
  # crop the data if 'crop' == true
  if (crop) {
    #all_data <- suppressWarnings(all_data[sf::st_intersects(all_data, bbox, sparse = FALSE), ])
    cli::cli_alert_info('Start cropping data ...')
    all_data <- suppressWarnings(sf::st_crop(all_data, bbox_proj))
    cli::cli_alert_success('Cropped building footprints using bbox')
  }

  # export data as sf polygons
  if(out_type == 'poly') {
    end_time <- Sys.time()
    process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    cli::cli_alert_success(paste0("Completed. Time taken: ", base::round(process_time), " seconds."))
    return(all_data)
  }

  # auto-generate raster outputs
  if (out_type %in% c("binary_rast", "graduated_rast", "rast", "all")) {
    if (isTRUE(mask) || out_type == 'binary_rast' || out_type == 'all') {
      cli::cli_alert_info('Generate binary raster ...')
      binary <- rasterize_binary(all_data, bbox, res=cell_size)
      cli::cli_alert_success('Generated binary building footprints raster')
    }

    if (isTRUE(mask) && out_type != "binary_rast") {
      cli::cli_alert_info('Mask building height raster ...')
      graduated <- rasterize_height(all_data, bbox, res=cell_size, mask=binary)
      cli::cli_alert_success('Generated masked building footprints raster')
    } else {
      cli::cli_alert_info('Generate binary raster ...')
      graduated <- rasterize_height(all_data, bbox, res=cell_size)
      cli::cli_alert_success('Generated building height raster')
    }

    if (out_type == "binary_rast") {
      end_time <- Sys.time()
      process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      cli::cli_alert_success(paste0("Completed. Time taken: ", base::round(process_time), " seconds."))
      return(binary)
    }
    if (out_type == "graduated_rast") {
      end_time <- Sys.time()
      process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      cli::cli_alert_success(paste0("Completed. Time taken: ", base::round(process_time), " seconds."))
      return(graduated)
    }
    if (out_type == "rast") {
      end_time <- Sys.time()
      process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      cli::cli_alert_success(paste0("Completed. Time taken: ", base::round(process_time), " seconds."))
      return(list(binary = binary, graduated = graduated))
    }
    if (out_type == "all") {
      end_time <- Sys.time()
      process_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      cli::cli_alert_success(paste0("Completed. Time taken: ", base::round(process_time), " seconds."))
      return(list(poly = all_data, binary = binary, graduated = graduated))
    }
  }

  stop("Invalid out_type specified.")
}
