# gloBFPr
Access and analyze the Building Footprint Datasets.

## Overview
The gloBFPr package allows R users to search, download, and process global 
building footprint tiles with associated height information, derived from the 
3D-GloBFP dataset published by Che et al. (2024, 2025). The data is hosted on 
Zenodo and covers global urban areas in shapefile format. The package will look 
to include more global Building Dataset in the future.

## Features
📦 Access tiled metadata of 3D-GloBFP dataset

🔍 Search tiles by bounding box (BBOX) or area of interest

⬇️ Download only the necessary files

🏗️ Retrieve building polygons and height attributes

🌍 Generate rasters of binary presence or graduated height

🗺️ Output spatial data in sf or terra raster format

## Installation
Install the development version:
```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub
devtools::install_github("billbillbilly/gloBFPr")
```

The package will be on CRAN soon.

## Usage
1. Load metadata

```r
library(gloBFPr)
metadata <- get_metadata()
```

2. Search and download data by bounding box

```r
bbox <- c(-83.065644,42.333792,-83.045217,42.346988)
buildings_list <- search_3dglobdf(bbox = bbox,
                                  metadata = metadata, 
                                  out_type = "all", 
                                  # mask = TRUE,
                                  cell_size = 1)
```
This will return a list containing:
- poly: an sf object of 3D building footprints
- binary: a binary raster of building presence
- graduated: a raster representing building height in meters

Specify `cell_size = 1` to generate raster layers with 1-meter resolution, 
ensuring detailed spatial representation of building geometries within 
the defined area of interest.

Setting `mask = TRUE` ensures the height raster is masked by the building footprints.

## Note
The downloading process may take some time, depending on the number and size
of building footprint tiles.

This implementation relies on the current structure of the dataset as hosted on Figshare.
It may break if the dataset owner changes the file organization or metadata format.

## Citation

If you use the 3D-GloBFP dataset, please cite:

Che Yangzi, Li Xuecao, Liu Xiaoping, Wang Yuhao, Liao Weilin, Zheng Xianwei,
Zhang Xucai, Xu Xiaocong, Shi Qian, Zhu Jiajun, Zhang Honghui, Yuan Hua, &
Dai Yongjiu (2024). 3D-GloBFP: the first global three-dimensional building
footprint dataset. Earth Syst. Sci. Data, 16, 5357-5374

Che, Y., Li, X., Liu, X., Wang, Y., Liao, W., Zheng, X., Zhang, X., Xu, X.,
Shi, Q., Zhu, J., Zhang, H., Yuan, H., & Dai, Y. (2025).
3D-GloBFP: the first global three-dimensional building footprint dataset [Data set].
Zenodo. https://doi.org/10.5281/zenodo.15487037
