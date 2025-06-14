# gloBFPr <a href="https://github.com/billbillbilly/gloBFPr/"><img src="logo.svg" align="right" height="150" alt="forestdata website" /></a>

<!-- badges: start -->
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/billbillbilly/gloBFPr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/billbillbilly/gloBFPr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Access and analyze the Building Footprint Datasets.

## Overview
The `gloBFPr` package allows R users to search, download, and process global 
building footprint tiles with associated height information, derived from the 
3D-GloBFP dataset published by Che et al. (2024, 2025). The data is hosted on 
Zenodo and covers global urban areas in shapefile format. The package will look 
to include more global Building Dataset in the future.

<img src="images/cover.png" align="center" width="90%">

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

#### Output examples:
<p align="center">
  <img src="images/BFshp.png?raw=true" width="45%">
&nbsp; &nbsp; &nbsp; &nbsp;
  <img src="images/BHshp.png?raw=true" width="45%">
</p>

<p align="center">
  <img src="images/BF.png?raw=true" width="30%">
&nbsp; &nbsp; &nbsp; &nbsp;
  <img src="images/BH.png?raw=true" width="30%">
&nbsp; &nbsp; &nbsp; &nbsp;
  <img src="images/croppedBH.png?raw=true" width="30%">
</p>

## Note
The downloading process may take some time, depending on the number and size
of building footprint tiles.

This implementation relies on the current structure of the dataset as hosted on Figshare.
It may break if the dataset owner changes the file organization or metadata format.

Please read the function documentation carefully. The dataset may require proper citation when used.

## Issues and bugs
If you discover a bug not associated with connection to the API that is
not already a [reported
issue](https://github.com/billbillbilly/gloBFPr/issues), please [open
a new issue](https://github.com/billbillbilly/gloBFPr/issues/new)
providing a reproducible example.
