pkgname <- "gloBFPr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('gloBFPr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("get_3d_world")
### * get_3d_world

flush(stderr()); flush(stdout())

### Name: get_3d_world
### Title: get_3d_world
### Aliases: get_3d_world

### ** Examples





cleanEx()
nameEx("get_era5_met")
### * get_era5_met

flush(stderr()); flush(stdout())

### Name: get_era5_met
### Title: Fetch ERA5 conditions for OpenFOAM boundary conditions
### Aliases: get_era5_met

### ** Examples

## Not run: 
##D Sys.setenv(CDS_API_KEY = "your-token-here")
##D 
##D # Summer evening in Detroit
##D met <- get_era5_met(
##D   lon               = -83.05,
##D   lat               =  42.34,
##D   datetime          = "2023-07-15 22:00"
##D )
##D 
##D ## ---- Wind simulation ------------------------------------------------
##D prepare_foam_case(
##D   case_dir       = "path/to/case",
##D   stl_file       = "path/to/buildings.stl",
##D   domain         = list(xmin = 0, xmax = 500, ymin = 0, ymax = 500,
##D                         zmin = 0, zmax = 200),
##D   inlet_velocity = met$inlet_velocity,
##D   z_ref          = met$z_ref,
##D   T_ref          = met$T_ref
##D )
## End(Not run)




cleanEx()
nameEx("get_fused_dsm")
### * get_fused_dsm

flush(stderr()); flush(stdout())

### Name: get_fused_dsm
### Title: get_fused_dsm
### Aliases: get_fused_dsm

### ** Examples





cleanEx()
nameEx("get_metadata")
### * get_metadata

flush(stderr()); flush(stdout())

### Name: get_metadata
### Title: get_metadata
### Aliases: get_metadata

### ** Examples

meta <- gloBFPr::get_metadata(test=TRUE)




cleanEx()
nameEx("get_metrics")
### * get_metrics

flush(stderr()); flush(stdout())

### Name: get_metrics
### Title: get_metrics
### Aliases: get_metrics get_morphology get_neighbors get_bgvi get_dng

### ** Examples

library(gloBFPr)
data(globfp_example)
result <- gloBFPr::get_morphology(globfp_example[c(1:3),], quiet = TRUE)

result <- gloBFPr::get_neighbors(globfp_example[c(1:3),], radius = 100)

result <- gloBFPr::get_dng(#globfp_example[c(1:3),],
                           datasource = "metachm",
                           unit = "m2")

# Measure along real road and path centre lines instead
result <- gloBFPr::get_dng(#globfp_example[c(1:3),],
                           datasource = "metachm",
                           unit = "m2",
                           network = "osm")



cleanEx()
nameEx("globfp_example")
### * globfp_example

flush(stderr()); flush(stdout())

### Name: globfp_example
### Title: Test 3D-GloBFP dataset
### Aliases: globfp_example
### Keywords: datasets

### ** Examples

data(globfp_example)
head(globfp_example)



cleanEx()
nameEx("globfp_example_canopy_height")
### * globfp_example_canopy_height

flush(stderr()); flush(stdout())

### Name: globfp_example_canopy_height
### Title: Example canopy height raster for the 3D-GloBFP sample
### Aliases: globfp_example_canopy_height
### Keywords: datasets

### ** Examples

data(globfp_example_canopy_height)
canopy_height <- terra::rast(globfp_example_canopy_height)
canopy_height



cleanEx()
nameEx("globfp_example_dem")
### * globfp_example_dem

flush(stderr()); flush(stdout())

### Name: globfp_example_dem
### Title: Example DEM for the 3D-GloBFP sample
### Aliases: globfp_example_dem
### Keywords: datasets

### ** Examples

data(globfp_example_dem)
dem <- terra::rast(globfp_example_dem)
dem



cleanEx()
nameEx("plot_bgvi_viewshed")
### * plot_bgvi_viewshed

flush(stderr()); flush(stdout())

### Name: plot_bgvi_viewshed
### Title: Visualize an Individual Building BGVI Viewshed
### Aliases: plot_bgvi_viewshed

### ** Examples





cleanEx()
nameEx("search_3dglobdf")
### * search_3dglobdf

flush(stderr()); flush(stdout())

### Name: search_3dglobdf
### Title: search_3dglobdf
### Aliases: search_3dglobdf

### ** Examples

## Not run: 
##D buildings <- gloBFPr::search_3dglobdf(bbox=c(-84.485519,45.636118,-84.462774,45.650639))
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
