library(zellkonverter)
library(SpatialExperiment)
library(jsonlite)
library(png)
library(ggplot2)
library(ggspavis)

# Load h5ad as SingleCellExperiment
adata <- readH5AD("Human_Prostate_Cancer/Human_Prostate_Cancer_Acinar_Cell_Carcinoma_FFPE_09102021_10x_Visium.h5ad",reader = "R")

# Convert to SpatialExperiment
spe <- as(adata, "SpatialExperiment")

# Read spatial metadata
spatial_dir <- "Human_Prostate_Cancer/spatial"

# Load coordinates
coords <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), header = FALSE)
coords
colnames(coords) <- c("barcode", "tissue", "row", "col", "pxl_col_in_fullres", "pxl_row_in_fullres")
coords
spatialCoords(spe) <- as.matrix(coords[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])

# Add scalefactors and images
scale_factors <- fromJSON(file.path(spatial_dir, "scalefactors_json.json"))
scale_factors
hires_img <- readPNG(file.path(spatial_dir, "tissue_hires_image.png"))
hires_img
lowres_img <- readPNG(file.path(spatial_dir, "tissue_lowres_image.png"))
lowres_img
imgData(spe)
assay(spe)
colnames(colData(spe))
head(colData(spe))

image_path <- "Human_Prostate_Cancer/spatial/tissue_hires_image.png"
spe <- addImg(
  spe,
  sample_id = "sample01",     # must match sample_id in colData(spe)
  image_id = "H&E_hig_res",   # unique identifier for this image
  imageSource = image_path,  # path or URL to the image
  scaleFactor = 0.0728014,        # scaling factor relative to full resolution
  load = TRUE                # TRUE loads the image into memory
)

image_path <- "Human_Prostate_Cancer/spatial/tissue_lowres_image.png"
spe <- addImg(
  spe,
  sample_id = "sample01",     # must match sample_id in colData(spe)
  image_id = "H&E_lowres",   # unique identifier for this image
  imageSource = image_path,  # path or URL to the image
  scaleFactor = 0.02184042,        # scaling factor relative to full resolution
  load = TRUE                # TRUE loads the image into memory
)

#removal for imgdata
#spe <- rmvImg(spe, sample_id = "sample01", image_id = "H&E_hig_res")

# View summary
spe
imgData(spe)

## Check number of features/genes (rows) and spots (columns)
dim(spe)

## Check names of 'assay' tables
assayNames(spe)
assay(spe)

#Counts table and gene metadata
## Have a look at the counts table
assay(spe)[1:6,1:4]
assay(spe)[20:40, 2000:2010]
assay(spe)[33488:33508, 2000:2010]
## Have a look at the genes metadata
head(rowData(spe))
##Coordinates table and spot metadata
## Check the spatial coordinates
head(spatialCoords(spe))
## spot-level metadata
head(colData(spe))
#Image metadata
## Have a look at the image metadata
imgData(spe)

## retrieve the image
spi <- getImg(spe)
dev.off()  # closes the current graphics device
X11()
library(grid)
grid.newpage()
grid.raster(imgRaster(spi))

library(magrittr) 
## Extract the spot locations
spot_coords <- spatialCoords(spe) %>% as.data.frame
## Scale by low-res factor
lowres_scale <- imgData(spe)[imgData(spe)$image_id == 'H&E_lowres', 'scaleFactor']
spot_coords$x_axis <- spot_coords$pxl_col_in_fullres * lowres_scale
spot_coords$y_axis <- spot_coords$pxl_row_in_fullres * lowres_scale
## lowres image is 600x600 pixels
dim(imgRaster(spi))
## flip the Y axis
grid.raster(imgRaster(spi))
points(x=spot_coords$x_axis, y=spot_coords$y_axis)

#using ggplot2
ggplot(mapping = aes(1:600, 1:600)) +
  annotation_raster(imgRaster(spi), xmin = 1, xmax = 600, ymin = 1, ymax = 600) +
  geom_point(data=spot_coords, aes(x=x_axis, y=y_axis), alpha=0.2) + xlim(1, 600) + ylim(1, 600) +
  coord_fixed() + 
  theme_void()

## Add the annotation to the coordinate data frame
spot_coords$on_tissue <- as.logical(colData(spe)$in_tissue)

ggplot(mapping = aes(1:600, 1:600)) +
  annotation_raster(imgRaster(spi), xmin = 1, xmax = 600, ymin = 1, ymax = 600) +
  geom_point(data=spot_coords, aes(x=x_axis, y=y_axis, colour=on_tissue), alpha=0.2) + xlim(1, 600) + ylim(1, 600) +
  coord_fixed() + 
  theme_void()
plotSpots(spe, in_tissue = NULL, annotate='in_tissue', size=0.5)
