## Load packages {-}
library(SpatialExperiment)
library(STexampleData)
library(ggspavis)
library(ggplot2)
library(scater)
library(scran)
library(igraph)
library(pheatmap)
library(ggExtra)

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
head(rowData(spe))
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

#Spot-level Quality Control
## Plot spatial coordinates without annotations
plotSpots(spe)

#Calculating QC metrics
## Dataset dimensions before the filtering
dim(spe)
## Subset to keep only on-tissue spots
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)
## Classify genes as "mitochondrial" (is_mito == TRUE) 
## or not (is_mito == FALSE)
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
rowData(spe)$gene_name[is_mito==TRUE]
## Calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe))

colnames(colData(spe))

#Number of cells per spot
## Density and histogram of the number of cells in each spot
ggplot(data = as.data.frame(colData(spe)),
       aes(x = cell_count)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 1,
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.5,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Number of cells per spot") + 
  ylab("Density") + 
  theme_classic()
library_sizes <- Matrix::colSums(assay(spe, "raw_count"))
colData(spe)$library_size <- library_sizes


#Library size threshold plot
## Density and histogram of library sizes
ggplot(data = as.data.frame(colData(spe)),
       aes(x = library_size)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Library size") + 
  ylab("Density") + 
  theme_classic()
colnames(colData(spe))
## Scatter plot, library size against number of cells per spot
plotSpotQC(spe, plot_type = "scatter", x_metric ="n_genes_by_counts", y_metric = "library_size", y_threshold = 700)

help("plotSpotQC")

head(assay(spe))
p = ggplot(as.data.frame(colData(spe)), aes(x=n_genes_by_counts, y=library_size)) +
  geom_point(size=0.5) + 
  geom_smooth(se=FALSE) +
  geom_hline(yintercept = 700, colour='red') + 
  theme_minimal()
ggMarginal(p, type='histogram', margins = 'both')
p

## Select library size threshold
qc_lib_size <- colData(spe)$library_size < 5000
## Check how many spots are filtered out
table(qc_lib_size)
## Add threshold in colData
colData(spe)$qc_lib_size <- qc_lib_size

## Check putative spatial patterns of removed spots
plotSpotQC(spe, plot_type = "spot", 
       discard = "qc_lib_size")

#Number of expressed genes
## Density and histogram of expressed genes
ggplot(data = as.data.frame(colData(spe)),
       aes(x = n_genes_by_counts)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Genes expressed in each spot") + 
  ylab("Density") + 
  theme_classic()

## Select expressed genes threshold
qc_detected <- colData(spe)$n_genes_by_counts < 3000
## Check how many spots are filtered out
table(qc_detected)

#Percentage of mitochondrial expression
## Density and histogram of percentage of mitochondrial expression
ggplot(data = as.data.frame(colData(spe)),
       aes(x = subsets_mito_percent)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Percentage of mitochondrial expression") + 
  ylab("Density") + 
  theme_classic()

## Density and histogram of the number of cells in each spot
ggplot(data = as.data.frame(colData(spe)),
       aes(x = n_genes_by_counts)) +
  
  geom_density(alpha = 0.5,
               adjust = 1.5,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Number of cells per spot") + 
  ylab("Density") + 
  theme_classic()

#Normalisation of counts
## Calculate library size factors
assays(spe)[["counts"]] <- assays(spe)[["raw_count"]]
spe <- computeLibraryFactors(spe)

## Have a look at the size factors
summary(sizeFactors(spe))

## Density and histogram of library sizes
ggplot(data = data.frame(sFact = sizeFactors(spe)), 
       aes(x = sFact)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Library size") + 
  ylab("Density") + 
  theme_classic()

## Calculate logcounts and store in the spe object
spe <- logNormCounts(spe)

## Check that a new assay has been added
assayNames(spe)

#Selecting genes
#Highly Variable Genes
## Fit mean-variance relationship
dec <- modelGeneVar(spe)
## Visualize mean-variance relationship
fit <- metadata(dec)
fit_df <- data.frame(mean = fit$mean,
                     var = fit$var,
                     trend = fit$trend(fit$mean))

ggplot(data = fit_df, 
       aes(x = mean, y = var)) + 
  geom_point() + 
  geom_line(aes(y = trend), colour = "dodgerblue", linewidth = 1.5) + 
  labs(x = "mean of log-expression",
       y = "variance of log-expression") + 
  theme_classic()

## Select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)

## How many HVGs?
length(top_hvgs)

#Spatially variable genes (SVGs)
#Dimensionality reduction
## Set seed
set.seed(987)
## Compute PCA
spe <- runPCA(spe, subset_row = top_hvgs)
## Check correctness - names
reducedDimNames(spe)
## Check correctness - dimensions
dim(reducedDim(spe, "PCA"))

#UMAP: Uniform Manifold Approximation and Projection
## Set seed
set.seed(987)
## Compute UMAP on top 50 PCs
spe <- runUMAP(spe, dimred = "PCA")
## Check correctness - names
reducedDimNames(spe)
## Update column names for easier plotting
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
