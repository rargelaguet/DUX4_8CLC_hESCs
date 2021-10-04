library(SingleCellExperiment)
library(scran)
library(scater)
library(shiny)
library(iSEE)

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# Define I/O

# Define options
opts$samples <- c(
  "HNES1_DUX4_overexpression_L001",
  "HNES1_wildtype_L001"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%opts$samples]
table(sample_metadata$sample)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.10]

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]
dim(sce_filt)

##############################
## Dimensionality reduction ##
##############################

data <- scale(t(logcounts(sce_filt)), center = T, scale = F)

# PCA
npcs <- 15
reducedDim(sce_filt, "PCA") <- irlba::prcomp_irlba(data, n=npcs)$x#[,1:npcs]

# Filter PCA solution
# reducedDim(sce_filt, "PCA") <- reducedDim(sce_filt, "PCA")[,-3]

# UMAP
set.seed(42)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 30, min_dist = 0.3)

#######################
## Define color maps ##
#######################

# sce$celltype <- sce$celltype.mapped_mnn
# 
# celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(sce_filt$celltype)]
# all(unique(sce_filt$celltype) %in% names(celltype.colors))
# all(names(celltype.colors)%in%unique(sce_filt$celltype))
# sce_filt$celltype <- factor(sce_filt$celltype, levels=names(celltype.colors))
# 
# celltype_color_fun <- function(n){
#   return(celltype.colors)
# }

categorical_color_fun <- function(n){
  return(RColorBrewer::brewer.pal(n, "Set2"))
}

# Define color maps
ecm <- ExperimentColorMap(
  # List of colormaps for assays.
  # assays = list(
  #   counts = viridis::viridis,
  #   cufflinks_fpkm = fpkm_color_fun
  # ),
  # colData = list(
  #   celltype = celltype_color_fun
  # ),
  # Colormaps applied to all undefined continuous assays
  all_continuous = list(
    assays = viridis::viridis
  ),
  # Colormaps applied to all undefined categorical assays
  all_discrete = list(
    assays = categorical_color_fun
  )
  # Colormap applied to all undefined categorical covariates.
  # global_discrete <- list()
  # Colormap applied to all undefined continuous covariates.
  # global_continuous <- list()
)

#########################
## Define iSEE options ##
#########################

# sce_filt <- registerAppOptions(sce_filt, color.maxlevels=40)
# getAppOption("color.maxlevels", sce_filt)

##############
## Run iSEE ##
##############

app <- iSEE(sce_filt, colormap = ecm)
runApp(app)
