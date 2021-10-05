library(MOFA2)
library(SingleCellExperiment)
library(Seurat)

reticulate::use_python("/Users/argelagr/opt/anaconda3/envs/main/bin/python", required = T)
# setwd("/Users/argelagr/DUX4_hESCs_multiome")
here::i_am("rna/mofa/run_mofa.R")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outfile <- file.path(io$basedir,"results/rna/mofa/mofa.rds")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  # .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[pass_rnaQC==TRUE]
table(sample_metadata$sample)

###############
## Load data ##
###############

seurat <- load_Seurat(io$rna.seurat, normalise=T, cells=sample_metadata$cell)
dim(seurat)

#######################
## Feature selection ##
#######################

seurat <- FindVariableFeatures(seurat, selection.method = "vst")

#######################
# Create MOFA object ##
#######################

MOFAobject <- create_mofa(seurat)

####################
## Define options ##
####################

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42
# train_opts$maxiter <- 5

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(MOFAobject,
  model_options = model_opts,
  training_options = train_opts
)

#####################
## Train the model ##
#####################

model <- run_mofa(MOFAobject)


plot_weights(model, factor = 8, view="ATAC_chromVAR", nfeatures = 10, text_size = 3)
plot_weights(model, factor = 2, view="RNA", nfeatures = 15, text_size = 4)


plot_factor(model, factors = 2, color_by = "TRIM49", group_by = "sample", add_violin = T, add_boxplot = T, add_dots = F, legend=F) +
  # scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factors(model, factors = c(1,2), color_by = "sample") +
  # scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none"
  )
