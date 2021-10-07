# library(muscat)

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Seurat))

here::i_am("rna/pseudobulk/pseudobulk_rna.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",    help='SingleCellExperiment object')
p$add_argument('--seurat',         type="character",    help='Seurat object')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--group_by',    type="character",    help='Metadata column to group cells by')
p$add_argument('--normalisation_method',    type="character",    help='Metadata column to group cells by')
p$add_argument('--outdir',      type="character",    help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- "/Users/argelagr/data/DUX4_hESCs_multiome/results/rna/clustering/sample_metadata_after_clustering.txt.gz"
# args$sce <- "/Users/argelagr/data/DUX4_hESCs_multiome/processed/rna/SingleCellExperiment.rds"
# args$group_by <- "cluster"
# args$normalisation_method <- "log"
# args$outdir <- "/Users/argelagr/data/DUX4_hESCs_multiome/results/rna/pseudobulk"
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%opts$samples]
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Load SingleCellExperiment
seurat <- load_Seurat(io$rna.seurat, cells=sample_metadata$cell)
seurat@meta.data <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell")

################
## Pseudobulk ##
################

# assays(sce)$cpm <- edgeR::cpm(assay(sce), normalized.lib.sizes = FALSE, log = FALSE)

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = args$group_by,
  fun = c("sum"),
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

###################
## Normalisation ##
###################

if (args$normalisation_method=="deseq2") {
  suppressPackageStartupMessages(library(DESeq2))
	dds <- DESeqDataSet(sce_pseudobulk, design=~1)
	dds <- varianceStabilizingTransformation(dds)
	logcounts(sce_pseudobulk) <- assay(dds)

} else if (args$normalisation_method=="log") {

	logcounts(sce_pseudobulk) <- log(counts(sce_pseudobulk)+1)
}

# Create Seurat object from the SingleCellExperiment
sce_pseudobulk@colData$sample <- rownames(sce_pseudobulk@colData)  # at least one metadata column is needed
seurat_pseudobulk <- as.Seurat(sce_pseudobulk)

seurat_pseudobulk <- RenameAssays(seurat_pseudobulk, originalexp="RNA")

##########
## Save ##
##########

saveRDS(sce_pseudobulk, file.path(args$outdir,sprintf("SingleCellExperiment_pseudobulk_%s.rds",args$group_by)))
saveRDS(seurat_pseudobulk, file.path(args$outdir,sprintf("Seurat_pseudobulk_%s.rds",args$group_by)))
