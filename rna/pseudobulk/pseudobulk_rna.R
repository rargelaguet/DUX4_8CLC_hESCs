library(muscat)
library(DESeq2)

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outfile <- file.path(io$basedir,"results/rna/pseudobulk/SingleCellExperiment.rds")

# Options
opts$group.by <- "eight_cell_like_ricard"

###############
## Load data ##
###############

# Load cell metadata
# io$metadata <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/doublets/sample_metadata_after_doublets.txt.gz"
sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%opts$samples]

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

###################################
## Aggregate counts per celltype ##
###################################

# assays(sce)$cpm <- edgeR::cpm(assay(sce), normalized.lib.sizes = FALSE, log = FALSE)

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = opts$group.by,
  fun = c("sum"),
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

###########################
## Normalise with DESeq2 ##
###########################

# create DESeq object
dds <- DESeqDataSet(sce_pseudobulk, design=~1)

# This function calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) 
# and then transforms the count data (normalized by division by the size factors or normalization factors), 
# yielding a matrix of values which are now approximately homoskedastic 
dds <- varianceStabilizingTransformation(dds)

logcounts(sce_pseudobulk) <- assay(dds)

################################################
## Normalise with a simple log transformation ##
################################################

logcounts(sce_pseudobulk) <- log(counts(sce_pseudobulk)+1)

###################
## Sanity checks ##
###################

# cor(
#   colMeans(logcounts(sce_pseudobulk)),
#   metadata(sce_pseudobulk)$n_cells
# )

# hist(logcounts(sce_pseudobulk))

# logcounts(sce_pseudobulk)[1:10,1:2]
# counts(sce_pseudobulk)[1:10,1:2]

##########
## Save ##
##########

saveRDS(sce_pseudobulk, io$outfile)
