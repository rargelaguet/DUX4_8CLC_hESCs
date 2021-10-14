# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/processing/2_archR_pseudobulk_celltypes.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outdir',     type="character",    help='Output directory')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--matrices_to_pseudobulk',     type="character",       nargs="+",   help='Matrices to pseudobulk')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
# args$group_by <- "cluster"
# args$matrices_to_pseudobulk <- c("PeakMatrix", "GeneScoreMatrix_distal", "GeneScoreMatrix_TSS")
# args$threads <- 1
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & sample%in%opts$samples]
stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads = args$threads)

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(all(rownames(sample_metadata.to.archr) == rownames(getCellColData(ArchRProject.filt))))
ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = sample_metadata.to.archr[[args$group_by]],
  name = args$group_by,
  cells = rownames(sample_metadata.to.archr),
  force = TRUE
)

# print cell numbers
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,args$group_by)[[1]])

###################################################
## Pseudobulk into a SummarizedExperiment object ##
###################################################

# if (is.null(args$matrices_to_pseudobulk)) {
# 	args$matrices_to_pseudobulk <- getAvailableMatrices(ArchRProject.filt)
# }

stopifnot(args$matrices_to_pseudobulk%in%getAvailableMatrices(ArchRProject.filt))

se_list <- list()
for (i in args$matrices_to_pseudobulk) {
  
  # summarise
  se_list[[i]] <- getGroupSE(ArchRProject.filt, groupBy = args$group_by, useMatrix = i, divideN = TRUE)
  
  # save
  outfile <- sprintf("%s/pseudobulk_%s_summarized_experiment.rds",args$outdir,i)
  saveRDS(se_list[[i]], outfile)
}
