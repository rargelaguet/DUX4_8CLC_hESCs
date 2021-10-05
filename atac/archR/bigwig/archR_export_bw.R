# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

here::i_am("/Users/argelagr/DUX4_hESCs_multiome/atac/archR/bigwig/archR_export_bw.R")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
# io$metadata <- file.path(io$basedir,"/results/rna/dimensionality_reduction/sample_metadata_after_clustering.txt.gz")

# Options
opts$group.by <- "eight_cell_like_ricard"
opts$ncores <- 2

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & !is.na(eight_cell_like_ricard)] %>%
  .[sample%in%opts$samples]

stopifnot(sample_metadata$cell %in% rownames(ArchRProject))

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads = opts$ncores)

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])

###################
## Export bigwig ##
###################

# This function will group, summarize and export a bigwig for each group in an ArchRProject.
getGroupBW(
  ArchRProj = ArchRProject.filt,
  groupBy = opts$group.by,
  # groupBy = "Sample",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 100,
  ceiling = 4
)
