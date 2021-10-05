# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

here::i_am("atac/archR/bigwig/archR_export_bw.R")

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

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(all(rownames(sample_metadata.to.archr) == rownames(getCellColData(ArchRProject.filt))))
ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = sample_metadata.to.archr[[opts$group.by]],
  name = opts$group.by,
  cells = rownames(sample_metadata.to.archr),
  force = TRUE
)

stopifnot(opts$group.by %in% names(ArchRProject.filt@projectMetadata$GroupCoverages))

# print cell numbers
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,opts$group.by)[[1]])


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
