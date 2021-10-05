# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

here::i_am("atac/archR/pseudobulk/1_archR_add_GroupCoverage.R")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$metadata <- file.path(io$basedir,"/results/rna/dimensionality_reduction/sample_metadata_after_clustering.txt.gz")

# Options
opts$group.by <- "eight_cell_like_ricard"
opts$ncores <- 2

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads = opts$ncores)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & !is.na(eight_cell_like_ricard) & sample%in%opts$samples]

##################
## Subset ArchR ##
##################

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
  data = sample_metadata.to.archr[[opts$group.by]],
  name = opts$group.by,
  cells = rownames(sample_metadata.to.archr),
  force = TRUE
)

# print cell numbers
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,opts$group.by)[[1]])

#########################
## Add Group Coverages ##
#########################

# Check if group Coverages already exist
ArchRProject.filt@projectMetadata$GroupCoverages

# This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates 
# and then merge these replicates into a single insertion coverage file.
# Output: creates files in archR/GroupCoverages/celltype: [X]._.Rep[Y].insertions.coverage.h5
ArchRProject.filt <- addGroupCoverages(ArchRProject.filt, 
  groupBy = opts$group.by, 
  force = TRUE, 
  minCells = 50,
  maxCells = 200,
)

##########
## Save ##
##########

io$archR.projectMetadata <- paste0(io$archR.directory,"/projectMetadata.rds")
saveRDS(ArchRProject.filt@projectMetadata, io$archR.projectMetadata)

# saveArchRProject(ArchRProject.filt)