# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/pseudobulk/1_archR_add_GroupCoverage.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
# p$add_argument('--outdir',     type="character",    help='Output directory')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--min_cells',     type="integer",    default=25,   help='Minimum number of cells')
p$add_argument('--max_cells',     type="integer",    default=1000,   help='Maximum number of cells')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
# args$group_by <- "cluster"
# args$min_cells <- 25
# args$max_cells <- 100
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

#########################
## Add Group Coverages ##
#########################

# Check if group Coverages already exist
# ArchRProject.filt@projectMetadata$GroupCoverages

# This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates 
# and then merge these replicates into a single insertion coverage file.
# Output: creates files in archR/GroupCoverages/celltype: [X]._.Rep[Y].insertions.coverage.h5
ArchRProject.filt <- addGroupCoverages(ArchRProject.filt, 
  groupBy = args$group_by, 
  minCells = args$min_cells,
  maxCells = args$max_cells,
  force = TRUE
)

##########
## Save ##
##########

io$archR.projectMetadata <- paste0(io$archR.directory,"/projectMetadata.rds")
saveRDS(ArchRProject.filt@projectMetadata, io$archR.projectMetadata)

# saveArchRProject(ArchRProject.filt)