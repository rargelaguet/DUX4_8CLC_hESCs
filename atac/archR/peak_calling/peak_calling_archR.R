# https://www.ArchRProject.com/bookdown/calling-peaks-with-archr.html
# Note: this requires the creation of pseudobulk replicates with 'addGroupCoverages'
# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html
# see /.../atac/archR/pseudobulk/1_archR_add_GroupCoverage.R

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/peak_calling/peak_calling_archR.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outdir',     type="character",    help='Output directory')
p$add_argument('--pathToMacs2',     type="character",    help='Path to MACS2 software')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--pvalue_cutoff',     type="double",   help='MACS2 p-value cutoff')
p$add_argument('--extend_summits',     type="integer",   help='Number of bp to extend peak summits')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/results/atac/archR/test/qc/sample_metadata_after_qc.txt.gz"
# args$pathToMacs2 <- "/bi/group/reik/ricard/software/miniconda3/envs/main/bin/macs2"
# args$group_by <- "cluster"
# args$pvalue_cutoff <- 1e-3
# args$extend_summits <- 300
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


##################
## Peak calling ##
##################

ArchRProject.filt <- addReproduciblePeakSet(
  ArchRProj = ArchRProject.filt, 
  groupBy = args$group_by, 
  peakMethod = "Macs2",
  excludeChr = c("chrM", "chrY"),
  pathToMacs2 = args$pathToMacs2,
  cutOff = args$pvalue_cutoff,
  extendSummits = args$extend_summits,
  plot = FALSE,
  force = TRUE
)

################
## Save peaks ##
################

# Save PeakSet
# NOTE THAT THIS HAS TO GO TO ARCHR DIRECTORY
saveRDS(ArchRProject.filt@peakSet, file.path(args$outdir,"/PeakSet.rds"))

# fetch peaks in data.table format
dt <- getPeakSet(ArchRProject.filt) %>% as.data.table() %>% setnames("seqnames","chr")

# Save peak metadata
fwrite(dt, file.path(args$outdir,"peak_metadata.tsv.gz"), sep="\t")

# save peaks in bed format
fwrite(dt[,c("chr","start","end")], file.path(args$outdir,"peaks_archR_macs2.bed.gz"), sep="\t", col.names = F)

#####################
## Add peak matrix ##
#####################

ArchRProject@peakSet <- ArchRProject.filt@peakSet
ArchRProject <- addPeakMatrix(ArchRProject, binarize = FALSE)
saveArchRProject(ArchRProject)