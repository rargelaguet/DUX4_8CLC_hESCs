# https://www.ArchRProject.com/bookdown/calling-peaks-with-archr.html
# Note: this requires the creation of pseudobulk replicates with 'addGroupCoverages'
# see /.../atac/archR/pseudobulk/archR_pseudobulk_celltypes.R

here::i_am("atac/archR/peak_calling/peak_calling_archR.R")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$metadata <- paste0(io$basedir,"/results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/peak_calling")

# Options
opts$pvalue.cutoff <- 0.01
opts$group.by <- "eight_cell_like_ricard"

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & !is.na(eight_cell_like_ricard)] %>%
  .[sample%in%opts$samples]

stopifnot(sample_metadata$cell %in% rownames(ArchRProject))

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
## Sanity checks ##
###################


##################
## Peak calling ##
##################

# This function will get insertions from coverage files, call peaks, and merge peaks to get a "Union Reproducible Peak Set".

ArchRProject.filt <- addReproduciblePeakSet(
  ArchRProj = ArchRProject.filt, 
  groupBy = opts$group.by, 
  peakMethod = "Macs2",
  excludeChr = c("chrM", "chrY"),
  pathToMacs2 = io$pathToMacs2,
  cutOff = opts$pvalue.cutoff,
  extendSummits = 300,
  plot = FALSE,
  force = TRUE
)

##################
## Filter Peaks ##
##################

# NOTE: THIS DOESN'T UPDATE THE CONTENT THAT HAS BEEN STORED IN TEH ARROWFILE SUCH AS "sumy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/rowSums"))"
# opts$min.score <- 25
# ArchRProject.filt@peakSet <- ArchRProject.filt@peakSet[ArchRProject.filt@peakSet$score>=opts$min.score]

################
## Save peaks ##
################

# Save PeakSet
saveRDS(ArchRProject.filt@peakSet, paste0(io$archR.directory,"/PeakSet.rds"))

# fetch peaks in data.table format
dt <- getPeakSet(ArchRProject.filt) %>% as.data.table() %>% setnames("seqnames","chr")

# Save peak metadata
outfile <- file.path(io$archR.directory,"PeakCalls/peak_metadata.tsv.gz")
fwrite(dt, outfile, sep="\t")

# save peaks in bed format
to.save <- dt[,c("chr","start","end")]
fwrite(to.save, file.path(io$archR.directory,"PeakCalls/peaks_archR_macs2.bed.gz"), sep="\t", col.names = F)

#####################
## Add peak matrix ##
#####################

ArchRProject <- addPeakMatrix(ArchRProject, binarize = FALSE)
saveArchRProject(ArchRProject)