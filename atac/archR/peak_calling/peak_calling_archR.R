# https://www.ArchRProject.com/bookdown/calling-peaks-with-archr.html
# Note: this requires the creation of pseudobulk replicates with 'addGroupCoverages'
# see /.../atac/archR/pseudobulk/archR_pseudobulk_celltypes.R

########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
  pathToMacs2 <- "/Users/ricard/anaconda3/envs/base_new/bin/macs2"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
  pathToMacs2 <- "/nfs/research1/stegle/users/ricard/conda-envs/R4/bin/macs2"
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$metadata <- paste0(io$basedir,"/results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/peak_calling")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$pvalue.cutoff <- 0.01
opts$group.by <- "celltype.mapped"

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & !is.na(celltype.mapped)] %>%
  .[sample%in%opts$samples]

stopifnot(sample_metadata$cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])

###################
## Sanity checks ##
###################

stopifnot("celltype.mapped" %in% names(ArchRProject.filt@projectMetadata$GroupCoverages))

##################
## Peak calling ##
##################

# This function will get insertions from coverage files, call peaks, and merge peaks to get a "Union Reproducible Peak Set".
ArchRProject.filt <- addReproduciblePeakSet(
  ArchRProj = ArchRProject.filt, 
  groupBy = opts$group.by, 
  peakMethod = "Macs2",
  excludeChr = c("chrM", "chrY"),
  pathToMacs2 = pathToMacs2,
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
dt <- getPeakSet(ArchRProject.filt) %>% as.data.table() %>% setnames(c("seqnames"),c("chr"))# %>%
  # .[,c("strand","idx","nearestTSS","distToTSS","GroupReplicate","replicateScoreQuantile","groupScoreQuantile","Reproducibility"):=NULL]

# Save peak metadata
outfile <- paste0(io$archR.directory,"/PeakCalls/peak_metadata.tsv.gz")
fwrite(dt, outfile, sep="\t")

# save peaks in bed format
to.save <- dt[,c("chr","start","end")]
fwrite(to.save, paste0(io$archR.directory,"/PeakCalls/peaks_archR_macs2.bed.gz"), sep="\t", col.names = F)


#####################
## Add peak matrix ##
#####################

ArchRProject <- addPeakMatrix(ArchRProject, binarize = FALSE)
saveArchRProject(ArchRProject)

stop()

# ArchRProject.filt <- addPeakMatrix(ArchRProject.filt, binarize = FALSE)
# foo <- getMatrixFromProject(ArchRProject.filt, useMatrix = "PeakMatrix")
# bar <- getMatrixFromProject(ArchRProject, useMatrix = "PeakMatrix")
# getAvailableMatrices(ArchRProject.filt)

# Save matrix.mtx.gz
# mtx <- getMatrixFromProject(ArchRProject.filt, useMatrix = "PeakMatrix")@assays@data[[1]]
# outfile <- paste0(io$archR.directory,"/PeakMatrix/matrix.mtx.gz")
# Matrix::writeMM(mtx, file=outfile)

# Save barcodes.tsv.gz
# outfile <- paste0(io$archR.directory,"/PeakMatrix/barcodes.tsv.gz")
# fwrite(as.data.table(gsub("#","_", colnames(mtx))), outfile, col.names=F)

# save features.tsv.gz
# foo <- dt[,c("chr","start","end")] %>% .[,foo:=sprintf("%s_%s_%s",chr,start,end)] %>% .[,"foo"]
# outfile <- paste0(io$archR.directory,"/PeakMatrix/features.tsv.gz")
# fwrite(foo, outfile, sep="\t", col.names = F)

###################################################
## Calculate fraction of reads in peaks per cell ##
###################################################

peak.mtx <- getMatrixFromProject(ArchRProject.filt, useMatrix = "PeakMatrix", binarize = T)@assays@data[[1]]

peak_qc.dt <- data.table(
  cell = colnames(peak.mtx),
  fragments_in_peaks = Matrix::colSums(peak.mtx)
)

to.plot <- peak_qc.dt %>% 
  merge(sample_metadata) %>%
  .[,fraction_reads_in_tss:=ReadsInTSS_atac/nFrags_atac] %>%
  .[,fraction_fragments_in_peaks:=fragments_in_peaks/nFrags_atac]
  
ggscatter(to.plot[pass_atacQC==TRUE], x="fraction_fragments_in_peaks",y="nFrags_atac", color="pass_rnaQC", size=0.7) +
  yscale("log10", .format = TRUE) +
# ggscatter(to.plot[sample=="E8.5_rep2"], x="TSSEnrichment_atac",y="nFrags_atac", color="pass_rnaQC", size=0.7) +
  facet_wrap(~sample)

asd = to.plot[sample=="E8.5_rep2" & pass_atacQC==TRUE] %>%
  .[,c("fragments_in_peaks", "nFeature_RNA", "nCount_RNA", "pass_rnaQC", "TSSEnrichment_atac", "ReadsInTSS_atac", "nFrags_atac", "pass_atacQC", "fraction_reads_in_tss", "fraction_fragments_in_peaks")]
