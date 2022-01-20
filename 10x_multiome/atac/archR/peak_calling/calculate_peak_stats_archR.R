########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
io$outfile <- paste0(io$basedir, "/results/atac/archR/peak_calling/peak_stats_binarised.txt.gz")

# Options
opts$remove.small.celltypes <- TRUE

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples]

# subset celltypes with sufficient number of cells
if (opts$remove.small.celltypes) {
  opts$min.cells <- 25
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("celltype.mapped")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}
opts$celltypes <- unique(sample_metadata$celltype.mapped) %>% .[!is.na(.)]

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]

#######################
## Fetch Peak Matrix ##
#######################

atac.peakMatrix.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix", binarize = TRUE)
dim(atac.peakMatrix.se)

# Define peak names
# row.ranges.dt <- rowData(atac.peakMatrix.pseudobulk.se) %>% as.data.table %>% .[,idx:=sprintf("%s_%s_%s",seqnames,start,end)]
row.ranges.dt <- rowRanges(atac.peakMatrix.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,c("chr","start","end")] %>%
  .[,idx:=sprintf("%s_%s_%s",chr,start,end)]
rownames(atac.peakMatrix.se) <- row.ranges.dt$idx

###############################
## Load pseudobulk ATAC data ##
###############################

atac.peakMatrix.pseudobulk.se <- readRDS(io$archR.pseudobulk.peakMatrix.se)[,opts$celltypes]

opts$celltypes[!opts$celltypes %in% colnames(atac.peakMatrix.pseudobulk.se)]

# Define peak names
row.ranges.dt <- rowData(atac.peakMatrix.pseudobulk.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,c("chr","start","end")] %>%
  .[,idx:=sprintf("%s_%s_%s",chr,start,end)]
rownames(atac.peakMatrix.pseudobulk.se) <- row.ranges.dt$idx

##########################
## Calculate peak stats ##
##########################

peakStats.pseudobulk.dt <- data.table(
  peak = rownames(atac.peakMatrix.pseudobulk.se),
  var = assay(atac.peakMatrix.pseudobulk.se,"PeakMatrix") %>% sparseMatrixStats::rowVars(.) %>% round(5),
  mean = assay(atac.peakMatrix.pseudobulk.se,"PeakMatrix") %>% Matrix::rowMeans(.) %>% round(5)
)

peakStats.singlecell.dt <- data.table(
  peak = rownames(atac.peakMatrix.se),
  var = assay(atac.peakMatrix.se,"PeakMatrix") %>% sparseMatrixStats::rowVars(.) %>% round(5),
  mean = assay(atac.peakMatrix.se,"PeakMatrix") %>% Matrix::rowMeans(.) %>% round(5)
)

peakStats.dt <- merge(peakStats.singlecell.dt, peakStats.pseudobulk.dt, by="peak", suffixes=c("_singlecell","_pseudobulk"))

##########
## Save ##
##########

fwrite(peakStats.dt, io$outfile, sep="\t")
