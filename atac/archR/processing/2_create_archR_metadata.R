

#####################
## Define settings ##
#####################

here::i_am("atac/archR/processing/2_create_archR_metadata.R")
source(here::here("settings.R"))

# I/O
io$metadata <- paste0(io$basedir,"/results/rna/doublet_detection/sample_metadata_after_doublets.txt.gz")
io$metadata.out <- paste0(io$archR.directory,"/sample_metadata.txt.gz")

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[,c("cell", "sample", "barcode", "nFeature_RNA", "nCount_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA", "pass_rnaQC", "doublet_score", "doublet_call")]
       # "TSSEnrichment_atac", "ReadsInTSS_atac", "PromoterRatio_atac", "NucleosomeRatio_atac", "nFrags_atac", "BlacklistRatio_atac", "pass_atacQC")

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

######################
## Load ArchR stats ##
######################
  
# fetch archR's metadata
# note that QC is done later in the QC/qc.R script
archR_metadata <- getCellColData(ArchRProject) %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","cell") %>%
  .[,c("cell", "TSSEnrichment", "ReadsInTSS", "PromoterRatio", "NucleosomeRatio", "nFrags",  "BlacklistRatio")]# %>%
  # setnames("Sample","sample") %>%
  # .[,cell:=stringr::str_replace_all(cell,"#","_")] %>%
  # .[,sample:=strsplit(cell,"#") %>% map_chr(1)] %>%
  # .[,barcode:=strsplit(cell,"#") %>% map_chr(2)]

cols.to.rename <- c("TSSEnrichment","ReadsInTSS","PromoterRatio","NucleosomeRatio","nFrags","BlacklistRatio")
idx.cols.to.rename <- which(colnames(archR_metadata)%in%cols.to.rename)
colnames(archR_metadata)[idx.cols.to.rename] <- paste0(colnames(archR_metadata)[idx.cols.to.rename], "_atac")

###########
## Merge ##
###########

sample_metadata_all <- sample_metadata %>% 
  merge(archR_metadata,by="cell", all=TRUE) 

# Fill missing entries
sample_metadata_all %>%
  .[is.na(sample),sample:=strsplit(cell,"#") %>% map_chr(1)] %>%
  .[is.na(barcode),barcode:=strsplit(cell,"#") %>% map_chr(2)] %>%
  .[is.na(stage),stage:=strsplit(sample,"_") %>% map_chr(1)]

# sanity checks
table(sample_metadata$stage)
table(sample_metadata$sample)
stopifnot(all(!is.na(sample_metadata_all$sample)))
stopifnot(all(!is.na(sample_metadata_all$barcode)))
stopifnot(all(!is.na(sample_metadata_all$stage)))

#############################
## Update ArchR's metadata ##
#############################

metadata.to.archR <- sample_metadata_all %>% 
  .[cell%in%rownames(ArchRProject)] %>% setkey(cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(all(metadata.to.archR$TSSEnrichment_atac == getCellColData(ArchRProject, "TSSEnrichment")[[1]]))

for (i in colnames(metadata.to.archR)) {
  ArchRProject <- addCellColData(
    ArchRProject,
    data = metadata.to.archR[[i]], 
    name = i,
    cells = rownames(metadata.to.archR),
    force = TRUE
  )
}

head(getCellColData(ArchRProject))

##########
## Save ##
##########

fwrite(sample_metadata_all, io$metadata.out, sep="\t", na="NA", quote=F)

saveArchRProject(ArchRProject)
