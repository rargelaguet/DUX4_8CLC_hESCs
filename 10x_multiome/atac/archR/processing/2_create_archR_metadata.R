
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/processing/2_create_archR_metadata.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outfile',     type="character",    help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/results/rna/clustering/sample_metadata_after_clustering.txt.gz"
# args$outfile <- "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/processed/atac/archR/test/sample_metadata_after_archR.txt.gz"
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))

###################
## Load metadata ##
###################

# Note that this metadata file can be derived from the RNA pipeline
sample_metadata <- fread(args$metadata)# %>%
  # .[,c("cell", "sample", "barcode", "nFeature_RNA", "nCount_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA", "pass_rnaQC", "doublet_score", "doublet_call")]

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

######################
## Load ArchR stats ##
######################
  
# fetch archR's metadata
archR_metadata <- getCellColData(ArchRProject) %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","cell") %>%
  .[,c("cell", "TSSEnrichment", "ReadsInTSS", "PromoterRatio", "NucleosomeRatio", "nFrags",  "BlacklistRatio")]

cols.to.rename <- c("TSSEnrichment","ReadsInTSS","PromoterRatio","NucleosomeRatio","nFrags","BlacklistRatio")
idx.cols.to.rename <- which(colnames(archR_metadata)%in%cols.to.rename)
colnames(archR_metadata)[idx.cols.to.rename] <- paste0(colnames(archR_metadata)[idx.cols.to.rename], "_atac")

###########
## Merge ##
###########

sample_metadata_tosave <- sample_metadata %>% 
  merge(archR_metadata,by="cell", all=TRUE) 

# Fill missing entries
sample_metadata_tosave %>%
  .[is.na(sample),sample:=strsplit(cell,"#") %>% map_chr(1)] %>%
  .[is.na(barcode),barcode:=strsplit(cell,"#") %>% map_chr(2)]

# round
sample_metadata_tosave[,c("TSSEnrichment_atac","NucleosomeRatio_atac","PromoterRatio_atac","BlacklistRatio_atac"):=list(round(TSSEnrichment_atac,2),round(NucleosomeRatio_atac,2),round(PromoterRatio_atac,2),round(BlacklistRatio_atac,2))]
sample_metadata_tosave[,c("ribosomal_percent_RNA","mitochondrial_percent_RNA"):=list(round(ribosomal_percent_RNA,2),round(mitochondrial_percent_RNA,2))]

# sanity checks
table(sample_metadata$sample)
stopifnot(all(!is.na(sample_metadata_tosave$sample)))
stopifnot(all(!is.na(sample_metadata_tosave$barcode)))


#############################
## Update ArchR's metadata ##
#############################

metadata.to.archR <- sample_metadata_tosave %>% 
  .[cell%in%rownames(ArchRProject)] %>% setkey(cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

# stopifnot(all(metadata.to.archR$TSSEnrichment_atac == getCellColData(ArchRProject,"TSSEnrichment")[[1]]))

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

fwrite(sample_metadata_tosave, args$outfile, sep="\t", na="NA", quote=F)

saveArchRProject(ArchRProject)
