# suppressMessages(library(SingleCellExperiment))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(argparse))
# suppressMessages(library(Seurat))


#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/DUX4_hESCs_multiome"
  io$gene_metadata <- "/Users/ricard/data/ensembl/human/v93/bioMart/all_genes/Hsapiens_genes_BioMart.93.txt.gz"
  io$pathToMacs2 <- "/Users/ricard/anaconda3/envs/base_new/bin/macs2"
} else if (grepl("BI2404M",Sys.info()['nodename'])) {
    io$basedir <- "/Users/argelagr/data/DUX4_hESCs_multiome"
    io$gene_metadata <- "/Users/argelagr/data/ensembl/human/v93/bioMart/all_genes/Hsapiens_genes_BioMart.93.txt.gz"
    io$pathToMacs2 <- "/Users/argelagr/opt/anaconda3/envs/main/bin/macs2"
} else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
  if (grepl("DUX4_hESCs_multiome", Sys.info()['effective_user'])) {
    stop()
  } else if (grepl("argelag", Sys.info()['effective_user'])) {
    io$basedir <-'/bi/group/reik/ricard/data/DUX4_hESCs_multiome'
    io$gene_metadata <- "/bi/group/reik/ricard/data/ensembl/human/v93/bioMart/all_genes/Hsapiens_genes_BioMart.93.txt.gz"
    io$pathToMacs2 <- "/bi/group/reik/ricard/software/miniconda3/envs/main/bin/macs2"
  }
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")

# RNA
io$rna.seurat <- paste0(io$basedir,"/processed/rna/seurat.rds")
io$rna.sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
io$rna.pseudobulk.sce <- paste0(io$basedir,"/results/rna/pseudobulk/SingleCellExperiment.rds")

# ATAC archR
io$archR.directory <- paste0(io$basedir,"/processed/atac/archR")

#############
## Options ##
#############

opts <- list()

opts$samples <- c(
  "HNES1_DUX4_overexpression_L001",
  "HNES1_wildtype_L001"
)

opts$chr <- paste0("chr",c(1:22,"X","Y"))

opts$ZGA_genes <- c(
  "CCNA1","DUXA","KDM4E","KHDC1L","LEUTX","MBD3L2","MBD3L3","PRAMEF1","PRAMEF12", 
  "PRAMEF11","RFPL2","RFPL4A","RFPL4B","SLC34A2","TRIM43","TRIM43B","TRIM49","TRIM49B","ZNF296","ZSCAN4"
)

###################
## Load metadata ##
###################

# sample_metadata <- fread(io$metadata)
# .[pass_QC==T] %>% 
# .[batch%in%opts$batches] %>%
# .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ","_")] %>%
# .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/","_")] %>%
# .[,celltype.mapped:=factor(celltype.mapped, levels=names(opts$celltype.colors))]


###################
## Edit metadata ##
###################

# io$metadata <- "/Users/ricard/data/DUX4_hESCs_multiome/processed/atac/archR/sample_metadata_after_archR.txt.gz"
# sample_metadata <- fread(io$metadata)
# sample_metadata[,stage:=substr(sample,1,4)]
# table(sample_metadata$stage)
# fwrite(sample_metadata, io$metadata, sep="\t", quote=F, na="NA")


