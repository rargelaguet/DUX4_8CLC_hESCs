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
  io$archR.directory <- paste0(io$basedir,"/processed/atac/archR")
} else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
  if (grepl("DUX4_hESCs_multiome", Sys.info()['effective_user'])) {
    stop()
  } else if (grepl("argelag", Sys.info()['effective_user'])) {
    io$basedir <-'/bi/group/reik/ricard/data/DUX4_hESCs_multiome'
    io$gene_metadata <- "/bi/group/reik/ricard/data/ensembl/human/v93/bioMart/all_genes/Hsapiens_genes_BioMart.93.txt.gz"
  }
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")

# RNA
io$rna.seurat <- paste0(io$basedir,"/processed/rna/seurat.rds")
io$rna.sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
io$rna.pseudobulk.sce <- paste0(io$basedir,"/results/rna/pseudobulk/SingleCellExperiment.rds")


#############
## Options ##
#############

opts <- list()

opts$samples <- c(
  "HNES1_DUX4_overexpression_L001",
  "HNES1_wildtype_L001"
)

opts$chr <- paste0("chr",c(1:22,"X","Y"))

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


