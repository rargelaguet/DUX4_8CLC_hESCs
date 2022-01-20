
here::i_am("rna/differential/differential.R")

###########################################################
## Script to do differential expression between lineages ##
###########################################################

suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$groupA <- "FALSE"
args$groupB <- "TRUE"
args$outfile <- "/Users/argelagr/data/DUX4_hESCs_multiome/results/rna/differential/eight_cell_like_differential_rna.tsv.gz"
## END TEST

#########
## I/O ##
#########

source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("rna/differential/utils.R"))

# Sanity checks
# stopifnot(args$stages%in%opts$stages)
# stopifnot(args$groupA%in%opts$celltypes)
# stopifnot(args$groupB%in%opts$celltypes)

#############
## Options ##
#############

# Define groups
opts$groups <- c(args$groupA,args$groupB)

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1

# For a given gene, the minimum fraction of cells that must express it in at least one group
opts$min_detection_rate_per_group <- 0.40

###############
## Load data ##
###############

# Update cell metadata
# io$metadata <- paste0(io$basedir,"/results/rna/celltype_denoising/sample_metadata_after_celltype_denoising.txt.gz")
sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & eight_cell_like_ricard%in%opts$groups] %>%
  setnames("eight_cell_like_ricard","group") %>%
  .[,c("cell","group")]

# sample_metadata[,group:=as.factor( c("A","B")[as.numeric(stage_lineage%in%opts$groupB)+1] )]

# Sort cells so that groupA comes before groupB
sample_metadata[,group:=factor(group,levels=opts$groups)] %>% setorder(group)
table(sample_metadata$group)

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$rna.sce, cells = sample_metadata$cell, normalise = TRUE)
sce$group <- sample_metadata$group

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[symbol%in%rownames(sce)] %>%
  .[,c("symbol","ens_id")] %>%
  setnames("symbol","gene")

################
## Parse data ##
################

# calculate detection rate per gene
cdr.dt <- data.table(
  rownames(sce),
  rowMeans(logcounts(sce[,sce$group==opts$groups[1]])>0) %>% round(2),
  rowMeans(logcounts(sce[,sce$group==opts$groups[2]])>0) %>% round(2)
) %>% setnames(c("gene",sprintf("detection_rate_%s",opts$groups[1]),sprintf("detection_rate_%s",opts$groups[2])))
# .[,cdr_diff:=abs(out[,(sprintf("detection_rate_%s",opts$groups[1])),with=F][[1]] - out[,(sprintf("detection_rate_%s",opts$groups[2])),with=F][[1]])] %>%

# Filter genes
sce <- sce[rownames(sce)%in%gene_metadata$gene,]

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, opts$groups, opts$min_detection_rate_per_group) %>%
  # Add sample statistics
  .[,c("groupA_N","groupB_N"):=list(table(sample_metadata$group)[1],table(sample_metadata$group)[2])]%>% 
  # setnames(c("groupA_N","groupB_N"),c(sprintf("N_%s",opts$groups[1]),sprintf("N_%s",opts$groups[2]))) %>%
  # Add gene statistics
  merge(cdr.dt, all.y=T, by="gene") %>%
  merge(gene_metadata, all.y=T, by="gene") %>%
  # Calculate statistical significance
  .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC)] %>%
  .[is.na(sig),sig:=FALSE] %>%
  setorder(-sig, padj_fdr, na.last=T)

# Parse columns
out[,c("p.value","padj_fdr","logFC","log_padj_fdr"):=list(signif(p.value,digits=3), signif(padj_fdr,digits=3), round(logFC,3),round(log_padj_fdr,3))]

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
