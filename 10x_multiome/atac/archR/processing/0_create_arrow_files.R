suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/processing/0_create_arrow_files.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples',           type="character",  nargs='+',      help='Samples')
p$add_argument('--fragments_files',           type="character",  nargs='+',      help='ATAC Fragments files')
p$add_argument('--genome',           type="character", default="hg38",      help='Genome')
p$add_argument('--min_fragments',     type="integer",    default=1000,   help='Minimum number of ATAC fragments')
p$add_argument('--max_fragments',     type="integer",    default=1e7,    help='Maximum number of ATAC fragments')
p$add_argument('--min_tss_score',   type="double",     default=2.5,    help='Minimum TSS score threshold')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$fragments_files <- c(
#   "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/original/HNES1_wildtype_L001/atac_fragments.tsv.gz",
#   "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/original/HNES1_DUX4_overexpression_L001/atac_fragments.tsv.gz"
# )
# args$samples <- c("HNES1_wildtype_L001","HNES1_DUX4_overexpression_L001")
# args$genome <- "hg38"
# args$min_fragments <- 1000
# args$max_fragments <- 1e7
# args$min_tss_score <- 2.5
# args$threads <- 1
# args$outdir <- "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/processed/atac/archR/test"
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
setwd(args$outdir)


# ArchR options
addArchRThreads(threads=args$threads) 
addArchRGenome(args$genome)


########################
## create Arrow Files ##
########################

ArrowFiles <- createArrowFiles(
  inputFiles = args$fragments_files,
  sampleNames = args$samples,
  outputNames = args$samples,
  addTileMat = FALSE,
  addGeneScoreMat = FALSE,
  excludeChr = c("chrM", "chrY"),

  # QC metrics
  minFrags = args$min_fragments,  # The minimum number of fragments per cell
  maxFrags = args$max_fragments,  # The maximum number of fragments per cell
  minTSS = args$min_tss_score   # The minimum TSS enrichment score per cell
)
