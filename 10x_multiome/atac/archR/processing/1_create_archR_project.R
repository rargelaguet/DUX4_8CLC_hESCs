suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/processing/1_create_archR_project.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--arrow_files',           type="character",  nargs='+',      help='Arrow files')
p$add_argument('--genome',           type="character", default="hg38",      help='Genome')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$arrow_files <- c(
#   "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/processed/atac/archR/HNES1_DUX4_overexpression_L001.arrow",
#   "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/processed/atac/archR/HNES1_wildtype_L001.arrow"
# )
# args$genome <- "hg38"
# args$threads <- 1
# args$outdir <- "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/processed/atac/archR"
## END TEST ##

#####################
## Define settings ##
#####################

here::i_am("atac/archR/processing/1_create_archR_project.R")
source(here::here("settings.R"))

# ArchR options
addArchRGenome(args$genome)

############################
## create an ArchRProject ##
############################

ArchRProject <- ArchRProject(
  ArrowFiles = args$arrow_files, 
  outputDirectory = args$outdir,
  copyArrows = FALSE
)

##########
## Save ##
##########

saveArchRProject(ArchRProject)