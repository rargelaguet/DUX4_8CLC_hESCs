suppressPackageStartupMessages(library(ArchR))

#####################
## Define settings ##
#####################

here::i_am("atac/archR/processing/1_create_archR_project.R")
source(here::here("settings.R"))

# Options

# I/O
io$output.directory <- file.path(io$basedir,"processed/atac/archR")

io$arrowFiles <- c(
  "HNES1_wildtype_L001" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/HNES1_wildtype_L001.arrow"),
  "HNES1_DUX4_overexpression_L001" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/HNES1_DUX4_overexpression_L001.arrow")
)

# ArchR options
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

############################
## create an ArchRProject ##
############################

ArchRProject <- ArchRProject(
  ArrowFiles = io$arrowFiles, 
  outputDirectory = io$output.directory,
  copyArrows = FALSE
)

##########
## Save ##
##########

saveArchRProject(ArchRProject)
# saveArchRProject(ArchRProject, io$output.directory)