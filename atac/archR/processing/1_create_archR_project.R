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
  # "E7.5_rep1" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/E7.5_rep1.arrow"),
  # "E7.5_rep2" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/E7.5_rep2.arrow"),
  # "E8.0_rep1" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/E8.0_rep1.arrow"),
  # "E8.0_rep2" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/E8.0_rep2.arrow"),
  # "E8.5_rep1" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/E8.5_rep1.arrow"),
  # "E8.5_rep2" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/E8.5_rep2.arrow"),
  "E8.75_rep1" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/E8.75_rep1.arrow"),
  "E8.75_rep2" = file.path(io$basedir,"processed/atac/archR/ArrowFiles/E8.75_rep2.arrow")
)

# ArchR options
addArchRThreads(threads = 1) 
addArchRGenome("mm10")

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