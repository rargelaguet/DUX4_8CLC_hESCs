suppressPackageStartupMessages(library(ArchR))

#####################
## Define settings ##
#####################

here::i_am("atac/archR/processing/0_create_arrow_files.R")
source(here::here("settings.R"))

# Options
opts$samples <- c(
  # "E7.5_rep1",
  # "E7.5_rep2",
  # "E8.0_rep1",
  # "E8.0_rep2",
  # "E8.5_rep1",
  # "E8.5_rep2",
  "E8.75_rep1",
  "E8.75_rep2"
)

# I/O
io$output.directory <- file.path(io$basedir,"processed/atac/archR")
setwd(io$output.directory)

io$fragment_files <- opts$samples %>% 
  map_chr(~ sprintf("%s/original/%s/atac_fragments.tsv.gz",io$basedir,.)) %>%
  set_names(opts$samples)

# Options
opts$min.fragments <- 2500
opts$filterTSS.score <- 3

# ArchR options
addArchRThreads(threads = 2) 
addArchRGenome("mm10")


########################
## create Arrow Files ##
########################

# Steps:
# (1) Read accessible fragments from the provided input files.
# (2) Calculate quality control information for each cell (i.e. TSS enrichment scores and nucleosome info).
# (3) Filter cells based on quality control parameters.
# (4) Create a genome-wide TileMatrix using 500-bp bins.
# (5) Create a GeneScoreMatrix using the custom geneAnnotation that was defined when we called addArchRGenome().

ArrowFiles <- createArrowFiles(
  inputFiles = io$fragment_files,
  sampleNames = names(io$fragment_files),
  outputNames = names(io$fragment_files),
  addTileMat = FALSE,
  addGeneScoreMat = TRUE,
  excludeChr = c("chrM", "chrY"),

  # QC metrics
  minFrags = opts$min.fragments,  # The minimum number of fragments per cell
  minTSS = opts$filterTSS.score   # The minimum TSS enrichment score per cell
)
