suppressPackageStartupMessages(library(ArchR))

#####################
## Define settings ##
#####################

# if (grepl("ricard",Sys.info()['nodename'])) {
#   source("/Users/ricard/gastrulation_multiome_10x/settings.R")
#   source("/Users/ricard/gastrulation_multiome_10x/utils.R")
# } else if (grepl("ebi",Sys.info()['nodename'])) {
#   source("/homes/ricard/gastrulation_multiome_10x/settings.R")
#   source("/homes/ricard/gastrulation_multiome_10x/utils.R")
# } else {
#   stop("Computer not recognised")
# }
setwd(io$archR.directory)

################
## Define I/O ##
################

# io$outdir <- paste0(io$archR.directory,"/pdf")

####################
## Define options ##
####################

addArchRGenome("hg38")
addArchRThreads(threads = 1) 

########################
## Load ArchR project ##
########################

ArchRProject <- loadArchRProject(io$archR.directory)

# Load ArchR projectMetadata
# if (file.exists(io$archR.projectMetadata)) {
# 	ArchRProject@projectMetadata <- readRDS(io$archR.projectMetadata)
# }

# Load peaks
if (file.exists(io$archR.peakSet.granges)) {
	ArchRProject <- addPeakSet(ArchRProject, peakSet = readRDS(io$archR.peakSet.granges), force = TRUE)
}

# Load motif annotations over peaks
# io$archR.peakAnnotation <- sprintf("%s/Annotations/peakAnnotation.rds",io$archR.directory)
# if (file.exists(io$archR.peakSet.granges)) {
# 	ArchRProject@peakAnnotation <- readRDS(io$archR.peakAnnotation)
# }

# Add background peaks
io$archR.bgdPeaks <- file.path(io$archR.directory, "Background-Peaks.rds")
if (!"bgdPeaks" %in% metadata(getPeakSet(ArchRProject))$bgdPeaks) {
	if (file.exists(io$archR.bgdPeaks)) metadata(ArchRProject@peakSet)$bgdPeaks <- io$archR.bgdPeaks
}

##########
## TEST ##
##########

# ArchRProject@peakSet <- readRDS(io$archR.peakSet.granges)
# seqlevels(ArchRProject@peakSet) <- sort(seqlevels(ArchRProject@peakSet))
# ArchRProject@peakSet <- sort(ArchRProject@peakSet)


# getAvailableMatrices(ArchRProject)

# io$arrow.files <- opts$samples %>% 
#   # map_chr(~ sprintf("%s/%s.arrow",io$archR.directory,.))
#   map_chr(~ sprintf("%s.arrow",.))
# 
# ArchRProject <- ArchRProject(
#   ArrowFiles = io$arrow.files,
#   # outputDirectory = "ArchROutput",
#   outputDirectory = io$archR.directory,
#   copyArrows = FALSE
# )
# saveArchRProject(ArchRProject)
