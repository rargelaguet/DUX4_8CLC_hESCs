
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/differential/archr_differential_accessibility_peaks.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--matrix',    type="character",    help='Matrix to use, see getAvailableMatrices')
p$add_argument('--test',      type="character",    help='Statistical test')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--group_variable',    type="character",    help='group variable')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args <- list()
# args$metadata <- "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
# args$matrix <- "GeneScoreMatrix_TSS"
# args$test <- "wilcoxon"
# args$group_variable <- "cluster"
# args$groupA <- "cluster_1"
# args$groupB <- "cluster_2"
# args$outfile <- "/bi/group/reik/ricard/data/DUX4_hESCs_multiome/results/atac/archR/differential/differential_GeneScoreMatrix_TSS_cluster1_vs_cluster2.txt.gz"
# args$threads <- 1
## END TEST

# Sanity checks
stopifnot(args$test%in%c("binomial","ttest","wilcoxon"))

###################
## Load settings ##
###################

source(here::here("settings.R"))
source(here::here("utils.R"))

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & sample%in%opts$samples]

stopifnot(args$group_variable%in%colnames(sample_metadata))
sample_metadata <- sample_metadata %>% .[get(args$group_variable)%in%c(args$groupA,args$groupB)]
stopifnot(args$groupA %in% unique(sample_metadata[[args$group_variable]]))
stopifnot(args$groupB %in% unique(sample_metadata[[args$group_variable]]))

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R")) 

# Subset ArchR object
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

# Sanity checks
stopifnot(args$matrix %in% getAvailableMatrices(ArchRProject))

###########################
## Update ArchR metadata ##
###########################

foo <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")
stopifnot(all(rownames(foo) == rownames(getCellColData(ArchRProject.filt))))

ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = foo[[args$group_variable]], 
  name = args$group_variable,
  cells = rownames(foo),
  force = TRUE
)

stopifnot(c(args$groupA,args$groupB)%in%getCellColData(ArchRProject.filt)[[args$group_variable]])

#######################
## Differential test ##
#######################

markerTest <- getMarkerFeatures(
  ArchRProject.filt, 
  useMatrix = args$matrix,
  groupBy = args$group_variable,
  testMethod = args$test,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = args$groupA,
  bgdGroups = args$groupB
)

####################
## Prepare output ##
####################

dt.1 <- rowData(markerTest) %>% as.data.table %>%
  setnames("seqnames","chr") %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("chr","start","end","strand"):=NULL]

dt.2 <- data.table(
  # Log2FC = round(assay(markerTest,"Log2FC")[,1],3),
  MeanDiff = round(assay(markerTest,"MeanDiff")[,1],3),
  FDR = assay(markerTest,"FDR")[,1],
  AUC = round(assay(markerTest,"AUC")[,1],3)
)

dt <- cbind(dt.1,dt.2) %>% setorder(-AUC)

##########
## Save ##
##########

fwrite(dt, args$outfile, sep="\t")

################
## START TEST ##
################

# R.utils::sourceDirectory("/Users/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)

# ArchRProj = ArchRProject.filt
# groupBy = "celltype.predicted"
# useGroups = args$groupA
# bgdGroups = args$groupB
# useMatrix = "PeakMatrix"
# bias = c("TSSEnrichment","log10(nFrags)")
# normBy = NULL
# testMethod = "wilcoxon"
# maxCells = 500
# scaleTo = 10^4
# threads = 1
# k = 100
# bufferRatio = 0.8
# binarize = FALSE
# useSeqnames = NULL
# verbose = TRUE
# logFile = createLogFile("getMarkerFeatures")

##############
## END TEST ##
##############