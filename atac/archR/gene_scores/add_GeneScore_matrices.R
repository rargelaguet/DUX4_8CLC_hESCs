# https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html

#####################
## Define settings ##
#####################

here::i_am("atac/archR/gene_scores/add_GeneScore_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

#####################
## Define gene set ##
#####################

# getGeneAnnotation(ArchRProject)
genes.gr <- getGenes(ArchRProject)

getAvailableMatrices(ArchRProject)

#########################################
## Add Gene Scores using default model ##
#########################################

# Note that this will add the matrices to the arrowFiles
addGeneScoreMatrix(
  input = ArchRProject,
  genes = genes.gr,
  geneModel = "exp(-abs(x)/5000) + exp(-1)", # string should be a function of x, where x is the distance from the TSS.
  matrixName = "GeneScoreMatrix_distal",
  extendUpstream = c(1000, 1e+05),
  extendDownstream = c(1000, 1e+05),  # The minimum and maximum number of bp downstream of the transcription termination site to consider for gene activity score calculation.
  geneUpstream = 5000,                # Number of bp upstream the gene to extend the gene body.
  geneDownstream = 0,                 # Number of bp downstream the gene to extend the gene body.
  tileSize = 500,                     # The size of the tiles used for binning counts prior to gene activity score calculation.
  geneScaleFactor = 5,                # A numeric scaling factor to weight genes based on the inverse of their length 
  scaleTo = 10000,                    # Each column in the calculated gene score matrix will be normalized
  excludeChr = c("chrY", "chrM")
)

#############################################
## Add Gene Scores ignoring distal regions ##
#############################################


# TSS
addGeneScoreMatrix(
  input = ArchRProject,
  genes = genes.gr,
  useTSS = TRUE,
  extendTSS = TRUE,
  geneModel = "1",                    # string should be a function of x, where x is the distance from the TSS.
  matrixName = "GeneScoreMatrix_TSS",
  extendUpstream = c(0,0),
  extendDownstream = c(0,0),          # The minimum and maximum number of bp downstream of the transcription termination site to consider for gene activity score calculation.
  geneUpstream = 500,                # Number of bp upstream the gene to extend the gene body.
  geneDownstream = 100,                 # Number of bp downstream the gene to extend the gene body.
  tileSize = 100,                     # The size of the tiles used for binning counts prior to gene activity score calculation.
  geneScaleFactor = 1,                # A numeric scaling factor to weight genes based on the inverse of their length 
  scaleTo = 10000,                    # Each column in the calculated gene score matrix will be normalized
  excludeChr = c("chrY", "chrM"),
  force = TRUE
)


##########
## Test ##
##########

# atac.GeneScoreMatrix.se <- getMatrixFromProject(ArchRProject, binarize = FALSE, useMatrix = "GeneScoreMatrix_nodistal_v2")
# rownames(atac.GeneScoreMatrix.se) <- rowData(atac.GeneScoreMatrix.se)$name
# tmp <- rowMeans(assay(atac.GeneScoreMatrix.se)) %>% sort(decreasing = T)
# head(tmp)
# tmp[c("Auts2","Nav2")]
# tmp[c("Actb","Polr2f")]
