suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/processing/4_add_GeneScore_matrices.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ndims',           type="integer",    default=30,                  help='Number of LSI dimensions')
p$add_argument('--batch.variable',  type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--batch.method',    type="character",  default="MNN",               help='Batch correctin method ("Harmony" or "MNN")')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  default="celltype.predicted",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$pdfdir <- sprintf("%s/pdf",args$outdir); dir.create(io$pdfdir,showWarnings = F)

# Options
opts$lsi.iterations = 2
opts$lsi.cluster.resolution = 2

## START TEST ##
# args$metadata <- io$metadata
# args$samples <- opts$samples[1:2]
# args$nfeatures <- 15000
# args$matrix <- "PeakMatrix"
# args$ndims <- 30
# args$colour_by <- c("celltype.predicted","sample","log_nFrags_atac","doublet_call")
# args$batch.variable <- "stage"
# args$batch.method <- "Harmony"
# args$outdir <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction")
## END TEST ##

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  # .[pass_atacQC==TRUE & doublet_call==FALSE & sample%in%args$samples] %>%
  .[pass_atacQC==TRUE & sample%in%args$samples] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

table(sample_metadata$sample)
table(sample_metadata$celltype.predicted)

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))


###################
## Sanity checks ##
###################

stopifnot(args$colour_by %in% colnames(sample_metadata))

if (length(args$batch.variable)>0) {
  stopifnot(args$batch.variable%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch.variable]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch.variable))
    args$batch.variable <- NULL
  } else {
    library(batchelor)
  }
}

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]

ArchRProject.filt@sampleColData <- ArchRProject.filt@sampleColData[args$samples,,drop=F]

# Update ArchR metadata
stopifnot(sample_metadata$cell == rownames(getCellColData(ArchRProject.filt)))
ArchRProject.filt <- addCellColData(ArchRProject.filt,
  data = sample_metadata$stage, 
  name = "stage",
  cells = sample_metadata$cell,
  force = TRUE
)

#######################
## Feature selection ##
#######################

if (args$matrix=="PeakMatrix_filt") {
  
  # Load peak variability estimates
  peak.variability.dt <- fread(io$archR.peak.variability)
  
  # Define highly variable peaks
  peaks <- peak.variability.dt %>% 
    .[,peak:=sub("_",":",peak)] %>% .[,peak:=sub("_","-",peak)] %>%
    setorder(-variance_pseudobulk) %>% 
    head(n=args$nfeatures) %>% 
    .$peak
  
  # Subset peaks in the ArchR object
  names(ArchRProject.filt@peakSet) <- sprintf("%s:%s-%s",seqnames(ArchRProject.filt@peakSet), start(ArchRProject.filt@peakSet), end(ArchRProject.filt@peakSet))
  ArchRProject.filt@peakSet[peaks,] <- ArchRProject.filt@peakSet[peaks,]
  
  ArchRProject.filt <- addFeatureMatrix(
    input = ArchRProject.filt,
    features = ArchRProject.filt@peakSet,
    matrixName = "PeakMatrix_filt",
    binarize = TRUE
  )
  
  # foo <- getMatrixFromProject(ArchRProject.filt, useMatrix = "PeakMatrix_filt", binarize = TRUE)
  # foo <- getMatrixFromProject(ArchRProject.filt, useMatrix = "PeakMatrix", binarize = TRUE)
}

###########################
## Latent Semantic Index ##
###########################

# Iterative LSI: two iterations
ArchRProject.filt <- addIterativeLSI(
  ArchRProj = ArchRProject.filt,
  useMatrix = args$matrix, 
  name = "IterativeLSI", 
  firstSelection = "Top",
  depthCol = "nFrags",
  iterations = opts$lsi.iterations, 
  # clusterParams = list(
  #  resolution = opts$lsi.cluster.resolution, 
  #  sampleCells = 10000, 
  #  n.start = 10
  # ), 
  saveIterations = FALSE,
  varFeatures = args$nfeatures, 
  force = TRUE
)

# Correlation between latent dimensions and number of peaks
# cor(getReducedDims(ArchRProject.filt, "IterativeLSI"),ArchRProject.filt$nFrags)[,1] %>% abs %>% sort(decreasing = T)

############################
## LSI + Batch correction ##
############################

if (length(args$batch.variable)>0) {
  print(sprintf("Applying %s batch correction for variable: %s", args$batch.method, args$batch.variable))
  outfile <- sprintf("%s/%s_lsi_features%d_dims%d_%sbatchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$nfeatures, args$ndims, args$batch.method, paste(args$batch.variable,collapse="-"))
  
  # Harmony
  if (args$batch.method=="Harmony") {
    ArchRProject.filt <- addHarmony(
      ArchRProj = ArchRProject.filt,
      reducedDims = "IterativeLSI",
      name = "IterativeLSI_Harmony",
      groupBy = "stage",
      force = TRUE
    )
    lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI_Harmony") %>% round(3) %>% 
      as.data.table(keep.rownames = T) %>% setnames("rn","cell")
    
  # MNN  
  } else if  (args$batch.method=="MNN") {
    
    stop("Not implemented")
    
    lsi.corrected <- getReducedDims(ArchRProject.filt, "IterativeLSI") %>% as.data.table %>%
      split(ArchRProject.filt$sample) %>% map(as.matrix) %>% reducedMNN %>% .[["corrected"]]
    rownames(lsi.corrected) <- rownames(getReducedDims(ArchRProject.filt, "IterativeLSI") )
    ArchRProject.filt@reducedDims[["IterativeLSI_MNN"]] <- SimpleList(
      matDR = lsi.corrected,
      params = NA,
      date = Sys.time(),
      scaleDims = NA,
      corToDepth = NA
    )
    lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI_MNN") %>% round(3) %>% 
      as.data.table(keep.rownames = T) %>% setnames("rn","cell")
  } else {
    stop("Batch correction method not recognised")
  }
} else {
  outfile <- sprintf("%s/%s_lsi_features%d_ndims%d.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$nfeatures, args$ndims)
  lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI") %>% round(3) %>% 
    as.data.table(keep.rownames = T) %>% setnames("rn","cell")
}

# Save LSI coordinates
fwrite(lsi.dt, outfile)

##########
## UMAP ##
##########

for (i in args$n_neighbors) {
  for (j in args$min_dist) {
    
    # Define the latent space to run UMAP on
    if (length(opts$batch.correction)>0) {
      if (args$batch.method=="Harmony") {
        dimred <- "IterativeLSI_Harmony"
      } else if  (args$batch.method=="MNN") {
        dimred <- "IterativeLSI_MNN"
      } else {
        stop("Batch correction method not recognised")
      }
    } else {
      dimred <- "IterativeLSI"
    }
    
    # Run UMAP
    ArchRProject.filt <- addUMAP(
      ArchRProj = ArchRProject.filt, 
      reducedDims = dimred,
      name = "UMAP",
      metric = "cosine",
      nNeighbors = i, 
      minDist = j, 
      seed = args$seed,
      saveModel = FALSE,
      force = TRUE
    )
    
    # Fetch UMAP coordinates
    umap.dt <- getEmbedding(ArchRProject.filt,"UMAP") %>%
      as.data.table(keep.rownames = T) %>%
      setnames(c("cell","umap1","umap2"))
    
    # Plot
    to.plot <- umap.dt %>%
      merge(sample_metadata,by="cell")
    
    for (k in args$colour_by) {
      
      p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", fill=k)) +
        geom_point(size=1.5, shape=21, stroke=0.05) +
        # ggrastr::geom_point_rast(size=1.5, shape=21, stroke=0.05) +  # DOES NOT WORK IN THE CLUSTER
        theme_classic() +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
        )
      
      # if (k%in%c("celltype.mapped","celltype.predicted")) {
      #   p <- p + scale_fill_manual(values=opts$celltype.colors) +
      #     theme(
      #       legend.position="none",
      #       legend.title=element_blank()
      #     )
      # }
      
      # Save UMAP plot
      outfile <- sprintf("%s/%s_umap_nfeatures%d_ndims%d_neigh%d_dist%s_%s.pdf",io$pdfdir, paste(args$samples,collapse="-"), args$nfeatures, args$ndims, i, j, k)
      pdf(outfile, width=7, height=5)
      print(p)
      dev.off()
    }
    
    # Save UMAP coordinates
    outfile <- sprintf("%s/%s_umap_nfeatures%d_ndims%d_neigh%d_dist%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$nfeatures, args$ndims, i, j)
    fwrite(umap.dt, outfile)
  }
}

