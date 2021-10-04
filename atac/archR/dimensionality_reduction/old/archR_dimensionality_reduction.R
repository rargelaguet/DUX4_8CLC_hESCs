
#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
# io$metadata <- paste0(io$basedir,"/results/atac/archR/celltype_assignment/sample_metadata_after_archR.txt.gz")
io$outdir <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/test")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

# LSI parameters
# opts$lsi.iterations <- 1
# opts$lsi.cluster.resolution <- 2
opts$lsi.iterations <- 2
opts$lsi.cluster.resolution <- 2
opts$lsi.varFeatures <- 25000
opts$lsi.dims <- 25
opts$matrix <- "PeakMatrix"

# UMAP parameters
opts$umap.seed <- 42
opts$umap.neighbours <- 25
opts$umap.minDist <- 0.3

opts$batch.correction <- FALSE

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

###########################
## Latent Semantic Index ##
###########################

# Iterative LSI: two iterations
ArchRProject.filt <- addIterativeLSI(ArchRProject.filt,
  useMatrix = opts$matrix, 
  name = "IterativeLSI", 
  firstSelection = "Top", # "Top" or "Var"
  depthCol = "nFrags",
  iterations = opts$lsi.iterations, 
  clusterParams = list(
    resolution = opts$lsi.cluster.resolution, 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = opts$lsi.varFeatures, 
  force = TRUE
)

# Correlation between latent dimensions and number of peaks
cor(getReducedDims(ArchRProject.filt, "IterativeLSI"),ArchRProject.filt$nFrags)[,1] %>% abs %>% sort(decreasing = T)

#############################
## Batch effect correction ##
#############################

# Harmony
if (opts$batch.correction) {
  
  ArchRProject.filt@cellColData$stage <- substr(ArchRProject.filt@cellColData$sample,1,4)
  
  ArchRProject.filt <- addHarmony(
    ArchRProj = ArchRProject.filt,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "stage",
    force = TRUE
  )
  
}

# MNN
# Probably makes sense to first integrate E7.5 samples, and then E8.5 samples. sue "merge.order	"
if (opts$batch.correction) {
  stop()
  dimred <- "MNN" 
  library(batchelor)
  pca <- getReducedDims(ArchRProject.filt, "IterativeLSI") %>% as.data.table %>%
    split(ArchRProject.filt$sample) %>% map(as.matrix)
  pca.corrected <- reducedMNN(pca)$corrected
  ArchRProject.filt@reducedDims[["MNN"]] <- SimpleList(
    matDR = pca.corrected,
    params = NA,
    date = Sys.time(),
    scaleDims = NA,
    corToDepth = NA
  )
}

##########
## UMAP ##
##########

if (opts$batch.correction) {
  dimred <- "Harmony"
} else {
  dimred <- "IterativeLSI"
}

ArchRProject.filt <- addUMAP(ArchRProject.filt, 
  reducedDims = dimred,
  name = "UMAP",
  metric = "cosine",
  nNeighbors = opts$umap.neighbours, 
  minDist = opts$umap.minDist, 
  seed = opts$umap.seed,
  force = TRUE
)
head(ArchRProject.filt@embeddings$UMAP$df)

###########
## t-SNE ##
###########

# ArchRProject.filt <- addTSNE(ArchRProject.filt,
#   reducedDims = "IterativeLSI",
#   name = "TSNE",
#   perplexity = 30
# )
# head(ArchRProject.filt@embeddings$TSNE$df)

##########
## Plot ##
##########


to.plot <- getEmbedding(ArchRProject.filt,"UMAP") %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("cell","umap1","umap2")) %>% merge(sample_metadata,by="cell")
  
p1 <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype.mapped), size=2, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# pdf(sprintf("%s/archr_umap_celltype_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p1)
# dev.off()

p2 <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(color=log10(nFrags_atac)), size=1) +
  scale_color_gradient(low = "gray80", high = "red") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# pdf(sprintf("%s/archr_umap_nFragsATAC_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p2)
# dev.off()

p3 <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=sample), size=1.5, shape=21, color="black", stroke=0.05) +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# pdf(sprintf("%s/archr_umap_sample_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p3)
# dev.off()

##########
## Save ##
##########
  
# save umap coordinates
to.save <- getEmbedding(ArchRProject.filt,"UMAP") %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("archR_cell","V1","V2"))
# fwrite(to.save, sprintf("%s/archr_umap_coordinates_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.txt",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))

# save options
# options.to.save <- opts[c("lsi.iterations", "lsi.cluster.resolution", "lsi.varFeatures", "lsi.dims", "umap.seed", "umap.neighbours", "umap.minDist")]
# saveRDS(options.to.save, sprintf("%s/hyperparameters.rds",io$outdir))

#############
## Explore ##
#############

getAvailableMatrices(ArchRProject.filt)
deviations.se <- getMatrixFromProject(ArchRProject.filt, "DeviationMatrix")
