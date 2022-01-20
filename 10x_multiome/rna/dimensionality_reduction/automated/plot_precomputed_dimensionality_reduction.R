
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

# I/O
io$indir <- paste0(io$basedir,"/results/rna/dimensionality_reduction")
io$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction/test")

# Options
opts$samples <- list(
  # "E7.5_rep1" = "E7.5_rep1",
  # "E7.5_rep2" = "E7.5_rep2",
  # "E8.5_rep1" = "E8.5_rep1",
  # "E8.5_rep2" = "E8.5_rep2",
  # "E7.5" = c("E7.5_rep1","E7.5_rep2"),
  # "E8.5" = c("E8.5_rep1","E8.5_rep2")
  "all_cells" = c("E7.5_rep1","E7.5_rep2","E8.0_rep1","E8.0_rep2","E8.5_rep1","E8.5_rep2")
)

opts$npcs <- c(30)
opts$nfeatures <- c(2500)
# opts$batch.correction <- "sample"

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
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


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]

#############################
## Load PCA representation ##
#############################

pca.dt <- names(opts$samples) %>% map(function(i) {
  opts$nfeatures %>% map(function(j) {
    opts$npcs %>% map(function(k) {
      # file <- sprintf("%s/%s/%s_pca_features%d_pcs%d.txt.gz",io$indir, i, paste(opts$samples[[i]],collapse="-"), j, k)
      file <- sprintf("%s/%s/%s_pca_features%d_pcs%d_batchcorrectionbysample.txt.gz",io$indir, i, paste(opts$samples[[i]],collapse="-"), j, k)
      # file <- sprintf("%s/%s_pca_features%d_pcs%d_batchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs,paste(args$batch.correction,collapse="-"))
      if (file.exists(file)) {
        dt <- fread(file) %>% 
          merge(sample_metadata[,c("cell","celltype.mapped")], by="cell") %>%
          .[,c("nfeatures","npcs"):=list(i,j)]
      } else {
        print(sprintf("%s does not exist",file))
        return(NULL)
      }
    }) %>% rbindlist
  }) %>% rbindlist 
}) %>% rbindlist

##############
## Run UMAP ##
##############

##############################
## Load UMAP representation ##
##############################

opts$umap.neighbours <- 25
opts$umap.distance <- 0.3

umap.dt <- names(opts$samples) %>% map(function(i) {
  opts$nfeatures %>% map(function(j) {
    opts$npcs %>% map(function(k) {
      # file <- sprintf("%s/%s/%s_pca_features%d_pcs%d.txt.gz",io$indir, i, paste(opts$samples[[i]],collapse="-"), j, k)
      file <- sprintf("%s/%s/%s_umap_features%d_pcs%d_neigh%s_dist%s.txt.gz",io$indir, i, paste(opts$samples[[i]],collapse="-"), j, k, opts$umap.neighbours, opts$umap.distance)
      # file <- sprintf("%s/%s_pca_features%d_pcs%d_batchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs,paste(args$batch.correction,collapse="-"))
      if (file.exists(file)) {
        dt <- fread(file) %>% 
          merge(sample_metadata[,c("cell","celltype.mapped","stage")], by="cell") %>%
          .[,c("sample","nfeatures","npcs","neighbours","min_dist"):=list(i,j,k,opts$umap.neighbours,opts$umap.distance)]
      } else {
        print(sprintf("%s does not exist",file))
        return(NULL)
      }
    }) %>% rbindlist
  }) %>% rbindlist 
}) %>% rbindlist

##########
## Plot ##
##########

to.plot <- umap.dt[sample=="all_cells"]

p <- ggplot(to.plot, aes_string(x="UMAP1", y="UMAP2", fill="stage")) +
  geom_point(size=1.5, shape=21, stroke=0.05) +
  scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none",
    legend.title=element_blank()
  )

# outfile <- sprintf("%s/%s_umap_features%d_pcs%d_neigh%d_dist%s_%s.pdf",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs, i, j, k)
# pdf(outfile, width=7, height=5)
print(p)
# dev.off()
