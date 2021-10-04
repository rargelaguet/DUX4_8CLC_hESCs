
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

# I/O
io$indir <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix")
io$outdir <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix/test")

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

opts$ndims <- c(50)
opts$nfeatures <- c(50000)

# opts$matrix <- "PeakMatrix"

opts$celltypes = c(
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

opts$aggregated.celltypes <- c(
  # "Erythroid1" = "Erythroid",
  # "Erythroid2" = "Erythroid",
  # "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

opts$remove.ExE.celltypes <- FALSE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & celltype.predicted%in%opts$celltypes] %>%
  .[,celltype.predicted:=factor(celltype.predicted,levels=opts$celltypes)]

# Remove ExE celltypes
if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.predicted%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

# Rename cell types
sample_metadata %>%
  .[,aggregated_celltype:=stringr::str_replace_all(celltype.predicted,opts$aggregated.celltypes)]

opts$aggregated_celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(sample_metadata$aggregated_celltype)]
stopifnot(unique(sample_metadata$aggregated_celltype)%in%names(opts$aggregated_celltype.colors))
sample_metadata %>%
  .[,aggregated_celltype:=factor(aggregated_celltype, levels=names(opts$celltype.colors))]

#############################
## Load LSI representation ##
#############################

lsi.dt <- names(opts$samples) %>% map(function(i) {
  opts$nfeatures %>% map(function(j) {
    opts$ndims %>% map(function(k) {
      # file <- sprintf("%s/%s/%s_pca_features%d_pcs%d.txt.gz",io$indir, i, paste(opts$samples[[i]],collapse="-"), j, k)
      file <- sprintf("%s/%s/%s_lsi_features%d_ndims%d.txt.gz",io$indir, i, paste(opts$samples[[i]],collapse="-"), j, k)
      # file <- sprintf("%s/%s_pca_features%d_pcs%d_batchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$ndims,paste(args$batch.correction,collapse="-"))
      if (file.exists(file)) {
        dt <- fread(file)# %>% merge(sample_metadata[,c("cell","celltype.predicted")], by="cell") %>% .[,c("nfeatures","ndims"):=list(i,j)]
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

opts$umap.neighbours <- c(15,30,45)
opts$umap.distance <- c(0.15,0.3,0.45)

umap.dt <- names(opts$samples) %>% map(function(i) {
  opts$nfeatures %>% map(function(j) {
    opts$ndims %>% map(function(k) {
      opts$umap.neighbours %>% map(function(neigh) {
        opts$umap.distance %>% map(function(dist) {
          # file <- sprintf("%s/%s/%s_pca_features%d_pcs%d.txt.gz",io$indir, i, paste(opts$samples[[i]],collapse="-"), j, k)
          file <- sprintf("%s/%s/%s_umap_nfeatures%d_ndims%d_neigh%s_dist%s.txt.gz",io$indir, i, paste(opts$samples[[i]],collapse="-"), j, k, neigh, dist)
          # file <- sprintf("%s/%s_pca_features%d_pcs%d_batchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$ndims,paste(args$batch.correction,collapse="-"))
          if (file.exists(file)) {
            dt <- fread(file) %>% 
              merge(sample_metadata[,c("cell","celltype.predicted","aggregated_celltype","stage")], by="cell") %>%
              .[,c("sample","nfeatures","ndims","neighbours","min_dist"):=list(i,j,k,neigh,dist)]
          } else {
            print(sprintf("%s does not exist",file))
            return(NULL)
          }
        }) %>% rbindlist
      }) %>% rbindlist
    }) %>% rbindlist
  }) %>% rbindlist 
}) %>% rbindlist

##########
## Plot ##
##########

to.plot <- umap.dt[sample=="all_cells" & neighbours==30 & min_dist==0.45]

p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", fill="stage")) +
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

# outfile <- sprintf("%s/%s_umap_features%d_pcs%d_neigh%d_dist%s_%s.pdf",args$outdir, paste(args$samples,collapse="-"), args$features, args$ndims, i, j, k)
# pdf(outfile, width=7, height=5)
# print(p)
# dev.off()


p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", fill="aggregated_celltype")) +
  # facet_wrap(~neighbours+min_dist) +
  # geom_point(size=1.5, shape=21, stroke=0.05) +
  ggrastr::geom_point_rast(size=1.5, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(fill = guide_legend(override.aes = list(size=2.5))) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="right",
    legend.title=element_blank()
  )

pdf(paste0(io$outdir,"/umap_by_celltype.pdf"), width=10, height=5.5)
# pdf(paste0(io$outdir,"/umap_by_celltype.pdf"))
print(p)
dev.off()
