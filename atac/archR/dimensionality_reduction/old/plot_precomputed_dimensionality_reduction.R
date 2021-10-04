#####################
## Define settings ##
#####################

source(here::here("settings.R"))

# I/O
io$umap <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/dimensionality_reduction/umap.txt.gz"
io$outdir <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/precomputed"); dir.create(io$outdir, showWarnings = F)

# Options
opts$motif_annotation <- "Motif_cisbp_lenient"

opts$aggregated.celltypes <- c(
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

###############
## Load UMAP ##
###############

umap.dt <- fread(io$umap.atac)

################
## Parse data ##
################

to.plot <- sample_metadata %>%
  .[,aggregated_celltype:=stringr::str_replace_all(celltype.predicted,opts$aggregated.celltypes)] %>%
  .[,aggregated_celltype:=factor(aggregated_celltype, levels=names(opts$celltype.colors))] %>%
  .[,aggregated_celltype:=stringr::str_replace_all(aggregated_celltype,"_"," ")] %>%
  data.table::merge.data.table(umap.dt,by="cell") %>%
  droplevels
  
names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$aggregated_celltype)]

stopifnot(all(unique(to.plot$aggregated_celltype) %in% names(opts$celltype.colors)))
unique(to.plot$aggregated_celltype)[!unique(to.plot$aggregated_celltype) %in% names(opts$celltype.colors)]

##########
## Plot ##
##########

p <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  # geom_point(aes(colour=aggregated_celltype), size=0.1) +
  ggrastr::geom_point_rast(aes(colour=aggregated_celltype), size=0.05) +
  scale_color_manual(values=opts$celltype.colors) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.title = element_blank(),
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/umap.pdf"), width=4.5, height=4.5)
print(p)
dev.off()


# Plot each sample separately
for (i in unique(to.plot$sample)) {
  print(i)
  
  to.plot_i <- to.plot %>% copy %>%
    .[,alpha:=1.0] %>%
    .[sample!=i,c("celltype","alpha"):=list("None",0.25)]
  
  p <- ggplot() +
    ggrastr::geom_point_rast(aes(x=umap1, y=umap2), size=0.25, color="grey", alpha=0.25, data=to.plot_i[sample!=i]) +
    ggrastr::geom_point_rast(aes(x=umap1, y=umap2, fill=aggregated_celltype), size=0.75, stroke=0.1, shape=21, alpha=1.0, data=to.plot_i[sample==i]) +
    scale_fill_manual(values=opts$celltype.colors) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position="none"
    )
  
  pdf(sprintf("%s/umap_sample%s.pdf",io$outdir,i), width=5, height=5)
  print(p)
  dev.off()
}


##########################
## Plot chromVAR scores ##
##########################

# io$archr.chromvar.se <- sprintf("%s/chromVAR_deviations_summarized_experiment_%s_correlated_peaks_archr.rds",io$archr.chromvar.dir,opts$motif_annotation)
io$chromvar.se <- sprintf("%s/results/atac/archR/chromvar/chromVAR_deviations_Motif_cisbp_lenient_archr_chip.rds",io$basedir)
atac.chromvar.se <- readRDS(io$chromvar.se) %>% .[,colnames(.) %in%sample_metadata$cell]

TFs.to.plot <- c("FOXA2","TAL1","GATA1","FOXC2")

i <- "GATA1"
for (i in TFs.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(atac.chromvar.se),
    chromvar = assay(atac.chromvar.se[i,],"z")[1,]
  ) %>% merge(umap.dt, by="cell") %>% merge(sample_metadata[,c("cell","celltype.predicted")])
  
  ggplot() +
    ggrastr::geom_point_rast(aes(x=umap1, y=umap2), size=0.15, color="grey", alpha=0.25, data=to.plot[chromvar<1]) +
    ggrastr::geom_point_rast(aes(x=umap1, y=umap2, color=chromvar), size=0.30, alpha=1.0, data=to.plot[chromvar>1]) +
    # scale_fill_gradient(low = "gray80", high = "blue") +
    scale_color_gradient2(low = "gray50", mid="gray90", high = "red") +
    # scale_fill_gradientn(colours = terrain.colors(10)) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position="none"
    )
  
  pdf(sprintf("%s/umap_sample%s.pdf",io$outdir,i), width=5, height=5)
  print(p)
  dev.off()
  
  
  
}
