library(pheatmap)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/differential")

# Options
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  # "Mixed_mesoderm",
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
  # "Erythroid1",
  # "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)
# opts$celltypes <- c("Epiblast", "Primitive_Streak", "Caudal_epiblast")

# RNA
opts$rna.min_Log2FC <- 1

# ATAC
opts$atac.matrix <- "PeakMatrix"
opts$atac.min_MeanDiff <- 0.15

##########################################
## Load differential expression results ##
##########################################

rna.file <- paste0(io$outdir,"/rna_diff_precomputed.txt.gz")

if (file.exists(rna.file)) {
  rna.dt <- fread(rna.file) %>% .[celltypeA%in%opts$celltypes & celltypeB%in%opts$celltypes]
} else {
  rna.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
    file <- sprintf("%s/%s_vs_%s.txt.gz", io$rna.differential,i,j)
    if (file.exists(file)) {
      fread(file, select = c(1,2)) %>% 
        .[,.(nhits_positive=sum(logFC>=opts$rna.min_Log2FC,na.rm=T), 
             nhits_negative=sum(logFC<=-(opts$rna.min_Log2FC),na.rm=T))] %>%
        .[,nhits:=nhits_positive+nhits_negative] %>%
        .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
        return
    }
  }) %>% rbindlist }) %>% rbindlist %>%
    .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
    .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)]
  
  rna.dt <- rna.dt %>% rbind(
    rna.dt %>% copy %>% setnames(c("nhits_negative","nhits_positive","nhits","celltypeB","celltypeA")) %>% .[,c("nhits_positive","nhits_negative","nhits","celltypeA","celltypeB")]
  )
  fwrite(rna.dt, paste0(io$outdir,"/rna_diff_precomputed.txt.gz"))
}

#####################
## Plot DE results ##
#####################

to.plot <- rna.dt
celltype.order <- to.plot %>% .[,.(median(nhits)),by="celltypeA"] %>% setorder(-V1) %>% .$celltypeA
to.plot <- to.plot %>% .[,celltypeA:=factor(celltypeA,levels=celltype.order)]

p <- ggplot(to.plot, aes(x=factor(celltypeA), y=nhits)) +
  # geom_point(aes(fill = celltypeA), shape=21, size=1) +
  geom_boxplot(aes(fill = celltypeA), alpha=0.9, outlier.shape=NA, coef=1) +
  coord_flip(ylim = c(-5,1500)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(y="Number of DE genes", x="") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

pdf(sprintf("%s/rna_number_DE_genes_boxplots.pdf",io$outdir), width=5, height=6)
print(p)
dev.off()

#############################################
## Load differential accessibility results ##
#############################################

# NOTE: NEGATIVE MEANS MORE UPREGULATED IN CELLTYPE B
acc.file <- paste0(io$outdir,"/acc_diff_precomputed.txt.gz")

if (file.exists(acc.file)) {
  acc.dt <- fread(acc.file) %>% .[celltypeA%in%opts$celltypes & celltypeB%in%opts$celltypes]
} else {
  io$archR.diff.dir <- "/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/differential/PeakMatrix"
  acc.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
    file <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$archR.diff.dir,opts$atac.matrix,i,j)
    if (file.exists(file)) {
      fread(file, select = c(1,2)) %>% 
        .[,.(nhits_positive=sum(MeanDiff>=opts$atac.min_MeanDiff,na.rm=T), 
             nhits_negative=sum(MeanDiff<=-(opts$atac.min_MeanDiff),na.rm=T))] %>%
        .[,nhits:=nhits_positive+nhits_negative] %>%
        .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
        return
    }
  }) %>% rbindlist }) %>% rbindlist %>%
    .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
    .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)]
    # .[,sig:=abs(logFC)>=opts$min.Log2FC & padj_fdr<=opts$min.FDR] %>%
  
  acc.dt <- acc.dt %>% rbind(
    acc.dt %>% copy %>% setnames(c("nhits_negative","nhits_positive","nhits","celltypeB","celltypeA")) %>% .[,c("nhits_positive","nhits_negative","nhits","celltypeA","celltypeB")]
  )
  
  fwrite(acc.dt, paste0(io$outdir,"/acc_diff_precomputed.txt.gz"))
}

#####################
## Plot DA results ##
#####################

to.plot <- acc.dt
celltype.order <- to.plot %>% .[,.(median(nhits)),by="celltypeA"] %>% setorder(-V1) %>% .$celltypeA
to.plot <- to.plot %>% .[,celltypeA:=factor(celltypeA,levels=celltype.order)]

p <- ggplot(to.plot, aes(x=factor(celltypeA), y=nhits)) +
  # geom_point(aes(fill = celltypeA), shape=21, size=1) +
  geom_boxplot(aes(fill = celltypeA), alpha=0.9, outlier.shape=NA, coef=1) +
  coord_flip(ylim = c(-5,9500)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(y="Number of DA peaks", x="") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

pdf(sprintf("%s/acc_number_DA_peaks_boxplots.pdf",io$outdir), width=5, height=6)
print(p)
dev.off()

#########################################
## Scatterplot of DE genes vs DA peaks ##
#########################################

to.plot <- merge(rna.dt, acc.dt, by=c("celltypeA","celltypeB"), suffixes=c("_rna","_acc")) %>%
  .[,transition:=sprintf("%s_to_%s",celltypeA,celltypeB)] %>%
  .[,residuals:=lm(formula=nhits_rna~nhits_acc, data=.)[["residuals"]]]

to.plot.text <- to.plot[order(-abs(residuals))] %>% head(n=25)

to.plot[nhits_acc>7000,nhits_acc:=7000]

p <- ggscatter(to.plot, x="nhits_acc", y="nhits_rna", size=1.5, fill="gray50", shape=21,
          add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  # scale_fill_manual(values=opts$celltype.colors, drop=F) +
  stat_cor(method = "pearson") +
  coord_cartesian(ylim=c(0,750)) +
  # ggrepel::geom_text_repel(aes(label=transition), size=3, max.overlaps=Inf, data=to.plot.text) +
  labs(x="Number of differentially expressed genes", y="Number of differentially accessible peaks") +
  theme(
    legend.position = "none",
    axis.text = element_text(size=rel(0.5)),
    axis.title = element_text(size=rel(0.85))
  )

pdf(sprintf("%s/number_DE_genes_vs_number_DA_peaks_scatterplot.pdf",io$outdir), width=6, height=4.5)
print(p)
dev.off()

##########################
## Boxplot of residuals ##
##########################

celltype.order <- to.plot %>% .[,.(mean=mean(residuals)),by="celltypeA"] %>% setorder(-mean) %>% .$celltypeA
to.plot <- to.plot %>% .[,celltypeA:=factor(celltypeA,levels=celltype.order)]

p <- ggplot(to.plot, aes(x=factor(celltypeA), y=residuals)) +
  # geom_point(aes(fill = celltypeA), shape=21, size=1) +
  geom_boxplot(aes(fill = celltypeA), alpha=0.9, outlier.shape=NA, coef=1) +
  coord_flip(ylim=c(-600,600)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(y="Residuals", x="") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", size=rel(0.75)),
  )

pdf(sprintf("%s/rna_vs_acc_differential_residuals_boxplots.pdf",io$outdir), width=5, height=6)
print(p)
dev.off()

##########################
## Heatmap of residuals ##
##########################

# to.plot <- to.plot %>%
#   dcast(celltypeA~celltypeB, value.var="residuals", drop=FALSE) %>%
#   matrix.please 
# 
# # Fill NAs
# for (i in rownames(to.plot)) {
#   for (j in colnames(to.plot)) {
#     if (is.na(to.plot[i,j])) to.plot[i,j] <- to.plot[j,i] 
#     if (is.na(to.plot[i,j])) to.plot[i,j] <- to.plot[j,i] 
#   }
# }
# 
# pheatmap(
#   mat = to.plot, 
#   cluster_rows  = F, cluster_cols = F, 
#   width = 6, height = 6,
#   fontsize = 6, 
#   filename = paste0(io$outdir,"/acc_vs_rna_differential_residuals_heatmap.pdf")
# )

#############################
## Heatmap of RNA DE genes ##
#############################

to.plot.rna <- rna.dt %>%
  dcast(celltypeA~celltypeB, value.var="nhits", drop=FALSE) %>%
  matrix.please 

# Fill NAs
for (i in rownames(to.plot.rna)) {
  for (j in colnames(to.plot.rna)) {
    if (is.na(to.plot.rna[i,j])) to.plot.rna[i,j] <- to.plot.rna[j,i] 
  }
}

to.plot.rna <- to.plot.rna[rev(rownames(to.plot.rna)), rev(rownames(to.plot.rna))] 
to.plot.rna[to.plot.rna>1500] <- 1500

# Plot
pheatmap(
  mat = to.plot.rna, 
  cluster_rows  = F, cluster_cols = F, 
  width = 5, height = 6.5,
  fontsize = 7,
  legend = F,
  filename = paste0(io$outdir,"/rna_nDE_genes_heatmap.pdf")
)

# Plot legend
pheatmap(
  mat = to.plot.rna, 
  legend = T,
  filename = paste0(io$outdir,"/rna_nDE_peaks_heatmap_legend.pdf")
)

##############################
## Heatmap of ATAC DA genes ##
##############################

to.plot.acc <- acc.dt %>%
  dcast(celltypeA~celltypeB, value.var="nhits", drop=FALSE) %>%
  matrix.please 

# Fill NAs
for (i in rownames(to.plot.acc)) {
  for (j in colnames(to.plot.acc)) {
    if (is.na(to.plot.acc[i,j])) to.plot.acc[i,j] <- to.plot.acc[j,i] 
  }
}

to.plot.acc <- to.plot.acc[rev(rownames(to.plot.acc)), rev(rownames(to.plot.acc))] 
to.plot.acc[to.plot.acc>8000] <- 8000

# Plot
pheatmap(
  mat = to.plot.acc, 
  # color = colorRampPalette(c("gray80", "purple"))(100),
  cluster_rows  = F, cluster_cols = F, 
  width = 5, height = 6.5,
  fontsize = 7, 
  legend = F,
  filename = paste0(io$outdir,"/acc_nDA_peaks_heatmap.pdf")
)

# Plot legend
pheatmap(
  mat = to.plot.acc, 
  legend = T,
  filename = paste0(io$outdir,"/acc_nDA_peaks_heatmap_legend.pdf")
)

#########################################################
## number of DA peaks versus fraction of open DA peaks ##
#########################################################

to.plot <- acc.dt %>% copy %>% 
  .[,fraction_upregulated_hits:=nhits_positive/nhits] %>%
  .[,.(fraction_upregulated_hits=mean(fraction_upregulated_hits), nhits=mean(nhits)), by="celltypeA"]

ggplot(to.plot, aes(x=fraction_upregulated_hits, y=nhits)) +
  geom_point(aes(fill = celltypeA), shape=21, size=3) +
  coord_cartesian(xlim=c(0.2,0.8)) +
  # geom_hline(yintercept=min(to.plot$nhits), linetype="dashed", size=0.5) +
  geom_vline(xintercept=0.5, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(x="Average fraction of upregulated peaks", y="Average number of DA peaks") +
  theme(
    legend.position = "none",
    # axis.title.y = element_blank(),
    axis.text = element_text(color="black", size=rel(0.75)),
    axis.title = element_text(color="black", size=rel(1)),
  )

################################################################################
## fraction of upregulated RNA genes versus  fraction of upregulated DA peaks ##
################################################################################

foo <- rna.dt %>% copy %>% 
  .[,fraction_upregulated_hits:=nhits_positive/nhits] %>%
  .[,.(fraction_upregulated_hits=mean(fraction_upregulated_hits), nhits=mean(nhits)), by="celltypeA"] %>%
  .[,class:="RNA"]

bar <- acc.dt %>% copy %>% 
  .[,fraction_upregulated_hits:=nhits_positive/nhits] %>%
  .[,.(fraction_upregulated_hits=mean(fraction_upregulated_hits), nhits=mean(nhits)), by="celltypeA"] %>%
  .[,class:="ATAC"]

to.plot <- rbind(foo,bar) %>% dcast(celltypeA~class, value.var=c("fraction_upregulated_hits","nhits")) %>%
  .[,nhits_scaled:=6*minmax.normalisation(nhits_ATAC)]

p <- ggplot(to.plot, aes(x=fraction_upregulated_hits_ATAC, y=fraction_upregulated_hits_RNA)) +
  geom_point(aes(fill=celltypeA, size=nhits_scaled), shape=21) +
  # coord_cartesian(xlim=c(0.25,0.75), ylim=c(0.25,0.75)) +
  # geom_hline(yintercept=min(to.plot$nhits), linetype="dashed", size=0.5) +
  geom_vline(xintercept=0.5, linetype="dashed", size=0.5) +
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors) +
  scale_size_continuous(range = c(3,10)) +
  ggrepel::geom_text_repel(aes(label=celltypeA), size=3, max.overlaps=Inf, data=to.plot[abs(fraction_upregulated_hits_ATAC)>0.6 | abs(fraction_upregulated_hits_RNA)>0.6]) +
  theme_classic() +
  labs(x="Average fraction of upregulated peaks", y="Average fraction of upregulated genes") +
  theme(
    legend.position = "none",
    # axis.title.y = element_blank(),
    axis.text = element_text(color="black", size=rel(0.85)),
    axis.title = element_text(color="black", size=rel(1)),
  )

pdf(sprintf("%s/fraction_DEgenes_vs_fraction_DApeaks.pdf",io$outdir), width=6, height=5)
print(p)
dev.off()

################
## PAGA graph ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/load_paga_graph.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/load_paga_graph.R")
} else {
  stop("Computer not recognised")
}


# Define colors
rna.col.seq <- round(seq(0,1,0.1), 2)
atac.col.seq <- round(seq(0,1,0.1), 2)
rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))
atac.colors <- colorRampPalette(c("gray92", "darkblue"))(length(atac.col.seq)) 


p.net <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  node.size = 0,
  edge.size = 0.15,
  edge.color = "grey",
  label = FALSE,
  label.size = 2.3
)


to.plot.values <- to.plot[,c("celltypeA","nhits_RNA")] %>% .[,nhits_RNA:=minmax.normalisation(nhits_RNA)] %>% matrix.please %>% .[opts$celltypes,]
colors <- round(to.plot.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
p.rna <- p.net + geom_text(label = "\u25D0", aes(x=x, y=y), color=colors, size=20, family = "Arial Unicode MS",
                           data = p.net$data[,c("x","y")] %>% dplyr::mutate(expr=colors)) +
  scale_colour_manual(values=colors) + 
  labs(title="Number of DE genes") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

to.plot.values <- to.plot[,c("celltypeA","nhits_ATAC")] %>% .[,nhits_ATAC:=minmax.normalisation(nhits_ATAC)] %>% matrix.please %>% .[opts$celltypes,]
colors <- round(to.plot.values,1) %>% map(~ atac.colors[which(atac.col.seq == .)]) %>% unlist
p.acc <- p.net + geom_text(label = "\u25D1", aes(x=x, y=y), color=colors, size=20, family = "Arial Unicode MS",
                           data = p.net$data[,c("x","y")] %>% dplyr::mutate(acc=colors)) +
  scale_fill_manual(values=colors) + 
  labs(title="Number of DA peaks") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

p <- cowplot::plot_grid(plotlist=list(p.rna,p.acc), nrow=1, scale = 0.9)

outfile <- sprintf("%s/paga_DE_genes_DA_peaks.png",io$outdir)  
png(outfile, width = 600, height = 400)
print(p)
dev.off()

#################################################
## Number of DE genes versus cell type numbers ##
#################################################

# opts$celltypes  <- opts$celltypes %>% head(n=3)

foo <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s.txt.gz", io$rna.differential,i,j)
  if (file.exists(file)) {
    fread(file, select = c(1,2,3,6,7)) %>% 
      .[!is.na(logFC),.(nhits=sum(logFC>=1,na.rm=T), ncells=min(groupA_N,groupB_N))] %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      return
  } }) %>% rbindlist }) %>% rbindlist %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)]

foo <- foo %>% rbind(
  foo %>% copy %>% setnames(c("nhits","ncells","celltypeB","celltypeA")) %>% .[,c("nhits","ncells","celltypeA","celltypeB")]
)

bar <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s.txt.gz", io$rna.differential,i,j)
  if (file.exists(file)) {
    fread(file, select = c(1,2,3,6,7)) %>% 
      .[!is.na(logFC),.(nhits=sum(p.value<=0.01,na.rm=T), ncells=min(groupA_N,groupB_N))] %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      return
  } }) %>% rbindlist }) %>% rbindlist %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)]

bar <- bar %>% rbind(
  bar %>% copy %>% setnames(c("nhits","ncells","celltypeB","celltypeA")) %>% .[,c("nhits","ncells","celltypeA","celltypeB")]
)


ggscatter(foo, x="ncells", y="nhits", size=1.5,
          add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") 


