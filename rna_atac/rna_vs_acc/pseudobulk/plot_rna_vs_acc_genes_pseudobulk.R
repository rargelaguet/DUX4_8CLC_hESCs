#####################
## define settings ##
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
# io$pseudobulk.GeneScoreMatrix <- sprintf("%s/pseudobulk/pseudobulk_%s_summarized_experiment.rds",io$archR.directory,opts$assay)
io$archR.pseudobulk.GeneMatrix.se <- sprintf("%s/pseudobulk/pseudobulk_GeneScoreMatrix_nodistal_summarized_experiment.rds",io$archR.directory)
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk")
dir.create(paste0(io$outdir,"/per_gene"), showWarnings = F)
dir.create(paste0(io$outdir,"/per_celltype"), showWarnings = F)


# Options
opts$celltypes <- c(
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

opts$assay <- "GeneScoreMatrix_nodistal"
# opts$assay <- "GeneScoreMatrix"


#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$rna.atlas.marker_genes) %>%
  .[celltype%in%opts$celltypes]

##########################
## Load pseudobulk ATAC ##
##########################

# Load SummarizedExperiment
atac.gene_scores.se <- readRDS(io$archR.pseudobulk.GeneMatrix.se)[,opts$celltypes]

# Rename features
rownames(atac.gene_scores.se) <- rowData(atac.gene_scores.se)$name

# Prepare data.table
atac.dt <- assay(atac.gene_scores.se,opts$assay) %>% as.data.table %>% 
  .[,gene:=rowData(atac.gene_scores.se)[,"name"]] %>%
  melt(id.vars=c("gene"), variable.name="celltype") %>%
  # .[,celltype:=stringr::str_replace_all(celltype,opts$celltypes.to.rename)] %>%
  .[,.(atac=mean(value)), by=c("gene","celltype")]

#########################
## Load pseudobulk RNA ##
#########################

# Load SingleCellExperiment
sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

# Prepare data.table
rna.dt <- as.matrix(logcounts(sce.pseudobulk)) %>% t %>% 
  as.data.table(keep.rownames="celltype") %>% 
  melt(id.vars="celltype", value.name="expr", variable.name="gene")


###########
## Merge ##
###########

# Sanity checks
foo <- unique(atac.dt$celltype)
bar <- unique(rna.dt$celltype)
length(intersect(foo,bar))


foo <- unique(atac.dt$gene)
bar <- unique(rna.dt$gene)
length(intersect(foo,bar))


# Merge
dt <- merge(rna.dt, atac.dt, by = c("gene","celltype"))


##################
## Filter genes ##
##################

# dt <- dt[grep("Rik",gene,invert = T)]
# dt <- dt[gene!="Xist"]
# dt <- dt[!grepl("mt-",gene)]
# dt <- dt[!grepl("Rps|Rpl",gene)]
# dt <- dt[!grepl("Rik",gene)]
# dt <- dt[!gene%in%fread(io$gene_metadata)[chr=="chrY",symbol]]


###############################
## Scatterplot per cell type ##
###############################

opts$min.atac <- 0
opts$max.atac <- 5

opts$min.expr <- min(rna.dt$expr)
opts$max.expr <- 15

dt[expr>opts$max.expr,expr:=opts$max.expr]
dt[atac>opts$max.atac,atac:=opts$max.atac]


for (i in unique(dt$celltype)) {
  
  to.plot <- dt[celltype==i]# %>% head(n=1000)
  
  p <- ggplot(to.plot, aes(x=expr, y=atac)) +
    ggrastr::geom_point_rast(color="black", size=0.5) +
    geom_smooth(method="lm") + stat_cor(method = "pearson") +
    coord_cartesian(xlim=c(opts$min.expr,opts$max.expr+0.25), ylim=c(opts$min.atac,opts$max.atac+0.25)) +
    ggrepel::geom_text_repel(aes(label=gene), size=3, data=to.plot[expr>4 & atac>3]) +
    labs(x="RNA expression", y="Gene accessibility", title=i) +
    theme_classic() + 
    theme(
      plot.title = element_text(hjust=0.5),
      axis.text = element_text(color="black"),
      legend.position = "none"
    )
  
  pdf(sprintf("%s/per_celltype/%s_acc_vs_rna_pseudobulk.pdf",io$outdir,i), width = 8, height = 6)
  print(p)
  dev.off()
}


##########################
## Scatterplot per gene ##
##########################

genes.to.plot <- dt %>%
  .[,.(var_expr=var(expr), var_atac=var(atac)),by="gene"] %>% 
  setorder(-var_expr) %>% 
  head(n=100) %>% .$gene


genes.to.plot <- c("Ubb", "Top2a", "Actg1", "Gapdh")

for (i in genes.to.plot) {
  outfile <- sprintf("%s/per_gene/individual_genes/%s_rna_vs_acc_pseudobulk.pdf",io$outdir,i)
  if (!file.exists(outfile)) {
    to.plot <- dt[gene==i]
    
    p <- ggplot(to.plot, aes(x=atac, y=expr)) +
      ggrastr::geom_point_rast(aes(fill=celltype), color="black", shape=21, size=3, stroke=0.5) +
      geom_smooth(method="lm", color="black", fill="black", alpha=0.25) + stat_cor(method = "pearson") +
      scale_fill_manual(values=opts$celltype.colors) +
      # coord_cartesian(xlim=c(opts$min.expr,opts$max.expr+0.01), ylim=c(opts$min.acc,opts$max.acc+0.01)) +
      # ggrepel::geom_text_repel(aes(label=celltype), size=3, data=to.plot[expr>2 | atac>2]) +
      labs(x="Gene accessibility", y="RNA expression", title=i) +
      theme_classic() + 
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(color="black"),
        legend.position = "none"
      )
    
    pdf(outfile, width = 6, height = 5)
    # png(sprintf("%s/per_gene/%s_acc_vs_rna_pseudobulk.png",io$outdir,i), width = 500, height = 400)
    print(p)
    dev.off()
  }
}


##########################
## Correlation analysis ##
##########################

cor.dt <- dt %>% copy %>%
  .[,c("atac","expr"):=list(atac + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>%
  .[, .(V1 = unlist(cor.test(atac, expr)[c("estimate", "p.value")])), by = c("gene")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  .[, sig := padj_fdr<=0.10] %>% 
  setorder(padj_fdr, na.last = T)

# Save
fwrite(cor.dt, paste0(io$outdir,"/per_gene/cor_rna_vs_acc_pseudobulk.txt.gz"), sep="\t", quote=F)


# Plot

to.plot <- cor.dt

negative_hits <- to.plot[sig==TRUE & r<0,gene]
positive_hits <- to.plot[sig==TRUE & r>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$r), na.rm=T)
ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)

p <- ggplot(to.plot, aes(x=r, y=-log10(padj_fdr+1e-100))) +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.25) +
  ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & r>0],n=50),
      aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3, max.overlaps=100) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & r<0],n=10),
      aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3, max.overlaps=100) +
  scale_color_manual(values=c("black","red")) +
  scale_size_manual(values=c(0.75,1.25)) +
  scale_x_continuous(limits=c(-xlim-0.2,xlim+0.2)) +
  scale_y_continuous(limits=c(0,ylim+6)) +
  annotate("text", x=0, y=ylim+6, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.15, y=ylim+6, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.15, y=ylim+6, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

# pdf(sprintf("%s/volcano_plots/volcano_pearson_correlation.pdf",io$outdir), width = 9, height = 6)
png(sprintf("%s/per_gene/volcano_pearson_correlation.png",io$outdir), width = 800, height = 500)
print(p)
dev.off()
