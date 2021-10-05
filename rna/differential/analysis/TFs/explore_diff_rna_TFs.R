
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- sprintf("%s/TFs",io$rna.differential)


# Options
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
)# %>% head(n=3)

opts$min.Log2FC <- 1
opts$min.FDR <- 1e-2


#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_vs_%s.txt.gz", io$rna.differential,i,j)
  if (file.exists(file)) {
    fread(file, select = c(1,2,4)) %>% 
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>%
    return
  }
}) %>% rbindlist }) %>% rbindlist %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,sig:=abs(logFC)>=opts$min.Log2FC & padj_fdr<=opts$min.FDR] %>%
  .[is.na(sig),sig:=FALSE]

# dt[,mean(sig),by=c("celltypeA","celltypeB")] %>% View

#################
## Subset data ##
#################

# Select TFs
# opts$motif_annotation <- "Motif_JASPAR2020_human"
opts$motif_annotation <- "Motif_cisbp"

source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")

dt.filt <- dt %>%
  .[,gene:=toupper(gene)] %>% 
  .[gene%in%unique(motif2gene.dt$gene)]

# How many TFs are dynamic during gastrulation?
# Answer: 40% (it is a conservative estimate, because some differential expression changes might be undetected)
dt.filt %>% .[,.(mean(sig,na.rm=T)), by=c("gene")] %>% .[,mean(V1>0)]

#################################################
## Plot fraction of significant tests per gene ##
#################################################

to.plot <- dt.filt %>%
  .[,.(mean(sig,na.rm=T)), by=c("gene")] %>%
  setorder(-V1) %>%
  head(n=50) %>%
  .[,gene:=factor(gene,levels=rev(gene))]

p <- ggplot(to.plot, aes_string(x="gene", y="V1"), fill="gray70") +
  geom_point(size=2) +
  geom_segment(aes_string(xend="gene"), size=0.5, yend=0) +
  coord_flip() +
  labs(x="",y="Fraction of significant tests") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(sprintf("%s/rna_fraction_significant_tests_per_gene.pdf",io$outdir), width = 5, height = 9)
print(p)
dev.off()

#####################################################
## Plot fraction of significant tests per celltype ##
#####################################################

to.plot <- dt.filt %>%
  .[,.(mean(sig,na.rm=T)), by=c("celltypeA")] %>%
  setorder(-V1) %>% .[,celltypeA:=factor(celltypeA,levels=rev(celltypeA))]

p <- ggplot(to.plot, aes_string(x="celltypeA", y="V1"), fill="gray70") +
  geom_point(size=2) +
  geom_segment(aes_string(xend="celltypeA"), size=0.5, yend=0) +
  coord_flip() +
  labs(x="",y="Fraction of significant tests") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(sprintf("%s/rna_fraction_significant_tests_per_celltype.pdf",io$outdir), width = 5, height = 7)
print(p)
dev.off()


#############
## Heatmap ##
#############

to.plot <- dt %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltypeA","celltypeB")] %>%
  dcast(celltypeA~celltypeB, value.var="N", drop=FALSE) %>%
  matrix.please 

for (i in rownames(to.plot)) {
  for (j in colnames(to.plot)) {
    if (is.na(to.plot[i,j])) to.plot[i,j] <- to.plot[j,i] 
  }
}
diag(to.plot) <- 0
mean(is.na(to.plot))

pdf(sprintf("%s/heatmap_differential.pdf",io$outdir), width=8, height=6)
pheatmap(to.plot, cluster_rows  = F, cluster_cols = F, fontsize = 6)
dev.off()
