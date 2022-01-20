
#####################
## Define settings ##
#####################

# Load default settings
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
io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/peak_calling/motifmatchr"); dir.create(io$outdir, showWarnings = F)

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(io$archR.peak.metadata) %>%
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)]

#################################
## Load motifmatcher Positions ##
#################################

motifmatcher_positions.se <- readRDS(sprintf("%s/Annotations/Motif_cisbp-Positions-In-Peaks.rds",io$archR.directory))

# Rename TFs
names(motifmatcher_positions.se) <- names(motifmatcher_positions.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
motifmatcher.se <- motifmatcher.se[!duplicated(names(motifmatcher_positions.se))]

motifmatcher_positions.dt <- unlist(motifmatcher_positions.se) %>% as.data.table %>% 
  .[,motif:=as.factor(names(unlist(motifmatcher_positions.se)))] %>% 
  setnames("seqnames","chr") %>% 
  .[,location:=sprintf("%s:%s-%s",chr,start,end)]

#################################
## Load motifmatcher Matches ##
#################################

motifmatcher.se <- readRDS(sprintf("%s/Annotations/Motif_cisbp-Matches-In-Peaks.rds",io$archR.directory))

# Rename TFs
colnames(motifmatcher.se) <- colnames(motifmatcher.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
motifmatcher.se <- motifmatcher.se[,!duplicated(colnames(motifmatcher.se))]

# Rename peaks
tmp <- rowRanges(motifmatcher.se)
rownames(motifmatcher.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))

# Subset pekas
# motifmatcher.se <- motifmatcher.se[unique(cor_dt$peak),]

###############################################################
## Plot distribution of motifmatchr scores across all motifs ##
###############################################################

to.plot <- motifmatcher_positions.dt[chr%in%c("chr1","chr2")]

p <- ggdensity(to.plot, x="score", fill="strand") +
  labs(x="motifmarchr scores") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.8))
  )
pdf(paste0(io$outdir,"/histogram_motifmatchr_scores_allmotifs.pdf"), width = 6, height = 5)
print(p)
dev.off()

#######################################################
## Plot distribution of motifmatchr scores per motif ##
#######################################################

max.score <- max(motifmatcher_positions.dt$score)
min.score <- min(motifmatcher_positions.dt$score)

for (i in names(motifmatcher_positions.se)) {
  
  print(i)
  to.plot <- motifmatcher_positions.dt[motif==i]
  
  p <- ggdensity(to.plot, x="score", fill="strand") +
    coord_cartesian(xlim=c(min.score-0.1,max.score+0.1)) +
    labs(x=sprintf("%s motifmatchr scores (N=%d)",i,nrow(to.plot))) +
    theme(
      legend.title = element_blank(),
      axis.text = element_text(size=rel(0.8))
    )
  
  pdf(sprintf("%s/histograms_per_TF/%s_histogram_motifmatchr_scores.pdf",io$outdir,i), width = 6, height = 5)
  print(p)
  dev.off()
}


##########################################
## Plot number of binding events per TF ##
##########################################

to.plot <- data.table(
  TF = colnames(motifmatcher.se),
  N = colSums(assay(motifmatcher.se))
) %>% setorder(-N) %>% .[,TF:=factor(TF)] %>% 
  .[,log_N:=log10(N)]

p <- ggboxplot(to.plot, x="TF", y="log_N") +
  labs(x="Number of motif matches per TF") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.8))
  )
pdf(paste0(io$outdir,"/boxplot_number_of_motifs_per_TF.pdf"), width = 6, height = 5)
print(p)
dev.off()

p <- gghistogram(to.plot, x="log_N", fill="gray70") +
  labs(x="Number of motif matches per TF (log10)") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.8))
  )
pdf(paste0(io$outdir,"/histogram_number_of_motifs_per_TF.pdf"), width = 6, height = 5)
print(p)
dev.off()

##################################################
## Plot distribution of binding events per peak ##
##################################################

to.plot <- data.table(
  peak = rownames(motifmatcher.se),
  N = rowSums(assay(motifmatcher.se))
) %>% merge(peak_metadata.dt[,c("peak","peakType","GC","score")], by=c("peak"))

p <- ggdensity(to.plot, x="N", fill="peakType") +
  labs(x="Number of TF motifs per peak") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.8))
  )

pdf(paste0(io$outdir,"/histogram_number_of_motifs_per_peak.pdf"), width = 6, height = 5)
print(p)
dev.off()

##########
## Save ##
##########

to.save <- data.table(
  motif = colnames(motifmatcher.se),
  N = colSums(assay(motifmatcher.se))
) %>% 
  merge(unique(motifmatcher_positions.dt[,c("width","motif")]),by="motif") %>%
  setorder(-N)
fwrite(to.save, paste0(io$outdir,"/number_motifmatches_per_TF.txt.gz"), sep="\t", quote=F)

to.save <- data.table(
  peak = rownames(motifmatcher.se),
  N = rowSums(assay(motifmatcher.se))
) %>% merge(peak_metadata.dt[,c("peak","peakType","GC","score")], by=c("peak")) %>%
  setnames("score","peakScore")
fwrite(to.save, paste0(io$outdir,"/number_motifmatches_per_peak.txt.gz"), sep="\t", quote=F)

#############
## Explore ##
#############


to.plot <- seq(1,20,by=0.5) %>% 
  map( function(x) {
    motifmatcher_positions.dt[score>=x,.N,by="motif"] %>% .[,min_score:=x]
  }) %>% rbindlist


p <- ggscatter(to.plot[,.(N=sum(N)),by="min_score"], x="min_score", y="N") +
  yscale("log10", .format = TRUE) +
  labs(x="Minimum score", y="Total number of motif matches") +
  theme(
    axis.text = element_text(size=rel(0.5))
  )

pdf(paste0(io$outdir,"/scatterplot_score_vs_number_motifmatches.pdf"), width = 6, height = 5)
print(p)
dev.off()

# Which TFs are lost with a minimum_score
# target.motifs <- to.plot[min_score==5 & N>=10,motif] %>% as.character
# all.motifs <- names(motifmatcher_positions.se)
# all.motifs[!all.motifs%in%target.motifs]

##############################
## Inspect individual peaks ##
##############################

which(assay(motifmatcher.se["chr14:38320846-38321446",])[1,])

foo <- motifmatcher_positions.dt %>% 
  .[chr=="chr14"] %>% .[start>38320846 & end<38321446]

foo <- motifmatcher_positions.dt %>% .[chr=="chr1" & strand=="+"] %>%
  setorder(-score)
