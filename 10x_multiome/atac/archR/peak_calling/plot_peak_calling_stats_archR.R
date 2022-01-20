#####################
## Define settings ##
#####################

# load default setings
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
# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/peak_calling")

# Options
opts$min.score <- 20

##########################
## Load peak annotation ##
##########################

# peakSet.gr <- readRDS(io$archR.peakSet.granges)
# io$archR.peak.metadata <- paste0(io$archR.directory,"/PeakCalls/all_peaks/peak_metadata.tsv.gz")
peakSet.dt <- fread(io$archR.peak.metadata) %>%
  .[,peak:=sprintf("%s_%s_%s",chr,start,end)] %>%
  .[,peakType:=factor(peakType,levels=c("Promoter","Intronic","Exonic","Distal"))]

# Filter
# peakSet.gr.filt <- peakSet.gr[peakSet.gr$score>opts$min.score,]

# Load peak stats
io$archR.peak.stats <- paste0(io$basedir,"/results/atac/archR/peak_calling/peak_stats_binarised.txt.gz")
peakStats.dt <- fread(io$archR.peak.stats)

###############################
## Plot mean versus variance ##
###############################

to.plot <- peakSet.dt[chr=="chr1",c("peak","peakType")] %>% 
  merge(peakStats.dt,by="peak")

p <- ggscatter(to.plot, x="var_singlecell", y="var_pseudobulk", size=1) +
  labs(x="Variance (single-cell)", y="Variance (pseudobulk)") +
  theme(
    axis.text = element_text(size=rel(0.7))
  )

pdf(paste0(io$outdir,"/foo.pdf"), width = 6, height = 5)
print(p)
dev.off()

###############################################
## Plot minimum score versus number of peaks ##
###############################################

to.plot <- seq(1,250,by=1) %>% map( function(x) {
    data.table(N=peakSet.dt[score>=x,.N], min_score=x)
  }) %>% rbindlist

p <- ggscatter(to.plot, x="min_score", y="N", size=1) +
  geom_vline(xintercept=20, linetype="dashed", color="orange") +
  # yscale("log10", .format = TRUE) +
  labs(x="Minimum peak score", y="Total number of peaks") +
  theme(
    axis.text = element_text(size=rel(0.7))
  )

pdf(paste0(io$outdir,"/scatterplot_score_vs_number_peaks.pdf"), width = 6, height = 5)
print(p)
dev.off()

#####################################
## Plot fraction of peaks per type ##
#####################################

to.plot <- peakSet.dt %>%
  .[score>=opts$min.score] %>%
  .[,.N,by="peakType"]

p <- ggpie(to.plot, x="N", label="peakType", fill="peakType") +
  theme(
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/pieplot_peakType.pdf"), width = 5, height = 4)
print(p)
dev.off()

#############################################
## Plot average accessibility per peakType ##
#############################################

to.plot <- peakSet.dt %>% 
  .[score>=opts$min.score] %>%
  merge(peakStats.dt,by="peak")

p <- ggboxplot(to.plot, x="peakType", y="mean_singlecell", fill="peakType", outlier.shape = NA) +
  coord_cartesian(ylim=c(0,1.0)) +
  labs(x="", y="Average chromatin accessibility") +
  theme(
    axis.text = element_text(size=rel(0.75)),
    legend.position = "none"
  )

pdf(paste0(io$outdir,"/accessibility_boxplot_peakType.pdf"), width = 7, height = 5)
print(p)
dev.off()