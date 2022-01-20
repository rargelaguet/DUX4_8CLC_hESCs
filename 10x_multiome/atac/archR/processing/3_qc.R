
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/processing/3_qc.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outdir',     type="character",    help='Output directory')
p$add_argument('--min_tss_enrichment',     type="integer",    default=8,   help='Minimum TSS enrichment')
p$add_argument('--min_number_fragments',     type="integer",    default=3000,    help='Maximum number of ATAC fragments')
p$add_argument('--max_blacklist_ratio',     type="double",    default=0.05,    help='Maximum Blacklist Ratio')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- "/Users/argelagr/data/DUX4_hESCs_multiome/processed/atac/archR/sample_metadata_after_archR.txt.gz"
# args$min_tss_enrichment <- 8
# args$min_number_fragments <- 3000
# args$max_blacklist_ratio <- 0.05
# args$threads <- 2
# args$outdir <- "/Users/argelagr/data/DUX4_hESCs_multiome/results/atac/archR/test/qc"
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))

# Options
opts$chr <- paste0("chr",1:3)
opts$test <- TRUE

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata)

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads=args$threads) 

##################
## Subset ArchR ##
##################

if (opts$test) {
  cells.to.use <- split(ArchRProject$cellNames,ArchRProject$sample) %>% map(~ head(.,n=100)) %>% unlist
  ArchRProject <- ArchRProject[cells.to.use,]
}

# Subset chr for faster computations
tss.granges <- getTSS(ArchRProject)
tss.granges <- tss.granges[seqnames(tss.granges)%in%opts$chr]

#########################
## Plot TSS Enrichment ##
#########################

data_tss.dt <- opts$samples %>% map(function(i) {
  plotTSSEnrichment(
    ArchRProj = ArchRProject[ArchRProject$Sample==i,], 
    groupBy = "Sample", 
    returnDF = TRUE,
    TSS = tss.granges
  ) %>% as.data.table %>% return
}) %>% rbindlist %>% 
  setnames("group","sample") %>% 
  melt(id.vars=c("sample","x"))
fwrite(data_tss.dt, sprintf("%s/qc_TSSenrichment.txt.gz",args$outdir))

# data_tss.dt <- fread(sprintf("%s/qc_TSSenrichment.txt.gz",args$outdir)) %>% 
#   .[,.(value=mean(value)), by = c("sample","x","variable")]

to_plot_tss.dt <- data_tss.dt %>% 
  .[,.(value=mean(value)), by = c("sample","x","variable")] %>%
  .[variable=="normValue"]

p <- ggline(to_plot_tss.dt, x="x", y="value", plot_type="l") +
  facet_wrap(~sample, scales="fixed", nrow=1) +
  # scale_colour_manual(values=opts$sample.colors) +
  # scale_x_continuous(breaks=seq(-2000,2000,1000)) +
  labs(x="Distance from TSS (bp)", y="TSS enrichment (normalised)") +
  theme(
    axis.text.y = element_text(size=rel(0.5)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=rel(0.75)),
    legend.position = "none",
    legend.title = element_blank()
  )

pdf(sprintf("%s/qc_TSSenrichment.pdf",args$outdir), width=8, height=4)
print(p)
dev.off()

#####################################
## Plot Fragment size distribution ##
#####################################

# to.plot.fragmentsize <- plotFragmentSizes(ArchRProject, groupBy = "Sample", returnDF=T) %>% 
#   as.data.table %>% setnames("group","sample")
data_fragmentsize.dt <- opts$samples %>% map(function(i) {
  plotFragmentSizes(
    ArchRProj = ArchRProject[ArchRProject$Sample==i,], 
    groupBy = "Sample", 
    returnDF = TRUE
  ) %>% as.data.table %>% return
}) %>% rbindlist %>% 
  setnames("group","sample")

fwrite(data_fragmentsize.dt, sprintf("%s/qc_FragmentSizeDistribution.txt.gz",args$outdir))

# data_fragmentsize.dt <- fread(sprintf("%s/qc_FragmentSizeDistribution.txt.gz",args$outdir)) %>% 
#   .[,.(fragmentPercent=mean(fragmentPercent)), by = c("sample","fragmentSize")]

to_plot_fragmentsize.dt <- data_fragmentsize.dt %>% 
  .[,.(fragmentPercent=mean(fragmentPercent)), by = c("sample","fragmentSize")]

p <- ggline(to_plot_fragmentsize.dt, x="fragmentSize", y="fragmentPercent", plot_type="l") +
  facet_wrap(~sample, scales="fixed", nrow=1) +
  scale_x_continuous(breaks=seq(125,750,125)) +
  # scale_colour_manual(values=opts$sample.colors) +
  labs(x="Fragment Size (bp)", y="Percentage of fragments (%)") +
  theme(
    axis.text = element_text(size=rel(0.55)),
    axis.title = element_text(size=rel(0.75)),
    legend.position = "none",
    legend.title = element_blank()
  )

pdf(sprintf("%s/qc_FragmentSizeDistribution.pdf",args$outdir), width=8, height=4)
print(p)
dev.off()

##################################
## Plot histogram of QC metrics ##
##################################

to.plot <- sample_metadata %>%
  .[!is.na(nFrags_atac)] %>%
  # .[,log_nFrags:=log2(nFrags_atac)] %>%
  melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","nFrags_atac","BlacklistRatio_atac"))
  # melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","nFrags_atac"))

tmp <- data.table(
  variable = c("TSSEnrichment_atac", "nFrags_atac", "BlacklistRatio_atac"),
  value = c(args$min_tss_enrichment, args$min_number_fragments, args$max_blacklist_ratio)
)
# tmp <- data.table(
#   variable = c("TSSEnrichment_atac", "nFrags_atac"),
#   value = c(args$min_tss_enrichment, args$min_number_fragments)
# )
# p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
p <- gghistogram(to.plot, x="value", y="..density..", bins=70, fill="sample") +
  geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
  facet_wrap(~variable, scales="free") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size=rel(0.5))
  )
pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outdir), width=8, height=5)
print(p)
dev.off()


#############
## Call QC ##
#############

sample_metadata %>%
  .[,pass_atacQC:=TSSEnrichment_atac>=args$min_tss_enrichment & nFrags_atac>=args$min_number_fragments & BlacklistRatio_atac<=args$max_blacklist_ratio] %>%
  .[is.na(pass_atacQC),pass_atacQC:=FALSE]

print(sample_metadata[,mean(pass_atacQC,na.rm=T),by="sample"])
# print(sample_metadata[,mean(is.na(nFrags_atac)),by="sample"])

# Save
fwrite(sample_metadata, file.path(args$outdir,"sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")


###########################################
## Plot QC statistics after QC filtering ##
###########################################

# Barplot of the fraction of cells that pass QC for each sample

to.plot <- sample_metadata %>%
  .[,mean(pass_atacQC),by="sample"]

p <- ggbarplot(to.plot, x="sample", y="V1", fill="sample") +
  labs(x="", y="Fraction of cells that pass QC") +
  coord_cartesian(ylim=c(0,1)) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
    axis.text.y =  element_text(size=rel(0.8))
  )

# pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outdir), width=9, height=7)
pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outdir))
print(p)
dev.off()


# Boxplots of QC metrics
# to.plot <- sample_metadata %>%
#   .[pass_atacQC==TRUE & TSSEnrichment_atac<21] %>%
#   .[,log_nFrags:=log10(nFrags_atac)] %>%
#   # melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","log_nFrags","BlacklistRatio_atac"))
#   melt(id.vars=c("sample","cell","sample"), measure.vars=c("TSSEnrichment_atac","log_nFrags"))

# # Boxplots
# p <- ggboxplot(to.plot, x="sample", y="value", fill="sample", outlier.shape=NA) +
#   # scale_fill_manual(values=opts$sample.colors) +
#   facet_wrap(~variable, scales="free_y") +
#   theme(
#     # legend.position = "none",
#     legend.title = element_blank(),
#     # axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
#     axis.text.x = element_text(colour="black",size=rel(0.65)),  
#     axis.text.y = element_text(colour="black",size=rel(0.75)),  
#     axis.title.x = element_blank()
#   )
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outdir), width=10, height=6)
# print(p)
# dev.off()