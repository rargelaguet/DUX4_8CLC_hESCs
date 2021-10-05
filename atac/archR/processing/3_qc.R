########################
## Load ArchR project ##
########################

here::i_am("atac/archR/processing/3_qc.R")
source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- file.path(io$basedir,"results/atac/archR/qc")

# Options
opts$samples <- c(
  "HNES1_DUX4_overexpression_L001",
  "HNES1_wildtype_L001"
)

# QC thresholds
opts$min.TSSEnrichment <- 8
opts$min.log_nFrags <- 11.5
opts$max.BlacklistRatio <- 0.05
opts$chr <- paste0("chr",1:3)
opts$test <- FALSE


########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata)

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads = 4) 

# Subset
# ArchRProject <- ArchRProject[sample_metadata[pass_atacQC==TRUE,cell]

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

to.plot.tss <- opts$samples %>% map(function(i) {
  plotTSSEnrichment(
    ArchRProj = ArchRProject[ArchRProject$Sample==i,], 
    groupBy = "Sample", 
    returnDF = TRUE,
    TSS = tss.granges
  ) %>% as.data.table %>% return
}) %>% rbindlist %>% 
  setnames("group","sample") %>% 
  melt(id.vars=c("sample","x"))
fwrite(to.plot.tss, sprintf("%s/qc_TSSenrichment.txt.gz",io$outdir))

to.plot.tss <- fread(sprintf("%s/qc_TSSenrichment.txt.gz",io$outdir)) %>% 
  .[,.(value=mean(value)), by = c("sample","x","variable")]

p <- ggline(to.plot.tss[variable=="normValue"], x="x", y="value", plot_type="l") +
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

pdf(sprintf("%s/qc_TSSenrichment.pdf",io$outdir), width=8, height=4)
print(p)
dev.off()

#####################################
## Plot Fragment size distribution ##
#####################################

# to.plot.fragmentsize <- plotFragmentSizes(ArchRProject, groupBy = "Sample", returnDF=T) %>% 
#   as.data.table %>% setnames("group","sample")
to.plot.fragmentsize <- opts$samples %>% map(function(i) {
  plotFragmentSizes(
    ArchRProj = ArchRProject[ArchRProject$Sample==i,], 
    groupBy = "Sample", 
    returnDF = TRUE
  ) %>% as.data.table %>% return
}) %>% rbindlist %>% 
  setnames("group","sample")

fwrite(to.plot.fragmentsize, sprintf("%s/qc_FragmentSizeDistribution.txt.gz",io$outdir))

# to.plot <- to.plot.fragmentsize %>% .[variable=="fragmentPercent"] %>% .[,fragmentSize:=1:.N,by="sample"] %>% setnames("value","fragmentPercent") %>% .[,variable:=NULL]

to.plot.fragmentsize <- fread(sprintf("%s/qc_FragmentSizeDistribution.txt.gz",io$outdir)) %>% 
  .[,.(fragmentPercent=mean(fragmentPercent)), by = c("sample","fragmentSize")]

# to.plot.fragmentsize2 <- to.plot.fragmentsize %>% dcast(sample~variable, value.var="value")
p <- ggline(to.plot.fragmentsize, x="fragmentSize", y="fragmentPercent", plot_type="l") +
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

pdf(sprintf("%s/qc_FragmentSizeDistribution.pdf",io$outdir), width=8, height=4)
print(p)
dev.off()

##################################
## Plot histogram of QC metrics ##
##################################

to.plot <- sample_metadata %>%
  .[!is.na(nFrags_atac)] %>%
  .[,log_nFrags:=log2(nFrags_atac)] %>%
  # melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","log_nFrags","BlacklistRatio_atac"))
  melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","log_nFrags"))

# tmp <- data.table(
#   variable = c("TSSEnrichment_atac", "log_nFrags", "BlacklistRatio_atac"),
#   value = c(opts$min.TSSEnrichment, opts$min.log_nFrags, opts$max.BlacklistRatio)
# )
tmp <- data.table(
  variable = c("TSSEnrichment_atac", "log_nFrags"),
  value = c(opts$min.TSSEnrichment, opts$min.log_nFrags)
)
# p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
p <- gghistogram(to.plot, x="value", y="..density..", bins=70, fill="sample") +
  geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
  facet_wrap(~variable, scales="free") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size=rel(0.5))
  )
pdf(sprintf("%s/qc_metrics_histogram.pdf",io$outdir), width=8, height=5)
print(p)
dev.off()




#############
## Call QC ##
#############

sample_metadata %>%
  .[,pass_atacQC:=TSSEnrichment_atac>=opts$min.TSSEnrichment & log2(nFrags_atac)>=opts$min.log_nFrags & BlacklistRatio_atac<=opts$max.BlacklistRatio] %>%
  .[is.na(pass_atacQC),pass_atacQC:=FALSE]

print(sample_metadata[,mean(pass_atacQC,na.rm=T),by="sample"])
# print(sample_metadata[,mean(is.na(nFrags_atac)),by="sample"])

# Save
outfile <- paste0(io$outdir,"/sample_metadata_after_qc.txt.gz")
fwrite(sample_metadata, outfile, quote=F, na="NA", sep="\t")


###########################################
## Plot QC statistics after QC filtering ##
###########################################

# Barplot of the fraction of cells that pass QC for each sample

to.plot <- sample_metadata %>%
  .[,mean(pass_atacQC),by="sample"]

p <- ggbarplot(to.plot, x="sample", y="V1", fill="gray70") +
    labs(x="", y="Fraction of cells that pass QC") +
  coord_cartesian(ylim=c(0,1)) +
    theme(
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
    )

# pdf(sprintf("%s/qc_metrics_barplot.pdf",io$outdir), width=9, height=7)
pdf(sprintf("%s/qc_metrics_barplot.pdf",io$outdir))
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
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",io$outdir), width=10, height=6)
# print(p)
# dev.off()