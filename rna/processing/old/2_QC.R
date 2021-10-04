suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(argparse))

#####################
## Define arguments ##
#####################

p <- ArgumentParser(description='')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--nFeature_RNA',       type="integer",                    help='Minimum number of expressed genes')
p$add_argument('--nCount_RNA',       type="integer",                    help='Minimum number of reads')
p$add_argument('--mitochondrial_percent_RNA',       type="integer",                    help='Maximum percentage of mitochondrial reads')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE))


#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
    source("/Users/ricard//gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
    source("/homes/ricard//gastrulation_multiome_10x/settings.R")
}

## START TEST ##
# args <- list()
args$outputdir <- paste0(io$basedir,"/results/rna_soupX/qc")
args$samples <- opts$samples #he qcad(opts$samples,n=2)
args$nFeature_RNA <- 2500
# args$log_nCount_RNA <- log2(3000)
args$mitochondrial_percent_RNA <- 50
args$ribosomal_percent_RNA <- 25
## END TEST ##

# Sanity checks
stopifnot(args$samples%in%opts$samples)

###############
## Load data ##
###############

io$metadata <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/processed/rna_soupX/metadata.txt.gz"
metadata <- fread(io$metadata) %>% 
    .[sample%in%args$samples] %>%
    # .[,pass_rnaQC:=nFeature_RNA>args$nFeature_RNA & nCount_RNA>2**args$log_nCount_RNA & mitochondrial_percent_RNA<args$mitochondrial_percent_RNA]
    .[,pass_rnaQC:=nFeature_RNA>args$nFeature_RNA & mitochondrial_percent_RNA<args$mitochondrial_percent_RNA & ribosomal_percent_RNA<args$ribosomal_percent_RNA]

#####################
## Plot QC metrics ##
#####################

to.plot <- metadata %>% copy %>%
    # .[,log_nCount_RNA:=log2(nCount_RNA)] %>%
    melt(id.vars=c("sample","cell"), measure.vars=c("nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA"))

## Box plot 

p <- ggboxplot(to.plot, x="sample", y="value") +
    facet_wrap(~variable, scales="free_y", nrow=1) +
    theme(
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
        axis.title.x = element_blank()
    )

pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir), width=12, height=6)
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir))
print(p)
dev.off()

## histogram 

tmp <- data.table(
    variable = c("nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA"),
    value = c(args$nFeature_RNA, args$mitochondrial_percent_RNA, args$ribosomal_percent_RNA)
)

p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free", nrow=1) +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.5))
    )
    
pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outputdir), width=12, height=6)
# pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outputdir))
print(p)
dev.off()


########################################################
## Plot fraction of cells that pass QC for each sample ##
########################################################

to.plot <- metadata %>%
    .[,mean(pass_rnaQC),by="sample"]

p <- ggbarplot(to.plot, x="sample", y="V1", fill="gray70") +
    labs(x="", y="Fraction of cells that pass QC") +
    # facet_wrap(~stage)
    theme(
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
    )

pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outputdir))
print(p)
dev.off()

##########
## Save ##
##########

fwrite(metadata, paste0(args$outputdir,"/sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")

