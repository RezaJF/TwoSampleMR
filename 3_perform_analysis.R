# 2-Sample MR study
library(data.table)
library(devtools)
library(TwoSampleMR) # install_github("MRCIEU/TwoSampleMR")
library(MRInstruments)
library(ieugwasr)
library(MendelianRandomization)
library(MRPRESSO)
library(ggplot2)
library(dplyr)
library(tidyjson)
library(randomForest)
library(car)

options(stringsAsFactors=F)
library("data.table")
library("optparse")
library("GenABEL")

option_list <- list(
  make_option("--datafile", type="character", default="",
    help="datafile for the anaysis with exposure and outcome harmonized"),
  make_option("--method", type="logical", default="FALSE",
    help="methods for the Mendelian Randomization: TwoSampleMR::mr_method_list(), [default='mr_ivw, mr_egger_regression, mr_simple_median,  mr_weighted_median, mr_simple_mode, mr_weighted_mode']"),
  make_option("--MRPRESSO", type="logical", default=FALSE,
    help="MR-PRESSO method, [default='FALSE']"),
  make_option("--heterogeneity", type="logical", default=FALSE,
    help="Heterogeneity statistics: ivw, egger, [default='FALSE']"),
  make_option("--pleiotropy", type="logical", default=FALSE,
    help="Horizontal pleiotropy: egger, [default='FALSE']"),
  make_option("--singleSNP", type="logical", default=FALSE,
    help="Single SNP analysis: ivw, egger, [default='FALSE']"),
  make_option("--leaveoneout", type="logical", default=FALSE,
    help="Leave-one-out analysis: ivw, egger, [default='FALSE']"),
  make_option("--scatterplot", type="logical", default=FALSE,
    help="Scatter plot, [default='FALSE']"),
  make_option("--forestplot", type="logical", default=FALSE,
    help="Forest plot: singleSNP, [default='FALSE']"),
  make_option("--leaveoneoutplot", type="logical", default=FALSE,
    help="Leave-one-out plot: leaveoneout, [default='FALSE']"),
  make_option("--funnelplot", type="logical", default=FALSE,
    help="Funnel plot: singleSNP, [default='FALSE']"),
  make_option("--mrreport", type="logical", default=FALSE,
    help="MR report as Html file, [default='FALSE']"),
  make_option("--directionality", type="logical", default=FALSE,
    help="MR-Steiger directionality test, [default='FALSE']"),
  make_option("--mrsteiger", type="logical", default=FALSE,
    help="MR-Steiger sensitivity test for directionality, [default='FALSE']"),
  make_option("--mrwrapper", type="logical", default=FALSE,
    help="MR-results from all 11 methods, [default='FALSE']"),
  make_option("--mrmoe", type="logical", default=FALSE,
    help="MR-Mixture of Experts: AUC for all MR-methods, [default='FALSE']"),
  make_option("--outdir", type="character", default="./output",
    help="Output directory, [default='~/output']")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

outdir <- opt$outdir
if(!file.exists(outdir)) system(paste("mkdir -p",outdir))

# Get datafile with instruments for exposure and outcome
dat <- read.table(opt$datafile, sep = "\t", header = TRUE, stringsAsFactor = FALSE)
dat <- as.data.frame(dat, fix.empty.names=FALSE)

# Perform MR
res <- if(opt$method) {mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_simple_median",  "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode", "mr_raps"))}
write.table(res, paste0(outdir,"/", "MR_methods.txt"), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)

"performMR"

# Perform MR-PRESSO global method
res_presso <- if(opt$MRPRESSO) {run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)} # mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
if(opt$MRPRESSO) write.table(res_presso, paste0(outdir,"/", "MRPRESSO", ".txt"), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)

"MRPRESSO"

# Sensitivity analysis
# Heterogeneity statistics
res_heterogeneity <- if(opt$heterogeneity) {mr_heterogeneity(dat, method_list=c("mr_ivw", "mr_egger_regression"))}
if(opt$heterogeneity) write.table(res_heterogeneity, paste0(outdir,"/", "heterogeneity", ".txt"), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)

"Heterogeneity"

# Horizontal pleiotropy
res_pleiotropy <- if(opt$pleiotropy) {mr_pleiotropy_test(dat)}
if(opt$pleiotropy) write.table(res_pleiotropy, paste0(outdir,"/", "pleiotropy", ".txt"), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)

"Pleiotropy"

# Single SNP analysis
res_singleSNP <- if(opt$singleSNP) {mr_singlesnp(dat, all_method=c("mr_ivw", "mr_egger_regression"))}
if(opt$singleSNP) write.table(res_singleSNP, paste0(outdir,"/", "singleSNP", ".txt"), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)

"singleSNP"

# Leave-one-out analysis
res_loo <- if(opt$leaveoneout) {mr_leaveoneout(dat)}
if(opt$leaveoneout) write.table(res_loo, paste0(outdir,"/", "leaveoneout", ".txt"), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)

"LOO"

# Plotting
# Scatter plot
p1 <- if(opt$scatterplot) {mr_scatter_plot(res, dat)}
if(opt$scatterplot) ggsave(p1[[1]], file = paste0(outdir,"/", "scatterplot", ".png") ,  width=7, height=7, dpi = 300)

"scatterplot"

# Forest plot
p2 <- if(opt$forestplot) {mr_forest_plot(res_singleSNP)}
if(opt$forestplot) ggsave(p2[[1]], file = paste0(outdir,"/", "forestplot_ivw", ".png") ,  width=7, height=7, dpi = 300)

"forestplot"

# Leave-one-out plot
p3 <- if(opt$leaveoneoutplot) {mr_leaveoneout_plot(res_loo)}
if(opt$leaveoneoutplot) ggsave(p3[[1]], file = paste0(outdir,"/", "leaveoneoutplot_ivw", ".png") ,  width=7, height=7, dpi = 300)

"looplot"

# Funnel plot
p4 <- if(opt$funnelplot) mr_funnel_plot(res_singleSNP)
if(opt$funnelplot) ggsave(p4[[1]], file = paste0(outdir,"/", "funnelplot_ivw", ".png") ,  width=7, height=7, dpi = 300)

"funnelplot"

# Generate MR report
if(opt$mrreport) {mr_report(dat,output_path = opt$outdir,output_type = "html",author = "Laxmi",study = paste0(outdir,"/", "TwoSampleMR"),path = system.file("reports", package = "TwoSampleMR"))}

"MRreport"

# MR Steiger directionality test
res_directionality <- if(opt$directionality) {directionality_test(dat)}
if(opt$directionality) write.table(res_directionality, paste0(outdir,"/", "directionality", ".txt"), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)

"directionality"

# MR Steiger directionality test - sensitivity analysis
res_mrsteiger <- if(opt$mrsteiger) {mr_steiger(p_exp = dat$pval.exposure, p_out = dat$pval.outcome, n_exp = dat$samplesize.exposure, n_out = dat$samplesize.outcome, r_xxo = 1, r_yyo = 1, r_exp=0, r_out=0)}
if(opt$mrsteiger) write.table(res_mrsteiger[1:12], paste0(outdir,"/", "mrsteiger", ".txt"), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)
if(opt$mrsteiger) png(filename = paste0(outdir,"/", "mrsteiger_sensitivity_plot", ".png") ,  width=7, height=7, res = 300)
if(opt$mrsteiger) res_mrsteiger[13]
if(opt$mrsteiger) dev.off()

"MRsteiger"

# MR-MoE: Using a mixture of experts machine learning approach
load("/mnt/cargo/laxmi/rf.rdata") # Load the downloaded RData object. This loads the rf object

"loadrf"

res_wrapper <- if(opt$mrwrapper) {mr_wrapper(dat)} # Obtain estimates from all methods, and generate data metrics
if(opt$mrwrapper) save(res_wrapper, file = paste0(outdir,"/", "mrwrapper", ".RData"))

"MRwrapper"

res_moe <- if(opt$mrmoe) {mr_moe(res_wrapper, rf)} # MR-MoE - predict the performance of each method
if(opt$mrmoe) save(res_moe, file = paste0(outdir,"/", "mrmoe", ".RData"))

"MRmoe"

"done step 3"
