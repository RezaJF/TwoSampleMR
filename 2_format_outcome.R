# 2-Sample MR study
library(data.table)
library(devtools)
library(TwoSampleMR) # install_github("MRCIEU/TwoSampleMR")
library(MRInstruments)
library(ieugwasr)
library(MendelianRandomization)
library(MRPRESSO)
library(ggplot2)

options(stringsAsFactors=F)
library(data.table)
library(optparse)
library(GenABEL)

option_list <- list(
  make_option("--exposure", type="character", default="",
    help="instrument for the exposure"),
  make_option("--outcome", type="character", default="",
    help="instrument for the outcome"),
  make_option("--Phenotype_type", type="character", default="",
    help="type of exposure variable, [default='']"),
  make_option("--Phenotype", type="character", default="",
    help="exposure variable, [default='']"),
  make_option("--SNP", type="character", default="",
    help="RSID/rs number of the SNP, [default='']"),
  make_option("--beta", type="character", default="",
    help="instrument downloaded from the source, [default='']"),
  make_option("--se", type="character", default="",
    help="standard error of the estimate, [default='']"),
  make_option("--eaf", type="character", default="",
    help="effect/alternate allele frequency, [default='']"),
  make_option("--effect_allele", type="character", default="",
    help="effect/alternate allele, [default='']"),
  make_option("--other_allele", type="character", default="",
    help="reference allele, [default='']"),
  make_option("--pval", type="character", default="",
    help="p value of the estimate, [default='']"),
  make_option("--chr", type="character", default="",
    help="chromosome number, [default='']"),
  make_option("--position", type="character", default="",
    help="base pair position of SNP, [default='']"),
  make_option("--units", type="character", default="",
    help="measurement unit, [default='']"),
  make_option("--ncase", type="character", default="",
    help="numbers of cases, [default='']"),
  make_option("--ncontrol", type="character", default="",
    help="numbers of controls, [default='']"),
  make_option("--samplesize", type="character", default="",
    help="sample size of the study, [default='']"),
  make_option("--gene", type="character", default="",
    help="name of the gene, [default='']"),
  make_option("--id", type="character", default="",
    help="id of the variant, [default='']"),
  make_option("--outdir", type="character", default="./output",
    help="Output directory [default='~/output']")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

outdir <- opt$outdir
if(!file.exists(outdir)) system(paste("mkdir -p",outdir))

# Get instruments for exposure

expo <- fread(opt$exposure, sep="\t",header=TRUE, stringsAsFactor=FALSE)

expo <- as.data.frame(expo, fix.empty.names=FALSE)

dim(expo)

# Get instruments for outcome that matches to exposure
file <- fread(opt$outcome, sep="\t",header=TRUE, stringsAsFactor=FALSE)

dim(file)

if(!opt$Phenotype == "NA") file$Phenotype <- opt$Phenotype
if(!opt$SNP == "NA") colnames(file)[colnames(file) == opt$SNP] <- "SNP"
if(!opt$beta == "NA") colnames(file)[colnames(file) == opt$beta] <- "beta"
if(!opt$se == "NA") colnames(file)[colnames(file) == opt$se] <- "se"
if(!opt$eaf == "NA") colnames(file)[colnames(file) == opt$eaf] <- "eaf"
if(!opt$effect_allele == "NA") colnames(file)[colnames(file) == opt$effect_allele] <- "effect_allele"
if(!opt$other_allele == "NA") colnames(file)[colnames(file) == opt$other_allele] <- "other_allele"
if(!opt$pval == "NA") colnames(file)[colnames(file) == opt$pval] <- "pval"
if(!opt$chr == "NA") colnames(file)[colnames(file) == opt$chr] <- "chr"
if(!opt$position == "NA") colnames(file)[colnames(file) == opt$position] <- "position"
if(!opt$units == "NA") file$units <- opt$units
if(!opt$ncase == "NA") colnames(file)[colnames(file) == opt$ncase] <- "ncase"
if(!opt$ncontrol == "NA") colnames(file)[colnames(file) == opt$ncontrol] <- "ncontrol"
if(!opt$samplesize == "NA") colnames(file)[colnames(file) == opt$samplesize] <- "samplesize"
if(!opt$gene == "NA") colnames(file)[colnames(file) == opt$gene] <- "gene"
if(!opt$id == "NA") colnames(file)[colnames(file) == opt$id] <- "id"

file[file$SNP!="",] <- file
dim(file)

file <- as.data.frame(file, fix.empty.names=FALSE)

out <- format_data(file, type="outcome")

out <- setDT(out)[SNP %chin% expo$SNP]
dim(out)
head(out)

expo <- setDT(expo)[SNP %chin% out$SNP]
dim(expo)
head(expo)

write.table(out, paste0(outdir, "/", opt$Phenotype, "_matched.txt"), sep = "\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

# Harmonise the exposure and outcome data
dat <- harmonise_data(
    exposure_dat = expo,
    outcome_dat = out
)

# Drop duplicate exposure-outcome summary sets

#if(!opt$Phenotype_type == "NA") dat <-power_prune(dat,method=1,dist.outcome=opt$Phenotype_type)
#if(!opt$Phenotype_type == "NA") dat <-power_prune(dat,method=2,dist.outcome=opt$Phenotype_type)

write.table(dat, paste0(outdir,"/",basename(gsub("_clumped.txt$","_",opt$exposure)), opt$Phenotype, "_harmonized.txt"), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)

#write.table(dat, paste0(outdir,"/",basename(gsub("_clumped.txt$","_",opt$exposure)), basename(gsub(".txt.gz$","_harmonised.txt",opt$outcome))), sep = "\t", quote=F, col.names = TRUE, row.names = FALSE)

