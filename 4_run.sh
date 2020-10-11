#!/bin/bash

cd /mnt/scratch/laxmi/analysis/TwoSampleMR/

# Store results from MR

group=smoking_paternal_effect_birthweight
outdir=/mnt/scratch/laxmi/MRresult/${group}
mkdir -p ${outdir}

output=/mnt/scratch/laxmi/MRinstrument
mkdir -p ${output}

step=$1

# STEP 1: Format datafile of exposure 
if [[ $step == 1 ]] || [[ $step == "all" ]]
then

input=/mnt/cargo/laxmi/MRinstrument/smokingInitiation2019.txt.gz

Phenotype=smoking
SNP=RSID
beta=BETA
se=SE
eaf=AF
effect_allele=ALT
other_allele=REF
pval=PVALUE
chr=CHROM
position=POS
units=logical
ncase=NA
ncontrol=NA
samplesize=N
gene=NA
id=NA

# prepare_dataset_from_localFile
Rscript ./01_format_exposure_for_twoSampleMR.R \
        --input ${input} \
        --output ${output} \
        --Phenotype ${Phenotype} \
        --SNP ${SNP} \
        --beta ${beta} \
        --se ${se} \
        --eaf ${eaf} \
        --effect_allele ${effect_allele} \
        --other_allele ${other_allele} \
        --pval ${pval} \
        --chr ${chr} \
        --position ${position} \
        --units ${units} \
        --ncase ${ncase} \
        --ncontrol ${ncontrol} \
        --samplesize ${samplesize} \
        --gene ${gene} \
        --id ${id} \
        --outdir ${outdir}

wait

echo "done step 1"

fi

# STEP 2: Format datafile of outcome 

if [[ $step == 2 ]] || [[ $step == "all" ]]
then

exposure=${output}/smoking_clumped.txt
outcome=/mnt/cargo/laxmi/MRinstrument_short_version/Paternal_instrument_smoking.txt

Phenotype_type=continuous
Phenotype=birthweight
SNP=SNP
beta=est
se=SE
eaf=MAF
effect_allele=A1
other_allele=A2
pval=Pval_Estimate
chr=CHR
position=BP
units=grams
ncase=NA
ncontrol=NA
samplesize=N
gene=NA
id=NA

# prepare_dataset_from_localFile
Rscript ./02_format_outcome_for_twoSampleMR.R \
        --exposure ${exposure} \
        --outcome ${outcome} \
        --Phenotype_type ${Phenotype_type} \
        --Phenotype ${Phenotype} \
        --SNP ${SNP} \
        --beta ${beta} \
        --se ${se} \
        --eaf ${eaf} \
        --effect_allele ${effect_allele} \
        --other_allele ${other_allele} \
        --pval ${pval} \
        --chr ${chr} \
        --position ${position} \
        --units ${units} \
        --ncase ${ncase} \
        --ncontrol ${ncontrol} \
        --samplesize ${samplesize} \
        --gene ${gene} \
        --id ${id} \
        --outdir ${outdir}

wait

echo "done step 2"

fi

# STEP 3: perform_and_plot_TwoSampleMR
if [[ $step == 3 ]] || [[ $step == "all" ]]
then

datafile=${outdir}/smoking_birthweight_harmonized.txt

Rscript ./03_perform_and_plot_TwoSampleMR.R \
        --datafile ${datafile} \
        --method TRUE \
        --MRPRESSO TRUE \
        --heterogeneity TRUE \
        --pleiotropy TRUE \
        --singleSNP TRUE \
        --leaveoneout TRUE \
        --scatterplot TRUE \
        --forestplot TRUE \
        --leaveoneoutplot TRUE \
        --funnelplot TRUE \
        --mrreport TRUE \
        --directionality TRUE \
        --mrsteiger TRUE \
        --mrwrapper TRUE \
        --mrmoe TRUE \
        --outdir ${outdir}

wait

echo "done step 3"

fi
