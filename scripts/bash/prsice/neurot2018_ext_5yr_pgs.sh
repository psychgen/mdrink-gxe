#!/bin/bash
#
# Job name (this is what the cluster will show as running - keep short):
#SBATCH --job-name=neurot2018_ext_5yr
#
#Project:
#SBATCH --account=p471
#
#Wall clock limit (change based on resources required):
#SBATCH --time=02:00:00
#
#SBATCH --ntasks=1
#
#Output filename customization
#This specification is Jobname_User_JobID
#SBATCH --output=./reports/%x_%u_%j.out
#
# Max memory usage (change based on resources required):
#SBATCH --mem-per-cpu=32G

## Set up job enviroment:
source /cluster/bin/jobsetup
module purge
module load R/3.5.0

# Name used for output file location and filestems (keep short)
ANALYSIS_SNAME="neurot2018_ext_5yr"
# Filepath for inputs
INM="/cluster/p/p471/cluster"
# Filepath for outputs
OUTM="/cluster/p/p471/cluster/common/raw_prsice_output/"
# Where is prsice installed?
PRDIR="/cluster/p/p471/cluster/common/prsice/"

# Who is running the PRS creation?
USR="laurie"
# When?
DATE="2022-08-04"
# With which MAF/clumping params?
MAFSET="maf0.01"
CLUMPSET="clump250_1_0.1"

## Files for PRSice
#BASE
# Amend with filepath/name for your summary statistics file
BASE=${INM}/common/gwas_sumstats/neurot2018_ext_5yr

#TARGET
TARGET=${INM}/data/genetic_data/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc
TNAME = MoBaPsychGen_v1-ec-eur-batch-basic-qc

#Makes the output directories if they do not already exist
mkdir -p ${OUTM}${USR}/${TNAME}/${DATE}/${MAFSET}/${CLUMPSET}/${ANALYSIS_SNAME}

set -o errexit

#It is important for the user to take responsibility for ensuring
#settings are as you want them - the script may run either way

Rscript ${PRDIR}PRSice.R \
--prsice ${PRDIR}PRSice_linux \
--base $BASE \
--target $TARGET \
--out ${OUTM}${USR}/${TNAME}/${DATE}/${MAFSET}/${CLUMPSET}/${ANALYSIS_SNAME}/$ANALYSIS_SNAME \
--A1 effect_allele \
--A2 other_allele \
--all-score  \
--bar-levels 5e-08,5e-07,5e-06,5e-05,5e-04,0.001,0.01,0.005,0.1,0.5,1 \
--bp base_pair_location \
--chr chormosome \
--clump-kb 250 \
--clump-p 1.000000 \
--clump-r2 0.100000 \
--fastscore  \
--lower 5e-08 \
--maf 0.01 \
--pvalue P_ix \
--print-snp \
--snp SNP \
--stat beta \
--thread 1 \
--no-regress \
--no-default \
--x-range chr6:25000000-34000000 \
--beta \
--binary-target F
