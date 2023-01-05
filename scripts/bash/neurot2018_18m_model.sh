#!/bin/bash
#
# Job name (this is what the cluster will show as running - keep short):
#SBATCH --job-name=neurot2018_18m
#
#Project:
#SBATCH --account=p471
#
#Wall clock limit (change based on resources required):
#SBATCH --time=16:00:00 
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
module load PLINK/2.00a2.3_x86_64

INM="/cluster/p/p471/cluster/data/genetic_data/MoBaPsychGen_v1"
OUTM="/cluster/p/p471/cluster/projects/mdrink_gxe"

plink2 \
--bfile ${INM}/MoBaPsychGen_v1-ec-eur-batch-basic-qc \
--linear interaction no-x-sex \
--keep ${OUTM}/singleton_ids.tab \
--extract ${OUTM}/neurot2018_snps.tab \
--pheno ${OUTM}/phenos_18m.tab \
--covar ${OUTM}/covs_18m.tab \
--out ${OUTM}/neurot2018_18m 




