#03_run_SNP_models.R

#This script runs the SNP-wise interaction models from clumped SNP lists, saving
#the pvalues for re-ranking SNPs in creating the xPGS

# load packages
library(tidyverse)
library(genotools)
library(tictoc)

# get clumped SNP lists
snp_lists <- genotools::get_pgs_snps(c("neurot2018","ptsd2019","adhd","height2"),
                        maf = "0.01",
                        clump = "250_1_0.1")

map(snp_lists,length)

write.table(snp_lists$neurot2018$SNP , file = paste0("./data/neurot2018_snps.tab"), quote= F, sep=" ", row.names=F ,na="-9", col.names = F)
write.table(snp_lists$ptsd2019$SNP , file = paste0("./data/ptsd2019_snps.tab"), quote= F, sep=" ", row.names=F,na="-9", col.names = F)
write.table(snp_lists$adhd$SNP , file = paste0("./data/adhd_snps.tab"), quote= F, sep=" ", row.names=F,na="-9", col.names = F)
write.table(snp_lists$height2$SNP , file = paste0("./data/height2_snps.tab"), quote= F, sep=" ", row.names=F,na="-9", col.names = F)
write.table(snp_lists$height2$SNP %>% sample(size=10) , file = paste0("./data/test_10_snps.tab"), quote= F, sep=" ", row.names=F,na="-9", col.names = F)


# identify singletons, who will be the training sample for the xPGS

load( file = './data/processed_data.RData')

sibs <- alldata %>% 
  filter(duplicated(m_id))

sibpairs <- alldata %>% 
  filter(m_id %in% sibs$m_id)

singletons <- alldata %>% 
  filter(!m_id %in% sibs$m_id)

# save a singletons id list for inclusion

singleton_ids <- singletons %>%
  select(preg_id, BARN_NR) %>% 
  left_join(covs) %>% 
  select(FID,IID) %>% 
  drop_na(IID)

write.table(singleton_ids, file = paste0("./data/singleton_ids.tab"), quote= F, sep=" ", row.names=F,na="-9")



# read in covariate file for IDs, and create versions with mdrink, sex, parity, 10PCs, and batch

geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc"
covs_dir = "//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1"

covs <- readr::read_tsv(paste0(covs_dir,"/",geno_data,"-cov.txt"), col_types = readr::cols(.default = "c")) %>%
  dplyr::mutate_at(dplyr::vars(dplyr::matches("PC")), as.numeric) %>%
  dplyr::select(id_moba = dplyr::matches("ID_",ignore.case = FALSE), dplyr::everything()) %>%
  dplyr::mutate(preg_id_BARN_NR = ifelse(Role =="Child", id_moba,NA),
                f_id = ifelse(Role == "Father",id_moba,NA),
                m_id = ifelse(Role == "Mother",id_moba,NA)) %>%
  tidyr::separate(preg_id_BARN_NR, into=c("preg_id", "BARN_NR"), sep="_") %>% 
  select(preg_id, BARN_NR, f_id, m_id, FID,IID, paste0("PC",seq(1,10)), genotyping_batch, imputation_batch)




# use cov file stem to create cov and pheno files

for(wave in c("18m","3yr","5yr")){
  
  covfile <- singletons %>% 
    select(preg_id, BARN_NR, sex, parity, paste0("mdrink_",wave)) %>% 
    left_join(covs) %>% 
    select(FID,IID, paste0("PC",seq(1,10)), genotyping_batch,  sex, parity, paste0("mdrink_",wave)) %>% 
    drop_na(IID)

  phenofile <- singletons %>% 
    select(preg_id, BARN_NR,  sex, parity, paste0("cbcl_ext_c_",wave), paste0("cbcl_int_c_",wave)) %>% 
    left_join(covs) %>% 
    select(FID,IID, paste0("cbcl_ext_c_",wave), paste0("cbcl_int_c_",wave)) %>% 
    drop_na(IID)
  
  write.table(covfile, file = paste0("./data/covs_",wave,".tab"), quote= F, sep=" ", row.names=F,na="-9")
  write.table(phenofile, file = paste0("./data/phenos_",wave,".tab"), quote= F, sep=" ", row.names=F,na="-9")
          
}


# create submission scripts 

for(pgs in c("neurot2018","ptsd2019","adhd","height2")){
  for(wave in c("18m","3yr","5yr")){
   
    fileConn<-file(paste0("./scripts/bash/",pgs,"_",wave,"_model.sh"), "wb")
    writeLines(paste0(c(
'#!/bin/bash
#
# Job name (this is what the cluster will show as running - keep short):
#SBATCH --job-name='),pgs,c('_'),wave,c('
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

plink2 \\
--bfile ${INM}/MoBaPsychGen_v1-ec-eur-batch-basic-qc \\
--linear interaction no-x-sex \\
--keep ${OUTM}/singleton_ids.tab \\
--extract ${OUTM}/'),pgs,c('_snps.tab \\
--pheno ${OUTM}/phenos_'),wave,c('.tab \\
--covar ${OUTM}/covs_'),wave,c('.tab \\
--out ${OUTM}/'),pgs,c('_'),wave,c(' 



')), fileConn)
    
    close(fileConn)
  
   }
}


# make a meta submission script
fileConn2<-file(paste0("./scripts/bash/meta_sub.sh"), "wb")
writeLines(paste0(c(
  '#!/bin/bash
#
# Job name (this is what the cluster will show as running - keep short):
#SBATCH --job-name=meta_sub_pgs
#
#Project:
#SBATCH --account=p471
#
#Wall clock limit (change based on resources required):
#SBATCH --time=00:15:00
#
#SBATCH --ntasks=1
#
#Output filename customization
#This specification is Jobname_User_JobID
#SBATCH --output=./%x_%u_%j.out
#
# Max memory usage (change based on resources required):
#SBATCH --mem-per-cpu=2G

## Set up job enviroment:
module init

mkdir -p ./reports

for FILE in ./*model.sh; do
echo ${FILE}
sbatch ${FILE}
sleep 1 # pause to be kind to the scheduler
done')) , fileConn2)
close(fileConn2)

# Move files to cluster and run there, filtering results to only the relevant interaction term rows with:

# for FILE in ./*.glm.linear; do
# echo ${FILE}
# awk '{ if ($7 ~ /ADDxmdrink/) { print } }' ${FILE} > ${FILE}.filtered
# rm ${FILE}
# done


