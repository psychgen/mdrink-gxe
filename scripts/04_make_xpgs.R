#04_make_xPGS.R

#This script makes and processes the xPGS, ready for modelling

# Load required packages
library(genotools)
library(tidyverse)

## First, recreate sstats files with pvalues from SNP models

# Get snp lists
snp_lists <- genotools::get_pgs_snps(c("neurot2018","ptsd2019","adhd","height2"),
                                     maf = "0.01",
                                     clump = "250_1_0.1")

# Load threshold-PGS.PC corrs to see which threshold to use for each

load(file="./output/pgs_pca_threshold_rs.RData")

snp_lists <- snp_lists$neurot2018 %>% 
  mutate(pgs="neurot2018") %>% 
  bind_rows(snp_lists$ptsd2019 %>% 
              mutate(pgs="ptsd2019")) %>% 
  bind_rows(snp_lists$adhd %>% 
              mutate(pgs="adhd")) %>%
  bind_rows(snp_lists$height2 %>% 
              mutate(pgs="height2")) %>% 
  left_join(corrs %>% 
              group_by(pgs) %>% 
              slice_max(r) %>% 
              mutate(threshold=as.numeric(str_remove_all(threshold, "_p<|_res")),
                     pgs=str_remove_all(pgs,".pgs.pc"))) %>% 
  mutate(include= ifelse(P<threshold,1,0)) %>% 
  filter(include==1)


# Filter down to only post-clumping SNPs (only do this once - means that can work with files in R)

genotools::pgs_metadata %>% 
  filter(Pheno_shortname %in% c("neurot2018","ptsd2019","adhd","height2")) %>% 
  .$Sumstats_filename

#ADHD
if(! "adhd_filtered" %in% list.files("./data/sstats") ){
  adhd <- readr::read_tsv("./data/sstats/adhd_eur_jun2017")
  adhd_flt <- adhd %>% 
    filter(SNP %in% filter(snp_lists,pgs=="adhd")$SNP )
  rm(adhd)
  readr::write_tsv(adhd_flt, file="./data/sstats/adhd_filtered" )
} else {adhd_flt <-readr::read_tsv( file="./data/sstats/adhd_filtered" )}

#Neurot
if(! "neurot_filtered" %in% list.files("./data/sstats") ){
  neurot <- readr::read_tsv("./data/sstats/29255261-GCST005232-EFO_0007660-build37.f.tsv.gz")
  neurot_flt <- neurot %>% 
    filter(variant_id %in% filter(snp_lists,pgs=="neurot2018")$SNP)
  rm(neurot)
  readr::write_tsv(neurot_flt, file="./data/sstats/neurot_filtered" )
} else {neurot_flt <-readr::read_tsv( file="./data/sstats/neurot_filtered" )}

#Height
if(! "height_filtered" %in% list.files("./data/sstats") ){
  height <- readr::read_tsv("./data/sstats/Meta-analysis_Wood_et_al+UKBiobank_2018.txt")
  height_flt <- height %>% 
    filter(SNP %in% filter(snp_lists,pgs=="height2")$SNP)
  rm(height)
  readr::write_tsv(height_flt, file="./data/sstats/height_filtered" )
} else {height_flt <-readr::read_tsv( file="./data/sstats/height_filtered" )}

#PTSD
if(! "ptsd_filtered" %in% list.files("./data/sstats") ){
  ptsd <- readr::read_tsv("./data/sstats/pts_all_freeze2_overall.results.gz")
  ptsd_flt <- ptsd %>% 
    filter(SNP %in% filter(snp_lists,pgs=="ptsd2019")$SNP)
  rm(ptsd)
  readr::write_tsv(ptsd_flt, file="./data/sstats/ptsd_filtered" )
} else {ptsd_flt <-readr::read_tsv( file="./data/sstats/ptsd_filtered" )}

## Read in results from SNP models (pre-filtered so that only interaction terms are present)

filelist <- paste0("./output/snp_mods/",list.files("./output/snp_mods"))

#assuming tab separated values with a header    
datalist <- lapply(filelist, function(x)read_tsv(x, col_names = F)) 
datalist <- setNames(datalist, list.files("./output/snp_mods"))


snp_res <- bind_rows(datalist, .id = "model") %>% 
  `colnames<-`( c("model","CHROM","POS","ID","REF", "ALT","A1","TEST","OBS_CT","BETA","SE","T_STAT","P_ix")) %>% 
  filter(!is.na(P_ix)) %>% 
  separate(model, into=c("pgs","outcome","family","type","null"), sep="\\.") %>% 
  select(-family,-type,-null) %>% 
  mutate(pgs=str_sub(pgs,end=-5)) %>% 
  right_join(snp_lists %>% 
               select("CHROM"=CHR,"ID"=SNP,"POS"=BP,pgs))

save(snp_res, file="./output/single_snp_ix_model_results.RData")
load(file="./output/single_snp_ix_model_results.RData")
# Selectively join the filtered sstats files to the P_ix column, creating new sstats files for
# each pgs outcome combination, plus one in which the average pvalue is taken across waves

adhd_sstats_new <- adhd_flt %>% 
  right_join(snp_res %>% 
               filter(pgs=="adhd") %>% 
               select(outcome, SNP=ID, ref_ix = REF, alt_ix = ALT, a1_ix=A1, B_ix=BETA, P_ix )) %>% 
  mutate(P_ix = case_when(A1==a1_ix  ~ P_ix,
                          A2==a1_ix  ~ 1-P_ix, #IMPORTANT - make sure tested allele and p value correspond
                          TRUE ~ NA_real_),
         B_ix = case_when(A1==a1_ix  ~ B_ix,
                           A2==a1_ix  ~ -B_ix, #IMPORTANT - make sure tested allele and beta correspond
                           TRUE ~ NA_real_))

height_sstats_new <- height_flt %>% 
  right_join(snp_res %>% 
               filter(pgs=="height2") %>% 
               select(outcome, SNP=ID, ref_ix = REF, alt_ix = ALT, a1_ix=A1,B_ix=BETA, P_ix )) %>% 
  mutate(P_ix = case_when(Tested_Allele==a1_ix  ~ P_ix,
                          Other_Allele==a1_ix  ~ 1-P_ix, #IMPORTANT - make sure tested allele and p value correspond
                          TRUE ~ NA_real_),
         B_ix = case_when(Tested_Allele==a1_ix  ~ B_ix,
                          Other_Allele==a1_ix  ~ -B_ix, #IMPORTANT - make sure tested allele and beta correspond
                          TRUE ~ NA_real_))

neurot_sstats_new <- neurot_flt %>% 
  right_join(snp_res %>% 
               filter(pgs=="neurot2018") %>% 
               select(outcome, variant_id=ID, ref_ix = REF, alt_ix = ALT, a1_ix=A1,B_ix=BETA, P_ix )) %>% 
  mutate(P_ix = case_when(effect_allele==a1_ix  ~ P_ix,
                          other_allele==a1_ix  ~ 1-P_ix, #IMPORTANT - make sure tested allele and p value correspond
                          TRUE ~ NA_real_),
         B_ix = case_when(effect_allele==a1_ix  ~ B_ix,
                           other_allele==a1_ix  ~ -B_ix, #IMPORTANT - make sure tested allele and beta correspond
                           TRUE ~ NA_real_)) %>% 
  rename("SNP" = variant_id) ## Change in name for SNP column for neurot must be reflected in PRSice file

ptsd_sstats_new <- ptsd_flt %>% 
  right_join(snp_res %>% 
               filter(pgs=="ptsd2019") %>% 
               select(outcome, SNP=ID, ref_ix = REF, alt_ix = ALT, a1_ix=A1,B_ix=BETA, P_ix )) %>% 
  mutate(P_ix = case_when(A1==a1_ix  ~ P_ix,
                          A2==a1_ix  ~ 1-P_ix, #IMPORTANT - make sure tested allele and p value correspond
                          TRUE ~ NA_real_),
         B_ix = case_when(A1==a1_ix  ~ B_ix,
                          A2==a1_ix  ~ -B_ix, #IMPORTANT - make sure tested allele and beta correspond
                          TRUE ~ NA_real_))

all_sstats_new <- list(adhd_sstats_new, height_sstats_new, neurot_sstats_new, ptsd_sstats_new) %>% 
  setNames(c("adhd","height2","neurot2018","ptsd2019"))

# apply function to make the separate sumstats files and submission scripts

source("./scripts/04a_xpgs_funs.R")

make_sstats(names(all_sstats_new))

# additionally make overall sumstats files (excl height),
# if only to get N effective tests by clumping across pgs traits

all_pgs_ex_height <-snp_res %>% 
  filter(pgs!="height2") %>% 
  mutate(out= ifelse(str_detect(outcome, "ext"),"ext","int")) %>% 
  group_by(ID,out) %>% 
  summarise(P_ix = mean(P_ix, na.rm=T),
            B_ix = mean(BETA, na.rm=T),
            CHROM= CHROM,
            POS=POS,
            REF=REF,
            ALT=ALT,
            A1=A1) %>% 
  distinct() %>% 
  mutate(A2= ifelse(A1==ALT, REF, ALT))

int= all_pgs_ex_height %>% filter(out=="int")
ext= all_pgs_ex_height %>% filter(out=="ext")


readr::write_tsv(int, file= "./data/sstats/allexclheight_int")

readr::write_tsv(ext, file= "./data/sstats/allexclheight_ext")

for(sstats in c("allexclheight_int","allexclheight_ext" )){
  genotools::make_prsice( user= "laurie",
                          jobname = sstats,
                          sumstats_filename = sstats,
                          outcome_type = "continuous",
                          cpu_time="08:00:00",
                          memory="32G",
                          inputs_dir="/cluster/p/p471/cluster",
                          outputs_dir="/cluster/p/p471/cluster/common/raw_prsice_output/",
                          genotype_dir="data/genetic_data/MoBaPsychGen_v1",
                          genotype_data="MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                          prsice_dir="/cluster/p/p471/cluster/common/prsice/",
                          A1 = "A1",
                          A2 = "A2",
                          stat= "B_ix",
                          pvalue= "P_ix",
                          snp= "ID",
                          chr = "CHROM",
                          bp = "POS",
                          thresholds = "5e-08,5e-07,5e-06,5e-05,5e-04,0.001,0.01,0.005,0.1,0.5,1",
                          clump_kb = "250",
                          clump_p = "1.000000",
                          clump_r2 = "0.100000",
                          lower = "5e-08",
                          maf = "0.01",
                          mhc= "exclude",
                          add_metadata = FALSE)
  
}
  

### --- transfer sstats files and PGS to cluster and run PRSice
#-#-#-#-#-#
### --- transfer files resulting from PRSice run to d/d/common/pgs_directory...




new_pgs_sel <- genotools::available_pgs() 

geno <- fetch_pgs(new_pgs_sel,
                  maf = "0.01",
                  clump = "250_1_0.1")  

geno_procd <- process_pgs(geno)


geno_pcs_kids <-  geno_procd %>% 
  select(IID, preg_id, BARN_NR, matches("_res")) %>% 
  drop_na(preg_id) %>% 
  pgs_pca(indid = "IID",
          pgs_var_stem =new_pgs_sel) %>% 
  select( preg_id, BARN_NR, matches("pgs.pc")) %>% 
  rename_with(~paste0(.x,".child"), matches("pgs.pc")) %>% 
  mutate(BARN_NR=as.factor(BARN_NR))


geno_pcs_dads <-  geno_procd %>% 
  select(IID, f_id, matches("_res")) %>% 
  drop_na(f_id) %>% 
  pgs_pca(indid = "IID",
          pgs_var_stem = new_pgs_sel) %>% 
  select( f_id, matches("pgs.pc")) %>% 
  rename_with(~paste0(.x,".father"), matches("pgs.pc"))

geno_pcs_mums <-  geno_procd %>% 
  select(IID, m_id, matches("_res")) %>% 
  drop_na(m_id) %>% 
  pgs_pca(indid = "IID",
          pgs_var_stem = new_pgs_sel) %>% 
  select(IID, m_id, matches("pgs.pc")) %>% 
  rename_with(~paste0(.x,".mother"), matches("pgs.pc")) 

# Reattach father ids for join

load( file = './data/unscaled_data.RData')
fid_alldata <- alldata %>% 
  select(preg_id, BARN_NR, m_id, f_id) %>% 
  mutate(BARN_NR = as.factor(BARN_NR))
# Read in phenotypic data
load( file = './data/processed_data.RData')
alldata <- alldata %>% 
  left_join(fid_alldata)


alldata_new <- alldata %>% 
  left_join(geno_pcs_kids) %>% 
  left_join(geno_pcs_dads) %>% 
  left_join(geno_pcs_mums) %>% 
  select(-f_id)

#Ensure that paternal xPGS that vary between sibs are set to missing, as per the original PGS
alldata_new <- alldata_new %>% 
  mutate(across(.cols=matches("father"), ~ ifelse(is.na(ptsd2019.pgs.pc.father), NA, .x) ))

save(alldata_new, file = './data/processed_data_xpgs.RData')
