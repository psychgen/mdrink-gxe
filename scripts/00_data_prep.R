# 00_data_preparation.R


# load packages for data prep ####


library(tidyverse)
library(phenotools)
library(genotools)

# variables required:

# * Maternal at-risk drinking at 1.5, 3, 5 years (exposures):
# **  "EE608", "EE609",  "GG487", "GG488", "LL514", "LL515"
# * Child emotional and behavioural problems at 1.5, 3, and 5 years (outcomes):
# ** "cbcl_ext_c_18m", "cbcl_ext_c_5yr", "cbcl_ext_c_3yr", "cbcl_int_c_18m", "cbcl_int_c_5yr", "cbcl_int_c_3yr"
# * Child PGS (moderators):
# ** "neurot2018","ptsd2019","adhd","height2"
# * Prenatal at-risk drinking, sex, parity, maternal and paternal PGS (covariates):
# ** "AA1465","CC1165", "PARITET_5","KJONN"



# curate dataset ####

pheno <- curate_dataset(variables_required=
                          list(
                            moba=c(
                              "cbcl_ext_c_18m", "cbcl_int_c_18m",                      # CBCL 1.5yr
                              "cbcl_ext_c_3yr", "cbcl_int_c_3yr",                      # CBCL 3yr
                              "cbcl_ext_c_5yr", "cbcl_int_c_5yr",                      # CBCL 5yr
                              "EE608", "EE609",  "GG487", "GG488", "LL514", "LL515",   # exposures
                              "AA1465","CC1165", "PARITET_5","KJONN")),                # covariates
                        PDB = "2306",
                        moba_data_version = 12,
                        completion_threshold=0.5,
                        return_items=FALSE,
                        consistent_items=TRUE,
                        transformations=NULL,
                        exclusions=NULL,
                        out_format="merged_df")

save(pheno, file="./data/curated_dataset.RData")
load(file="./data/curated_dataset.RData")


pheno <-pheno %>% 
  mutate(mdrink_we_18m =  case_when(EE608_raw==1~3,
                                    EE608_raw==2~2,
                                    EE608_raw==3~1,
                                    EE608_raw==4~0.5,
                                    EE608_raw==5~0,
                                    EE608_raw==6~ NA_real_),
         mdrink_wd_18m =  case_when(EE609_raw==1~3,
                                    EE609_raw==2~3,
                                    EE609_raw==3~2,
                                    EE609_raw==4~1,
                                    EE609_raw==5~0,
                                    EE609_raw==6~ NA_real_),
         mdrink_we_3yr =  case_when(GG487_raw==1~3,
                                    GG487_raw==2~2,
                                    GG487_raw==3~1,
                                    GG487_raw==4~0.5,
                                    GG487_raw==5~0,
                                    GG487_raw==6~ NA_real_),
         mdrink_wd_3yr =  case_when(GG488_raw==1~3,
                                    GG488_raw==2~3,
                                    GG488_raw==3~2,
                                    GG488_raw==4~1,
                                    GG488_raw==5~0,
                                    GG488_raw==6~ NA_real_),
         mdrink_we_5yr =  case_when(LL514_raw==1~3,
                                    LL514_raw==2~2,
                                    LL514_raw==3~1,
                                    LL514_raw==4~0.5,
                                    LL514_raw==5~0,
                                    LL514_raw==6~ NA_real_),
         mdrink_wd_5yr =  case_when(LL515_raw==1~3,
                                   LL515_raw==2~3,
                                   LL515_raw==3~2,
                                   LL515_raw==4~1,
                                   LL515_raw==5~0,
                                   LL515_raw==6~ NA_real_),
         mdrink_pre_q1 =  case_when(AA1465_raw==1~3,
                                    AA1465_raw==2~3,
                                    AA1465_raw==3~2,
                                    AA1465_raw==4~1,
                                    AA1465_raw==5~0,
                                    AA1465_raw==6~ NA_real_),
         mdrink_pre_q3 =  case_when(CC1165_raw==1~3,
                                    CC1165_raw==2~3,
                                    CC1165_raw==3~2,
                                    CC1165_raw==4~1,
                                    CC1165_raw==5~0,
                                    CC1165_raw==6~ NA_real_)) %>% 
  rowwise() %>%
  mutate(mdrink_18m= round(mean(c(mdrink_we_18m,mdrink_wd_18m),na.rm=TRUE),0),
         mdrink_3yr=  round(mean(c(mdrink_we_3yr,mdrink_wd_3yr),na.rm=TRUE),0),
         mdrink_5yr=  round(mean(c(mdrink_we_5yr,mdrink_wd_5yr),na.rm=TRUE),0),
         mdrink_pre=  round(mean(c(mdrink_pre_q1,mdrink_pre_q3),na.rm=TRUE),0)) %>% 
  ungroup() %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  select(-matches("_wd_|_we_|_q1|_q3"))


geno <- fetch_pgs(c("neurot2018","ptsd2019","adhd","height2"),
                  maf = "0.01",
                  clump = "250_1_0.1")  

geno_procd <- process_pgs(geno)


geno_pcs_details <- geno_procd %>% 
  select(IID,  matches("_res")) %>% 
 pgs_pca(indid = "IID",
          pgs_var_stem = c("neurot2018","ptsd2019","adhd","height2"),
          return=c("pcs"))

corrs <- list(cor(geno_pcs_details %>% 
      select(matches("neurot2018"))) %>%  
  as.data.frame() %>% 
  rownames_to_column("pgs") %>% 
  filter(str_detect(pgs, "pgs.pc")) %>% 
  rename_all(.funs = list(~str_remove_all(.,"neurot2018"))),
  cor(geno_pcs_details %>% 
        select(matches("ptsd2019"))) %>%  
    as.data.frame() %>% 
    rownames_to_column("pgs") %>% 
    filter(str_detect(pgs, "pgs.pc")) %>% 
    rename_all(.funs = list(~str_remove_all(.,"ptsd2019"))),
  cor(geno_pcs_details %>% 
        select(matches("adhd"))) %>%  
    as.data.frame() %>% 
    rownames_to_column("pgs") %>% 
    filter(str_detect(pgs, "pgs.pc")) %>% 
    rename_all(.funs = list(~str_remove_all(.,"adhd"))),
  cor(geno_pcs_details %>% 
        select(matches("height2"))) %>%  
    as.data.frame() %>% 
    rownames_to_column("pgs") %>% 
    filter(str_detect(pgs, "pgs.pc")) %>% 
    rename_all(.funs = list(~str_remove_all(.,"height2"))))


highest_rcol <- function(d) {
  
  r <- d%>% unlist() %>% which.min %>% names() %>% as_tibble
  return(r)
  
}

#max r2 threshold
corrs <- corrs %>% 
  purrr::reduce(bind_rows) %>% 
  select(-.pgs.pc) %>% 
  pivot_longer(matches("_p"), names_to="threshold",values_to="r")


save(corrs, file="./output/pgs_pca_threshold_rs.RData")


# NEUROT
# _p<0.001_res 
# 
# 
# PTSD
# _p<0.01_res 
# 
# 
# ADHD
# _p<0.005_res 
#  
# 
# HEIGHT
# _p<0.001_res 


geno_pcs_kids <-  geno_procd %>% 
  select(IID, preg_id, BARN_NR, matches("_res")) %>% 
  drop_na(preg_id) %>% 
  pgs_pca(indid = "IID",
          pgs_var_stem = c("neurot2018","ptsd2019","adhd","height2")) %>% 
  select( preg_id, BARN_NR, matches("pgs.pc")) %>% 
  rename_with(~paste0(.x,".child"), matches("pgs.pc")) %>% 
  mutate(BARN_NR=as.numeric(BARN_NR))



geno_pcs_dads <-  geno_procd %>% 
  select(IID, f_id, matches("_res")) %>% 
  drop_na(f_id) %>% 
  pgs_pca(indid = "IID",
          pgs_var_stem = c("neurot2018","ptsd2019","adhd","height2")) %>% 
  select( f_id, matches("pgs.pc")) %>% 
  rename_with(~paste0(.x,".father"), matches("pgs.pc"))

geno_pcs_mums <-  geno_procd %>% 
  select(IID, m_id, matches("_res")) %>% 
  drop_na(m_id) %>% 
  pgs_pca(indid = "IID",
          pgs_var_stem = c("neurot2018","ptsd2019","adhd","height2")) %>% 
  select(IID, m_id, matches("pgs.pc")) %>% 
  rename_with(~paste0(.x,".mother"), matches("pgs.pc")) 


alldata <- pheno %>% 
  left_join(geno_pcs_kids) %>% 
  left_join(geno_pcs_dads) %>% 
  left_join(geno_pcs_mums)

save(alldata, file = './data/unscaled_data.RData')

# Final tidy up, including setting sex coding to 0-1, scaling numeric variables and adding a flag variable for genetic data availablity

alldata <- alldata %>% 
  select(preg_id,BARN_NR,m_id,"sex"=KJONN_raw,"parity"=PARITET_5_raw, matches("cbcl"),matches("mdrink"),matches("pgs")) %>% 
  filter(sex %in% c(1,2)) %>% 
  mutate(sex=as.factor(sex),
         BARN_NR= as.factor(BARN_NR),
         across(where(is.numeric),scale),
         across(where(is.numeric),as.numeric),
         sex=as.numeric(sex)-1,
         include= ifelse((!is.na(adhd.pgs.pc.child)&(!is.na(cbcl_ext_c_18m)|!is.na(mdrink_18m)|
                                                       !is.na(cbcl_ext_c_3yr)|!is.na(mdrink_3yr)|
                                                       !is.na(cbcl_ext_c_5yr)|!is.na(mdrink_5yr))),
                         1,0))

# write out processed dataset ####

#save(alldata, file = './data/processed_data.RData')
#see below for new savepoint


### Father's PGS differ within sibships ###

#It was observed in the multi-level models that a) fathers' pgs could not be treated
#as a between variable, suggesting variation within a cluster, and b) that treating
#it as a within variable and specifying the variance term to bring it into the model and
#stop listwise deletion (as with other covariates) was resulting in large effects
#that are likely attributable to missingness (i.e., some kind of selection effect
#to do with parental separation). This effect was absent in lavaan models which handles
#missingness on covariates differently (without list-wise deletion or needing them to
#be modelled).

# To explore this:

sibs <- alldata %>% 
  filter(duplicated(m_id))

sibpairs <- alldata %>% 
  filter(m_id %in% sibs$m_id)

diff_pgsfs_within_sibs <- sibpairs %>% group_by(m_id) %>%
summarize(pgs_fs = paste0(unique(neurot2018.pgs.pc.father), collapse = ',')) %>% 
  tidyr::separate(pgs_fs, into = c("val1","val2","val3"), sep="," ) %>% 
  filter(!is.na(val2))

#Approx 2,700, of which...
not_due_to_missingness <- diff_pgsfs_within_sibs %>% 
  filter(!str_detect(val1,"NA")&!str_detect(val2,"NA"))

#only about 60 have truly varying paternal PGS (i.e., the father is another MoBa father)

#simplest solution is to set all these paternal PGS to missing, and treat as a between
#level variable as per maternal pgs

alldata <- alldata %>% 
  mutate(across(.cols=matches("father"), ~ ifelse(m_id %in% diff_pgsfs_within_sibs$m_id, NA, .x) ))

save(alldata, file = './data/processed_data.RData')

# Also do this for the unscaled data for the descriptives

load(file='./data/unscaled_data.Rdata')

alldata <- alldata %>% 
  mutate(across(.cols=matches("father"), ~ ifelse(m_id %in% diff_pgsfs_within_sibs$m_id, NA, .x) ))

save(alldata, file = './data/unscaled_data.RData')


