#02_run_multilevel_models.R

# This script runs multilevel models in Mplus, using any accepted constraints 
# from the linear modelling and only running with tier 1,2, and 3 covariates

# Load packages

library(tidyverse)
library(MplusAutomation)
library(tictoc)

# Read in data

load( file = './data/processed_data.RData')

# Prepare for mplus - prepareMplusData needs numeric or factor vars only

mplusdata <- alldata %>% 
  select(preg_id, BARN_NR, m_id, sex, parity, matches("mdrink"), matches("cbcl"),matches("pgs")) %>% 
  mutate(indid= factor(paste0(preg_id,"_",BARN_NR)),
         m_id=factor(m_id)) %>% 
  rename_with(~str_remove_all(.,"cbcl_|_c_")) %>% 
  rename_with(~str_replace_all(.,"18m","_1"))%>% 
  rename_with(~str_replace_all(.,"3yr","_2"))%>% 
  rename_with(~str_replace_all(.,"5yr","_3"))%>% 
  rename_with(~str_replace_all(.,"mdrink__","mdrink_"))

# Make mini version to test without long wait times:

mplusdata_mini <- mplusdata %>% 
  slice_sample(n=2000)

# Source functions script

source('./scripts/02a_multilevel_model_funs.R')

# Create data files and model syntax files for each pgs using sourced function

# -- Create a data frame with inputs 

args <- expand.grid(names(alldata) %>% 
                             .[str_detect(.,"pgs")] %>% 
                             str_remove_all(".pgs.pc.child|.pgs.pc.father|.pgs.pc.mother") %>% 
                             unique(),
                    c("unconstrained","effects","residuals")) %>% 
  `colnames<-`(c("pgs","constraints")) %>% 
  mutate_all(as.character) %>% 
  as_tibble()

# -- Map inputs to function from funs files

pmap(args, prep_model_syntax,  mplusdata=mplusdata)

# Run models using Mplus via MplusAutomation (expect several hours per model with full sample...)


filepath1 <- "./scripts/mplus"

tic()
runModels(filepath1, recursive =F,replaceOutfile="modifiedDate", Mplus_command = "C:/Program Files/Mplus/Mplus" )
toc()

# Move the results files to the output file and clean up the scripts folder
file.copy(from=paste0(filepath1,"/",list.files(filepath1)[str_detect(list.files(filepath1),".inp",negate=T)]),
          to="./output/mplus",
          overwrite = TRUE, recursive = F,
          copy.mode = TRUE)

junk <- dir(path=filepath1, pattern=".out|.dat") 
file.remove(paste0(filepath1,"/",junk))

# Read in output

mplusOutput <- readModels("./scripts/mplus/output/gmm/step1", recursive=FALSE)
mySummaries <- readModels("./scripts/mplus/output/gmm/step1", recursive=FALSE, what="summaries")


