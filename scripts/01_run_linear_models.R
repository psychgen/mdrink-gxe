#01_run_linear_models.R

# This script prepares and runs linear models, testing constraints and the effect
# of sequentially adding covariates

# Load packages

library(tidyverse)
library(lavaan)

# Read in data

load( file = './data/processed_data.RData')

# Source script with model functions

source('./scripts/01a_linear_model_funs.R')

# Create a data frame with inputs and to hold results

linear_mods <- expand.grid(names(alldata) %>% 
                .[str_detect(.,"cbcl")] %>% 
              str_sub(end=-5) %>% 
              unique() , 
              names(alldata) %>% 
                .[str_detect(.,"pgs")] %>% 
                str_remove_all(".pgs.pc.child|.pgs.pc.father|.pgs.pc.mother") %>% 
                unique()) %>% 
  `colnames<-`(c("out","mod")) %>% 
  as_tibble()

# Run the models across each combination of outcome and moderator, saving results in nested dataframe

linear_mods_res <- linear_mods %>% 
  mutate(results = pmap(., run_test_models, data=alldata),
         main_effs_tab = map(results, make_main_tab),
         int_effs_tab = map(results, make_int_tab),
         fit_comp_tab = map(results, make_fit_tab))

# Save the result

save(linear_mods_res, file="./output/linear_mods_res.RData")
