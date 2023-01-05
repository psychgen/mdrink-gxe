#06a_multilevel_funs.R

# This script prepares the multilevel models in Mplus, using any accepted constraints 
# from the linear modelling and only running with tier 1,2, and 3 covariates

library(tidyverse)

# Function to create a unique Mplus data file and input script for each pgs 

prep_model_syntax <- function(pgs, constraints="unconstrained",mplusdata){

   mplusdata_tmp <- mplusdata %>% 
      select(indid, m_id, sex, parity, matches("mdrink"), int_1, int_2, int_3, ext_1, ext_2, ext_3, matches(paste0("^",pgs,".pgs.pc")) ) %>% 
      rename_with(~str_remove_all(.,pgs)) %>% 
      rename_with(~str_replace_all(.,".pgs.pc.","pgs_")) %>% 
      rename_with(~str_replace_all(.,"child","c")) %>% 
      rename_with(~str_replace_all(.,"father","f")) %>% 
      rename_with(~str_replace_all(.,"mother","m")) 
   
   # Create the interaction variables in a three step process:
   
   # 1) Nest the relevant variables and pivot to wide (nested) data
   message("Creating interactions... \n")
   df_lists <- mplusdata_tmp %>% 
      select(mdrink_1, mdrink_2, mdrink_3,mdrink_pre,
             pgs_m, pgs_f,pgs_c,sex,parity) %>% 
      summarise_all(list) %>% 
      pivot_longer(cols = everything(), 
                   names_to = "var", 
                   values_to = "vector") %>% 
      print()
   
   # 2) Expand to make two-way combinations of all rows
   
   df_lists_comb <- expand(df_lists,
                           nesting(var, vector),
                           nesting(var2 = var, vector2 = vector)) %>% 
      print()
   
   # -- Filter out redundant combinations
   
   df_lists_comb <- df_lists_comb %>% 
      filter(var != var2) %>% 
      arrange(var, var2) %>% 
      mutate(vars = paste0(var, ".", var2)) %>% 
      select(contains("var"), everything()) %>% 
      print()
   
   c_sort_collapse <- function(...){
      c(...) %>% 
         sort() %>% 
         str_c(collapse = ".")
   }
   
   df_lists_comb_as <- df_lists_comb %>% 
      mutate(vars = map2_chr(.x = var, 
                             .y = var2, 
                             .f = c_sort_collapse)) %>%
      distinct(vars, .keep_all = TRUE)
   
   # -- Create interactions 
   
   ints <- df_lists_comb_as %>% 
      mutate(interaction = map2(vector, vector2, `*`))
   
   # 3) Pivot back and unnest, naming interactions as needed for Mplus syntax
   
   ints_long <- ints %>% 
      select(vars, interaction)  %>%  
      pivot_wider(values_from = interaction,
                  names_from = vars) %>%
      unnest(cols = everything()) %>% 
      select("icx11" = mdrink_1.mdrink_pre,
             "icx12" = mdrink_1.parity,
             "m1"   =  mdrink_1.pgs_c,
             "icx13" = mdrink_1.pgs_f,
             "icx14" = mdrink_1.pgs_m,
             "icx15" = mdrink_1.sex,
             "icx21" = mdrink_2.mdrink_pre,
             "icx22" = mdrink_2.parity,
             "m2"   = mdrink_2.pgs_c,
             "icx23" = mdrink_2.pgs_f,
             "icx24" = mdrink_2.pgs_m,
             "icx25" = mdrink_2.sex,
             "icx31" = mdrink_3.mdrink_pre,
             "icx32" = mdrink_3.parity,
             "m3"   = mdrink_3.pgs_c,
             "icx33" = mdrink_3.pgs_f,
             "icx34" = mdrink_3.pgs_m,
             "icx35" = mdrink_3.sex,
             "icm1" = mdrink_pre.pgs_c,
             "icm2" = parity.pgs_c,
             "icm3" = pgs_c.pgs_f,
             "icm4" = pgs_c.pgs_m,
             "icm5" = pgs_c.sex)
   
   # -- rejoin to main data
   mplusdata_tmp_w_ints <- mplusdata_tmp %>% 
      cbind(ints_long) 
   
   message("Preparing data... \n")
   prepareMplusData(mplusdata_tmp_w_ints, filename = paste0("./data/xpgs_sibs_",pgs,"_mplus.dat"),
                    inpfile = TRUE)
   
   # Read in the resulting input file, and use it as the basis of input files for
   # each model, with each outcome 
   
   stem <-readLines( paste0("./data/xpgs_sibs_",pgs,"_mplus.inp"))
   
   if(str_detect(pgs, "int")){
      outcomes = c("int")
   } else if(str_detect(pgs, "ext")){ 
      outcomes= c("ext")
   } else {
   outcomes= c("int","ext")
   }
   
   for(out in outcomes){
      
      stem1 <- str_replace(stem, "Your title goes here", paste(pgs,out, sep = "_")) %>% 
         str_replace("./data","//ess01/P471/data/durable/projects/mdrink_gxe/data\n")
      
      filecon <- paste0("./scripts/mplus/xpgs/", pgs, "_", out, "_",constraints,".inp")
      
      message(paste0("Preparing syntax for... ",out," \n") ) 
      
      writeLines(c(stem1,paste0("
USEVARIABLES= mdrink_1 mdrink_2 mdrink_3 ",out,"_1 ",out,"_2 ",out,"_3 
pgs_c sex parity mdrink_pre pgs_f pgs_m m1 m2 m3 
icx11 icx12 icx13 icx14 icx15 icm1 icm2 icm3 icm4 icm5
icx21 icx22 icx23 icx24 icx25 icx31 icx32 icx33 icx34 icx35;

CLUSTER= m_id;
WITHIN= sex parity pgs_c ;
BETWEEN= pgs_m pgs_f; 

DEFINE:
center sex (GRANDMEAN);
CENTER parity (GROUPMEAN);

ANALYSIS:
TYPE = TWOLEVEL ;
!MODEL = NOCOV;
H1ITERATIONS = 1;
PROCESSORS = 5;
MCONVERGENCE = 0.015; 
LOGCRITERION = 0.005;
MITERATIONS = 2000;


MODEL:
   %WITHIN%
   !Covariates
   ",out,"_1 ON parity sex mdrink_pre; 
   ",out,"_2 ON parity sex mdrink_pre; 
   ",out,"_3 ON parity sex mdrink_pre; 

   !Main effect, individual level
   ",out,"_1 ON mdrink_1  (b1_1);
   ",out,"_2 ON mdrink_2  (b1_2);
   ",out,"_3 ON mdrink_3  (b1_3);
   
   !Within-level moderation effect
   ",out,"_1 ON m1  (b3_1);
   ",out,"_2 ON m2  (b3_2);
   ",out,"_3 ON m3  (b3_3);
   
   !Moderator direct effect
   ",out,"_1 ON pgs_c (b2_1);
   ",out,"_2 ON pgs_c (b2_2);
   ",out,"_3 ON pgs_c (b2_3);

   !Covariate interactions, within level;
   ",out,"_1 ON icx11-icx15;
   ",out,"_2 ON icx21-icx25;
   ",out,"_3 ON icx31-icx35;
   ",out,"_1 ON icm1-icm5;
   ",out,"_2 ON icm1-icm5;
   ",out,"_3 ON icm1-icm5;

   !Variance maternal variable, individual level
   mdrink_1(ew_1);
   mdrink_2(ew_2);
   mdrink_3(ew_3);
   
   !Variance child outcome , individual level
   ",out,"_1 (ow_1);
   ",out,"_2 (ow_2);
   ",out,"_3 (ow_3);
   
      
   !Covariance, child outcomes, individual level
   ",out,"_1 WITH ",out,"_2;
   ",out,"_2 WITH ",out,"_3;
   ",out,"_1 WITH ",out,"_3;

   !Variance interaction terms, individual level
   m1 (mw_1);
   m2 (mw_2);
   m3 (mw_3);
   icx11-icx15;
   icx21-icx25;
   icx31-icx35;
   icm1-icm5;

   !Variance of the moderator and moderation effect to avoid listwise deletion;
   pgs_c m1;

   !Within level covariate variances (only when covariate is not ubiquitous to
   !avoid listwise deletion);
   mdrink_pre ;
    
   %BETWEEN%
   
   !EXPOSURE
   !Latent variance component at between level, exposure variable 
   i_m_1 BY mdrink_1@1;
   i_m_1 BY mdrink_2@1;
   i_m_1 BY mdrink_3@1;
   !Direct effect of between component on child outcome (estimated)
   i_m_1 BY ",out,"_1 (g1_1);
   i_m_1 BY ",out,"_2 (g1_2);
   i_m_1 BY ",out,"_3 (g1_3);
   !Indirect effect, between level
   ",out,"_1 ON mdrink_1  (b1_1);
   ",out,"_2 ON mdrink_2  (b1_2);
   ",out,"_3 ON mdrink_3  (b1_3);
   
   !COVARIATES
   !Between level covariate effects
   ",out,"_1 ON pgs_m pgs_f;
   ",out,"_2 ON pgs_m pgs_f;
   ",out,"_3 ON pgs_m pgs_f;
   
   !Between level covariate interaction effects
   ! These are currently omitted as the complexity was breaking things...
   ! Does this risk introducing bias, given that the interaction of interest
   ! is decomposed across both levels?


   !INTERACTION EFFECT
   !Latent variance component at between level, interaction effect 
   mb1 BY m1@1;
   mb1 BY m2@1;
   mb1 BY m3@1;
   !Direct effect of between component of interaction effect on child outcome (estimated)
   mb1 BY ",out,"_1 (g3_1);
   mb1 BY ",out,"_2 (g3_2);
   mb1 BY ",out,"_3 (g3_3);
      !Indirect effect, between level
   ",out,"_1 ON m1  (b3_1);
   ",out,"_2 ON m2  (b3_2);
   ",out,"_3 ON m3  (b3_3);

   !Latent variance component mean @ 0
   [ i_m_1@0 ];
   [ mb1@0 ];
   
   !Between level obs maternal variable variance @ 0 
   mdrink_1@0;
   mdrink_2@0;
   mdrink_3@0;
   !Between level obs int variance @ 0;
   m1@0;
   m2@0;
   m3@0;
   
   !Between level obs child variance @ 0 
   ",out,"_1 (ob_1);
   ",out,"_2 (ob_2);
   ",out,"_3 (ob_3);
      
   !Covariance, child outcomes, between level
   ",out,"_1 WITH ",out,"_2;
   ",out,"_2 WITH ",out,"_3;
   ",out,"_1 WITH ",out,"_3;
   
   !Variance latent variable maternal var
   i_m_1(eb_1);
   mb1 (mb_1);
   
   !Between level covariate variances;
   pgs_m pgs_f;
   
   !Estimate means at this level
   [ mdrink_1 mdrink_2 mdrink_3 m1 m2 m3];
   [ ",out,"_1 ] (b0_1);
   [ ",out,"_2 ] (b0_2);
   [ ",out,"_3 ] (b0_3);
   
MODEL CONSTRAINT:
",
if(constraints=="onelevel"){
c("!Check concordance with linear models if no between effects estimated
0= g3_1;  
0= g3_2;
0= g3_3;
0= g1_1;
0= g1_2;
0= g1_3;
")},
if(constraints=="effects"|constraints=="residuals"){
c("
!Test equaility of effects across time
0=b1_1-b1_2;
0=b1_2-b1_3;
0=b2_1-b2_2;
0=b2_2-b2_3;
0=b3_1-b3_2;
0=b3_2-b3_3;
0=g1_1-g1_2;
0=g1_2-g1_3;
0=g3_1-g3_2;
0=g3_2-g3_3;
")}, 
if(constraints=="residuals"){
c("!Test equality of residual variances across time
0=ow_1-ow_2;
0=ow_2-ow_3;
0=ob_1-ob_2;
0=ob_2-ob_3;
 ")},
"   
OUTPUT:
sampstat;
CINTERVAL;
stdyx;
sval;

savedata: FORMAT IS f10.5;
results are ",pgs,"_",out,"_",constraints,"_decimals.dat;   
")), con =filecon, useBytes=T)
      
   }
   
}

