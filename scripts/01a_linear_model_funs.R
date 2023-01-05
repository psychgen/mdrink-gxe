#01a_linear_model_funs.R

# The purpose of this script is to make functions for running models in 
# 01_run_linear_models.R

library(tictoc)

run_test_models <- function(outname, #the root of the outcome name, before the wave-identifying suffix
                            modname, #the root of the moderating PGS, before .pgs.pc and the role-identifying suffix
                            data){ #the dataset
  
  
  crit = 0.05/(2*4) 
    
  tic()
  message(paste0("Running models for ",outname,", with moderator ",modname))
  ## TIER 1 COVS
  
  #First, specify models including only the tier 1 covs (parity and child sex) and their ix terms with the exposure and moderator as per the pre-reg 
  model_tier1 <-
    paste0('# Regression 18m
           ',outname,'_18m ~ b1_1*mdrink_18m +  b2_1*',
           modname,'.pgs.pc.child + b3_1*mdrink_18m:',
           modname,'.pgs.pc.child + parity + b5_1*mdrink_18m:parity + b6_1*',modname,'.pgs.pc.child:parity + ',
           'sex + b8_1*mdrink_18m:sex  + b9_1*',modname,'.pgs.pc.child:sex 
           # Regression 3yr
           ',outname,'_3yr ~ b1_2*mdrink_3yr + b2_2*',
           modname,'.pgs.pc.child + b3_2*mdrink_3yr:',
           modname,'.pgs.pc.child + parity + b5_2*mdrink_3yr:parity + b6_2*',modname,'.pgs.pc.child:parity + ',
           'sex + b8_2*mdrink_3yr:sex  + b9_2*',modname,'.pgs.pc.child:sex
           # Regression 5yr
           ',outname,'_5yr ~ b1_3*mdrink_5yr +  b2_3*',
           modname,'.pgs.pc.child + b3_3*mdrink_5yr:',
           modname,'.pgs.pc.child + parity + b5_3*mdrink_5yr:parity + b6_3*',modname,'.pgs.pc.child:parity + ',
           'sex + b8_3*mdrink_5yr:sex  + b9_3*',modname,'.pgs.pc.child:sex
           
           # Residual variances
           ',outname,'_18m ~~ outres1*',outname,'_18m
           ',outname,'_3yr ~~ outres2*',outname,'_3yr
           ',outname,'_5yr ~~ outres3*',outname,'_5yr
           
           # Covariances')
 
  #Fit the model with the tier 1 covariates
  
  fit_tier1 <-lavaan::sem(model_tier1,
                          data=data,
                          estimator="MLR",
                          missing="fiml.x",
                          se="robust",
                          cluster= "m_id")
  
  #Fit the model with the tier 1 covariates, with the first set of constraints (main and ix effects invariant across time - interpreting this
  #as all main and ix effects at present [i.e.,. inc covs], though perhaps this is overly restrictive and we could just test b1/b2/b3?)
  
  fit_tier1_constrained1 <-lavaan::sem(paste0(model_tier1,'
                                              # Constraints
                                              b1_1==b1_2
                                              b1_2==b1_3
                                              b2_1==b2_2
                                              b2_2==b2_3
                                              b3_1==b3_2
                                              b3_2==b3_3'),
                                       data=data,
                                       estimator="MLR",
                                       missing="fiml.x",
                                       se="robust",
                                       cluster= "m_id") 
  
  #Fit the model with the tier 1 covariates, with the first set of constraints (main and ix effects AND residuals invariant across time)
  fit_tier1_constrained2 <-lavaan::sem(paste0(model_tier1,'
                                              # Constraints
                                              b1_1==b1_2
                                              b1_2==b1_3
                                              b2_1==b2_2
                                              b2_2==b2_3
                                              b3_1==b3_2
                                              b3_2==b3_3
                                              
                                              outres1==outres2
                                              outres2==outres3'),
                                       data=data,
                                       estimator="MLR",
                                       missing="fiml.x",
                                       se="robust",
                                       cluster= "m_id")
  
  #Perform the LRT fit comparisons to evaluate the acceptability of the contstraints
  
  test_constraints_tier1 <- rbind(anova(fit_tier1,fit_tier1_constrained1),
            anova(fit_tier1,fit_tier1_constrained2)[2,]) %>% 
    rownames_to_column("model") %>%  
    as_tibble()
  
  #Identify the best model on the basis of the LRTs
  
  if(!any(test_constraints_tier1$`Pr(>Chisq)`>crit, na.rm=T)){
    best_tier1 = list("t1_unconstrained"= fit_tier1)
  }else if(all(test_constraints_tier1$`Pr(>Chisq)`>crit, na.rm=T)){
    best_tier1 = list("t1_betas+resids"= fit_tier1_constrained2)
  }else {
    best_tier1= list("t1_betas"= fit_tier1_constrained1)
  }
  
  t1_mods <-list("t1_unconstrained"= fit_tier1,"t1_betas"= fit_tier1_constrained1,"t1_betas+resids"= fit_tier1_constrained2)
  
  message("Model-fitting with Tier 1 covariates complete.")
  
  ## TIER 1+2 COVS
  #Repeat for models with tier 2 covs (prenatal exposure to at-risk drinking) and all ix terms added
  model_tier2 <-
    paste0('# Regression 18m
           ',outname,'_18m ~ b1_1*mdrink_18m +  b2_1*',
           modname,'.pgs.pc.child + b3_1*mdrink_18m:',
           modname,'.pgs.pc.child + parity + b5_1*mdrink_18m:parity + b6_1*',modname,'.pgs.pc.child:parity + ',
           'sex + b8_1*mdrink_18m:sex  + b9_1*',modname,'.pgs.pc.child:sex + ',
           'mdrink_pre + b11_1*mdrink_18m:mdrink_pre + b12_1*',modname,'.pgs.pc.child:mdrink_pre
           # Regression 3yr
           ',outname,'_3yr ~ b1_2*mdrink_3yr + b2_2*',
           modname,'.pgs.pc.child + b3_2*mdrink_3yr:',
           modname,'.pgs.pc.child + parity + b5_2*mdrink_3yr:parity + b6_2*',modname,'.pgs.pc.child:parity + ',
           'sex + b8_2*mdrink_3yr:sex  + b9_2*',modname,'.pgs.pc.child:sex + ',
           'mdrink_pre + b11_2*mdrink_3yr:mdrink_pre + b12_2*',modname,'.pgs.pc.child:mdrink_pre
           # Regression 5yr
           ',outname,'_5yr ~ b1_3*mdrink_5yr +  b2_3*',
           modname,'.pgs.pc.child + b3_3*mdrink_5yr:',
           modname,'.pgs.pc.child + parity + b5_3*mdrink_5yr:parity + b6_3*',modname,'.pgs.pc.child:parity + ',
           'sex + b8_3*mdrink_5yr:sex  + b9_3*',modname,'.pgs.pc.child:sex + ',
           'mdrink_pre + b11_3*mdrink_5yr:mdrink_pre + b12_3*',modname,'.pgs.pc.child:mdrink_pre
           
           # Residual variances
           ',outname,'_18m ~~ outres1*',outname,'_18m
           ',outname,'_3yr ~~ outres2*',outname,'_3yr
           ',outname,'_5yr ~~ outres3*',outname,'_5yr
           
           # Covariances')
  
  fit_tier2 <-lavaan::sem(model_tier2,
                          data=data,
                          estimator="MLR",
                          missing="fiml.x",
                          se="robust",
                          cluster= "m_id")
  
  fit_tier2_constrained1 <-lavaan::sem(paste0(model_tier2,'
                                              # Constraints
                                              b1_1==b1_2
                                              b1_2==b1_3
                                              b2_1==b2_2
                                              b2_2==b2_3
                                              b3_1==b3_2
                                              b3_2==b3_3'),
                                       data=data,
                                       estimator="MLR",
                                       missing="fiml.x",
                                       se="robust",
                                       cluster= "m_id") 
  fit_tier2_constrained2 <-lavaan::sem(paste0(model_tier2,'
                                              # Constraints
                                              b1_1==b1_2
                                              b1_2==b1_3
                                              b2_1==b2_2
                                              b2_2==b2_3
                                              b3_1==b3_2
                                              b3_2==b3_3
                                              
                                              outres1==outres2
                                              outres2==outres3'),
                                       data=data,
                                       estimator="MLR",
                                       missing="fiml.x",
                                       se="robust",
                                       cluster= "m_id")
  
  test_constraints_tier2 <- rbind(anova(fit_tier2,fit_tier2_constrained1),
                            anova(fit_tier2,fit_tier2_constrained2)[2,])  %>% 
    rownames_to_column("model") %>%  
    as_tibble()
  
  if(!any(test_constraints_tier2$`Pr(>Chisq)`>crit, na.rm=T)){
    best_tier2 = list("t2_unconstrained" =fit_tier2)
  }else if(all(test_constraints_tier2$`Pr(>Chisq)`>crit, na.rm=T)){
    best_tier2 = list("t2_betas+resids"=fit_tier2_constrained2)
  }else {
    best_tier2= list("t2_betas"=fit_tier2_constrained1)
  }
  
  t2_mods <-list("t2_unconstrained"= fit_tier2,"t2_betas"= fit_tier2_constrained1,"t2_betas+resids"= fit_tier2_constrained2)
  
  message("Model-fitting with Tier 1+2 covariates complete.")
  ## TIER 1+2+3 COVS 
  # Repeat the process with tier 3 covs (maternal and paternal PGS) and all ix terms added
  model_tier3 <-
    paste0('# Regression 18m
           ',outname,'_18m ~ b1_1*mdrink_18m +  b2_1*',
           modname,'.pgs.pc.child + b3_1*mdrink_18m:',
           modname,'.pgs.pc.child + parity + b5_1*mdrink_18m:parity + b6_1*',modname,'.pgs.pc.child:parity + ',
           'sex + b8_1*mdrink_18m:sex  + b9_1*',modname,'.pgs.pc.child:sex + ',
           'mdrink_pre + b11_1*mdrink_18m:mdrink_pre + b12_1*',modname,'.pgs.pc.child:mdrink_pre + ',
           modname,'.pgs.pc.mother + b14_1*mdrink_18m:',modname,'.pgs.pc.mother + b15_1*',modname,'.pgs.pc.child:',modname,'.pgs.pc.mother +',
           modname,'.pgs.pc.father + b17_1*mdrink_18m:',modname,'.pgs.pc.father + b18_1*',modname,'.pgs.pc.child:',modname,'.pgs.pc.father 
           # Regression 3yr
           ',outname,'_3yr ~ b1_2*mdrink_3yr +  b2_2*',
           modname,'.pgs.pc.child + b3_2*mdrink_3yr:',
           modname,'.pgs.pc.child + parity + b5_2*mdrink_3yr:parity + b6_2*',modname,'.pgs.pc.child:parity + ',
           'sex + b8_2*mdrink_3yr:sex  + b9_2*',modname,'.pgs.pc.child:sex + ',
           'mdrink_pre + b11_2*mdrink_3yr:mdrink_pre + b12_2*',modname,'.pgs.pc.child:mdrink_pre + ',
           modname,'.pgs.pc.mother + b14_2*mdrink_3yr:',modname,'.pgs.pc.mother + b15_2*',modname,'.pgs.pc.child:',modname,'.pgs.pc.mother +',
           modname,'.pgs.pc.father + b17_2*mdrink_3yr:',modname,'.pgs.pc.father + b18_2*',modname,'.pgs.pc.child:',modname,'.pgs.pc.father 
           # Regression 5yr
           ',outname,'_5yr ~ b1_3*mdrink_5yr +  b2_3*',
           modname,'.pgs.pc.child + b3_3*mdrink_5yr:',
           modname,'.pgs.pc.child + parity + b5_3*mdrink_5yr:parity + b6_3*',modname,'.pgs.pc.child:parity + ',
           'sex + b8_3*mdrink_5yr:sex  + b9_3*',modname,'.pgs.pc.child:sex + ',
           'mdrink_pre + b11_3*mdrink_5yr:mdrink_pre + b12_3*',modname,'.pgs.pc.child:mdrink_pre + ',
           modname,'.pgs.pc.mother + b14_3*mdrink_5yr:',modname,'.pgs.pc.mother + b15_3*',modname,'.pgs.pc.child:',modname,'.pgs.pc.mother +',
           modname,'.pgs.pc.father + b17_3*mdrink_5yr:',modname,'.pgs.pc.father + b18_3*',modname,'.pgs.pc.child:',modname,'.pgs.pc.father 
           
           # Residual variances
           ',outname,'_18m ~~ outres1*',outname,'_18m
           ',outname,'_3yr ~~ outres2*',outname,'_3yr
           ',outname,'_5yr ~~ outres3*',outname,'_5yr
           
           # Covariances')
  
  fit_tier3 <-lavaan::sem(model_tier3,
                          data=data,
                          estimator="MLR",
                          missing="fiml.x",
                          se="robust",
                          cluster= "m_id")
  
  fit_tier3_constrained1 <-lavaan::sem(paste0(model_tier3,'
                                              # Constraints
                                               b1_1==b1_2
                                              b1_2==b1_3
                                              b2_1==b2_2
                                              b2_2==b2_3
                                              b3_1==b3_2
                                              b3_2==b3_3'),
                                       data=data,
                                       estimator="MLR",
                                       missing="fiml.x",
                                       se="robust",
                                       cluster= "m_id") 
  fit_tier3_constrained2 <-lavaan::sem(paste0(model_tier3,'
                                              # Constraints
                                               b1_1==b1_2
                                              b1_2==b1_3
                                              b2_1==b2_2
                                              b2_2==b2_3
                                              b3_1==b3_2
                                              b3_2==b3_3
                                              
                                              outres1==outres2
                                              outres2==outres3'),
                                       data=data,
                                       estimator="MLR",
                                       missing="fiml.x",
                                       se="robust",
                                       cluster= "m_id")
  
  test_constraints_tier3 <- rbind(anova(fit_tier3,fit_tier3_constrained1),
                            anova(fit_tier3,fit_tier3_constrained2)[2,])  %>% 
    rownames_to_column("model") %>%  
    as_tibble()
  
  if(!any(test_constraints_tier3$`Pr(>Chisq)`>crit, na.rm=T)){
    best_tier3 = list("t3_unconstrained"=fit_tier3)
  }else if(all(test_constraints_tier3$`Pr(>Chisq)`>crit, na.rm=T)){
    best_tier3 = list("t3_betas+resids"=fit_tier3_constrained2)
  }else {
    best_tier3= list("t3_betas"=fit_tier3_constrained1)
  }

  t3_mods <-list("t3_unconstrained"= fit_tier3,"t3_betas"= fit_tier3_constrained1,"t3_betas+resids"= fit_tier3_constrained2)
  
  message("All model-fitting for this outcome-moderator combination is complete.")
  #Create a function to extract standardized estimates in the desired format
  get_ests <- function(x){
    
    model_est <- standardizedsolution(x[[1]], type="std.nox") %>% 
      filter (op == "~" ) %>% 
      as_tibble() %>% 
      mutate(model=names(x)[[1]]) %>% 
      rbind(standardizedsolution(x[[2]], type="std.nox") %>% 
              filter (op == "~" ) %>% 
              as_tibble() %>% 
              mutate(model=names(x)[[2]]) ) %>% 
      rbind(standardizedsolution(x[[3]], type="std.nox") %>% 
              filter (op == "~" ) %>% 
              as_tibble() %>% 
              mutate(model=names(x)[[3]]))
    return(model_est)
  
     
  }
  toc() 
  #Return results in a nested list, with each level 1 element labelled according to the most recent tier
  #of covariates added, and containing the fit stats and estimates for the best model in that tier, as well
  #as the LRT results for the model fit comparisons
  
  return(list("tier1"=list("fit"=best_tier1, "comparison"=test_constraints_tier1, "res" = get_ests(t1_mods)  ),
              "tier2"=list("fit"=best_tier2, "comparison"=test_constraints_tier2, "res" = get_ests(t2_mods)  ),
              "tier3"=list("fit"=best_tier3,  "comparison"=test_constraints_tier3,"res" = get_ests(t3_mods)  )
  ))
  
 
}


#Function to make a main effects results table from the output

make_main_tab <- function(res){
  map(res, function(x){
    x$res %>% 
      filter(str_detect(label, "b1_|b2_"))
  }) %>% 
    bind_rows() %>% 
    select(outcome=lhs,model,predictor=rhs,est.std,se,pvalue,lci=ci.lower,uci=ci.upper)
}

#Function to make an interactions effects results table from the output

make_int_tab <- function(res){
  map(res, function(x){
    x$res %>% 
      filter(str_detect(label, "b3_"))
  }) %>% 
    bind_rows() %>% 
    select(outcome=lhs,model,predictor=rhs,est.std,se,pvalue,lci=ci.lower,uci=ci.upper)
}

#Function to make an overall fit table from the output

make_fit_tab <- function(res){
  map(res, function(x){
    x$comparison 
  }) %>% 
    bind_rows() 
}
