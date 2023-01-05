#07a_collate_funs.R

get_multilevel_output <- function(outname,modname,pgs_type="standard",...){
  
  modlist <- stringr::str_replace_all(paste(modname, str_remove_all(outname, "cbcl_|_c|"), c("effects","residuals","unconstrained"), sep= "_")," ","")
  
  mains<- data.frame()
  ints <- data.frame()
  fits <- data.frame()
  bws<- data.frame()
  
  for(res_file in modlist){  
    
    if(pgs_type=="standard"){
      stem="./output/mplus/"
    }else if(pgs_type=="xpgs"){
      stem="./output/mplus/xpgs/other/"
    }
    
    if(file.exists( paste0(stem,res_file,".out"))){
      mplusOutput_raw <- read_lines(paste0(stem,res_file,".out"))
    }else{ mplusOutput_raw=NULL}
    
    
    
    if(length(grep("MODEL RESULTS",mplusOutput_raw) )>0 ){
      
      #Replicate main_effs_tab structure
      maineffs <- mplusOutput_raw[ (min(grep("MODEL RESULTS",mplusOutput_raw))+10):
                                     min(grep("STANDARDIZED MODEL RESULTS",mplusOutput_raw))-3] %>% 
        as.data.frame()%>% 
        `colnames<-`(c("var")) %>%
        mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
        filter(str_detect(var, "MDRINK_1|MDRINK_2|MDRINK_3|PGS_C")) %>% 
        slice_head(n=6) %>% 
        separate(var, into= c("predictor", "est.std", "se", "est_se", "pvalue"), sep=" ") %>% 
        bind_cols(mplusOutput_raw[ (min(grep("CONFIDENCE INTERVALS OF MODEL RESULTS",mplusOutput_raw))+10):
                                     min(grep("MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES",mplusOutput_raw))-3] %>% 
                    as.data.frame()%>% 
                    `colnames<-`(c("var")) %>%
                    mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
                    filter(str_detect(var, "MDRINK_1|MDRINK_2|MDRINK_3|PGS_C")) %>% 
                    slice_head(n=6) %>% 
                    separate(var, into= c("predictor","lowest", "lci", "low", "est.std", "high", "uci","highest"), sep=" ") %>% 
                    select(  lci,uci)) %>% 
        as_tibble() %>% 
        mutate(predictor= case_when(predictor=="PGS_C" ~ paste0(modname,".pgs.pcs.child"),
                                    predictor=="MDRINK_1" ~ "mdrink_18m",
                                    predictor=="MDRINK_2" ~ "mdrink_3yr",
                                    predictor=="MDRINK_3" ~ "mdrink_5yr"),
               outcome = paste0(outname, rep(c("_18m","_3yr","_5yr"), each=2 )),
               model=paste0("ml_", unlist(str_split(res_file,"_"))[[ length( unlist(str_split(res_file,"_")) )]]),
               pvalue=ifelse(pvalue=="0.000", "0.0004", pvalue))  %>%
        mutate(across(.cols=c(est.std,se,est_se,pvalue,lci,uci), as.numeric)) %>% 
        select(outcome, model, predictor, est.std, se, pvalue, lci, uci)
      
      #Replicate int_effs_tab structure
      inteffs <- mplusOutput_raw[ (min(grep("MODEL RESULTS",mplusOutput_raw))+10):
                                    min(grep("STANDARDIZED MODEL RESULTS",mplusOutput_raw))-3] %>% 
        as.data.frame()%>% 
        `colnames<-`(c("var")) %>%
        mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
        filter(str_detect(var, "^M1|^M2|^M3")) %>% 
        slice_head(n=3) %>% 
        separate(var, into= c("predictor", "est.std", "se", "est_se", "pvalue"), sep=" ") %>% 
        bind_cols(mplusOutput_raw[ (min(grep("CONFIDENCE INTERVALS OF MODEL RESULTS",mplusOutput_raw))+10):
                                     min(grep("MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES",mplusOutput_raw))-3]  %>% 
                    as.data.frame()%>% 
                    `colnames<-`(c("var")) %>%
                    mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
                    filter(str_detect(var, "^M1|^M2|^M3")) %>% 
                    slice_head(n=3) %>% 
                    separate(var, into= c("predictor","lowest", "lci", "low", "est.std", "high", "uci","highest"), sep=" ") %>% 
                    select(  lci,uci)) %>% 
        as_tibble() %>% 
        mutate(predictor= case_when(predictor=="M1" ~ paste0("mdrink_18m:",modname,".pgs.pc.child"),
                                    predictor=="M2" ~ paste0("mdrink_3yr:",modname,".pgs.pc.child"),
                                    predictor=="M3" ~ paste0("mdrink_5yr:",modname,".pgs.pc.child")),
               outcome = paste0(outname, rep(c("_18m","_3yr","_5yr"), each=1 )),
               model=paste0("ml_", unlist(str_split(res_file,"_"))[[ length( unlist(str_split(res_file,"_")) )]]),
               pvalue=ifelse(pvalue=="0.000", "0.0004", pvalue))  %>%
        mutate(across(.cols=c(est.std,se,est_se,pvalue,lci,uci), as.numeric)) %>% 
        select(outcome, model, predictor, est.std, se, pvalue, lci, uci)
      
      #Replicate fit_comps_tab structure
      
      fitstats <-  tibble( model=paste0("ml_", unlist(str_split(res_file,"_"))[[ length( unlist(str_split(res_file,"_")) )]]),
                           Nfree = mplusOutput_raw[ (min(grep("Number of Free Parameters",mplusOutput_raw)))] %>% 
                             str_remove("Number of Free Parameters") %>% str_trim(),
                           AIC = mplusOutput_raw[ (min(grep("Akaike",mplusOutput_raw)))] %>% 
                             str_remove("Akaike \\(AIC\\)") %>% str_trim(),
                           BIC = mplusOutput_raw[ (min(grep("Bayesian",mplusOutput_raw)))] %>% 
                             str_remove("Bayesian \\(BIC\\)") %>% str_trim(),
                           LL = mplusOutput_raw[ (min(grep("H0 Value",mplusOutput_raw)))] %>% 
                             str_remove("H0 Value") %>% str_trim(),
                           LL_corr =mplusOutput_raw[ (min(grep("H0 Scaling Correction Factor",mplusOutput_raw)))] %>% 
                             str_remove("H0 Scaling Correction Factor") %>% str_trim()) %>% 
        mutate(across(.cols=c(Nfree,AIC,BIC,LL,LL_corr), as.numeric))
      
      
      #wbweffs
      bweffs <- mplusOutput_raw[ (grep("Between Level",mplusOutput_raw)[[1]]):
                                   (min(grep("STANDARDIZED MODEL RESULTS",mplusOutput_raw))-3)] %>% 
        as.data.frame()%>% 
        `colnames<-`(c("var")) %>%
        mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
        filter(str_detect(var, "INT|EXT")) %>% 
        slice_head(n=6) %>% 
        separate(var, into= c("outcome", "est.std", "se", "est_se", "pvalue"), sep=" ") %>% 
        bind_cols(mplusOutput_raw[ (grep("Between Level",mplusOutput_raw)[[4]]+2):
                                     min(grep("MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES",mplusOutput_raw))-3]  %>% 
                    as.data.frame()%>% 
                    `colnames<-`(c("var")) %>%
                    mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
                    filter(str_detect(var, "INT|EXT")) %>% 
                    slice_head(n=6) %>% 
                    separate(var, into= c("outcome","lowest", "lci", "low", "est.std", "high", "uci","highest"), sep=" ") %>% 
                    select(  lci,uci)) %>% 
        as_tibble()%>% 
        mutate(predictor= c(paste0("BTWmdrink", c("_18m","_3yr","_5yr")), 
                            paste0("BTW_exp_", modname,"_ix", c("_18m","_3yr","_5yr") )),
               outcome = rep(paste0(outname, c("_18m","_3yr","_5yr")), 2 ),
               model=paste0("ml_", unlist(str_split(res_file,"_"))[[ length( unlist(str_split(res_file,"_")) )]]),
               pvalue=ifelse(pvalue=="0.000", "0.0004", pvalue))  %>%
        mutate(across(.cols=c(est.std,se,est_se,pvalue,lci,uci), as.numeric)) %>% 
        select(outcome, model, predictor, est.std, se, pvalue, lci, uci)
      
      mains <- rbind(mains, maineffs)
      ints <- rbind(ints, inteffs)
      fits <- rbind(fits, fitstats)
      bws <- rbind(bws,bweffs)
    } else {
      maineffs <- tibble(outcome=character(), model=character(), predictor=character(), 
                         est.std=numeric(), se=numeric(), pvalue=numeric(), lci=numeric(), uci=numeric()) %>% 
        add_row( outcome = paste0(outname, rep(c("_18m","_3yr","_5yr"), each=2 )),
                 model=paste0("ml_", unlist(str_split(res_file,"_"))[[ length( unlist(str_split(res_file,"_")) )]]))
      inteffs <- maineffs %>% distinct() %>% 
        add_row( outcome = paste0(outname, rep(c("_18m","_3yr","_5yr"), each=1 )),
                 model=paste0("ml_", unlist(str_split(res_file,"_"))[[ length( unlist(str_split(res_file,"_")) )]]))
      fitstats <- tibble(model=character(),Nfree=numeric(),AIC=numeric(),BIC=numeric(),LL=numeric(),LL_corr=numeric()) %>% 
        add_row( model=paste0("ml_", unlist(str_split(res_file,"_"))[[ length( unlist(str_split(res_file,"_")) )]]))
      bweffs=maineffs
      
      mains <- rbind(mains, maineffs)
      ints <- rbind(ints, inteffs)
      fits <- rbind(fits, fitstats)
      bws <- rbind(bws,bweffs)
    }
    
  }
  
  #Manual LRT, where L1 is the baseline model (code from statmodel.com/chidiff.shtml)
  L1 = fits %>% filter(model=="ml_unconstrained") %>% .$LL
  c1 = fits %>% filter(model=="ml_unconstrained") %>% .$LL_corr
  p1 = fits %>% filter(model=="ml_unconstrained") %>% .$Nfree
  L0 = fits %>% filter(model=="ml_effects") %>% .$LL
  c0 = fits %>% filter(model=="ml_effects") %>% .$LL_corr
  p0 = fits %>% filter(model=="ml_effects") %>% .$Nfree
  L02 = fits %>% filter(model=="ml_residuals") %>% .$LL
  c02 = fits %>% filter(model=="ml_residuals") %>% .$LL_corr
  p02 = fits %>% filter(model=="ml_residuals") %>% .$Nfree
  
  cd01 = (p0 * c0 - p1*c1)/(p0 - p1)
  TRd01 = -2*(L0 - L1)/cd01
  pch1 = 1- pchisq(TRd01,p1-p0)
  
  
  if(is.na(L1)&!is.na(L0)){
    L1=L0
    c1=c0
    p1=p0
  }
  
  cd02 = (p02 * c02 - p1*c1)/(p02 - p1)
  TRd02 = -2*(L02 - L1)/cd02
  pch2 = 1- pchisq(TRd02,p1-p02)
  
  finalfits <- fits %>% 
    mutate(model= factor(model, levels= c("ml_unconstrained",
                                          "ml_effects",
                                          "ml_residuals"))) %>% 
    arrange(model) %>% 
    mutate(Chisq=NA,
           Df=NA,
           `Chisq diff` = c(NA,TRd01,TRd02),
           `Df diff` = c(NA,p1-p0,p1-p02 ),
           `Pr(>Chisq)`=c(NA,pch1,pch2)) %>% 
    select(model,Df,AIC,BIC, Chisq,`Chisq diff`,`Df diff`, `Pr(>Chisq)`)
  
  return( list("mains"=mains, "ints"= ints, "fitcomps"= finalfits, "bws"=bws))
}



best_fit_p <- function(pval){
  crit = 0.05/(2*4) 
  min <- min(pval,na.rm=T)
  max <- max(pval,na.rm=T)
  if(min<crit & max <crit){
    return(NA)
  }else if(min<crit & max >crit){
    return(max)
  }else{
    return(pval[[3]])
  }
  
}

get_multilevel_output_xpgs <- function(outname,modname,...){
  
  modlist <- paste(modname, str_remove_all(outname, "cbcl_|_c"), c("effects","residuals","unconstrained"), sep= "_")
  
  mains<- data.frame()
  ints <- data.frame()
  fits <- data.frame()
  
  for(res_file in modlist){  
    
    if(file.exists( paste0("./output/mplus/xpgs/",res_file,".out"))){
      mplusOutput_raw <- read_lines(paste0("./output/mplus/xpgs/",res_file,".out"))
    }else{ mplusOutput_raw=NULL}
    
    
    
    if(length(grep("STANDARDIZED MODEL RESULTS",mplusOutput_raw) )>0 ){
      #Replicate main_effs_tab structure
      maineffs <- mplusOutput_raw[ (min(grep("STANDARDIZED MODEL RESULTS",mplusOutput_raw))+15):
                                     min(grep("R-SQUARE",mplusOutput_raw))-3] %>% 
        as.data.frame()%>% 
        `colnames<-`(c("var")) %>%
        mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
        filter(str_detect(var, "MDRINK_1|MDRINK_2|MDRINK_3|PGS_C")) %>% 
        slice_head(n=6) %>% 
        separate(var, into= c("predictor", "est.std", "se", "est_se", "pvalue"), sep=" ") %>% 
        bind_cols(mplusOutput_raw[ (min(grep("CONFIDENCE INTERVALS OF STANDARDIZED MODEL RESULTS",mplusOutput_raw))+15):
                                     min(grep("MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES",mplusOutput_raw))-3] %>% 
                    as.data.frame()%>% 
                    `colnames<-`(c("var")) %>%
                    mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
                    filter(str_detect(var, "MDRINK_1|MDRINK_2|MDRINK_3|PGS_C")) %>% 
                    slice_head(n=6) %>% 
                    separate(var, into= c("predictor","lowest", "lci", "low", "est.std", "high", "uci","highest"), sep=" ") %>% 
                    select(  lci,uci)) %>% 
        as_tibble() %>% 
        mutate(predictor= case_when(predictor=="PGS_C" ~ paste0(modname,".pgs.pcs.child"),
                                    predictor=="MDRINK_1" ~ "mdrink_18m",
                                    predictor=="MDRINK_2" ~ "mdrink_3yr",
                                    predictor=="MDRINK_3" ~ "mdrink_5yr"),
               outcome = paste0(outname, rep(c("_18m","_3yr","_5yr"), each=2 )),
               model=paste0("ml_",unlist(str_split(res_file,"_"))[[length(unlist(str_split(res_file,"_")))]]),
               pvalue=ifelse(pvalue=="0.000", "0.0004", pvalue))  %>%
        mutate(across(.cols=c(est.std,se,est_se,pvalue,lci,uci), as.numeric)) %>% 
        select(outcome, model, predictor, est.std, se, pvalue, lci, uci)
      
      #Replicate int_effs_tab structure
      inteffs <- mplusOutput_raw[ (min(grep("STANDARDIZED MODEL RESULTS",mplusOutput_raw))+15):
                                    min(grep("R-SQUARE",mplusOutput_raw))-3] %>% 
        as.data.frame()%>% 
        `colnames<-`(c("var")) %>%
        mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
        filter(str_detect(var, "^M1|^M2|^M3")) %>% 
        slice_head(n=3) %>% 
        separate(var, into= c("predictor", "est.std", "se", "est_se", "pvalue"), sep=" ") %>% 
        bind_cols(mplusOutput_raw[ (min(grep("CONFIDENCE INTERVALS OF STANDARDIZED MODEL RESULTS",mplusOutput_raw))+15):
                                     min(grep("MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES",mplusOutput_raw))-3] %>% 
                    as.data.frame()%>% 
                    `colnames<-`(c("var")) %>%
                    mutate(var = gsub("\\s+", " ", str_trim(var))) %>%
                    filter(str_detect(var, "^M1|^M2|^M3")) %>% 
                    slice_head(n=3) %>% 
                    separate(var, into= c("predictor","lowest", "lci", "low", "est.std", "high", "uci","highest"), sep=" ") %>% 
                    select(  lci,uci)) %>% 
        as_tibble() %>% 
        mutate(predictor= case_when(predictor=="M1" ~ paste0("mdrink_18m:",modname,".pgs.pc.child"),
                                    predictor=="M2" ~ paste0("mdrink_3yr:",modname,".pgs.pc.child"),
                                    predictor=="M3" ~ paste0("mdrink_5yr:",modname,".pgs.pc.child")),
               outcome = paste0(outname, rep(c("_18m","_3yr","_5yr"), each=1 )),
               model=paste0("ml_",unlist(str_split(res_file,"_"))[[length(unlist(str_split(res_file,"_")))]]),
               pvalue=ifelse(pvalue=="0.000", "0.0004", pvalue))  %>%
        mutate(across(.cols=c(est.std,se,est_se,pvalue,lci,uci), as.numeric)) %>% 
        select(outcome, model, predictor, est.std, se, pvalue, lci, uci)
      
      #Replicate fit_comps_tab structure
      
      fitstats <-  tibble(model=paste0("ml_",unlist(str_split(res_file,"_"))[[length(unlist(str_split(res_file,"_")))]]),
                          Nfree = mplusOutput_raw[ (min(grep("Number of Free Parameters",mplusOutput_raw)))] %>% 
                            str_remove("Number of Free Parameters") %>% str_trim(),
                          AIC = mplusOutput_raw[ (min(grep("Akaike",mplusOutput_raw)))] %>% 
                            str_remove("Akaike \\(AIC\\)") %>% str_trim(),
                          BIC = mplusOutput_raw[ (min(grep("Bayesian",mplusOutput_raw)))] %>% 
                            str_remove("Bayesian \\(BIC\\)") %>% str_trim(),
                          LL = mplusOutput_raw[ (min(grep("H0 Value",mplusOutput_raw)))] %>% 
                            str_remove("H0 Value") %>% str_trim(),
                          LL_corr =mplusOutput_raw[ (min(grep("H0 Scaling Correction Factor",mplusOutput_raw)))] %>% 
                            str_remove("H0 Scaling Correction Factor") %>% str_trim()) %>% 
        mutate(across(.cols=c(Nfree,AIC,BIC,LL,LL_corr), as.numeric))
      
      mains <- rbind(mains, maineffs)
      ints <- rbind(ints, inteffs)
      fits <- rbind(fits, fitstats)
    } else {
      maineffs <- tibble(outcome=character(), model=character(), predictor=character(), 
                         est.std=numeric(), se=numeric(), pvalue=numeric(), lci=numeric(), uci=numeric()) %>% 
        add_row( outcome = paste0(outname, rep(c("_18m","_3yr","_5yr"), each=2 )),
                 model=paste0("ml_",unlist(str_split(res_file,"_"))[[length(unlist(str_split(res_file,"_")))]]) )
      inteffs <- maineffs %>% distinct() %>% 
        add_row( outcome = paste0(outname, rep(c("_18m","_3yr","_5yr"), each=1 )),
                 model=paste0("ml_",unlist(str_split(res_file,"_"))[[length(unlist(str_split(res_file,"_")))]]))
      fitstats <- tibble(model=character(),Nfree=numeric(),AIC=numeric(),BIC=numeric(),LL=numeric(),LL_corr=numeric()) %>% 
        add_row(model=paste0("ml_",unlist(str_split(res_file,"_"))[[length(unlist(str_split(res_file,"_")))]]))
      
      mains <- rbind(mains, maineffs)
      ints <- rbind(ints, inteffs)
      fits <- rbind(fits, fitstats)
    }
    
  }
  
  #Manual LRT, where L1 is the baseline model (code from statmodel.com/chidiff.shtml)
  L1 = fits %>% filter(model=="ml_unconstrained") %>% .$LL
  c1 = fits %>% filter(model=="ml_unconstrained") %>% .$LL_corr
  p1 = fits %>% filter(model=="ml_unconstrained") %>% .$Nfree
  L0 = fits %>% filter(model=="ml_effects") %>% .$LL
  c0 = fits %>% filter(model=="ml_effects") %>% .$LL_corr
  p0 = fits %>% filter(model=="ml_effects") %>% .$Nfree
  L02 = fits %>% filter(model=="ml_residuals") %>% .$LL
  c02 = fits %>% filter(model=="ml_residuals") %>% .$LL_corr
  p02 = fits %>% filter(model=="ml_residuals") %>% .$Nfree
  
  cd01 = (p0 * c0 - p1*c1)/(p0 - p1)
  TRd01 = -2*(L0 - L1)/cd01
  pch1 = 1- pchisq(TRd01,p1-p0)
  
  cd02 = (p02 * c02 - p1*c1)/(p02 - p1)
  TRd02 = -2*(L02 - L1)/cd02
  pch2 = 1- pchisq(TRd02,p1-p02)
  
  finalfits <- fits %>% 
    mutate(model= factor(model, levels= c("ml_unconstrained",
                                          "ml_effects",
                                          "ml_residuals"))) %>% 
    arrange(model) %>% 
    mutate(Chisq=NA,
           Df=NA,
           `Chisq diff` = c(NA,TRd01,TRd02),
           `Df diff` = c(NA,p1-p0,p1-p02 ),
           `Pr(>Chisq)`=c(NA,pch1,pch2)) %>% 
    select(model,Df,AIC,BIC, Chisq,`Chisq diff`,`Df diff`, `Pr(>Chisq)`)
  
  return( list("mains"=mains, "ints"= ints, "fitcomps"= finalfits))
}

