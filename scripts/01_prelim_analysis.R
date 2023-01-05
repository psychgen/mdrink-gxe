#01_prelim_analysis.R

require(tidyverse)
require(lavaan)

load(file = './data/processed_data.RData')


names(alldata)

alldata <- alldata %>% 
  select(preg_id,BARN_NR,m_id,"sex"=KJONN_raw,"parity"=PARITET_5_raw, matches("cbcl"),matches("mdrink"),matches("pgs")) %>% 
  filter(sex %in% c(1,2)) %>% 
  mutate(sex=sex-1) %>%   #Males0 Females1
  mutate(sex=as.factor(sex),
         BARN_NR= as.factor(BARN_NR),
         across(where(is.numeric),scale),
         across(where(is.numeric),as.numeric),
         sex=as.numeric(sex),
         include= ifelse((!is.na(adhd.pgs.pc.child)&(!is.na(cbcl_ext_c_18m)|!is.na(mdrink_18m)|
                                                       !is.na(cbcl_ext_c_3yr)|!is.na(mdrink_3yr)|
                                                       !is.na(cbcl_ext_c_5yr)|!is.na(mdrink_5yr))),
                         1,0))

# Restrict to those with


# Note that Centering variables is essential 

# Creating interaction terms prior to sem running

allres_lm <- data.frame()
main_effs <- data.frame()
allres_lav <- data.frame()

for(pgs_pheno in c("neurot2018","ptsd2019","adhd","height2")){
  for(out in c("int","ext")){

alldata$pgs_child <- alldata[,paste0(pgs_pheno,".pgs.pc.child")] %>% unlist()
alldata$pgs_mother <- alldata[,paste0(pgs_pheno,".pgs.pc.mother")]%>% unlist()
alldata$pgs_father <- alldata[,paste0(pgs_pheno,".pgs.pc.father")]%>% unlist()
alldata$out18m <- alldata[,paste0("cbcl_",out,"_c_18m")]%>% unlist()
alldata$out3yr <- alldata[,paste0("cbcl_",out,"_c_3yr")]%>% unlist()
alldata$out5yr <- alldata[,paste0("cbcl_",out,"_c_5yr")]%>% unlist()


model18m_unadj <- lm(out18m ~ (mdrink_18m+pgs_child)^2, data=alldata,  )
model18m_tier1 <- lm(out18m ~ (mdrink_18m+pgs_child+sex+parity)^2, data=alldata )
model18m_tier2 <- lm(out18m ~ (mdrink_18m+pgs_child+sex+parity+pgs_mother+pgs_father)^2, data=alldata )
model3yr_unadj <- lm(out3yr ~ (mdrink_3yr+pgs_child)^2, data=alldata )
model3yr_tier1 <- lm(out3yr ~ (mdrink_3yr+pgs_child+sex+parity)^2, data=alldata )
model3yr_tier2 <- lm(out3yr ~ (mdrink_3yr+pgs_child+sex+parity+pgs_mother+pgs_father)^2, data=alldata )
model5yr_unadj <- lm(out5yr ~ (mdrink_5yr+pgs_child)^2, data=alldata )
model5yr_tier1 <- lm(out5yr ~ (mdrink_5yr+pgs_child+sex+parity)^2, data=alldata )
model5yr_tier2 <- lm(out5yr ~ (mdrink_5yr+pgs_child+sex+parity+pgs_mother+pgs_father)^2, data=alldata )


allmods <- list(model18m_unadj,model18m_tier1,model18m_tier2,
             model3yr_unadj,model3yr_tier1,model3yr_tier2,
             model5yr_unadj,model5yr_tier1,model5yr_tier2)

extract_main_effs <- function(x){
  broom::tidy(x) %>% 
    filter(!str_detect(term,":"),
           str_detect(term, "mdrink|pgs_child")) 
}
extract_ints <- function(x){
  broom::tidy(x) %>% 
    filter(str_detect(term,":")) %>% 
    .[1,]
}

res <- map(allmods,extract_ints) %>% 
  purrr::reduce(bind_rows) %>% 
  mutate(adj = rep(c("unadj","t1covs","t2covs"),3),
         pgs=pgs_pheno,
         outcome=out)
allres_lm <- rbind(allres_lm,res)

maineffstmp <- map(allmods,extract_main_effs) %>% 
  purrr::reduce(bind_rows) %>% 
  mutate(adj = rep(rep(c("unadj","t1covs","t2covs"),each=2),3),
         pgs=pgs_pheno,
         outcome=out)
main_effs <- rbind(main_effs, maineffstmp)


  }
}

allres_lm <- allres_lm %>% 
  mutate(fdr_p = p.adjust(p.value, method="fdr"))

main_effs <- main_effs %>% 
  mutate(fdr_p = p.adjust(p.value, method="fdr"))








alldata$mdrink_18m_BY_pgs_child <- alldata$mdrink_18m * alldata$pgs_child
alldata$mdrink_3yr_BY_pgs_child <- alldata$mdrink_3yr * alldata$pgs_child
alldata$mdrink_5yr_BY_pgs_child <- alldata$mdrink_5yr * alldata$pgs_child

mod_unadj <- '
 
  out18m ~ b11*mdrink_18m
  out18m ~ b21*pgs_child
  out18m ~ b3*mdrink_18m_BY_pgs_child
  
  out3yr ~ b11*mdrink_3yr
  out3yr ~ b21*pgs_child
  out3yr ~ b3*mdrink_3yr_BY_pgs_child
  
  out5yr ~ b11*mdrink_5yr
  out5yr ~ b21*pgs_child
  out5yr ~ b3*mdrink_5yr_BY_pgs_child
  
  '
alldata$mdrink_18m_BY_pgs_mother <- alldata$mdrink_18m * alldata$pgs_mother
alldata$mdrink_3yr_BY_pgs_mother <- alldata$mdrink_3yr * alldata$pgs_mother
alldata$mdrink_5yr_BY_pgs_mother <- alldata$mdrink_5yr * alldata$pgs_mother



  }
}






summary(model18m_unadj)

alldata <- alldata %>% 
  mutate(pgs_child_BY_mdrink_18m = myData$X * myData$M)

moderation_model <- '
 
  out18m ~ b1*mdrink_18m
  out18m ~ b2*pgs_child
  out18m ~ b3*mdrink_18m_BY_pgs_child
  
  out3yr ~ b1*mdrink_3yr
  out3yr ~ b2*pgs_child
  out3yr ~ b3*mdrink_3yr_BY_pgs_child
  
  out5yr ~ b1*mdrink_5yr
  out5yr ~ b2*pgs_child
  out5yr ~ b3*mdrink_5yr_BY_pgs_child
  
  '

model <- moderation_model 
sem_fit <- lavaan::sem( model, std.lv = TRUE, data = alldata, fixed.x = FALSE, cluster= "m_id"   )
summary( sem_fit, standardized = TRUE,  rsquare = TRUE )
lavaan::standardizedSolution( sem_fit, type = "std.all" )
