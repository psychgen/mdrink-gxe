#07_collate_plot_res.R

#The purpose of this script is to collate and plot results from all models,
#across covariate tiers, PGS, outcomes, and filtering by accepted constraints

library(ggthemes)
library(tidyverse)
library(MplusAutomation)
library(tictoc)
library(patchwork)


#STANDARD PGS
#Linear model results
load(file="./output/linear_mods_res.RData")
linear_mods_res <- as_tibble(linear_mods_res)

#Replicate structure for multilevel model results

#Raw text approach far faster than MplusAutomation readModels to will use that

source("./scripts/07a_collate_funs.R")

all_mods_res <- linear_mods_res %>% 
  mutate( mlresults = pmap(., get_multilevel_output, pgs_type="standard"),
          mlmain_effs_tab = map(mlresults, function(x) x[["mains"]]),
          mlint_effs_tab = map(mlresults, function(x) x[["ints"]]),
          mlfit_comp_tab = map(mlresults, function(x) x[["fitcomps"]]),
          mlbtw_effs_tab = map(mlresults, function(x) x[["bws"]]))

save(all_mods_res, file="./output/all_mods_res.RData")

# Combine fit_comps_tab tables and use to create "most parsimonious" and "highest
# adjustment, most parsimonious" lkp tables for main effects and interacitons respectively

load("./output/all_mods_res.RData")

fit_tabs <- do.call(rbind, all_mods_res$fit_comp_tab) %>% 
  bind_rows(do.call(rbind, all_mods_res$mlfit_comp_tab) ) %>% 
  as_tibble() %>%
  mutate(out= c(rep(all_mods_res$out, each= nrow(all_mods_res$fit_comp_tab[[1]]) ),
                rep(all_mods_res$out, each= nrow(all_mods_res$mlfit_comp_tab[[1]]))),
         mod= c(rep(all_mods_res$mod, each= nrow(all_mods_res$fit_comp_tab[[1]]) ),
                rep(all_mods_res$mod, each= nrow(all_mods_res$mlfit_comp_tab[[1]]))),
         model=case_when(model=="fit_tier1"  ~ "t1_unconstrained",
                         model=="fit_tier1_constrained1"  ~ "t1_betas",
                         model=="fit_tier1_constrained2"  ~ "t1_betas+resids",
                         model=="fit_tier2"  ~ "t2_unconstrained",
                         model=="fit_tier2_constrained1"  ~ "t2_betas",
                         model=="fit_tier2_constrained2"  ~ "t2_betas+resids",
                         model=="fit_tier3"  ~ "t3_unconstrained",
                         model=="fit_tier3_constrained1"  ~ "t3_betas",
                         model=="fit_tier3_constrained2"  ~ "t3_betas+resids",
                         model=="ml_effects"  ~ "ml_betas",
                         model=="ml_residuals"  ~ "ml_betas+resids",
                         TRUE~ as.character(model)
         )) %>% 
  separate(model, into =c("control","constraints"),sep="_")

save(fit_tabs,file="./output/all_fit_tabs.RData")

most_pars <- fit_tabs %>%  
  group_by(out,mod,control) %>% 
  summarise(p_fit = best_fit_p(`Pr(>Chisq)`)) %>% 
  ungroup() %>% 
  left_join(fit_tabs %>% 
              select(out, mod, control,constraints,p_fit=`Pr(>Chisq)`)) %>% 
  filter(!(out=="cbcl_int_c"&mod=="neurot2018"&control=="ml"&constraints=="betas") ) #Manually filter out the model that did not converge

most_pars_most_adj <- most_pars %>%
  filter(control=="ml")


## Attach results to lookup tables for figures 2/3/4

# Figure 2: Attenuation of main effects (all adjustments, most parsimonious)

all_mains <- do.call(rbind, all_mods_res$main_effs_tab) %>% 
  bind_rows(do.call(rbind, all_mods_res$mlmain_effs_tab) ) %>% 
  as_tibble()%>% 
  mutate(out= c(rep(all_mods_res$out, each= nrow(all_mods_res$main_effs_tab[[1]]) ),
                rep(all_mods_res$out, each= nrow(all_mods_res$mlmain_effs_tab[[1]]))),
         mod= c(rep(all_mods_res$mod, each= nrow(all_mods_res$main_effs_tab[[1]]) ),
                rep(all_mods_res$mod, each= nrow(all_mods_res$mlmain_effs_tab[[1]]))),
         wave = str_remove_all(outcome, paste0(paste0(unique(most_pars$out),"_"),collapse="|")),
         model=case_when(model=="ml_effects"  ~ "ml_betas",
                         model=="ml_residuals"  ~ "ml_betas+resids",
                         TRUE~ as.character(model)
         )) %>% 
  separate(model, into =c("control","constraints"),sep="_") 

most_pars_mains <- most_pars %>% 
  left_join(all_mains) %>% 
  group_by(out, mod) %>% 
  filter(if(all(constraints!="unconstrained")) str_detect(outcome, "18m") else TRUE) %>% 
  ungroup() %>% 
  mutate(wave=ifelse(constraints=="unconstrained", wave, "all")) %>% 
  mutate(control=factor(control, levels=c("t1","t2","t3","ml"), labels=c("T1","T2","T3","ML"))) %>% 
  mutate(Predictor= factor(case_when(str_detect(predictor, "pgs")&mod=="adhd" ~ "ADHD[PGS]",
                                     str_detect(predictor, "pgs")&mod=="neurot2018" ~ "Neuroticism[PGS]",
                                     str_detect(predictor, "pgs")&mod=="ptsd2019" ~ "PTSD[PGS]",
                                     str_detect(predictor, "pgs")&mod=="height2" ~ "Height[PGS]",
                                     str_detect(predictor, "mdrink")&mod=="height2" ~ "Mat.drinking"),
                           levels= c("Mat.drinking", "ADHD[PGS]","Neuroticism[PGS]","PTSD[PGS]","Height[PGS]")))%>% 
  mutate(mod=factor(mod,levels=c("adhd","neurot2018", "ptsd2019","height2"))) %>% 
  drop_na(Predictor) %>% 
  droplevels()

save(all_mains,most_pars_mains, file="./output/all_main_effects.RData")

# New facet label names PGS variable
mod_labs <- c("Neuroticism", "PTSD","ADHD", "Height")
names(mod_labs) <- c("neurot2018", "ptsd2019","adhd","height2")
# New facet label names out variable
out_labs <- c("Emotional", "Behavioural")
names(out_labs) <- c("cbcl_int_c", "cbcl_ext_c")
pd=0.6

ggplot(most_pars_mains,aes(x=control,y=est.std, fill=wave,colour=wave))+
  geom_hline(aes(yintercept=0),linetype=2,size=0.6,colour="grey30")+
  geom_errorbar(aes(ymin=lci,ymax=uci), size=1, width=0, alpha=0.6, position =position_dodge(pd))+
  geom_point(size=3, position =position_dodge(pd),shape=21, alpha=1,stroke=1,colour="grey30")+
  facet_grid(out~Predictor, 
             labeller = labeller(out = out_labs, Predictor = label_parsed))+
  scale_fill_viridis_d("Wave",option="G",  end=0.8) +
  scale_colour_viridis_d("Wave",option="G",  end=0.8) +
  theme_bw()+
  theme(text=element_text(size=12),
        axis.title = element_text(),
        legend.key.width = unit(1.5, 'cm'),
        legend.text = element_text(size=10),
        legend.title =element_text(size=10),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -0.75)
        #plot.background = element_rect(fill="white")
  )+
  ylab("Main effect in linear regression (standardised)")+
  xlab("Adjustment for confounding \n(least to most adjusted)")

ggsave("./output/plots/fig2.tiff", device="tiff", width=18,height=14,units="cm",dpi=320,bg="white")


# Figure 3a: Interaction effects (most adjusted, most parsimonious)

all_ints <-do.call(rbind, all_mods_res$int_effs_tab) %>% 
  bind_rows(do.call(rbind, all_mods_res$mlint_effs_tab) ) %>% 
  as_tibble()%>% 
  distinct() %>% 
  mutate(out= c(rep(all_mods_res$out, each= nrow(all_mods_res$int_effs_tab[[1]]) ),
                rep(all_mods_res$out, each= nrow(all_mods_res$mlint_effs_tab[[1]]))),
         mod= c(rep(all_mods_res$mod, each= nrow(all_mods_res$int_effs_tab[[1]]) ),
                rep(all_mods_res$mod, each= nrow(all_mods_res$mlint_effs_tab[[1]]))),
         wave = str_remove_all(outcome, paste0(paste0(unique(most_pars$out),"_"),collapse="|")),
         model=case_when(model=="ml_effects"  ~ "ml_betas",
                         model=="ml_residuals"  ~ "ml_betas+resids",
                         TRUE~ as.character(model)
         )) %>% 
  separate(model, into =c("control","constraints"),sep="_") 

ints_final <- most_pars_most_adj %>%
  left_join(all_ints) %>% 
  group_by(out, mod) %>% 
  filter(if(all(constraints!="unconstrained")) str_detect(outcome, "18m") else TRUE) %>% 
  ungroup() %>% 
  mutate(wave=ifelse(constraints=="unconstrained", wave, "all"))

save(all_ints,ints_final, file= "./output/all_int_effects.RData")

inc= 0.05
plot_ints <- tibble( out = ints_final$out %>% rep(each=4/inc+1),
                     mod= ints_final$mod  %>% rep(each=4/inc+1),
                     wave= ints_final$wave %>% rep(each=4/inc+1),
                     pgs= rep(seq(-2,2, by=inc ),nrow(ints_final)),
                     b3= ints_final$est.std %>% rep(each=4/inc+1),
                     b3uci= ints_final$uci %>% rep(each=4/inc+1),
                     b3lci=ints_final$lci %>% rep(each=4/inc+1),
                     b1= ints_final %>% 
                       left_join(all_mains %>% 
                                   filter(control=="ml", str_detect(predictor,"mdrink")) %>% 
                                   mutate(wave=ifelse(constraints!="unconstrained"&wave=="18m","all",wave)) %>%  
                                   select(out,mod,control, constraints, wave, "b1"=est.std)) %>% 
                       .$b1 %>% rep(each=4/inc+1),
                     y= b1 + b3 * pgs,
                     y_lci= b1 + b3lci * pgs,
                     y_uci= b1 + b3uci * pgs,
                     y_sesoi_hi= b1 + 0 * pgs,
                     y_sesoi_lo= b1 + -0 * pgs) %>% 
  mutate(mod=factor(mod,levels=c("adhd","neurot2018", "ptsd2019","height2"),
                    labels=c("ADHD","Neuroticism", "PTSD","Height")))


linedat = plot_ints %>% 
  select(pgs,y_sesoi_hi,wave,out,mod) %>% 
  mutate(lt="hypothetical") %>% 
  add_row(lt="estimated", wave="all", out="cbcl_ext_c",mod="ADHD") %>%
  mutate(lt=factor(lt, levels=c("estimated","hypothetical")),
         mod= factor(mod,levels=c("ADHD","Neuroticism", "PTSD","Height"),
                     labels=c("ADHD","Neuroticism", "PTSD","Height") ))

p1 <- ggplot(plot_ints, aes(x=pgs,y=y,colour=wave,fill=wave)) + 
  geom_ribbon(aes(ymin=y_lci, ymax=y_uci, group=wave), alpha=0.3, colour="transparent" ) + 
  geom_line(data=linedat,aes( y=y_sesoi_hi,linetype=lt), size=0.7) + 
  scale_linetype_manual("Change in exposure effect\nacross PGS scores with...",values=c(1,2) , 
                        breaks= levels(linedat$lt),
                        labels=c(expression(paste(beta[3]," as estimated")),
                                 expression(paste(beta[3] == 0,"    "))))+
  geom_line(size=0.8, alpha=0.6)+
  facet_grid(out~mod, 
             labeller = labeller(out = out_labs))+
  coord_cartesian(ylim=c(-0.05,0.1))+
  scale_color_viridis_d("Wave", option="G", end=0.8) +
  scale_fill_viridis_d("Wave",option="G", end=0.8) +
  theme_bw()+
  theme(text=element_text(size=12),
        axis.title = element_text(),
        legend.key.width = unit(1.5, 'cm'),
        legend.text = element_text(size=10),
        legend.title =element_text(size=10),
        #legend.position="right",
        #legend.direction="vertical",
        #plot.background = element_rect(fill="white")
  )+
  xlab("Polygenic score (standardised)")+
  ylab("Standardised effect of matrenal \nat-risk drinking on outcome")

p1
ggsave("./output/plots/fig3.tiff", device="tiff", width=18,height=14,units="cm",dpi=320,bg="white")


# Table 2: within and between effects from best ML models
# Need the best, most adjusted model results from each mains, ints, and mlbtw

all_btweffs <- do.call(rbind, all_mods_res$mlbtw_effs_tab) %>% 
  as_tibble()%>% 
  mutate(out= rep(all_mods_res$out, each= nrow(all_mods_res$mlbtw_effs_tab[[1]]) ),
         mod= rep(all_mods_res$mod, each= nrow(all_mods_res$mlbtw_effs_tab[[1]]) ),
         wave = str_remove_all(outcome, paste0(paste0(unique(most_pars$out),"_"),collapse="|")),
         model=case_when(model=="ml_effects"  ~ "ml_betas",
                         model=="ml_residuals"  ~ "ml_betas+resids",
                         TRUE~ as.character(model)
         )) %>% 
  separate(model, into =c("control","constraints"),sep="_") 

slctd_maineffs <- most_pars_most_adj %>%
  left_join(all_mains) %>% 
  group_by(out, mod) %>% 
  filter(if(all(constraints!="unconstrained")) str_detect(outcome, "18m") else TRUE) %>% 
  ungroup() %>% 
  mutate(wave=ifelse(constraints=="unconstrained", wave, "all")) %>% 
  filter(str_detect(predictor, "mdrink"))

slctd_btweffs<- most_pars_most_adj %>%
  left_join(all_btweffs) %>% 
  group_by(out, mod) %>% 
  filter(if(all(constraints!="unconstrained")) str_detect(outcome, "18m") else TRUE) %>% 
  ungroup() %>% 
  mutate(wave=ifelse(constraints=="unconstrained", wave, "all"),
         effect=ifelse(str_detect(predictor, "mdrink"),"Main", "Interaction")) %>% 
  group_by(effect) %>% 
  mutate(fdr_p= p.adjust(pvalue, method="fdr")) %>% ungroup() %>% 
  arrange(effect,fdr_p)

w_btw_effs <- slctd_maineffs %>% 
  mutate(effect="Main",
         level="Within",
         fdr_p= p.adjust(pvalue, method="fdr")) %>% 
  bind_rows(ints_final %>% 
              mutate(effect="Interaction",
                     level="Within",
                     fdr_p= p.adjust(pvalue, method="fdr"))) %>% 
  bind_rows(slctd_btweffs %>% 
              mutate(level="Between")) %>% 
  arrange(out,mod,wave,effect,level) %>% 
  mutate(mod=factor(mod,levels=c("adhd","neurot2018", "ptsd2019","height2"),
                    labels=c("ADHD","Neuroticism", "PTSD","Height")),
         out=factor(out,levels=c("cbcl_int_c","cbcl_ext_c"),labels=c("Emotional","Behavioural"))) %>% 
  select(out,mod,wave,effect,level, est.std,fdr_p,lci,uci)

save(w_btw_effs, file="./output/within_between_effects.RData")


# Single SNP mod results
load(file="./output/single_snp_ix_model_results.RData")




snp_lists <- genotools::get_pgs_snps(c("allexclheight_ext"),
                                     maf = "0.01",
                                     clump = "250_1_0.1")

# How many tests?
n_tests <-nrow(snp_lists$allexclheight_ext)

top_snps<- snp_res %>% 
  #  mutate(out_short= ifelse(str_detect(outcome,"_int_"), "int","ext")) %>% 
  #  group_by(pgs, out_short, ID) %>% 
  #  summarise(B_ix = mean(BETA, na.rm=T),
  #            P_ix = mean(P_ix, na.rm=T)) %>% 
  #  left_join(snp_res %>% 
  select(pgs,outcome,"B_ix"=BETA,P_ix, ID,CHROM,POS,A1) %>% 
  #              distinct()) %>% 
  filter(pgs!="height2") %>% 
  arrange(P_ix) %>% head(15) 



save(top_snps, file="./output/top_snps.RData")



#XPGS
#Linear model results
load(file="./output/linear_mods_res_xpgs.RData")
linear_mods_res_xpgs <- as_tibble(linear_mods_res_xpgs) %>% 
  filter(str_detect(mod,"allexcl"))

source("./scripts/07a_collate_funs.R")
#Replicate structure for multilevel model results

all_mods_res_xpgs <- linear_mods_res_xpgs %>% 
  mutate( mlresults = pmap(., get_multilevel_output, pgs_type="xpgs" ),
          mlmain_effs_tab = map(mlresults, function(x) x[["mains"]]),
          mlint_effs_tab = map(mlresults, function(x) x[["ints"]]),
          mlfit_comp_tab = map(mlresults, function(x) x[["fitcomps"]]),
          mlbtw_effs_tab = map(mlresults, function(x) x[["bws"]])) %>% 
  as_tibble()

save(all_mods_res_xpgs, file="./output/all_mods_res_xpgs.RData")

# Combine fit_comps_tab tables and use to create "most parsimonious" and "highest
# adjustment, most parsimonious" lkp tables for main effects and interacitons respectively


load("./output/all_mods_res_xpgs.RData")

fit_tabs_xpgs <- do.call(rbind, all_mods_res_xpgs$fit_comp_tab) %>% 
  bind_rows(do.call(rbind, all_mods_res_xpgs$mlfit_comp_tab) ) %>% 
  as_tibble() %>%
  mutate(out= c(rep(all_mods_res_xpgs$out, each= nrow(all_mods_res_xpgs$fit_comp_tab[[1]]) ),
                rep(all_mods_res_xpgs$out, each= nrow(all_mods_res_xpgs$mlfit_comp_tab[[1]]))),
         mod= c(rep(all_mods_res_xpgs$mod, each= nrow(all_mods_res_xpgs$fit_comp_tab[[1]]) ),
                rep(all_mods_res_xpgs$mod, each= nrow(all_mods_res_xpgs$mlfit_comp_tab[[1]]))),
         model=case_when(model=="fit_tier1"  ~ "t1_unconstrained",
                         model=="fit_tier1_constrained1"  ~ "t1_betas",
                         model=="fit_tier1_constrained2"  ~ "t1_betas+resids",
                         model=="fit_tier2"  ~ "t2_unconstrained",
                         model=="fit_tier2_constrained1"  ~ "t2_betas",
                         model=="fit_tier2_constrained2"  ~ "t2_betas+resids",
                         model=="fit_tier3"  ~ "t3_unconstrained",
                         model=="fit_tier3_constrained1"  ~ "t3_betas",
                         model=="fit_tier3_constrained2"  ~ "t3_betas+resids",
                         model=="ml_effects"  ~ "ml_betas",
                         model=="ml_residuals"  ~ "ml_betas+resids",
                         TRUE~ as.character(model)
         )) %>% 
  separate(model, into =c("control","constraints"),sep="_")

save(fit_tabs_xpgs,file="./output/all_fit_tabs_xpgs.RData")

most_pars_xpgs <- fit_tabs_xpgs %>%  
  group_by(out,mod,control) %>% 
  summarise(p_fit = best_fit_p(`Pr(>Chisq)`)) %>% 
  ungroup() %>% 
  left_join(fit_tabs_xpgs %>% 
              select(out, mod, control,constraints,p_fit=`Pr(>Chisq)`))

most_pars_most_adj_xpgs <- most_pars_xpgs %>%
  filter(control=="ml")
# Figure 4: xPGS interaction effects (most adjusted, most parsimonious)

all_mains_xpgs <- do.call(rbind, all_mods_res_xpgs$main_effs_tab) %>% 
  bind_rows(do.call(rbind, all_mods_res_xpgs$mlmain_effs_tab) ) %>% 
  as_tibble()%>% 
  mutate(out= c(rep(all_mods_res_xpgs$out, each= nrow(all_mods_res_xpgs$main_effs_tab[[1]]) ),
                rep(all_mods_res_xpgs$out, each= nrow(all_mods_res_xpgs$mlmain_effs_tab[[1]]))),
         mod= c(rep(all_mods_res_xpgs$mod, each= nrow(all_mods_res_xpgs$main_effs_tab[[1]]) ),
                rep(all_mods_res_xpgs$mod, each= nrow(all_mods_res_xpgs$mlmain_effs_tab[[1]]))),
         wave = str_remove_all(outcome, paste0(paste0(unique(most_pars_xpgs$out),"_"),collapse="|")),
         model=case_when(model=="ml_effects"  ~ "ml_betas",
                         model=="ml_residuals"  ~ "ml_betas+resids",
                         TRUE~ as.character(model)
         )) %>% 
  separate(model, into =c("control","constraints"),sep="_")



ml_ints_xpgs <-do.call(rbind, map(all_mods_res_xpgs$mlint_effs_tab,distinct))  %>% 
  as_tibble()%>% 
  mutate(out= rep(all_mods_res_xpgs$out, each= nrow(unique(all_mods_res_xpgs$mlint_effs_tab[[1]])) ),
         mod= rep(all_mods_res_xpgs$mod, each= nrow(unique(all_mods_res_xpgs$mlint_effs_tab[[1]])) ),
         wave = str_remove_all(outcome, paste0(paste0(unique(most_pars_xpgs$out),"_"),collapse="|")),
         model=case_when(model=="ml_effects"  ~ "ml_betas",
                         model=="ml_residuals"  ~ "ml_betas+resids",
                         TRUE~ as.character(model)
         )) %>% 
  separate(model, into =c("control","constraints"),sep="_") 

ints_final_xpgs <- most_pars_most_adj_xpgs %>%
  left_join(ml_ints_xpgs) %>% 
  group_by(out, mod) %>% 
  filter(if(all(constraints!="unconstrained")) str_detect(outcome, "18m") else TRUE) %>% 
  ungroup() %>% 
  mutate(wave=ifelse(constraints=="unconstrained", wave, "all"),
         fdr_p=p.adjust(pvalue, method="fdr")) %>% 
  select(out,mod,wave,constraints,est.std,fdr_p,pvalue, lci,uci) %>% 
  drop_na(est.std)

save(ints_final_xpgs, file= "./output/all_int_effects_xpgs.RData")
# New facet label names PGS variable
mod_labs <- c("All SNPs (Beh. probs.)", "All SNPs (Emo. probs.)")
names(mod_labs) <- c("allexclheight_ext", "allexclheight_int")
# New facet label names out variable
out_labs <- c("Emotional", "Behavioural")
names(out_labs) <- c("cbcl_int_c", "cbcl_ext_c")
pd=0.6
inc= 0.05
plot_ints_xpgs <- tibble( out = ints_final_xpgs$out %>% rep(each=4/inc+1),
                          mod= ints_final_xpgs$mod  %>% rep(each=4/inc+1),
                          wave= ints_final_xpgs$wave %>% rep(each=4/inc+1),
                          pgs= rep(seq(-2,2, by=inc ),nrow(ints_final_xpgs)),
                          b3= ints_final_xpgs$est.std %>% rep(each=4/inc+1),
                          b3uci= ints_final_xpgs$uci %>% rep(each=4/inc+1),
                          b3lci=ints_final_xpgs$lci %>% rep(each=4/inc+1),
                          b1= ints_final_xpgs %>% 
                            left_join(all_mains_xpgs %>% 
                                        filter(control=="ml", str_detect(predictor,"mdrink")) %>% 
                                        mutate(wave=ifelse(constraints!="unconstrained"&wave=="18m","all",wave)) %>%  
                                        select(out,mod,control, constraints, wave, "b1"=est.std)) %>% 
                            .$b1 %>% rep(each=4/inc+1),
                          y= b1 + b3 * pgs,
                          y_lci= b1 + b3lci * pgs,
                          y_uci= b1 + b3uci * pgs,
                          y_sesoi_hi= b1 + 0 * pgs,
                          y_sesoi_lo= b1 + 0 * pgs) %>% 
  mutate(mod=factor(mod,levels=names(mod_labs),
                    labels=mod_labs))





linedat_xpgs = plot_ints_xpgs %>% 
  select(pgs,y_sesoi_hi,wave,out,mod) %>% 
  mutate(lt="hypothetical") %>% 
  add_row(lt="estimated", wave="all", out="cbcl_ext_c",mod="All SNPs (Beh. probs.)") %>%
  mutate(lt=factor(lt, levels=c("estimated","hypothetical")),
         mod=factor(mod,levels=mod_labs,
                    labels=mod_labs))

p1 <- ggplot(plot_ints_xpgs, aes(x=pgs,y=y,colour=wave,fill=wave)) + 
  geom_ribbon(aes(ymin=y_lci, ymax=y_uci, group=wave), alpha=0.3, colour="transparent" ) + 
  geom_line(data=linedat_xpgs,aes( y=y_sesoi_hi,linetype=lt), size=0.7) + 
  scale_linetype_manual("Change in exposure effect\nacross xPGS scores with...",values=c(1,2) , 
                        breaks= levels(linedat_xpgs$lt),
                        labels=c(expression(paste(beta[3]," as estimated")),
                                 expression(paste(beta[3] == 0,"    "))))+
  geom_line(size=0.8)+
  facet_grid(.~mod, 
             labeller = labeller(out = out_labs))+
  coord_cartesian(ylim=c(-0.25,0.25))+
  scale_color_viridis_d("Wave", option="G", end=0.8) +
  scale_fill_viridis_d("Wave",option="G", end=0.8) +
  theme_fivethirtyeight()+
  theme(text=element_text(size=12),
        axis.title = element_text(),
        legend.key.width = unit(1.5, 'cm'),
        legend.text = element_text(size=10),
        legend.title =element_text(size=10),
        #legend.position="right",
        #legend.direction="vertical",
        #plot.background = element_rect(fill="white")
  )+
  xlab("Interaction polygenic score (xPGS; standardised)")+
  ylab("Standardised effect of matrenal \nat-risk drinking on outcome")

p1
ggsave("./output/plots/figXxPGS.tiff", device="tiff", width=18,height=14,units="cm",dpi=320,bg="white")

