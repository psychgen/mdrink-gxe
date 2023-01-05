#04a_xpgs_funs.R

make_sstats <- function(pgs_names){
  
  for(pgs_name in pgs_names){
    
    overall_sstats <- all_sstats_new[[pgs_name]]
    
    if(!dir.exists(paste0("./data/sstats/",pgs_name))){ dir.create(paste0("./data/sstats/",pgs_name), recursive = T) } 
    
    for(out in c("int","ext")) {
      
      all<-  overall_sstats %>% 
        filter(str_detect(outcome, out))
      
      w1 <- all %>%  
        filter(str_detect(outcome, "18m"))
      
      w2 <- all %>%  
        filter(str_detect(outcome, "3yr"))
      
      w3 <- all %>%  
        filter(str_detect(outcome, "5yr"))
      
      all_avgd <- all %>% 
        group_by(SNP) %>% 
        summarise(P_ix = mean(P_ix, na.rm=T),
                  B_ix = mean(B_ix, na.rm=T)) %>% 
        left_join(all %>% 
                    select(-outcome,-P_ix,-B_ix)) %>% 
        distinct()
      
      all_mind <- all %>% 
        group_by(SNP) %>% 
        summarise(P_ix = min(P_ix, na.rm=T),
                  B_ix = ifelse(B_ix>0,max(B_ix, na.rm=T), min(B_ix,na.rm=T))) %>% 
        left_join(all %>% 
                    select(-outcome,-P_ix,-B_ix)) %>% 
        distinct()
      
      
      readr::write_tsv(w1, file= paste0("./data/sstats/",pgs_name,"/",pgs_name,"_",out,"_18m" ))
      readr::write_tsv(w2, file= paste0("./data/sstats/",pgs_name,"/",pgs_name,"_",out,"_3yr" ))
      readr::write_tsv(w3, file= paste0("./data/sstats/",pgs_name,"/",pgs_name,"_",out,"_5yr" ))
      readr::write_tsv(all_avgd, file= paste0("./data/sstats/",pgs_name,"/",pgs_name,"_",out,"_avg" ))
      readr::write_tsv(all_mind, file= paste0("./data/sstats/",pgs_name,"/",pgs_name,"_",out,"_min" ))
      
      
    }
    
   # if(!dir.exists(paste0("./scripts/bash/prsice/",pgs_name))){ dir.create(paste0("./scripts/bash/prsice/",pgs_name), recursive = T) } 
    
    input <- genotools::pgs_metadata %>% filter(Pheno_shortname==pgs_name) %>% 
      mutate(SNP = "SNP" ) #Accounts for renaming the neurot SNP column above
    
    # Use genotools to create PRSice scripts to make PGS based on these pvalues
    
    for(sstats in list.files(paste0("./data/sstats/",pgs_name) )){
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
                              A1 = input$A1,
                              A2 = input$A2,
                              stat= "B_ix",
                              pvalue= "P_ix",
                              snp= input$SNP,
                              chr = input$Chr,
                              bp = input$BP,
                              extract = input$extract,
                              thresholds = "5e-08,5e-07,5e-06,5e-05,5e-04,0.001,0.01,0.005,0.1,0.5,1",
                              clump_kb = "250",
                              clump_p = "1.000000",
                              clump_r2 = "0.100000",
                              lower = "5e-08",
                              maf = "0.01",
                              mhc= "exclude",
                              add_metadata = FALSE)
      
    }
  }
  

}



