library(dplyr)
library(pbapply)
##GSE174498_CEM_L-asparaginase 没有输出结果
# sample_meta <- readRDS("~/Drug_splicing/data/all_sample_meta.rds")
# sample_meta <- sample_meta %>% filter(unid != "GSE174498_CEM_L-asparaginase")
rmats_counts <- readRDS("~/Drug_splicing/data/rmats_filter.rds")
rmats_counts <- rmats_counts %>% filter(sample != "GSE174498_CEM_L-asparaginase")
all_samples <- unique(rmats_counts$sample)
dones <- list.files("~/Drug_splicing/scripts/Shiny/data/",pattern = "binding_pep.rds")
dones_meta <- list.files("~/Drug_splicing/scripts/Shiny/data/",pattern = "binding_pep_meta.rds")
dones <- gsub("_binding_pep.rds","",dones)
dones_meta <- gsub("_binding_pep_meta.rds","",dones_meta)
dones <- intersect(dones, dones_meta)
all_samples <- all_samples[which(!(all_samples %in% dones))]
samples <- readRDS("~/Drug_splicing/data/samples_with_HLA.rds")
all_samples <- intersect(all_samples, samples$unid)

for (i in 1:length(all_samples)){
  all_files <- list.files(paste0("~/GEO_data/SRA/gse_out/",all_samples[i],
                                 "/jcast_out/binding_out/res/"),pattern = "tsv",
                          full.names = T)
  if (length(all_files) > 20){
    cl <- parallel::makeCluster(20)
    dt <- pbsapply(all_files, 
                   function(x){
                     tmp <- data.table::fread(x,data.table = F)
                     colnames(tmp)[5] <- "score"
                     tmp$binding_type <- gsub("_.+","", gsub(".+//","",x))
                     return(tmp)
                   }, cl = cl, simplify = F)
    parallel::stopCluster(cl)
  }else{
    dt <- sapply(all_files,
                 function(x){
                   tmp <- data.table::fread(x,data.table = F)
                   colnames(tmp)[5] <- "score"
                   tmp$binding_type <- gsub(".+//","",x) %>% gsub("_.+","",.)
                   return(tmp)
                 },simplify = F)
  }

  dt <- bind_rows(dt) %>% distinct_all()
  dt_filter <- dt %>% 
    select(allele, peptide, percentile, binding_type) %>% 
    filter(percentile < 2) %>% distinct_all()
  ###peptide binding
  dt1 <- dt_filter %>% 
    tidyr::pivot_wider(names_from = c(allele, binding_type), values_from = percentile)
  dt_all <- dt %>% 
    select(allele, peptide, percentile, binding_type) %>% 
    filter(peptide %in% dt1$peptide) %>% 
    tidyr::pivot_wider(names_from = c(allele, binding_type), values_from = percentile)
  bc <- apply(dt_all[,2:ncol(dt_all)],1,function(x){sum(x < 2, na.rm = T)})
  dt_all$binding_counts <- bc
  dt_all <- dt_all %>% select(peptide, binding_counts, everything()) %>% 
    arrange(desc(binding_counts))
  
  all_mer <- readRDS(paste0("~/GEO_data/SRA/gse_out/",all_samples[i],
                            "/jcast_out/binding_out/all_mers.rds"))
  ###peptide metadata
  dt2 <- all_mer %>% filter(seq %in% dt_all$peptide)
  dt2 <- dt2 %>% 
    tidyr::separate_wider_delim(seq_id, delim = "|",
                                names = c("sp","UniProt_ID","UniProt_symbol","Gene_ID",
                                          "rMATS_type","rMTS_ID","Chr","Anch_Exon_SE",
                                          "Alt_Exon_SE","strand_phase","Msjc","Tier")) %>% 
    select(-sp) %>% rename(pep_len = type) %>% 
    mutate(Msjc = as.numeric(gsub("r","",Msjc)))
  
  dt3 <- rmats_counts %>% filter(sample %in% all_samples[i])
  dt_meta <- dt2 %>% distinct_all() %>% 
    mutate(rMATS_type2 = sub("[0-9]$", "", rMATS_type)) %>% 
    left_join(., dt3 %>% rename(rMTS_ID = ID, rMATS_type2 = type) %>% 
                mutate(rMTS_ID = as.character(rMTS_ID))) %>% 
    select(-GeneID) %>% 
    select(seq, pep_len, rMTS_ID, sample, everything())
  saveRDS(dt_meta, 
          paste0("~/Drug_splicing/scripts/Shiny/data/",
                 all_samples[i], "_binding_pep_meta.rds"))
  saveRDS(dt_all, 
          paste0("~/Drug_splicing/scripts/Shiny/data/",
                 all_samples[i], "_binding_pep.rds"))
  message("Complete ",i)
}


