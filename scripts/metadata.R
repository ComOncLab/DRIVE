library(dplyr)
all_meta <- list.files("~/GEO_data/metadata/",full.names = T)
meta_info <- lapply(all_meta,
                    function(x){
                      tt <- data.table::fread(x,data.table = F)
                      tt <- tt %>% select(run_accession,study_accession,organism_name,library_layout,
                                          study_geo_accession,experiment_geo_accession, experiment_title)
                      return(tt)
                    })
meta_info <- bind_rows(meta_info)
meta_info <- meta_info %>% distinct_all()
meta_info <- meta_info %>% filter(!is.na(experiment_geo_accession)) %>%
  filter(nchar(experiment_geo_accession) > 0)
meta_info$organism_name <- gsub(" ","_",meta_info$organism_name)
saveRDS(meta_info,"data/meta_info.rds")

meta_info <- readRDS("~/Drug_splicing/data/meta_info.rds")

dt <- meta_info %>% 
  group_by(experiment_geo_accession) %>% 
  summarise(srrs = paste(run_accession,collapse = ","),
            organ = unique(organism_name), pair = unique(library_layout)) %>% ungroup()
write.table(dt, "~/GEO_data/SRA/meta_info.txt",quote = F,row.names = F,col.names = T,sep = " ")

# currents <- list.files("~/GEO_data/SRA/out/")
# mouse <- dt %>% filter(organ == "Mus_musculus")
# mouse_curr <- currents[which(currents %in% mouse$experiment_geo_accession)]
# write.table(data.frame(tt=mouse_curr),"~/GEO_data/SRA/need_rm",row.names = F,col.names = F,quote = F)
# 
# currents <- list.files("~/GEO_data/SRA/out/",recursive = T, pattern = "gene_counts.rds")
# currents <- gsub("/salmon_res.+","",currents)
# all_curr <- list.files("~/GEO_data/SRA/out/")
# all_curr <- all_curr[which(!(all_curr %in% currents))]
# write.table(data.frame(tt=all_curr),"~/GEO_data/SRA/need_rm",row.names = F,col.names = F,quote = F)
# 
# dt_rem <- dt %>% filter(!(experiment_geo_accession %in% currents))
# write.table(dt_rem, "~/GEO_data/SRA/split_run/meta_info.txt",
#             quote = F,row.names = F,col.names = F,sep = " ")

###########
all_samples <- list.files("~/GEO_data/GSM_deepseek/",full.names = T)
res <- vector("list",length(all_samples))
for (i in seq_along(res)){
  tt <- data.table::fread(all_samples[i],data.table = F,sep = ":",header = F)
  rownames(tt) <- tt$V1
  tt$V1 <- NULL
  tt <- t(tt) %>% as.data.frame()
  tt$sample <- gsub(".+//","",all_samples[i]) %>% gsub("_out.+","",.)
  tt <- tt %>% select(sample,everything())
  res[[i]] <- tt
}

res <- bind_rows(res)
rownames(res) <- NULL
colnames(res) <- c("sample_name","cell_line","treat","treat_con","treat_time")
write.table(res, "~/GEO_data/deepseek_r1.txt",quote = F,row.names = F,sep = "\t")

###r2
all_samples <- list.files("~/GEO_data/GSM_deepseek_correct/",full.names = T)
res <- vector("list",length(all_samples))
for (i in seq_along(res)){
  tt <- readr::read_delim(all_samples[i],delim = ": ",
                          col_names = c("V1","V2"),show_col_types = FALSE) %>% as.data.frame()
  rownames(tt) <- tt$V1
  tt$V1 <- NULL
  tt <- t(tt) %>% as.data.frame()
  tt$sample <- gsub(".+//","",all_samples[i]) %>% gsub("_out.+","",.)
  tt <- tt %>% select(sample,everything())
  tt$`Treatment concentration` <- as.character(tt$`Treatment concentration`)
  tt$`Treatment time` <- as.character(tt$`Treatment time`)
  res[[i]] <- tt
}

res <- bind_rows(res)
rownames(res) <- NULL
colnames(res) <- c("sample_name","cell_line","treat","treat_con","treat_time")
res <- res[,1:5]
saveRDS(res,"data/drug_treat_info.rds")

drug_treat_info <- readRDS("~/Drug_splicing/data/drug_treat_info.rds")
gsm2gse <- lapply(list.files("~/GEO_data/GSM2GSE/",full.names = T),
                  function(x){
                    tt <- data.table::fread(x,data.table = F)
                    tt$gsm <- gsub(".+//","",x) %>% gsub(".txt","",.)
                    return(tt)
                  })
gsm2gse <- bind_rows(gsm2gse)
gsm2gse <- gsm2gse %>% distinct_all()
drug_treat_info <- left_join(drug_treat_info,
                             gsm2gse %>% rename(sample_name = gsm))
drug_treat_info <- drug_treat_info %>% arrange(desc(study_alias))

saveRDS(drug_treat_info,"data/drug_treat_info.rds")
xlsx::write.xlsx(drug_treat_info,"data/drug_treat_info.xlsx")

####人工检查
dt <- read.table("~/GEO_data/SRA/meta_info.txt",header = T)
test2000 <- readxl::read_xlsx("data/drug_treat_info_2000.xlsx")
test2000$treat_time <- stringr::str_to_sentence(test2000$treat_time)
test2000$treat_time <- gsub("days \\(aza day2/4\\)","d",test2000$treat_time)
test2000$treat_time <- gsub("days","d",test2000$treat_time)
test2000$treat_time <- gsub("day","d",test2000$treat_time)
test2000$treat_time <- gsub("hour after treatment","h",test2000$treat_time)
test2000$treat_time <- gsub("hours after treatment","h",test2000$treat_time)
test2000$treat_time <- gsub("hours","h",test2000$treat_time)
test2000$treat_time <- gsub("hour","h",test2000$treat_time)
test2000$treat_time <- gsub("hrs","h",test2000$treat_time)
test2000$treat_time <- gsub("hr","h",test2000$treat_time)
test2000$treat_time <- gsub(" ","",test2000$treat_time)
test2000$cell_line <- ifelse(test2000$cell_line == "KYSE-71","KYSE-70",test2000$cell_line)

test2000 <- left_join(test2000,dt %>% rename(sample_name = 1))
test2000 <- test2000 %>% filter(organ == "Homo_sapiens")
sum_test <- test2000 %>% 
  mutate(Control = ifelse(is.na(Control),"N","Y")) %>% 
  group_by(study_alias,cell_line,treat,treat_time) %>% 
  summarise(samples = paste0(sample_name,collapse = ","),
            is_control = unique(Control),
            is_pair = unique(pair)) %>% ungroup()
sum_test <- sum_test %>% filter(!grepl("mevalonate",treat))
sum_test <- sum_test %>% 
  rowwise() %>% 
  mutate(sample_counts = strsplit(samples,",")[[1]] %>% length()) %>% 
  ungroup()
sum_test <- sum_test %>% filter(sample_counts > 1)
gses <- unique(sum_test$study_alias)
res <- vector("list")
for (i in 1:length(gses)){
  tmp_gse <- sum_test %>% filter(study_alias == gses[i]) 
  cell_lines <- unique(tmp_gse$cell_line)
  for (j in 1:length(cell_lines)){
    tmp_gse_cell <- tmp_gse %>% filter(cell_line == cell_lines[j]) 
    control_cell <- tmp_gse_cell %>% filter(is_control == "Y")
    treat_cell <- tmp_gse_cell %>% filter(is_control == "N")
    control_time <- unique(control_cell$treat_time)
    if (length(control_time) > 1){
      for (k in 1:length(control_time)){
        control_tmp <- control_cell %>% filter(treat_time == control_time[k])
        treat_tmp <- treat_cell %>% 
          filter(treat_time == control_time[k])%>% 
          rowwise() %>% 
          mutate(ids = paste(study_alias, cell_line, treat, treat_time,sep = "_")) %>% 
          ungroup() %>% 
          mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
          select(ids, samples) %>% 
          rename(treat_sample = samples) %>% 
          mutate(control_sample = control_tmp$samples)
        res <- c(res,list(treat_tmp))
      }
    }else{
      if (length(unique(treat_cell$treat_time)) > 1){
        treat_tmp <- treat_cell %>% 
          rowwise() %>% 
          mutate(ids = paste(study_alias, cell_line, treat, treat_time, sep = "_")) %>% 
          ungroup() %>% 
          mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
          select(ids, samples) %>% 
          rename(treat_sample = samples) %>% 
          mutate(control_sample = control_cell$samples)
      }else{
        treat_tmp <- treat_cell %>% 
          rowwise() %>% 
          mutate(ids = paste(study_alias, cell_line, treat,sep = "_")) %>% 
          ungroup() %>% 
          mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
          select(ids, samples) %>% 
          rename(treat_sample = samples) %>% 
          mutate(control_sample = control_cell$samples)
      }
      res <- c(res,list(treat_tmp))
    }
  }
}
res <- bind_rows(res)
write.table(res, "~/GEO_data/SRA/split_run/run_gse/gse_2000.txt",quote = F,
            row.names = F, col.names = F, sep = " ")
all_controls <- res$control_sample %>% unique()
write.table(data.frame(x=all_controls),
            "~/GEO_data/SRA/split_run/run_gse/control_2000.txt",quote = F,
            row.names = F, col.names = F, sep = " ")
write.table(dt,
            "~/GEO_data/SRA/split_run/run_gse/meta_info.txt",quote = F,
            row.names = F, col.names = F, sep = " ")

############
drug_treat_info <- readRDS("~/Drug_splicing/data/drug_treat_info.rds")
meta_info <- readRDS("~/Drug_splicing/data/meta_info.rds")
drug_treat_info <- inner_join(drug_treat_info, 
                              meta_info %>% select(experiment_geo_accession,organism_name) %>% 
                                rename(sample_name = experiment_geo_accession) %>% distinct_all())
drug_treat_info <- drug_treat_info %>% filter(organism_name == "Homo sapiens")
test2000 <- readxl::read_xlsx("data/drug_treat_info_2000.xlsx")
drug_treat_info <- drug_treat_info %>% filter(!(sample_name %in% test2000$sample_name))
xlsx::write.xlsx(drug_treat_info, "data/drug_treat_info_3000.xlsx")

######继续检查
dt <- read.table("~/GEO_data/SRA/meta_info.txt",header = T)
test3000 <- readxl::read_xlsx("data/drug_treat_info_3000_1.xlsx")
test3000 <- test3000[1:996,]
test3000$treat <- stringr::str_to_sentence(test3000$treat)
test3000$treat_time <- stringr::str_to_sentence(test3000$treat_time)
test3000$treat_time <- gsub("days","d",test3000$treat_time)
test3000$treat_time <- gsub("day","d",test3000$treat_time)
test3000$treat_time <- gsub("hours","h",test3000$treat_time)
test3000$treat_time <- gsub("hour","h",test3000$treat_time)
test3000$treat_time <- gsub("hrs","h",test3000$treat_time)
test3000$treat_time <- gsub("hr","h",test3000$treat_time)
test3000$treat_time <- gsub(" ","",test3000$treat_time)

test3000 <- left_join(test3000,dt %>% rename(sample_name = 1))
test3000 <- test3000 %>% filter(organ == "Homo_sapiens")
sum_test <- test3000 %>% 
  mutate(Control = ifelse(is.na(control),"N","Y")) %>% 
  group_by(study_alias,cell_line,treat,treat_time) %>% 
  summarise(samples = paste0(sample_name,collapse = ","),
            is_control = unique(Control),
            is_pair = unique(pair)) %>% ungroup()
sum_test <- sum_test %>% 
  rowwise() %>% 
  mutate(sample_counts = strsplit(samples,",")[[1]] %>% length()) %>% 
  ungroup()

sum_test <- sum_test %>% filter(sample_counts > 1)
gses <- unique(sum_test$study_alias)
res <- vector("list")
for (i in 1:length(gses)){
  tmp_gse <- sum_test %>% filter(study_alias == gses[i]) 
  cell_lines <- unique(tmp_gse$cell_line)
  for (j in 1:length(cell_lines)){
    tmp_gse_cell <- tmp_gse %>% filter(cell_line == cell_lines[j]) 
    control_cell <- tmp_gse_cell %>% filter(is_control == "Y")
    treat_cell <- tmp_gse_cell %>% filter(is_control == "N")
    control_time <- unique(control_cell$treat_time)
    if (length(control_time) > 1){
      for (k in 1:length(control_time)){
        control_tmp <- control_cell %>% filter(treat_time == control_time[k])
        treat_tmp <- treat_cell %>% 
          filter(treat_time == control_time[k])%>% 
          rowwise() %>% 
          mutate(ids = paste(study_alias, cell_line, treat, treat_time,sep = "_")) %>% 
          ungroup() %>% 
          mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
          select(ids, samples) %>% 
          rename(treat_sample = samples) %>% 
          mutate(control_sample = control_tmp$samples)
        res <- c(res,list(treat_tmp))
      }
    }else{
      if (length(unique(treat_cell$treat_time)) > 1){
        treat_tmp <- treat_cell %>% 
          rowwise() %>% 
          mutate(ids = paste(study_alias, cell_line, treat, treat_time, sep = "_")) %>% 
          ungroup() %>% 
          mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
          select(ids, samples) %>% 
          rename(treat_sample = samples) %>% 
          mutate(control_sample = control_cell$samples)
      }else{
        treat_tmp <- treat_cell %>% 
          rowwise() %>% 
          mutate(ids = paste(study_alias, cell_line, treat,sep = "_")) %>% 
          ungroup() %>% 
          mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
          select(ids, samples) %>% 
          rename(treat_sample = samples) %>% 
          mutate(control_sample = control_cell$samples)
      }
      res <- c(res,list(treat_tmp))
    }
  }
}
res <- bind_rows(res)

write.table(res, "~/GEO_data/SRA/split_run/run_gse/gse_3000_1.txt",quote = F,
            row.names = F, col.names = F, sep = " ")
all_controls <- res$control_sample %>% unique()
write.table(data.frame(x=all_controls),
            "~/GEO_data/SRA/split_run/run_gse/control_3000_1.txt",quote = F,
            row.names = F, col.names = F, sep = " ")

####
done <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "A3SS.MATS.JCEC.txt")
done <- gsub("/rmats_out.+","",done)
all_gse <- read.table("~/GEO_data/SRA/split_run/run_gse/all_gse.txt")
all_gse <- all_gse %>% 
  filter(!(V1 %in% done))
all_con <- read.table("~/GEO_data/SRA/split_run/run_gse/all_con.txt")
all_con <- all_con %>% filter(V1 %in% all_gse$V3)

write.table(all_con, "~/GEO_data/SRA/split_run/run_gse/rem_con",quote = F,
            row.names = F, col.names = F, sep = " ")
write.table(all_gse, "~/GEO_data/SRA/split_run/run_gse/rem_gse",quote = F,
            row.names = F, col.names = F, sep = " ")

#####last
dt <- read.table("~/GEO_data/SRA/meta_info.txt",header = T)
# write.table(dt,
#             "~/GEO_data/SRA/split_run/run_gse/meta_info.txt",quote = F,
#             row.names = F, col.names = F, sep = " ")

test3000 <- readxl::read_xlsx("data/drug_treat_info_3000_2.xlsx")
test2000 <- readxl::read_xlsx("data/drug_treat_info_2000.xlsx")
test2000$treat_time <- stringr::str_to_sentence(test2000$treat_time)
test2000$treat_time <- gsub("days \\(aza day2/4\\)","d",test2000$treat_time)
test2000$treat_time <- gsub("days","d",test2000$treat_time)
test2000$treat_time <- gsub("day","d",test2000$treat_time)
test2000$treat_time <- gsub("hour after treatment","h",test2000$treat_time)
test2000$treat_time <- gsub("hours after treatment","h",test2000$treat_time)
test2000$treat_time <- gsub("hours","h",test2000$treat_time)
test2000$treat_time <- gsub("hour","h",test2000$treat_time)
test2000$treat_time <- gsub("hrs","h",test2000$treat_time)
test2000$treat_time <- gsub("hr","h",test2000$treat_time)
test2000$treat_time <- gsub(" ","",test2000$treat_time)
test2000$cell_line <- ifelse(test2000$cell_line == "KYSE-71","KYSE-70",test2000$cell_line)
test2000$treat_con <- gsub(" ","",test2000$treat_con)
test2000$treat_con <- gsub("µ","u",test2000$treat_con)
test2000$treat_con <- gsub("μ","u",test2000$treat_con)
test2000$treat_con <- gsub("/","_",test2000$treat_con)
test2000$treat_con <- gsub("notspecified","Notspecified",test2000$treat_con)
test2000$treat_con <- gsub(",","_",test2000$treat_con)

test3000$treat <- stringr::str_to_sentence(test3000$treat)
test3000$treat_time <- stringr::str_to_sentence(test3000$treat_time)
test3000$treat_time <- gsub("days","d",test3000$treat_time)
test3000$treat_time <- gsub("day","d",test3000$treat_time)
test3000$treat_time <- gsub("hours","h",test3000$treat_time)
test3000$treat_time <- gsub("hour","h",test3000$treat_time)
test3000$treat_time <- gsub("hrs","h",test3000$treat_time)
test3000$treat_time <- gsub("hr","h",test3000$treat_time)
test3000$treat_time <- gsub(" ","",test3000$treat_time)
test3000$treat_con <- gsub(" ","",test3000$treat_con)
test3000$treat_con <- gsub("µ","u",test3000$treat_con)
test3000$treat_con <- gsub("μ","u",test3000$treat_con)
test2000$treat_con <- gsub("/","_",test2000$treat_con)
test3000$treat_con <- gsub("notspecified","Notspecified",test3000$treat_con)

test2000 <- test2000 %>% select(-Correct)
colnames(test2000) <- colnames(test3000)
all_dt <- bind_rows(test2000, test3000)
all_dt <- left_join(all_dt,dt %>% rename(sample_name = 1))
all_dt <- all_dt %>% filter(organ == "Homo_sapiens")

###如果control哟pair single，只要single
all_dt <- all_dt %>% 
  group_by(study_alias,cell_line,treat,treat_time) %>% 
  mutate(is_rm = ifelse((length(unique(pair)) > 1) & (pair == "SINGLE"),"yes","no")) %>% 
  ungroup() %>% 
  filter(is_rm != "yes") %>% select(-is_rm)
all_dt$treat_con[909] <- "0.5nMtrabectedinand5uMolaparib"

sum_test <- all_dt %>% 
  mutate(Control = ifelse(is.na(control),"N","Y")) %>% 
  group_by(study_alias,cell_line,treat,treat_time,treat_con) %>% 
  summarise(samples = paste0(sample_name,collapse = ","),
            is_control = unique(Control),
            is_pair = unique(pair)) %>% ungroup()
sum_test <- sum_test %>% filter(!grepl("mevalonate",treat))
sum_test <- sum_test %>% 
  rowwise() %>% 
  mutate(sample_counts = strsplit(samples,",")[[1]] %>% length()) %>% 
  ungroup()

sum_test <- sum_test %>% filter(sample_counts > 1)
gses <- unique(sum_test$study_alias)
res <- vector("list")
for (i in 1:length(gses)){
  tmp_gse <- sum_test %>% filter(study_alias == gses[i]) 
  cell_lines <- unique(tmp_gse$cell_line)
  for (j in 1:length(cell_lines)){
    tmp_gse_cell <- tmp_gse %>% filter(cell_line == cell_lines[j]) 
    control_cell <- tmp_gse_cell %>% filter(is_control == "Y")
    treat_cell <- tmp_gse_cell %>% filter(is_control == "N")
    control_time <- unique(control_cell$treat_time)
    if (length(control_time) > 1){
      for (k in 1:length(control_time)){
        control_tmp <- control_cell %>% filter(treat_time == control_time[k])
        treat_tmp <- treat_cell %>% 
          filter(treat_time == control_time[k])
        if (length(unique(treat_tmp$treat_con)) > 1){
          treat_tmp <- treat_tmp %>% 
            rowwise() %>% 
            mutate(ids = paste(study_alias, cell_line, treat, treat_time, treat_con, sep = "_")) %>% 
            ungroup() %>% 
            mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
            select(ids, samples) %>% 
            rename(treat_sample = samples) %>% 
            mutate(control_sample = control_tmp$samples)
        }else{
          treat_tmp <- treat_tmp %>% 
            rowwise() %>% 
            mutate(ids = paste(study_alias, cell_line, treat, treat_time,sep = "_")) %>% 
            ungroup() %>% 
            mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
            select(ids, samples) %>% 
            rename(treat_sample = samples) %>% 
            mutate(control_sample = control_tmp$samples)
        }
        res <- c(res,list(treat_tmp))
      }
    }else{
      if (length(unique(treat_cell$treat_time)) > 1){
        if (length(unique(treat_cell$treat_con)) > 1){
          treat_tmp <- treat_cell %>% 
            rowwise() %>% 
            mutate(ids = paste(study_alias, cell_line, treat, treat_time, treat_con, sep = "_")) %>% 
            ungroup() %>% 
            mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
            select(ids, samples) %>% 
            rename(treat_sample = samples) %>% 
            mutate(control_sample = control_cell$samples)
        }else{
          treat_tmp <- treat_cell %>% 
            rowwise() %>% 
            mutate(ids = paste(study_alias, cell_line, treat, treat_time,sep = "_")) %>% 
            ungroup() %>% 
            mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
            select(ids, samples) %>% 
            rename(treat_sample = samples) %>% 
            mutate(control_sample = control_cell$samples)
        }
      }else{
        if (length(unique(treat_cell$treat_con)) > 1){
          treat_tmp <- treat_cell %>% 
            rowwise() %>% 
            mutate(ids = paste(study_alias, cell_line, treat, treat_con, sep = "_")) %>% 
            ungroup() %>% 
            mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
            select(ids, samples) %>% 
            rename(treat_sample = samples) %>% 
            mutate(control_sample = control_cell$samples)
        }else{
          treat_tmp <- treat_cell %>% 
            rowwise() %>% 
            mutate(ids = paste(study_alias, cell_line, treat,sep = "_")) %>% 
            ungroup() %>% 
            mutate(ids = gsub(" ","_",ids) %>% gsub("/","-",.)) %>% 
            select(ids, samples) %>% 
            rename(treat_sample = samples) %>% 
            mutate(control_sample = control_cell$samples)
        }
      }
      res <- c(res,list(treat_tmp))
    }
  }
}
res <- bind_rows(res)

###已有的
done <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "A3SS.MATS.JCEC.txt")
done <- gsub("/rmats_out.+","",done)
done[which(!(done %in% res$ids))] -> tt1
dt <- data.frame(unid = tt1)
dt <- dt %>% 
  rowwise() %>% 
  mutate(match_counts = grep(unid,res$ids,fixed = T) %>% length()) %>% ungroup() %>% 
  filter(match_counts == 1) %>% 
  rowwise() %>% 
  mutate(ids = res$ids[grep(unid,res$ids,fixed = T)]) %>% ungroup()

res <- left_join(res, dt)
res <- res %>% 
  mutate(unid = ifelse(ids %in% done, ids, unid))
res <- res %>% filter(!is.na(unid))
saveRDS(res, "data/all_samples.rds")

tt <- left_join(all_samples, 
                sum_test %>% rename(treat_sample = samples) %>% 
                  select(-is_control,-is_pair, -sample_counts))
saveRDS(tt,"data/all_sample_meta.rds")
xlsx::write.xlsx(tt %>% select(-match_counts,-ids),"data/all_sample_meta.xlsx")

####
meta_dt <- read.table("~/GEO_data/SRA/meta_info.txt",header = T)
meta_dt <- meta_dt %>% filter(organ == "Homo_sapiens")
curr_dt <- readRDS("data/all_samples.rds")
curr_samples <- strsplit(paste(c(paste(curr_dt$treat_sample,collapse = ","),
                                 paste(curr_dt$control_sample,collapse = ",")),collapse = ","),
                         ",")[[1]] %>% unique()

gsm2gse <- lapply(list.files("~/GEO_data/GSM2GSE/",full.names = T),
                  function(x){
                    tt <- data.table::fread(x,data.table = F)
                    tt$gsm <- gsub(".+//","",x) %>% gsub(".txt","",.)
                    return(tt)
                  })
gsm2gse <- bind_rows(gsm2gse)
gsm2gse <- gsm2gse %>% distinct_all()
meta_dt <- left_join(meta_dt %>% rename(gsm = experiment_geo_accession),gsm2gse)
meta_dt <- meta_dt %>% 
  rowwise() %>% 
  mutate(is_done = ifelse(gsm %in% curr_samples,"yes","no")) %>% ungroup()

tt <- meta_dt %>% 
  group_by(study_alias) %>% summarise(c_n =  length(unique(is_done))) %>% ungroup() %>% 
  filter(c_n > 1)
tt1 <- meta_dt %>% filter(study_alias %in% tt$study_alias) %>% 
  filter(is_done == "no")
####没有的样本是不需要的，基本是敲除或者耐药的，不是WT的细胞系背景
meta_dt <- meta_dt %>% filter(is_done == "yes")
all_srr <- strsplit(paste(meta_dt$srrs,collapse = ","),split = ",")[[1]]

test3000 <- readxl::read_xlsx("data/drug_treat_info_3000_2.xlsx")
test2000 <- readxl::read_xlsx("data/drug_treat_info_2000.xlsx")
dt <- bind_rows(test2000 %>% select(1,2), test3000 %>% select(1,2))

#############
all_sample_meta <- readxl::read_xlsx("data/all_sample_meta.xlsx")
drug_info <- all_sample_meta %>% tidyr::separate_longer_delim(cols = "treat",delim = "+") %>% 
  select(treat, DrugType, SMILE) %>% distinct_all()
tt <- drug_info %>% filter(is.na(DrugType))
tt1 <- drug_info %>% filter(!is.na(DrugType))
tt <- left_join(tt %>% select(treat), tt1)
drug_info <- bind_rows(tt, tt1) %>% distinct_all()
xlsx::write.xlsx(drug_info,"data/drug_meta.xlsx")

############
drug_info <- readxl::read_xlsx("data/drug_meta.xlsx")
drug_class <- drug_info %>% select(DrugType) %>% distinct_all() %>% 
  filter(!is.na(DrugType))
xlsx::write.xlsx(drug_class,"data/drug_class.xlsx")
drug_class <- readxl::read_xlsx("data/drug_class.xlsx")

###细胞类型
all_sample_meta <- readxl::read_xlsx("data/all_sample_meta.xlsx")
####每个细胞系的癌症类型，在细胞系cellosaurus中搜索
cells <- readr::read_file("~/Genome_data/cellosaurus.txt")
cells <- cells %>% stringr::str_split("//\n") %>% .[[1]]
cells <- cells[-167128]
library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  40, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = cells,
  .packages = c("dplyr","readr")
) %dopar% {
  cell_dt <- readr::read_delim(i,delim = "   ",col_names = c("attr","value"))
  if ("DI" %in% cell_dt$attr){
    return(
      data.frame(
        ID = cell_dt$value[which(cell_dt$attr == "ID")],
        AC = cell_dt$value[which(cell_dt$attr == "AC")],
        SY = ifelse("SY" %in% cell_dt$attr,cell_dt$value[which(cell_dt$attr == "SY")],NA),
        DI_NCIt = ifelse(length(grep("NCIt",cell_dt$value)) != 0, 
                         cell_dt$value[grep("NCIt",cell_dt$value)],NA),
        DI_ORDO = ifelse(length(grep("ORDO",cell_dt$value)) != 0, 
                         cell_dt$value[grep("ORDO",cell_dt$value)],NA)
      )
    )
  }else{
    return(NA)
  }
}
parallel::stopCluster(cl = my.cluster)
res <- res[lengths(res) > 1]
res <- bind_rows(res)
res <- res %>% 
  rowwise() %>% 
  mutate(ID2 = paste(ID,SY,sep = "; ")) %>% ungroup()
split_cells <- sapply(res$ID2,
                      function(x){
                        strsplit(x,"; ")[[1]]
                      })
get_match_row <- function(cell_name){
  tmp <- lapply(split_cells,function(x){cell_name %in% x}) %>% unlist %>% unname
  if (length(tmp) == 0){
    return(NA)
  }else{
    return(which(tmp)[1])
  }
}

cell_lines <- all_sample_meta %>% select(cell_line) %>% distinct_all()
cell_lines <- cell_lines %>% rowwise() %>% 
  mutate(match_row = get_match_row(cell_line)) %>% 
  mutate(ID = ifelse(!is.na(match_row),res$ID[match_row],NA),
         ID2 = ifelse(!is.na(match_row),res$SY[match_row],NA),
         DI = ifelse(!is.na(match_row),res$DI_NCIt[match_row],NA)) %>% 
  ungroup()
cell_lines$DI <- gsub(".+; ","",cell_lines$DI)
xlsx::write.xlsx(cell_lines,"data/cell_line_meta.xlsx")

##########



