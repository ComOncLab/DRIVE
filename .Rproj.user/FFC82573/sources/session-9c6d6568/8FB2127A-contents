library(dplyr)
all_dir <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "A3SS.MATS.JCEC.txt",
                      full.names = T)
all_dir <- gsub("/rmats_out/.+","",all_dir)
as_type <- c("A3SS","A5SS","MXE","RI","SE")
for (i in all_dir){
  #unlink(paste0(i,"/rmats_filter_out"), recursive = TRUE)
  #dir.create(paste0(i,"/rmats_filter_out"))
  for (j in 1:length(as_type)){
    file_name <- paste0(i,"/rmats_out/",as_type[j],".MATS.JC.txt")
    dt <- data.table::fread(file_name,data.table = F,check.names = T,
                            colClasses=list(character=c("IJC_SAMPLE_1","IJC_SAMPLE_2",
                                                        "SJC_SAMPLE_1","SJC_SAMPLE_2",
                                                        "IncLevel1","IncLevel2"),
                                            numeric=c("PValue","FDR","IncLevelDifference")))
    dt$PValue <- as.numeric(dt$PValue)
    dt$FDR <- as.numeric(dt$FDR)
    dt$IncLevelDifference <- as.numeric(dt$IncLevelDifference)
    dt <- dt %>% filter(FDR < 0.1)
    save_name <- paste0(i,"/rmats_filter_out/",as_type[j],".MATS.JC.txt")
    write.table(dt,save_name,sep = "\t",quote = F,row.names = F)
  }
}
#####之前有问题的需要重新跑jcast
all_dir <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "A3SS.MATS.JCEC.txt",
                      full.names = T)
all_dir <- gsub("/rmats_out/.+","",all_dir)
pros <- c()
for (i in all_dir){
  file_name <- paste0(i,"/rmats_out/A3SS.MATS.JC.txt")
  dt <- data.table::fread(file_name,data.table = F,check.names = T)
  if (is.character(dt$FDR)){
    pros <- append(pros, gsub(".+//","",i))
  }
}
write.table(data.frame(x=pros),"~/GEO_data/SRA/split_run/run_jcast/rems.txt",
            quote = F,row.names = F,col.names = F,sep = " ")

########有些样本筛 FDR<0.1之后就没有差异可变剪切了
all_files <- list.files("~/GEO_data/SRA/gse_out/")
rmats_filter <- vector("list",length = length(all_files))
for (i in seq_along(rmats_filter)){
  in_files <- list.files(paste0("~/GEO_data/SRA/gse_out/",all_files[i],
                                "/rmats_filter_out/"),full.names = T)
  dt <- lapply(in_files,
               function(x){
                 tt <- data.table::fread(x, data.table = F)
                 data.frame(type = x %>% gsub(".+//","",.) %>% gsub(".MATS.JC.txt","",.),
                            counts = nrow(tt))
               }) %>% bind_rows()
  dt$ids <- all_files[i]
  rmats_filter[[i]] <- dt
}
rmats_filter <- bind_rows(rmats_filter)
saveRDS(rmats_filter,"~/Drug_splicing/data/rmats_fdr01.rds")

####
all_samples <- read.table("~/GEO_data/SRA/split_run/run_jcast/all_samples.txt")
done <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "psq_")
done <- data.frame(name = done) %>% 
  mutate(ids = gsub("/jcast_out.+","",name)) %>% 
  mutate(psq = gsub(".+/","",name))
rems <- all_samples %>% filter(!(V1 %in% done$ids))
rms_filter <- readRDS("~/Drug_splicing/data/rmats_fdr01.rds")
rms_filter <- rms_filter %>% 
  group_by(ids) %>% summarise(t_counts = sum(counts)) %>% ungroup()
rems <- left_join(rems %>% rename(ids = 1), rms_filter)
###最后一行
all_files <- paste0("~/GEO_data/SRA/gse_out/",rems$ids,"/jcast_out/")
last_line <- c()
for (i in 1:length(all_files)){
  tt <- list.files(all_files[i],pattern = "log",recursive = T,full.names = T)
  tt1 <- system(paste0("grep -v '^\\s*$' ",tt," | tail -n 1"), intern = TRUE)
  last_line <- append(last_line, tt1)
}
####没有跑完的要么是FDR筛选后没有的，要么是Sequence discarded due to low coverage.

####拼起来
done <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "psq_")
done <- data.frame(name = done) %>% 
  mutate(ids = gsub("/jcast_out.+","",name)) %>% 
  mutate(psq = gsub(".+/","",name))
done_sum <- done %>% 
  group_by(ids) %>% 
  summarise(
    c_c = sum(psq=="psq_canonical.fasta"),
    t1_c = sum(psq=="psq_T1.fasta"),
    t2_c = sum(psq=="psq_T2.fasta"),
    t3_c = sum(psq=="psq_T3.fasta"),
    t4_c = sum(psq=="psq_T4.fasta")) %>% ungroup()
tie_names <- c("psq_T1.fasta","psq_T2.fasta",
               "psq_T3.fasta","psq_T4.fasta")
for (i in 1:nrow(done_sum)){
  psq_file <- list.files(paste0("~/GEO_data/SRA/gse_out/",done_sum$ids[i],"/jcast_out/"),
                         pattern = "jcast",full.names = T)
  psq_file <- Hmisc::escapeRegex(psq_file)
  psq_file <- gsub(";","\\;",psq_file,fixed = T)
  if (!file.exists(gsub("//jcast_.+","/all_seq.fasta",psq_file))){
    system(paste0("cat ",
                  paste(paste0(psq_file,"/",tie_names[as.logical(done_sum[i,3:6])]),collapse = " "),
                  " > ",gsub("//jcast_.+","/all_seq.fasta",psq_file))) 
  }
}
saveRDS(done_sum,"data/jcast_done_samples.rds")

####读取可变剪切过滤后的结果
all_files <- readRDS("~/Drug_splicing/data/rmats_fdr01.rds")
all_files <- all_files %>% 
  group_by(ids) %>% summarise(t_counts = sum(counts)) %>% ungroup() %>% 
  filter(t_counts > 0)

rmats_counts <- vector("list",length = length(all_files$ids))
for (i in seq_along(rmats_counts)){
  in_files <- list.files(paste0("~/GEO_data/SRA/gse_out/",all_files$ids[i],
                                "/rmats_filter_out/"),full.names = T)
  dt <- lapply(in_files,
               function(x){
                 tt <- data.table::fread(x,data.table = F,check.names = T,
                                         colClasses=list(character=c("IJC_SAMPLE_1","IJC_SAMPLE_2",
                                                                     "SJC_SAMPLE_1","SJC_SAMPLE_2",
                                                                     "IncLevel1","IncLevel2"),
                                                         numeric=c("PValue","FDR","IncLevelDifference")))
                 tt$PValue <- as.numeric(tt$PValue)
                 tt$FDR <- as.numeric(tt$FDR)
                 tt$IncLevelDifference <- as.numeric(tt$IncLevelDifference)
                 if (nrow(tt) == 0){
                   return(NA)
                 }
                 tt$type <- x %>% gsub(".+//","",.) %>% gsub(".MATS.JC.txt","",.)
                 return(tt)
               })
  dt <- dt[lengths(dt) > 1]
  dt_counts <- lapply(dt,
                      function(x){
                        tt <- x %>% 
                          select(ID,GeneID,geneSymbol,IJC_SAMPLE_1,IJC_SAMPLE_2,
                                 SJC_SAMPLE_1,SJC_SAMPLE_2,PValue,FDR,type)
                        }) %>% bind_rows() %>% 
    mutate(sample = all_files$ids[i])
  rmats_counts[[i]] <- dt_counts
}
rmats_counts <- bind_rows(rmats_counts)
saveRDS(rmats_counts,"data/rmats_filter.rds")

####筛选需要的肽，sample1 是treat，sample2 是control
##如果 IJC_SAMPLE_1 > IJC_SAMPLE_2，选择1
###如果 SJC_SAMPLE_1 > SJC_SAMPLE_2，选择2
all_fasta <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "all_seq.fasta",
                        full.names = T)

library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  40, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

seq_df <- foreach(
  i = 1:length(all_fasta),
  .export = c("all_fasta"),
  .packages = c("dplyr","seqinr")
) %dopar% {
  tt <- seqinr::read.fasta(all_fasta[i],seqtype = "AA",as.string = T)
  type_ids <- lapply(names(tt),
                     function(x){
                       tmp <- strsplit(x,"\\|")[[1]][c(5,6)]
                       data.frame(
                         type = tmp[1],
                         ID = tmp[2]
                       )
                     }) %>% bind_rows()
  type_ids$sample <- gsub(".+gse_out//","",all_fasta[i]) %>% gsub("/jcast.+","",.)
  return(type_ids)
}
parallel::stopCluster(cl = my.cluster)
seq_df <- bind_rows(seq_df)

####
rmats_counts <- readRDS("data/rmats_filter.rds")
rmats_counts <- rmats_counts %>% 
  mutate(com_id = paste0(sample,"_",type,"_",ID))
seq_df <- seq_df %>% 
  mutate(type2 = stringr::str_sub(type, 1, -2)) %>% 
  mutate(com_id = paste0(sample,"_",type2,"_",ID))
rmats_counts <- rmats_counts %>% 
  filter(com_id %in% seq_df$com_id)

choose_12 <- function(i1,i2,s1,s2){
  i1 <- strsplit(i1,",")[[1]] %>% as.numeric() %>% mean()
  i2 <- strsplit(i2,",")[[1]] %>% as.numeric() %>% mean()
  s1 <- strsplit(s1,",")[[1]] %>% as.numeric() %>% mean()
  s2 <- strsplit(s2,",")[[1]] %>% as.numeric() %>% mean()
  if (i1 > i2){
    if (s1 > s2){
      return("both")
    }else{
      return("one")
    }
  }else {
    if (s1 > s2){
      return("two")
    }else{
      return("no")
    }
  }
}

rmats_counts <- rmats_counts %>% 
  rowwise() %>% 
  mutate(sele = choose_12(IJC_SAMPLE_1,IJC_SAMPLE_2,SJC_SAMPLE_1,SJC_SAMPLE_2)) %>% 
  ungroup()
rmats_counts <- rmats_counts %>% 
  filter(sele != "no")
seq_df <- inner_join(seq_df, rmats_counts %>% select(com_id, sele) %>% distinct_all())  
seq_df <- seq_df %>% 
  group_by(com_id) %>% 
  mutate(is_keep = case_when(
    sele == "one" ~ (type == paste0(type2,"1")),
    sele == "two" ~ (type == paste0(type2,"2")),
    sele == "both" ~ TRUE,
  )) %>% ungroup()
seq_df <- seq_df %>% filter(is_keep)
saveRDS(seq_df,"data/keep_var2.rds")  

####保留需要的肽
seq_df <- readRDS("~/Drug_splicing/data/keep_var2.rds")
rems <- read.table("~/GEO_data/SRA/split_run/run_jcast/rems.txt")
seq_df <- seq_df %>% filter(sample %in% rems$V1)
seq_df <- seq_df %>% 
  mutate(com_id = paste0(sample,"_",type,"_",ID))
all_fasta <- paste0("/home/wt/GEO_data/SRA/gse_out/",unique(seq_df$sample),
                    "/jcast_out/all_seq.fasta")

library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  40, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

seq_df <- foreach(
  i = 1:length(all_fasta),
  .export = c("all_fasta","seq_df"),
  .packages = c("dplyr","seqinr")
) %dopar% {
  tt <- seqinr::read.fasta(all_fasta[i],seqtype = "AA",as.string = T)
  type_ids <- lapply(names(tt),
                     function(x){
                       tmp <- strsplit(x,"\\|")[[1]][c(5,6)]
                       data.frame(
                         type = tmp[1],
                         ID = tmp[2]
                       )
                     }) %>% bind_rows()
  type_ids <- type_ids %>% 
    mutate(com_id = paste0(gsub(".+gse_out/","",all_fasta[i]) %>% gsub("/j.+","",.),
                           "_",type,"_",ID)) %>% 
    mutate(ori_ids = 1:nrow(.))
  type_ids <- type_ids %>% 
    filter(com_id %in% seq_df$com_id)
  tt <- tt[type_ids$ori_ids]
  seqinr::write.fasta(tt, names = names(tt), file.out = all_fasta[i])
}
parallel::stopCluster(cl = my.cluster)

#####切割肽
generate_kmers_with_step <- function(sequence, k, step = 1) {
  seq_length <- nchar(sequence)
  if (seq_length < k) {
    return(NA)
  }
  # 生成起始位置序列
  start_positions <- seq(1, seq_length - k + 1, by = step)
  # 提取k-mer
  kmers <- sapply(start_positions, function(i) {
    substr(sequence, i, i + k - 1)
  })
  return(kmers)
}

standard_aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
all_fasta <- list.files("/home/wt/GEO_data/SRA/gse_out/",recursive = T,
                        pattern = "all_seq.fasta",full.names = T)
fa_ids <- gsub(".+gse_out//","",all_fasta) %>% gsub("/.+","",.)
all_fasta <- all_fasta[which(fa_ids %in% rems$V1)]

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
  i = 1:length(all_fasta),
  .export = c("all_fasta","standard_aa","generate_kmers_with_step"),
  .packages = c("dplyr","seqinr")
) %dopar% {
  ###保留标准氨基酸以及长度大于8的
  seqs <- seqinr::read.fasta(all_fasta[i],seqtype = "AA",as.string = T)
  res <- c()
  for (j in seq_along(seqs)){
    tmp_seq <- seqs[[j]] %>% as.character() %>% strsplit(.,"") %>% `[[`(1)
    non_stand <- which(!(tmp_seq %in% standard_aa))
    if (length(non_stand) != 0){
      res <- append(res, j)
    }
  }
  if (length(res) != 0){
    seqs <- seqs[-res]
  }
  seqs <- seqs[which(nchar(seqs) >= 8)]
  dir.create(gsub("all_seq.fasta","binding_out",all_fasta[i]))
  seqinr::write.fasta(seqs,names = names(seqs),
                      file.out = paste0(gsub("all_seq.fasta","binding_out",all_fasta[i]),
                                        "/filter.fasta"))
  ####切割肽
  res <- vector("list",length(seqs))
  for (j in seq_along(res)){
    tmp_seq <- seqs[[j]] %>% as.character()
    mer_8 <- generate_kmers_with_step(tmp_seq, k=8)
    mer_9 <- generate_kmers_with_step(tmp_seq, k=9)
    mer_10 <- generate_kmers_with_step(tmp_seq, k=10)
    mer_11 <- generate_kmers_with_step(tmp_seq, k=11)
    tmp_res <- data.frame(
      seq_id = names(seqs[j]),
      type = c(rep("8_mer",length(mer_8)),rep("9_mer",length(mer_9)),
               rep("10_mer",length(mer_10)),rep("11_mer",length(mer_11))),
      seq = c(mer_8, mer_9, mer_10, mer_11)
    )
    res[[j]] <- tmp_res
  }
  res <- bind_rows(res)
  res <- res %>% filter(!is.na(seq))
  saveRDS(res, paste0(gsub("all_seq.fasta","binding_out",all_fasta[i]),"/all_mers.rds"))
  ##用tmux进行并行
  dir.create(paste0(gsub("all_seq.fasta","binding_out",all_fasta[i]),"/split_fa"))
  n <- length(seqs)
  chunk_size <- 100
  group_indices <- rep(1:ceiling(n / chunk_size), each = chunk_size, length.out = n)
  sp_seqs <- split(seqs, group_indices)
  for (j in seq_along(sp_seqs)){
    seqinr::write.fasta(sp_seqs[[j]],names = names(sp_seqs[[j]]),
                        file.out = paste0(gsub("all_seq.fasta","binding_out",all_fasta[i]),
                                          "/split_fa/split_",j,".fasta"))
  }
  dir.create(paste0(gsub("all_seq.fasta","binding_out",all_fasta[i]),"/res"))
}
parallel::stopCluster(cl = my.cluster)

####每个细胞系的HLA分型，在细胞系cellosaurus中搜索
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
  hla_lines <- grep("HLA typing:",cell_dt$value)
  if (length(hla_lines) != 0){
    return(
      data.frame(
        ID = cell_dt$value[which(cell_dt$attr == "ID")],
        AC = cell_dt$value[which(cell_dt$attr == "AC")],
        SY = ifelse("SY" %in% cell_dt$attr,cell_dt$value[which(cell_dt$attr == "SY")],NA),
        HLA = gsub("HLA typing: ","",cell_dt$value[hla_lines]) %>% 
          gsub("\\(.+","",.) %>% paste(.,collapse = "-")
      )
    )
  }else{
    return(NA)
  }
}
parallel::stopCluster(cl = my.cluster)
res <- res[lengths(res) > 1]
res <- bind_rows(res)

all_sample_meta <- readxl::read_xlsx("data/all_sample_meta.xlsx")
all_fasta <- list.files("/home/wt/GEO_data/SRA/gse_out/",recursive = T,
                        pattern = "all_seq.fasta") %>% gsub("/jcast_out/.+","",.)
all_sample_meta <- all_sample_meta %>%
  filter(unid %in% all_fasta)
cell_lines <- all_sample_meta %>% select(cell_line) %>% distinct_all()
cell_lines <- cell_lines %>% rowwise() %>% 
  mutate(match_row = ifelse(length(grep(cell_line, res$ID))!=0,
                            grep(cell_line, res$ID),grep(cell_line, res$SY))) %>% 
  mutate(ID = ifelse(!is.na(match_row),res$ID[match_row],NA),
         ID2 = ifelse(!is.na(match_row),res$SY[match_row],NA),
         HLA = ifelse(!is.na(match_row),res$HLA[match_row],NA)) %>% 
  ungroup()
###人工检查选择HLA类型
xlsx::write.xlsx(cell_lines, "data/cell_line_HLA.xlsx")

cell_lines <- readxl::read_xlsx("data/cell_line_HLA.xlsx") %>% select(-1)
####没有的再用从RNA数据中call出的
cell_samples <- all_sample_meta %>%
  group_by(cell_line) %>%
  summarise(samples = strsplit(paste(c(paste(treat_sample,collapse = ","),
                                       paste(control_sample,collapse = ",")),
                                     collapse = ","),split = ",")[[1]] %>% unique() %>%
              paste(.,collapse = ",")) %>% ungroup()
cell_samples <- cell_samples %>%
  mutate(A1 = NA, A2 = NA,
         B1 = NA, B2 = NA,
         C1 = NA, C2 = NA)

# all_samples <- strsplit(paste(cell_samples$samples,collapse = ","),",")[[1]] %>% unique()
# for (i in 1:length(all_samples)){
#   file.copy(paste0("~/GEO_data/SRA/out/",all_samples[i],"/",all_samples[i],".json"),
#             paste0("~/GEO_data/SRA/HLA_typing/",all_samples[i],".genotype.json"),overwrite = T)
# }
all_hla <- data.table::fread("~/GEO_data/SRA/HLA_typing/genotypes.tsv",data.table = F)
all_hla <- all_hla[,1:7]
all_hla <- apply(all_hla,1,function(x){ifelse(nchar(x)==0,NA,x)}) %>% t() %>% as.data.frame()
###
all_hla[,2:7] <- apply(all_hla[,2:7],1,
                 function(x){
                   ifelse(grepl("[D-Z]",x),NA,x)
                   }) %>% t()
all_hla[,2:7] <- apply(all_hla[,2:7],1,
                       function(x){
                         stringr::str_replace(x, "^([^:]*:[^:]*):.*$", "\\1")
                       }) %>% t()

for (i in 1:nrow(cell_samples)){
  tmp <- all_hla %>% filter(subject %in% strsplit(cell_samples$samples[i],",")[[1]])
  if (nrow(tmp) != 0){
    cell_samples$A1[i] <- ifelse(all(is.na(tmp$A1)),NA,sort(table(tmp$A1),decreasing = T) %>% names() %>% `[`(1))
    cell_samples$A2[i] <- ifelse(all(is.na(tmp$A2)),NA,sort(table(tmp$A2),decreasing = T) %>% names() %>% `[`(1))
    cell_samples$B1[i] <- ifelse(all(is.na(tmp$B1)),NA,sort(table(tmp$B1),decreasing = T) %>% names() %>% `[`(1))
    cell_samples$B2[i] <- ifelse(all(is.na(tmp$B2)),NA,sort(table(tmp$B2),decreasing = T) %>% names() %>% `[`(1))
    cell_samples$C1[i] <- ifelse(all(is.na(tmp$C1)),NA,sort(table(tmp$C1),decreasing = T) %>% names() %>% `[`(1))
    cell_samples$C2[i] <- ifelse(all(is.na(tmp$C2)),NA,sort(table(tmp$C2),decreasing = T) %>% names() %>% `[`(1))
  }
}
###人工检查完成
cell_lines <- readxl::read_xlsx("data/cell_line_HLA.xlsx") %>% select(-1)
all_sample_meta <- readRDS("data/all_sample_meta.rds")
all_fasta <- list.files("/home/wt/GEO_data/SRA/gse_out/",recursive = T,
                        pattern = "all_seq.fasta") %>% gsub("/jcast_out/.+","",.)
all_sample_meta <- left_join(all_sample_meta, cell_lines %>% select(cell_line, HLA))
all_sample_meta <- all_sample_meta %>% filter(!is.na(HLA))
all_sample_meta <- all_sample_meta %>% 
  filter(unid %in% all_fasta)

split_hla <- function(hla_str){
  split1 <- strsplit(hla_str,"; ")[[1]]
  split1 <- strsplit(split1,",")
  cat_hla <- lapply(split1,
                    function(x){
                      if (length(x) > 1){
                        paste0("HLA-",c(x[1],paste0(substr(x[1],1,1),"*",x[2:length(x)])))
                      }else{
                        paste0("HLA-",x)
                      }
                    }) %>% unlist() %>% unique() %>% paste(.,collapse = ",")
  return(cat_hla)
}

all_sample_meta <- all_sample_meta %>% 
  rowwise() %>% 
  mutate(HLA2 = split_hla(HLA)) %>% ungroup() %>% 
  select(unid, cell_line, HLA2) %>% rename(HLA = HLA2)
saveRDS(all_sample_meta, "data/samples_with_HLA.rds")  
  
#####开始运行
samples <- readRDS("~/Drug_splicing/data/samples_with_HLA.rds")
#rems <- read.table("~/GEO_data/SRA/split_run/run_jcast/rems.txt")
#samples <- samples %>% filter(unid %in% rems$V1)
samples$unid2 <- Hmisc::escapeRegex(samples$unid)
samples$unid2 <- gsub(";","\\;",samples$unid2,fixed = T)
ba_json <- jsonlite::read_json("~/GEO_data/SRA/split_run/run_binding/ba.json")
el_json <- jsonlite::read_json("~/GEO_data/SRA/split_run/run_binding/el.json")
pre_ba <- readLines("~/GEO_data/SRA/split_run/run_binding/run_batch_ba.sh")
pre_el <- readLines("~/GEO_data/SRA/split_run/run_binding/run_batch_el.sh")
for (i in 1:length(samples$unid)){
  main_path <- paste0("/home/wt/GEO_data/SRA/gse_out/",samples$unid[i],"/jcast_out/binding_out/")
  main_path2 <- paste0("/home/wt/GEO_data/SRA/gse_out/",samples$unid2[i],"/jcast_out/binding_out/")
  split_fa <- list.files(paste0(main_path,"/split_fa"))
  write.table(data.frame(x=split_fa),paste0(main_path,"/all_splits"),
              row.names = F, col.names = F, quote = F, sep = " ")
  
  tmp_ba_json <- ba_json
  tmp_ba_json$alleles <- samples$HLA[i]
  jsonlite::write_json(tmp_ba_json, paste0(main_path,"/ba.json"),auto_unbox = TRUE, pretty = TRUE)
  tmp_el_json <- el_json
  tmp_el_json$alleles <- samples$HLA[i]
  jsonlite::write_json(tmp_el_json, paste0(main_path,"/el.json"),auto_unbox = TRUE, pretty = TRUE)
  
  pre_ba_tmp <- pre_ba
  pre_ba_tmp[6] <- paste0("\t\tcd ",main_path2)
  pre_ba_tmp <- c(paste0("cd ",main_path2),pre_ba_tmp)
  pre_ba_tmp[11] <- gsub("ba_",paste0(samples$unid[i],"_ba_"),pre_ba_tmp[11])
  pre_ba_tmp[5] <- paste0("\ttmux new-session -d -s ",samples$unid2[i],
                          "_ba_${split_name} \"")
  writeLines(pre_ba_tmp,paste0(main_path,"/run_batch_ba.sh"))
  pre_el_tmp <- pre_el
  pre_el_tmp[6] <- paste0("\t\tcd ",main_path2)
  pre_el_tmp <- c(paste0("cd ",main_path2),pre_el_tmp)
  pre_el_tmp[11] <- gsub("el_",paste0(samples$unid[i],"_el_"),pre_el_tmp[11])
  pre_el_tmp[5] <- paste0("\ttmux new-session -d -s ",samples$unid2[i],
                          "_el_${split_name} \"")
  writeLines(pre_el_tmp,paste0(main_path,"/run_batch_el.sh"))
}

all_ba <- paste0("bash /home/wt/GEO_data/SRA/gse_out/",samples$unid2,
                 "/jcast_out/binding_out/run_batch_ba.sh")
writeLines(all_ba,"~/GEO_data/SRA/split_run/run_binding/split_ba/run_ba_all.sh")
all_el <- paste0("bash /home/wt/GEO_data/SRA/gse_out/",samples$unid2,
                 "/jcast_out/binding_out/run_batch_el.sh")
writeLines(all_el,"~/GEO_data/SRA/split_run/run_binding/split_el/run_el_all.sh")

###########查看完成情况
files <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "all_split")
all_split <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,
                           pattern = "tsv")
all_ba_done <- all_split[grepl("ba_split",all_split)]

all_counts <- sapply(files,
                     function(x){
                       tt <- read.table(paste0("~/GEO_data/SRA/gse_out/",x),header = F)
                       data.frame(unid = gsub("/.+","",x),
                                  counts = nrow(tt))
                     },simplify = F)
all_counts <- bind_rows(all_counts)
saveRDS(all_counts,"~/Drug_splicing/data/split_counts.rds")

all_ba_done <- data.frame(ids = all_ba_done) %>% 
  tidyr::separate_wider_delim(cols = ids, delim = "/jcast_out/binding_out/res/", 
                              names = c("unid","sp"))
all_ba_done <- all_ba_done %>% group_by(unid) %>% summarise(done_counts = n()) 
rems_ba <- left_join(all_counts,all_ba_done) %>% 
  filter((done_counts != counts) | (is.na(done_counts)))
rems_ba$unid <- Hmisc::escapeRegex(rems_ba$unid)
rems_ba$unid <- gsub(";","\\;",rems_ba$unid,fixed = T)

all_ba_commd <- read.table("~/GEO_data/SRA/split_run/run_binding/split_ba/run_ba_all.sh")
all_ba_commd <- all_ba_commd %>% 
  mutate(unid = gsub(".+gse_out/","",V2) %>% gsub("/jcast_out/.+","",.))
rems_ba_commd <- all_ba_commd %>% filter(unid %in% rems_ba$unid)
write.table(rems_ba_commd %>% select(1,2) %>% rowwise() %>% 
              mutate(comm = paste(V1,V2)) %>% ungroup() %>% select(comm),
            "~/GEO_data/SRA/split_run/run_binding/split_ba/rems_run",
            quote = F, row.names = F,col.names = F,sep = " ")

split_counts <- readRDS("~/Drug_splicing/data/split_counts.rds")
split_counts$unid <- Hmisc::escapeRegex(split_counts$unid)
split_counts$unid <- gsub(";","\\;",split_counts$unid,fixed = T)
rems_ba_counts <- split_counts %>% filter(unid %in% rems_ba_commd$unid)
rems_ba_counts$unid <- factor(rems_ba_counts$unid,levels = rems_ba_commd$unid)
rems_ba_counts <- rems_ba_counts %>% arrange(unid)
saveRDS(rems_ba_counts, "~/GEO_data/SRA/split_run/run_binding/split_ba/ba_rems.rds")

###el
files <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "all_split")
all_split <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,
                        pattern = "tsv")
all_el_done <- all_split[grepl("el_split",all_split)]

all_el_done <- data.frame(ids = all_el_done) %>% 
  tidyr::separate_wider_delim(cols = ids, delim = "/jcast_out/binding_out/res/", 
                              names = c("unid","sp"))
all_el_done <- all_el_done %>% group_by(unid) %>% summarise(done_counts = n()) 
all_counts <- readRDS("~/Drug_splicing/data/split_counts.rds")
rems_el <- left_join(all_counts,all_el_done) %>% 
  filter((done_counts != counts) | (is.na(done_counts)))
rems_el$unid <- Hmisc::escapeRegex(rems_el$unid)
rems_el$unid <- gsub(";","\\;",rems_el$unid,fixed = T)

all_el_commd <- read.table("~/GEO_data/SRA/split_run/run_binding/split_el/run_el_all.sh")
all_el_commd <- all_el_commd %>% 
  mutate(unid = gsub(".+gse_out/","",V2) %>% gsub("/jcast_out/.+","",.))
rems_el_commd <- all_el_commd %>% filter(unid %in% rems_el$unid)
write.table(rems_el_commd %>% select(1,2) %>% rowwise() %>% 
              mutate(comm = paste(V1,V2)) %>% ungroup() %>% select(comm),
            "~/GEO_data/SRA/split_run/run_binding/split_el/rems_run",
            quote = F, row.names = F,col.names = F,sep = " ")

split_counts <- readRDS("~/Drug_splicing/data/split_counts.rds")
split_counts$unid <- Hmisc::escapeRegex(split_counts$unid)
split_counts$unid <- gsub(";","\\;",split_counts$unid,fixed = T)
rems_el_counts <- split_counts %>% filter(unid %in% rems_el_commd$unid)
rems_el_counts$unid <- factor(rems_el_counts$unid,levels = rems_el_commd$unid)
rems_el_counts <- rems_el_counts %>% arrange(unid)
saveRDS(rems_el_counts, "~/GEO_data/SRA/split_run/run_binding/split_el/el_rems.rds")

