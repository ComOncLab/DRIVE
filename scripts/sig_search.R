###改造WebCmap包，输入两个signature，可以算多个相似性值
###构建database，FC 或者 stat
all_deseq <- list.files("~/GEO_data/SRA/gse_out/",pattern = "deseq2_diff.rds",recursive = T)
mapping <- readRDS("~/Drug_splicing/data/human_id_mapping.rds")
mapping <- mapping %>% select(gene_id_version, symbol) %>% distinct_all()

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
  i = all_deseq,
  .export = c("mapping"),
  .packages = c("dplyr")
) %dopar% {
  tmp <- readRDS(paste0("~/GEO_data/SRA/gse_out/",i))
  tmp <- left_join(tmp,mapping %>% 
                     rename(gene_id = gene_id_version) %>% 
                     select(gene_id, symbol))
  tmp_fc <- tmp %>% select(symbol, stat) %>% 
    group_by(symbol) %>% 
    slice_max(stat, with_ties = F) %>% ungroup() %>% 
    filter(nchar(symbol) > 0) %>% as.data.frame()
  colnames(tmp_fc)[2] <- gsub("/deseq2_diff.rds","",i)
  rownames(tmp_fc) <- tmp_fc$symbol
  tmp_fc$symbol <- NULL
  return(tmp_fc)
}
parallel::stopCluster(cl = my.cluster)
all_deseq_dt <- Reduce(bind_cols,res)
saveRDS(all_deseq_dt,"data/deseq_stat.rds")

###########precalculated permutation results
deseq_fc <- readRDS("~/Drug_splicing/data/deseq_fc.rds")
##随机选100个Up 100 Down 
cal_Xsum <- function(gene_sig, up_gene_idx, down_gene_idx){
  ##gene_sig 是数值 vector
  sum(gene_sig[up_gene_idx], na.rm = TRUE) - sum(gene_sig[down_gene_idx], na.rm = TRUE)
}

all_samples <- colnames(deseq_fc)
per_res <- vector("list",length(all_samples))
for (samples in 1:length(all_samples)){
  dt <- deseq_fc %>% select(all_samples[samples])
  dt <- dt[,1] %>% na.omit()
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
    i = 1:10000,
    .export = c("cal_Xsum","dt"),
    .packages = c("dplyr")
  ) %dopar% {
    sampled_up_down <- sample(1:length(dt),200,replace = F)
    sampled_up <- sampled_up_down[1:100]
    sampled_down <- sampled_up_down[101:200]
    sampled_res <- cal_Xsum(dt, sampled_up, sampled_down)
    return(data.frame(sampled_id = i, sampled_xsum = sampled_res))
  }
  parallel::stopCluster(cl = my.cluster)
  res <- bind_rows(res)
  res <- res %>% select(-1)
  colnames(res) <- all_samples[samples]
  per_res[[samples]] <- res
  message("Complete ",samples)
}

########
cal_zhang <- function(gene_sig, up_gene_name, down_gene_name){
  ##gene_sig is dataframe, first column is gene name, second is FC (log2FC)
  #1 Convert the gene signature vector to ranked list
  refSort <- gene_sig[order(abs(gene_sig$FC), decreasing =TRUE, na.last = NA), ]
  refSort$refRank <- rank(abs(refSort$FC)) * sign(refSort$FC)
  #2 Compute the maximal theoretical score
  queryVector <- c(rep(1, length(up_gene_name)), rep(-1, length(down_gene_name)))
  names(queryVector) <- c(up_gene_name, down_gene_name)
  maxTheoreticalScore <- sum(abs(refSort$refRank)[1:length(queryVector)] * abs(queryVector))
  #3 The final score
  Connection_strength_score <- sum(queryVector * refSort$refRank[match(names(queryVector), refSort$gene)], 
                                   na.rm = TRUE) / maxTheoreticalScore
  return(Connection_strength_score)
}

all_samples <- colnames(deseq_fc)
per_res <- vector("list",length(all_samples))
for (samples in 1:length(all_samples)){
  dt <- deseq_fc %>% select(all_samples[samples])
  dt$gene <- rownames(dt)
  colnames(dt)[1] <- "FC"
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
    i = 1:10000,
    .export = c("cal_zhang","dt"),
    .packages = c("dplyr")
  ) %dopar% {
    sampled_up_down <- sample(1:nrow(dt),200,replace = F)
    sampled_up <- sampled_up_down[1:100]
    sampled_down <- sampled_up_down[101:200]
    sampled_res <- cal_zhang(dt, dt$gene[sampled_up], dt$gene[sampled_down])
    return(data.frame(sampled_id = i, sampled_zhang = sampled_res))
  }
  parallel::stopCluster(cl = my.cluster)
  res <- bind_rows(res)
  res <- res %>% select(-1)
  colnames(res) <- all_samples[samples]
  per_res[[samples]] <- res
  message("Complete ",samples)
}

#####




