library(dplyr)
deseq_fc <- readRDS("~/Drug_splicing/data/deseq_fc.rds")
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
per_res <- bind_cols(per_res)
saveRDS(per_res, "~/Drug_splicing/data/Zhang_permutation.rds")




