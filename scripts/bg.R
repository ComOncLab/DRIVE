# ####差异分析
# rm(list = ls())
# library(doParallel)
# library(foreach)
# library(dplyr)
# sample_meta <- readxl::read_xlsx("~/Drug_splicing/data/all_sample_meta.xlsx")
# do_diff_exp <- function(treat_group, control_group, exp_dt){
#   treat_group <- strsplit(treat_group,",")[[1]]
#   control_group <- strsplit(control_group,",")[[1]]
#   exp_dt <- exp_dt %>% 
#     tidyr::pivot_longer(cols = 2:ncol(exp_dt),names_to = "sample",values_to = "exp")
#   exp_dt <- left_join(
#     exp_dt,
#     data.frame(
#       sample = c(treat_group, control_group),
#       sample_type = c(rep("Treat",length(treat_group)),rep("Control",length(control_group)))
#     )
#   )
#   ###并行计算所有基因
#   all_genes <- unique(exp_dt$gene)
#   my.cluster <- parallel::makeCluster(
#     50, 
#     type = "PSOCK"
#   )
#   #register it to be used by %dopar%
#   doParallel::registerDoParallel(cl = my.cluster)
#   res <- foreach(
#     gene_i = all_genes,
#     .export = c("exp_dt"),
#     .packages = c("dplyr")
#   ) %dopar% {
#     tmp_exp <- exp_dt %>% filter(gene == gene_i)
#     tmp_res <- wilcox.test(exp~sample_type,data = tmp_exp)
#     return(
#       data.frame(
#         gene = gene_i,
#         p = tmp_res$p.value,
#         diff = median(tmp_exp$exp[tmp_exp$sample_type == "Treat"]) - median(tmp_exp$exp[tmp_exp$sample_type == "Control"])
#       )
#     )
#   }
#   parallel::stopCluster(cl = my.cluster)
#   res <- bind_rows(res)
#   return(res)
# }
# 
# for (i in 1:nrow(sample_meta)){
#   dt <- readRDS(paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/tpm.rds"))
#   diff_res <- do_diff_exp(sample_meta$treat_sample[i],
#                           sample_meta$control_sample[i],
#                           dt)
#   saveRDS(diff_res, paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/tpm_diff.rds"))
#   message("Complete ",i)
# }

# ####相关性
# library(doParallel)
# library(foreach)
# library(dplyr)
# cal_ken_par <- function(var_dt, ncores){
#   all_com <- combn(colnames(var_dt), 2, simplify = F)
#   all_com <- sapply(all_com, paste, collapse = "-")
#   ###对列进行并行
#   my.cluster <- parallel::makeCluster(
#     ncores, 
#     type = "PSOCK"
#   )
#   #register it to be used by %dopar%
#   doParallel::registerDoParallel(cl = my.cluster)
#   
#   res <- foreach(
#     i = all_com,
#     .export = c("var_dt"),
#     .packages = c("dplyr")
#   ) %dopar% {
#     com_tmp <- strsplit(i,"-")[[1]]
#     cor_tmp <- kendallknight::kendall_cor(var_dt[[com_tmp[1]]],var_dt[[com_tmp[2]]])
#     return(data.frame(coms = i, cor = cor_tmp))
#   }
#   parallel::stopCluster(cl = my.cluster)
#   res <- bind_rows(res)
#   return(res)
# }
# rank_res <- readRDS("~/Drug_splicing/data/deseq2_stat.rds")
# cor_res <- cal_ken_par(rank_res,ncores = 40)
# saveRDS(cor_res,"~/Drug_splicing/data/deseq2_cor.rds")

###########precalculated permutation results
library(dplyr)
deseq_fc <- readRDS("~/Drug_splicing/data/deseq_fc.rds")
##随机选100个Up 100 Down
cal_Xsum <- function(gene_sig, up_gene_name, down_gene_name, K=100){
  ##gene_sig is dataframe, first column is gene name, second is FC (log2FC)
  gene_sig <- gene_sig[order(gene_sig$FC,decreasing = TRUE,na.last = NA), ]
  sig_up_gene <- head(gene_sig$gene, K)
  sig_down_gene <- tail(gene_sig$gene, K)
  gene_sig_gene <- c(sig_up_gene, sig_down_gene)
  XUpInDisease <- intersect(up_gene_name, gene_sig_gene)
  XDownInDisease <- intersect(down_gene_name, gene_sig_gene)
  sum(gene_sig$FC[which(is.element(gene_sig$gene, XUpInDisease))], na.rm = TRUE) - sum(gene_sig$FC[which(is.element(gene_sig$gene, XDownInDisease))], na.rm = TRUE)
}

all_samples <- colnames(deseq_fc)
per_res <- vector("list",length(all_samples))
for (samples in 1:length(all_samples)){
  dt <- deseq_fc %>% select(all_samples[samples]) %>% rename(FC = 1)
  dt$gene <- rownames(dt)
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
    sampled_up_down <- sample(1:nrow(dt),200,replace = F)
    sampled_up <- sampled_up_down[1:100]
    sampled_down <- sampled_up_down[101:200]
    sampled_res <- cal_Xsum(dt, dt$gene[sampled_up], dt$gene[sampled_down])
    return(data.frame(sampled_id = i, sampled_xsum = sampled_res))
  }
  parallel::stopCluster(cl = my.cluster)
  res <- bind_rows(res)
  res <- res %>% select(-1)
  colnames(res) <- all_samples[samples]
  per_res[[samples]] <- res
  message("Complete ",samples)
}
per_res <- bind_cols(per_res)
saveRDS(per_res, "~/Drug_splicing/data/Xsum_permutation.rds")


