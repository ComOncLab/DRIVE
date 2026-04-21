library(dplyr)
library(tximport)
####先转化成基因表达
sample_meta <- readRDS("data/all_sample_meta.rds")
dt <- readRDS("~/Drug_splicing/data/human_id_mapping.rds")
all_sf <- list.files("~/GEO_data/SRA/out/",recursive = T,pattern = "quant.sf")
all_samples <- gsub("/.+","",all_sf)
library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(sample_meta),
  .export = c("sample_meta","dt","all_samples"),
  .packages = c("dplyr","tximport")
) %dopar% {
  gene_exp_dir <- paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/gene_exp")
  dir.create(gene_exp_dir)
  tmp_treat_sample <- strsplit(sample_meta$treat_sample[i],",")[[1]] %>% 
    intersect(.,all_samples)
  tmp_control_sample <- strsplit(sample_meta$control_sample[i],",")[[1]] %>% 
    intersect(.,all_samples)
  for (j in c(tmp_control_sample, tmp_treat_sample)){
    file.copy(from = paste0("~/GEO_data/SRA/out/",j,"/salmon_res/quant.sf"),
              to = paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/gene_exp/",j,".sf"))
  }
  salmon_files <- list.files(gene_exp_dir,full.names = T)
  names(salmon_files) <- gsub(".+/","",salmon_files) %>% gsub(".sf","",.)
  txi <- tximport(salmon_files, type = "salmon", tx2gene = dt %>% select(1,2))
  saveRDS(txi, paste0(gene_exp_dir,"/txi_out.rds"))
}
parallel::stopCluster(cl = my.cluster)

# ###所有样本的基因表达 
# sample_meta <- readxl::read_xlsx("data/all_sample_meta.xlsx")
# all_samples <- strsplit(paste(c(paste(sample_meta$treat_sample,collapse = ","),
#                                 paste(sample_meta$control_sample,collapse = ",")),
#                               collapse = ","),split = ",")[[1]] %>% unique()
# tt <- list.files("~/GEO_data/SRA/out/",recursive = T,pattern = "gene_counts.rds")
# tt <- gsub("/.+","",tt)
# ###其他的14个样本跑不出来，Fragment incompatibility prior below threshold.
# ###每行进行差异分析
# library(doParallel)
# library(foreach)
# #create the cluster
# my.cluster <- parallel::makeCluster(
#   40, 
#   type = "PSOCK"
# )
# #register it to be used by %dopar%
# doParallel::registerDoParallel(cl = my.cluster)
# 
# res <- foreach(
#   i = 1:nrow(sample_meta),
#   .export = c("sample_meta","tt"),
#   .packages = c("dplyr","IOBR")
# ) %dopar% {
#   tmp_treat_sample <- strsplit(sample_meta$treat_sample[i],",")[[1]] %>% 
#     intersect(., tt)
#   tmp_control_sample <- strsplit(sample_meta$control_sample[i],",")[[1]] %>% 
#     intersect(., tt)
#   counts_res <- sapply(c(tmp_control_sample, tmp_treat_sample),
#                        function(x){
#                          dt <- readRDS(paste0("/home/wt/GEO_data/SRA/out/",x,
#                                               "/salmon_res/gene_counts.rds"))
#                          colnames(dt)[2] <- x
#                          return(dt)
#                        },simplify = F)
#   counts_res <- Reduce(function(x, y) inner_join(x, y, by = "gene_id"), counts_res)
#   saveRDS(counts_res,paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/counts.rds"))
#   # rownames(counts_res) <- counts_res$gene_id
#   # counts_res$gene_id <- NULL
#   # counts_res <- as.matrix(counts_res) %>% round()
#   # tpm_res <- count2tpm(countMat = counts_res)
#   # tpm_res <- log2(as.matrix(tpm_res) + 1) %>% as.data.frame()
#   # tpm_res$gene <- rownames(tpm_res)
#   # tpm_res <- tpm_res %>% select(gene,everything())
#   # saveRDS(tpm_res, paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/tpm.rds"))
# }
# parallel::stopCluster(cl = my.cluster)

####差异分析，DESeq2
rm(list = ls())
library(dplyr)
sample_meta <- readRDS("data/all_sample_meta.rds")

do_diff_exp <- function(treat_group, control_group, txi_out){
  ###用txi 的输出进行差异分析
  txi_out <- readRDS(txi_out)
  treat_group <- strsplit(treat_group,",")[[1]] %>% 
    intersect(., colnames(txi_out$counts))
  control_group <- strsplit(control_group,",")[[1]] %>% 
    intersect(., colnames(txi_out$counts))
  coldata <- data.frame(ids = c(control_group, treat_group),
                        condition = c(rep("untreated",length(control_group)),
                                      rep("treated",length(treat_group))))
  sampleTable <- data.frame(ids = colnames(txi_out$counts)) %>% 
    left_join(., coldata)
  sampleTable$condition <- factor(sampleTable$condition)
  dds <- DESeq2::DESeqDataSetFromTximport(txi = txi_out,
                                          colData = sampleTable,
                                          design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = "untreated")
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds) %>% as.data.frame()
  res$gene_id <- rownames(res)
  return(res)
}

library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(sample_meta),
  .export = c("sample_meta"),
  .packages = c("dplyr")
) %dopar% {
  txi_file <- paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/gene_exp/txi_out.rds")
  diff_res <- do_diff_exp(treat_group = sample_meta$treat_sample[i],
                          control_group = sample_meta$control_sample[i],
                          txi_out = txi_file)
  saveRDS(diff_res, paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/deseq2_diff.rds"))
}
parallel::stopCluster(cl = my.cluster)

###合并全部
all_deseq <- list.files("~/GEO_data/SRA/gse_out/",pattern = "deseq2_diff.rds",recursive = T)
mapping <- readRDS("~/Drug_splicing/data/human_id_mapping.rds")
mapping <- mapping %>% select(gene_id_version, symbol) %>% distinct_all()
all_deseq_dt <- sapply(all_deseq,
                       function(x){
                         tmp <- readRDS(paste0("~/GEO_data/SRA/gse_out/",x))
                         tmp <- left_join(tmp,mapping %>% 
                                            rename(gene_id = gene_id_version) %>% 
                                            select(gene_id, symbol))
                         saveRDS(tmp,paste0("~/Drug_splicing/scripts/Shiny/data/",
                                            gsub("/.+","",x),"_deseq.rds"))
                         },
                       simplify = F)
###提取FC
all_deseq <- list.files("scripts/Shiny/data/",pattern = "deseq.rds",recursive = T)
library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:length(all_deseq),
  .export = c("all_deseq"),
  .packages = c("dplyr")
) %dopar% {
  tmp <- readRDS(paste0("~/Drug_splicing/scripts/Shiny/data/",all_deseq[i]))
  tmp <- tmp %>% 
    select(symbol, log2FoldChange) %>% 
    group_by(symbol) %>% 
    slice_max(abs(log2FoldChange), with_ties = F) %>% ungroup() %>% 
    filter(nchar(symbol) > 0) %>% as.data.frame()
  rownames(tmp) <- tmp$symbol
  tmp$symbol <- NULL
  colnames(tmp) <- gsub("_deseq.rds","",all_deseq[i])
  return(tmp)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_cols(res)
saveRDS(res,"~/Drug_splicing/data/deseq_fc.rds")

###差异通路富集分析
sample_meta <- readRDS("data/all_sample_meta.rds")
mapping <- readRDS("~/Drug_splicing/data/human_id_mapping.rds")
mapping <- mapping %>% select(gene_id_version, symbol) %>% distinct_all()
go_pathways <- fgsea::gmtPathways("data/c5.go.v2026.1.Hs.symbols.gmt")

library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(sample_meta),
  .export = c("sample_meta","mapping","go_pathways"),
  .packages = c("dplyr","fgsea")
) %dopar% {
  tmp <- readRDS(paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/deseq2_diff.rds"))
  tmp <- tmp %>%
    filter(!is.na(stat)) %>% 
    left_join(.,mapping %>% rename(gene_id = gene_id_version) %>% 
                select(gene_id, symbol)) %>% 
    group_by(symbol) %>% 
    slice_max(abs(stat), with_ties = F) %>% ungroup() %>% 
    filter(nchar(symbol) > 0)
  ranks <- tmp$stat
  names(ranks) <- tmp$symbol
  gsea_res <- fgsea::fgsea(pathways = go_pathways, stats = ranks, nperm = 1000)
  saveRDS(gsea_res, paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/diff_gsea.rds"))
}
parallel::stopCluster(cl = my.cluster)

#####hallmarks genes
# hallmarks <- data.table::fread("data/Menyhart_JPA_CancerHallmarks_core.txt",data.table = F)
# rownames(hallmarks) <- hallmarks$V1
# hallmarks$V1 <- NULL
# hallmarks <- as.data.frame(t(hallmarks))
# hallmarks <- lapply(hallmarks, function(x){x[which(nchar(x) > 0)]})
apm <- fgsea::gmtPathways("data/KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION.v2026.1.Hs.gmt")
# hallmarks <- c(hallmarks,immune)
# metasta_genes <- readxl::read_xlsx("data/met_genes.xlsx",sheet = 2,col_names = F) %>% 
#   rename(genes = 1)
# meta_path <- list(meta = metasta_genes$genes)
escape_up <- list(UP_Immune_Escape = c("TGFB1", "IL10", "VEGFA", "CD274", "PVR","PDCD1LG2", 
                                       "NECTIN2", "NECTIN3", "CD80", "CD86", 
                                       "CD47", "TNFSF4", "TNFSF9", "CCL2","TNFRSF14",
                                       "CSF1", "CEACAM1", "LGALS9",
                                       ##LIN_TUMOR_ESCAPE_FROM_IMMUNE_ATTACK
                                       "ACKR3","C14orf39", "CHD7", "DYNAP", "GJB4", "GPR149", 
                                       "HTRA1", "IL2RG", "JAG1", "KHDRBS3", "MMD", "NPPB", 
                                       "NRP2", "PLA2G7", "RGS16", "ST3GAL6", "SYCP3", 
                                       "TNNT2", "VCAM1"))
# escape_down <- readxl::read_xlsx("data/Table_S23_In vivo screen results.xlsx",sheet = 3)
# escape_down <- escape_down %>% filter(library.annotation == "core") %>% 
#   filter(inVitro.class == "Suppressor")
# library(nichenetr)
# escape_down <- nichenetr::convert_mouse_to_human_symbols(escape_down$Gene) %>% unname()
# escape_down <- escape_down[!is.na(escape_down)]
# saveRDS(escape_down,"data/nature_escape_gene.rds")
escape_down <- readRDS("data/nature_escape_gene.rds")
escape_down <- list(DOWN_Immune_Escape = unique(c(escape_down, apm[[1]])))
immune_pathway <- c(apm, escape_down, escape_up)

sample_meta <- readRDS("data/all_sample_meta.rds")
mapping <- readRDS("~/Drug_splicing/data/human_id_mapping.rds")
mapping <- mapping %>% select(gene_id_version, symbol) %>% distinct_all()

library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(sample_meta),
  .export = c("sample_meta","mapping","immune_pathway"),
  .packages = c("dplyr","fgsea")
) %dopar% {
  tmp <- readRDS(paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/deseq2_diff.rds"))
  tmp <- tmp %>%
    filter(!is.na(stat)) %>% 
    left_join(.,mapping %>% rename(gene_id = gene_id_version) %>% 
                select(gene_id, symbol)) %>% 
    group_by(symbol) %>% 
    slice_max(abs(stat), with_ties = F) %>% ungroup() %>% 
    filter(nchar(symbol) > 0)
  ranks <- tmp$stat
  names(ranks) <- tmp$symbol
  gsea_res <- fgsea::fgsea(pathways = immune_pathway, stats = ranks, nperm = 1000)
  gsea_res <- gsea_res %>% as.data.frame() %>% mutate(id = sample_meta$unid[i])
  return(gsea_res)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_rows(res) 
saveRDS(res,"data/immune_path_gsea.rds")



