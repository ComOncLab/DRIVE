library(dplyr)
##GSE174498_CEM_L-asparaginase 没有输出结果
# sample_meta <- readRDS("~/DRIVE/data/all_sample_meta.rds")
# sample_meta <- sample_meta %>% filter(unid != "GSE174498_CEM_L-asparaginase")
rmats_counts <- readRDS("~/DRIVE/data/rmats_filter.rds")
rmats_counts <- rmats_counts %>% filter(sample != "GSE174498_CEM_L-asparaginase")
rmats_counts <- rmats_counts %>% filter((abs(IncLevelDifference) > 0.1) & (FDR <0.05))
all_samples <- unique(rmats_counts$sample)
samples <- readRDS("~/DRIVE/data/samples_with_HLA.rds")
all_samples <- intersect(all_samples, samples$unid)
# done_files <- list.files("scripts/Shiny/data/",pattern = "pep_meta.rds",full.names = T)
# done_meta_info <- file.info(done_files)
# done_meta_info <- done_meta_info %>% filter(grepl("2026-04-23",mtime))
# done_files <- list.files("scripts/Shiny/data/",pattern = "pep.rds",full.names = T)
# done_info <- file.info(done_files)
# done_info <- done_info %>% filter(grepl("2026-04-22",mtime))
# done_samples <- intersect(gsub(".+//","",rownames(done_meta_info)) %>% 
#                             gsub("_binding_pep_meta.rds","",.),
#                           gsub(".+//","",rownames(done_info)) %>%
#                             gsub("_binding_pep.rds","",.))
# all_samples <- intersect(all_samples, samples$unid)
# all_samples <- all_samples[which(!(all_samples %in% done_samples))]

library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  10, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:length(all_samples),
  .packages = c("dplyr")
) %dopar% {
  all_files <- list.files(paste0("~/GEO_data/SRA/gse_out/",all_samples[i],
                                 "/jcast_out/binding_out/res/"),pattern = "tsv",
                          full.names = T)
  dt <- sapply(all_files,
               function(x){
                 tmp <- data.table::fread(x,data.table = F)
                 colnames(tmp)[5] <- "score"
                 tmp$binding_type <- gsub(".+//","",x) %>% gsub("_.+","",.)
                 return(tmp)
               },simplify = F)
  dt <- bind_rows(dt) %>% distinct_all()
  dt_filter <- dt %>% 
    select(allele, peptide, percentile, binding_type) %>% 
    filter(percentile < 2) %>% distinct_all()
  dt_filter_summ <- dt_filter %>% 
    group_by(peptide, allele) %>% 
    summarise(both_sb = n()) %>% ungroup() %>% filter(both_sb == 2) ###ba and el < 2
  dt_filter <- dt_filter %>% filter(peptide %in% dt_filter_summ$peptide)
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
    inner_join(., dt3 %>% rename(rMTS_ID = ID, rMATS_type2 = type) %>% 
                mutate(rMTS_ID = as.character(rMTS_ID))) %>% ##FDR < 0.05 psi > 0.1 
    select(-GeneID) %>% 
    select(seq, pep_len, rMTS_ID, sample, everything())
  dt_all <- dt_all %>% filter(peptide %in% dt_meta$seq)
  saveRDS(dt_meta, 
          paste0("~/DRIVE/scripts/Shiny/data/",
                 all_samples[i], "_binding_pep_meta.rds"))
  saveRDS(dt_all, 
          paste0("~/DRIVE/scripts/Shiny/data/",
                 all_samples[i], "_binding_pep.rds"))
}
parallel::stopCluster(cl = my.cluster)


######统计
all_files <- list.files("scripts/Shiny/data/",pattern = "_pep.rds")
sample_meta <- readRDS("~/DRIVE/data/all_sample_meta.rds")
sample_meta <- sample_meta %>% 
  filter(unid %in% gsub("_binding_pep.rds","",all_files))

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
  i = sample_meta$unid,
  .packages = c("dplyr")
) %dopar% {
  dt <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",i,"_binding_pep.rds"))
  dt <- dt %>%
    select(-binding_counts) %>% 
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "allel", values_to = "rank")
  dt_filter <- dt %>% filter(rank < 0.5)
  dt_filter_summ <- dt_filter %>% 
    mutate(allel = gsub("_.+","",allel)) %>% 
    group_by(peptide, allel) %>% 
    summarise(both_sb = n()) %>% ungroup() %>% filter(both_sb == 2) 
  dt_filter <- dt_filter %>% filter(peptide %in% dt_filter_summ$peptide)
  data.frame(
    unid = i,
    sb_counts = length(unique(dt_filter$peptide))
  )
}
parallel::stopCluster(cl = my.cluster)

res <- bind_rows(res)
sample_meta <- left_join(sample_meta, res)

spres <- readRDS("~/DRIVE/data/rmats_fdr01.rds")
spres <- spres %>% 
  group_by(ids) %>% summarise(t_counts = sum(counts)) %>% ungroup()
sample_meta <- left_join(sample_meta, spres %>% rename(unid = ids))
saveRDS(sample_meta,"data/sample_neo.rds")

###
sample_meta <- readRDS("data/sample_neo.rds")
sample_meta <- sample_meta %>% mutate(log2counts = log2(sb_counts))
ord <- sample_meta %>% group_by(Tissue_Source2) %>% 
  summarise(med = median(log2counts)) %>% 
  ungroup() %>% arrange(desc(med))
p1 <- ggboxplot(sample_meta,x="Tissue_Source2",y="log2counts",add = "jitter",
                color = "Tissue_Source2",order = ord$Tissue_Source2,xlab = F)+
  theme(legend.position = "none")+
  rotate_x_text(45)

dt <- sample_meta %>% filter(!is.na(DrugClass))
ord <- dt %>% group_by(DrugClass) %>% summarise(med = median(log2counts)) %>% 
  ungroup() %>% arrange(desc(med))
p2 <- ggboxplot(dt,x="DrugClass",y="log2counts",add = "jitter",
                color = "DrugClass",order = ord$DrugClass,xlab = F)+
  theme(legend.position = "none")+
  rotate_x_text(90)
library(patchwork)
p1 / p2
ggsave("Figs/cell_drug_neocounts.pdf",width = 7,height = 8)

###药物排名，至少三个细胞系
dt_summ <- sample_meta %>% group_by(Drug) %>% 
  summarise(counts = length(unique(cell_line))) %>% ungroup() %>% 
  filter(counts >= 3)
dt <- sample_meta %>% filter(Drug %in% dt_summ$Drug)
drug_summ <- dt %>% group_by(Drug) %>% 
  summarise(median_c = median(log2counts)) %>% ungroup() %>% 
  arrange(desc(median_c)) %>% 
  slice_head(n=10)
dt_filter <- dt %>% filter(Drug %in% drug_summ$Drug)

p1 <- ggboxplot(dt_filter,x="Drug",y="log2counts",add = "jitter",
          color = "Drug",order = drug_summ$Drug, ylab = "Log2Counts")+
  theme(legend.position = "none")+
  rotate_x_text(45)

dt <- dt %>% 
  mutate(norm_counts = sb_counts / t_counts)
drug_summ <- dt %>% group_by(Drug) %>% 
  summarise(median_c = median(norm_counts)) %>% ungroup() %>% 
  arrange(desc(median_c)) %>% 
  slice_head(n=10)
dt_filter <- dt %>% filter(Drug %in% drug_summ$Drug)
p2 <- ggboxplot(dt_filter,x="Drug",y="norm_counts",add = "jitter",
          color = "Drug",order = drug_summ$Drug, 
          ylab = "Normalized Counts (Neoantigen / Differential AS)")+
  theme(legend.position = "none")+
  rotate_x_text(45)
p1 / p2
ggsave("Figs/top10_drug_neocounts.pdf",width = 7,height = 8)

sample_meta <- sample_meta %>% mutate(log2AScounts = log2(t_counts))
ggscatter(sample_meta, x = "log2counts", y = "log2AScounts",
          color = "black", shape = 21, size = 3,
          add = "reg.line",  
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.coeff.args = list(method = "pearson", label.x = 10, label.sep = "\n"),
          xlab = "Neoantigen Counts", ylab = "Differential AS Counts"
)
ggsave("Figs/cor_neo_as.pdf",width = 6, height = 5)

#########不同类型的可变剪切导致的新抗原数量
all_files <- list.files("scripts/Shiny/data/",pattern = "_pep_meta.rds")
sample_meta <- readRDS("~/DRIVE/data/all_sample_meta.rds")
sample_meta <- sample_meta %>% 
  filter(unid %in% gsub("_binding_pep_meta.rds","",all_files))

my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = sample_meta$unid,
  .packages = c("dplyr")
) %dopar% {
  dt_meta <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",i,"_binding_pep_meta.rds"))
  dt <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",i,"_binding_pep.rds"))
  dt <- dt %>%
    select(-binding_counts) %>% 
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "allel", values_to = "rank")
  dt_filter <- dt %>% filter(rank < 0.5)
  dt_filter_summ <- dt_filter %>% 
    mutate(allel = gsub("_.+","",allel)) %>% 
    group_by(peptide, allel) %>% 
    summarise(both_sb = n()) %>% ungroup() %>% filter(both_sb == 2) 
  dt_filter <- dt_filter %>% filter(peptide %in% dt_filter_summ$peptide)
  
  dt_meta <- dt_meta %>% filter(seq %in% dt_filter$peptide)
  dt_summ <- dt_meta %>% 
    group_by(rMATS_type2) %>% 
    summarise(counts = length(unique(seq))) %>% ungroup()
  dt_summ$unid <- i
  return(dt_summ)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_rows(res)
saveRDS(res, "data/AS_type_neo.rds")

###
res <- readRDS("data/AS_type_neo.rds")
res <- res %>% mutate(log2counts = log2(counts))
res_summ <- res %>% group_by(rMATS_type2) %>% 
  summarise(median_c = median(counts)) %>% ungroup() %>% 
  arrange(desc(median_c))
my_c <- list(c("SE","A3SS"),c("SE","RI"),c("SE","MXE"),c("SE","A5SS"),
             c("A3SS","RI"),c("A3SS","MXE"),c("A3SS","A5SS"),
             c("RI","MXE"),c("RI","A5SS"),
             c("MXE","A5SS"))
ggviolin(res,x="rMATS_type2",y="log2counts",add = "boxplot",
          color = "rMATS_type2",order = res_summ$rMATS_type2, 
          xlab = F)+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = my_c)
ggsave("Figs/as_type_neo_diff.pdf",width = 6,height = 6)  
  
####桑吉图
rmats_counts <- readRDS("~/DRIVE/data/rmats_filter.rds")
rmats_counts <- rmats_counts %>% filter(sample != "GSE174498_CEM_L-asparaginase")
rmats_counts <- rmats_counts %>% filter((abs(IncLevelDifference) > 0.1) & (FDR <0.05))

all_files <- list.files("scripts/Shiny/data/",pattern = "_pep_meta.rds")
sample_meta <- readRDS("~/DRIVE/data/all_sample_meta.rds")
sample_meta <- sample_meta %>% 
  filter(unid %in% gsub("_binding_pep_meta.rds","",all_files))

my.cluster <- parallel::makeCluster(
  10, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = sample_meta$unid,
  .packages = c("dplyr")
) %dopar% {
  dt_meta <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",i,"_binding_pep_meta.rds"))
  dt <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",i,"_binding_pep.rds"))
  dt <- dt %>%
    select(-binding_counts) %>% 
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "allel", values_to = "rank")
  dt_filter <- dt %>% filter(rank < 0.5)
  dt_filter_summ <- dt_filter %>% 
    mutate(allel = gsub("_.+","",allel)) %>% 
    group_by(peptide, allel) %>% 
    summarise(both_sb = n()) %>% ungroup() %>% filter(both_sb == 2) 
  dt_filter <- dt_filter %>% filter(peptide %in% dt_filter_summ$peptide)
  
  dt_meta_filter <- dt_meta %>% filter(seq %in% dt_filter$peptide)
  dt_binding <- dt_meta_filter %>% 
    group_by(rMATS_type2, rMTS_ID) %>% 
    summarise(binding = length(unique(seq))) %>% ungroup()
  
  dt_mer <- readRDS(paste0("~/GEO_data/SRA/gse_out/",i,"/jcast_out/binding_out/all_mers.rds"))
  dt_mer <- dt_mer %>% filter(!(seq %in% dt_meta$seq)) ###ba 或者 el 至少有一个 大于 2
  dt_mer <- dt_mer %>% 
    tidyr::separate_wider_delim(seq_id, delim = "|",
                                names = c("sp","UniProt_ID","UniProt_symbol","Gene_ID",
                                          "rMATS_type","rMTS_ID","Chr","Anch_Exon_SE",
                                          "Alt_Exon_SE","strand_phase","Msjc","Tier")) %>% 
    select(-sp) %>% rename(pep_len = type) %>% 
    mutate(rMATS_type2 = sub("[0-9]$", "", rMATS_type))
  dt3 <- rmats_counts %>% filter(sample %in% i)
  dt_mer <- dt_mer %>% distinct_all() %>% 
    inner_join(., dt3 %>% rename(rMTS_ID = ID, rMATS_type2 = type) %>% 
                 mutate(rMTS_ID = as.character(rMTS_ID))) ##FDR < 0.05 psi > 0.1 
  
  dt_nobind <- dt_mer %>% 
    group_by(rMATS_type2, rMTS_ID) %>% 
    summarise(no_binding = length(unique(seq))) %>% ungroup()
  dt_all <- full_join(dt_binding, dt_nobind)
  dt_all$unid <- i
  return(dt_all)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_rows(res)
saveRDS(res,"data/neo_as_counts.rds")

####
res <- readRDS("data/neo_as_counts.rds")
res <- res %>% 
  tidyr::pivot_longer(cols = c(binding, no_binding), 
                      names_to = "type", values_to = "counts")
res <- res %>% mutate(counts = ifelse(is.na(counts),0,counts))
res_summ <- res %>% 
  group_by(rMATS_type2, type) %>% 
  summarise(total_counts = sum(counts,na.rm = T)) %>% 
  ungroup()

library(ggalluvial)
ggplot(data = res_summ,
       aes(axis1 = rMATS_type2, axis2 = type,
           y = total_counts)) +
  scale_x_discrete(limits = c("rMATS_type2", "type"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = type)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()+
  scale_fill_npg()+
  guides(fill = "none")
ggsave("Figs/as_neo_Alluvial.pdf",width = 6,height = 7)

res_summ %>% 
  group_by(rMATS_type2) %>%
  summarise(binding_per = total_counts[type == "binding"]/total_counts[type == "no_binding"])
# A tibble: 5 × 2
# rMATS_type2 binding_per
# <chr>             <dbl>
#   1 A3SS             0.0206
# 2 A5SS             0.0203
# 3 MXE              0.0208
# 4 RI               0.0212
# 5 SE               0.0209

##############人类常见HLA，考虑基因表达
samples <- readRDS("~/DRIVE/data/samples_with_HLA.rds")
top_hlas <- samples %>% 
  filter(grepl("HLA-A*02:01",HLA,fixed = T) | grepl("HLA-A*24:02",HLA,fixed = T) | grepl("HLA-C*04:01",HLA,fixed = T))
all_samples <- unique(top_hlas$unid)

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
  i = all_samples,
  .packages = c("dplyr")
) %dopar% {
  dt <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",i,"_binding_pep.rds"))
  dt <- dt %>%
    select(-binding_counts) %>% 
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "allel", values_to = "rank")
  dt_filter <- dt %>% 
    filter(gsub("_.+","",allel) %in% c("HLA-A*02:01","HLA-A*24:02","HLA-C*04:01"))
  dt_filter <- dt_filter %>% filter(rank < 0.5)
  dt_filter_summ <- dt_filter %>% 
    mutate(allel = gsub("_.+","",allel)) %>% 
    group_by(peptide, allel) %>% 
    summarise(both_sb = n()) %>% ungroup() %>% filter(both_sb == 2) 
  dt_filter <- dt_filter %>% filter(peptide %in% dt_filter_summ$peptide)
  
  dt_meta <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",i,"_binding_pep_meta.rds"))
  dt_meta <- dt_meta %>% filter(seq %in% dt_filter$peptide)
  
  gene_exp <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",i,"_deseq.rds"))
  gene_exp <- gene_exp %>% 
    select(symbol, log2FoldChange) %>% 
    group_by(symbol) %>% 
    slice_max(abs(log2FoldChange), with_ties = F, na_rm = T) %>% ungroup() %>% 
    filter(nchar(symbol) > 0) %>% as.data.frame()
  dt_meta <- left_join(dt_meta, gene_exp %>% rename(geneSymbol = symbol))
  dt_meta <- dt_meta %>% group_by(seq) %>% 
    summarise(max_fc = max(log2FoldChange,na.rm = T),
              max_fc_gene = geneSymbol[which.max(log2FoldChange)]) %>% ungroup()
  
  dt_filter <- left_join(dt_filter, dt_meta %>% rename(peptide = seq))
  dt_filter <- dt_filter %>% filter(max_fc > 1)
  tmp_res <- dt_filter %>% 
    mutate(allel = gsub("_.+","",allel)) %>% 
    group_by(allel) %>% summarise(neo_counts = length(unique(peptide)))
  tmp_res$unid <- i
  return(tmp_res)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_rows(res)
saveRDS(res,"data/topHLA_neo.rds")

###
res <- readRDS("data/topHLA_neo.rds")
sample_meta <- readRDS("~/DRIVE/data/all_sample_meta.rds")
res <- left_join(res,sample_meta)

###至少有三个细胞系的药物
res_summ <- res %>% group_by(Drug) %>% summarise(cell_c = length(unique(cell_line))) %>% 
  ungroup() %>% filter(cell_c >= 3)
res_filter <- res %>% filter(Drug %in% res_summ$Drug)
dt <- res_filter %>% 
  group_by(allel, Drug) %>% 
  summarise(mean_neo_c = median(log2((neo_counts)))) %>% ungroup()
dt <- dt %>% 
  tidyr::pivot_wider(names_from = Drug, values_from = mean_neo_c) %>% as.data.frame()
rownames(dt) <- dt$allel
dt$allel <- NULL
dt <- dt %>% t() %>% as.data.frame() %>% 
  arrange(desc(`HLA-A*02:01`)) %>% t() %>% as.data.frame()
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(3, 5, 12), c("#FEE5C1", "#FC8C59", "#B70503"))
p <- Heatmap(dt,cluster_rows = F,cluster_columns = F,col = col_fun,
        name = "Median of log2(#Neo)",border="black",
        rect_gp = gpar(col = "white", lwd = 2),column_names_max_height = unit(8, "cm"))
pdf("Figs/top3hla_drug.pdf",width = 13,height = 4)
draw(p)
dev.off()

############onvansertib 药物分析
sample_meta <- readRDS("~/DRIVE/data/all_sample_meta.rds")
samples <- readRDS("~/DRIVE/data/samples_with_HLA.rds")
onv_sample <- sample_meta %>% filter(Drug == "onvansertib")
onv_hla <- samples %>% filter(unid %in% onv_sample$unid)
indi_sample <- sample_meta %>% filter(Drug == "Indisulam")
indi_hla <- samples %>% filter(unid %in% indi_sample$unid)
overlap_hla <- intersect(paste0(indi_hla$HLA,collapse = ",") %>% strsplit(.,",") %>% `[[`(1),
                         paste0(onv_hla$HLA,collapse = ",") %>% strsplit(.,",") %>% `[[`(1))
all_samples <- c(onv_hla$unid, indi_hla$unid) %>% unique()
rmats_counts <- readRDS("~/DRIVE/data/rmats_filter.rds")
rmats_counts <- rmats_counts %>% filter((abs(IncLevelDifference) > 0.1) & (FDR <0.05))
dt_res <- rmats_counts %>% 
  filter(sample %in% all_samples)

res <- vector("list",length(all_samples))
for (i in 1:length(all_samples)){
  dt <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",all_samples[i],"_binding_pep.rds"))
  dt <- dt %>%
    select(-binding_counts) %>% 
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "allel", values_to = "rank")
  dt_filter <- dt %>% 
    filter(gsub("_.+","",allel) %in% overlap_hla)
  dt_filter <- dt_filter %>% filter(rank < 0.5)
  dt_filter_summ <- dt_filter %>% 
    mutate(allel = gsub("_.+","",allel)) %>% 
    group_by(peptide, allel) %>% 
    summarise(both_sb = n()) %>% ungroup() %>% filter(both_sb == 2) 
  dt_filter <- dt_filter %>% filter(peptide %in% dt_filter_summ$peptide)
  
  dt_meta <- readRDS(paste0("~/DRIVE/scripts/Shiny/data/",all_samples[i],
                            "_binding_pep_meta.rds"))
  dt_meta <- dt_meta %>% filter(seq %in% dt_filter$peptide)
  dt_meta <- dt_meta %>% mutate(ids = paste0(rMATS_type2,"_",rMTS_ID))
  
  dt_mer <- readRDS(paste0("~/GEO_data/SRA/gse_out/",all_samples[i],
                           "/jcast_out/binding_out/all_mers.rds"))
  dt_mer <- dt_mer %>% 
    tidyr::separate_wider_delim(seq_id, delim = "|",
                                names = c("sp","UniProt_ID","UniProt_symbol","Gene_ID",
                                          "rMATS_type","rMTS_ID","Chr","Anch_Exon_SE",
                                          "Alt_Exon_SE","strand_phase","Msjc","Tier")) %>% 
    select(-sp) %>% rename(pep_len = type) %>% 
    mutate(rMATS_type2 = sub("[0-9]$", "", rMATS_type))
  dt_mer <- dt_mer %>% mutate(ids = paste0(rMATS_type2,"_",rMTS_ID))
  
  dt3 <- rmats_counts %>% filter(sample %in% all_samples[i]) %>% 
    mutate(ids = paste0(type,"_",ID))
  dt_mer <- dt_mer %>% filter(ids %in% dt3$ids)
  dt_mer <- dt_mer %>% filter(ids %in% dt_meta$ids)
  
  dt_meta$dev_events <- length(unique(dt_mer$ids))
  res[[i]] <- dt_meta
}

res <- bind_rows(res)
saveRDS(res,"data/Indisulam_onvansertib_neoas.rds")
###
res <- readRDS("data/Indisulam_onvansertib_neoas.rds")
res_summ <- res %>% group_by(sample) %>% 
  summarise(neo_counts = length(unique(seq)),
            as_counts = unique(dev_events)) %>% ungroup() %>% 
  mutate(effc = neo_counts/as_counts) %>% 
  mutate(Drug = gsub(".+_","",sample))
ggboxplot(res_summ,x="Drug",y="effc",add = "jitter",color = "Drug",
          ylab = "#(Binding Peptides) / #(Differential AS)")+
  stat_compare_means()+
  scale_color_npg()+
  guides(color = 'none')
ggsave("Figs/Indisulam_onvansertib_neoas_compare.pdf",width = 4,height = 4)

##来自哪些基因
dt1 <- res %>% 
  filter(grepl("onvansertib",sample))
onv_genes <- unique(dt1$geneSymbol)
dt2 <- res %>% 
  filter(grepl("indisulam",sample))
ind_genes <- unique(dt2$geneSymbol)




