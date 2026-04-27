library(dplyr)
library(ggplot2)
library(ggpubr)
###Figure1
###测序数据统计
dt <- read.table("~/GEO_data/SRA/meta_info.txt",header = T)
all_sample_meta <- readxl::read_xlsx("data/all_sample_meta.xlsx")
all_samples <- strsplit(paste(c(paste(all_sample_meta$treat_sample,collapse = ","),
                                paste(all_sample_meta$control_sample,collapse = ",")),
                              collapse = ","),
                        ",")[[1]] %>% unique()
dt <- dt %>% filter(experiment_geo_accession %in% all_samples)
dt_summ <- dt %>% group_by(pair) %>% summarise(counts = n()) %>% ungroup()
ggbarplot(dt_summ,x="pair",y="counts",xlab = "Sequencing Strategies",
          ylab = "Sample Counts", palette = c("#707DA6","#CCAD9D"),fill = "pair")+
  theme(legend.position = "none")

###数据量
all_html <- list.files("~/GEO_data/SRA/out/",recursive = T,pattern = ".html")
for (i in 1:length(all_html)){
  file.copy(paste0("~/GEO_data/SRA/out/",all_html[i]),
            paste0("~/GEO_data/SRA/all_fq_html/",gsub(".+/","",all_html[i])))
}
all_srrs <- strsplit(paste(dt$srrs,collapse = ","),",")[[1]] %>% unique() 
###
all_html <- list.files("~/GEO_data/SRA/all_fq_html/",pattern = ".txt")
res <- sapply(all_html,
              function(x){
                tt <- read.table(paste0("~/GEO_data/SRA/all_fq_html/",x),sep = "\t")
                data.frame(srr = gsub(".txt","",x),
                           read_counts = tt$V1[1],
                           read_len = tt$V1[2])
              },simplify = F)
res <- bind_rows(res)
res <- res %>% 
  rowwise() %>% 
  mutate(uni = substr(read_counts, nchar(read_counts), nchar(read_counts))) %>% 
  mutate(read_counts2 = ifelse(
    uni == "M", as.numeric(gsub(" M","",read_counts)),
    as.numeric(gsub(" K","",read_counts)) / 1000
  )) %>% ungroup()
res$read_len <- as.numeric(res$read_len)
res <- res %>% filter(srr %in% all_srrs)
saveRDS(res,"data/srr_info.rds")

#######细胞系癌症类型分布，药物类型分布
library(ggpie)
cell_line_meta <- readRDS("~/DRIVE/data/cell_line_meta.rds")
#cell_line_meta$cell_line[35] <- "92-1"
#saveRDS(cell_line_meta,"data/cell_line_meta.rds")
p1 <- ggpie(data = cell_line_meta, group_key = "Tissue_Source2", count_type = "full",
      label_info = c("group","count"), label_type = "horizon",
      label_split = "\n",
      label_size = 4, label_pos = "out" )+
  theme(legend.position = "none")

drug_info <- readxl::read_xlsx("data/drug_class.xlsx")
p2 <- ggpie(data = drug_info, group_key = "DrugClass", count_type = "full",
      label_info = c("group","count"), label_type = "horizon",
      label_split = "\n",
      label_size = 4, label_pos = "out" )+
  theme(legend.position = "none")
library(patchwork)
p1 + p2
ggsave("Figs/cell_drug_info.pdf",width = 14, height = 8)

####细胞系乘以药物热图
drug_info <- readxl::read_xlsx("data/drug_class.xlsx")
cell_line_meta <- readRDS("~/DRIVE/data/cell_line_meta.rds")
all_sample_meta <- readxl::read_xlsx("data/all_sample_meta.xlsx")
drug_meta <- readxl::read_xlsx("data/drug_meta.xlsx")
all_sample_meta <- left_join(all_sample_meta %>% select(-DrugType,-SMILE) %>% 
                               rename(Drug = treat),
                             drug_meta %>% distinct(Drug,.keep_all = T)) %>% 
  left_join(., drug_info) %>% 
  left_join(., cell_line_meta %>% select(cell_line, Tissue_Source2))
all_sample_meta <- all_sample_meta %>% 
  rowwise() %>% 
  mutate(DrugClass = ifelse(grepl("\\+",Drug),"Combination Drugs",DrugClass)) %>% 
  ungroup()
#saveRDS(all_sample_meta, "data/all_sample_meta.rds")
dt <- all_sample_meta %>% 
  mutate(DrugClass = ifelse(is.na(DrugClass),"Other_Therapeutic",DrugClass)) %>% 
  group_by(DrugClass,Tissue_Source2) %>% 
  summarise(counts = n()) %>% ungroup()

ggplot(dt,aes(DrugClass, Tissue_Source2)) +
  geom_point(
    aes(fill = log2(counts+1), size = log2(counts+1)), 
    color = 'black',
    shape = 21
  ) +
  theme_minimal(
    base_size = 16
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(
    y = "Tissue",
    fill = 'log2(#Cell Lines)'
  ) + guides(size = "none")+
  scale_fill_gradient2(
    low = '#FEE5C1',
    high = '#B70503',
    mid = '#FC8C59',
    limits = c(1, 7),
    midpoint = 4
  ) +
  scale_size_area(
    limits = c(1, 7),
    max_size = 7
  ) +
  coord_cartesian() + 
  theme(legend.position = 'top')
ggsave("Figs/cell_drug_heatmap.pdf",height = 8,width = 6)

###########处理组，对照组样本数量
all_sample_meta <- readxl::read_xlsx("data/all_sample_meta.xlsx")
dt <- all_sample_meta %>% 
  select(unid, treat_sample, control_sample) %>% 
  rowwise() %>% 
  mutate(treat_sample_counts = length(strsplit(treat_sample,",")[[1]]),
         control_sample_counts = length(strsplit(control_sample,",")[[1]])) %>% 
  ungroup()
treat_summ <- dt %>% group_by(treat_sample_counts) %>% 
  summarise(counts = n()) %>% ungroup() %>% mutate(log2sample = log2(counts)) %>% 
  mutate(treat_sample_counts = as.character(treat_sample_counts))
control_summ <- dt %>% group_by(control_sample_counts) %>% 
  summarise(counts = n()) %>% ungroup() %>% mutate(log2sample = log2(counts)) %>% 
  mutate(control_sample_counts = as.character(control_sample_counts))
p1 <- ggbarplot(treat_summ, x="treat_sample_counts",y="log2sample",
          xlab = "Sample Counts of Treatment group",
          ylab = "Log2(# Cell Line - Treatment)",fill = "#D9AE2C")+
  geom_text(aes(label = counts), vjust = -0.5, size = 3)
p2 <- ggbarplot(control_summ, x="control_sample_counts",y="log2sample",
                xlab = "Sample Counts of Control group",
                ylab = "Log2(# Cell Line - Treatment)",fill = "#D88C27")+
  geom_text(aes(label = counts), vjust = -0.5, size = 3)

p1/p2
ggsave("Figs/treat_control_counts.pdf",width = 5,height = 7)

######Figure2########
###差异基因数量
sample_meta <- readRDS("data/all_sample_meta.rds")
diff_res <- vector("list",nrow(sample_meta))
for (i in seq_along(diff_res)){
  tmp <- readRDS(paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/deseq2_diff.rds"))
  diff_res[[i]] <- data.frame(ids = sample_meta$unid[i], 
                              diff_counts = sum((tmp$padj < 0.05) & (abs(tmp$log2FoldChange) > 1), na.rm = T))
}
diff_res <- bind_rows(diff_res)

# dt <- diff_res %>%
#   mutate(
#     group = cut(diff_counts, 
#                 breaks = c(0, 200, 400, 600, 1000, 2000, 8500),
#                 labels = c("[0-200)", "[200-400)", "[400-600)", "[600-1000)",
#                            "[1000-2000)","[2000-8500)"),
#                 include.lowest = T, right = F)
#   ) %>%
#   count(group) %>%
#   arrange(group)
gghistogram(diff_res,x="diff_counts",add = "median", rug = TRUE,
            xlab = "Number of DEGs",ylab = "Number of Treatments",fill = "#D9AE2C")+
  annotate("text", x = median(diff_res$diff_counts)+400, y = 160, 
           label = paste0("Median: ",median(diff_res$diff_counts)), 
           size = 5, color = "red")
ggsave("Figs/num_deg.pdf",width = 6,height = 4)

#####差异通路数量
sample_meta <- readRDS("data/all_sample_meta.rds")
diff_res <- vector("list",nrow(sample_meta))
for (i in seq_along(diff_res)){
  tmp <- readRDS(paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/diff_gsea.rds"))
  diff_res[[i]] <- data.frame(ids = sample_meta$unid[i], 
                              diff_counts = sum((tmp$padj < 0.1), na.rm = T))
}
diff_res <- bind_rows(diff_res)
# dt <- diff_res %>%
#   mutate(
#     group = cut(diff_counts, 
#                 breaks = c(0, 400, 800, 1200, 1600, 2000, 3600),
#                 labels = c("[0-400)", "[400-800)", "[800-1200)", "[1200-1600)",
#                            "[1600-2000)","[2000-3600)"),
#                 include.lowest = T, right = F)
#   ) %>%
#   count(group) %>%
#   arrange(group)
gghistogram(diff_res,x="diff_counts",add = "median", rug = TRUE,
            xlab = "Number of Enriched Pathways",ylab = "Number of Treatments",fill = "#D9AE2C")+
  annotate("text", x = median(diff_res$diff_counts)+300, y = 80, 
           label = paste0("Median: ",median(diff_res$diff_counts)), 
           size = 5, color = "red")

ggsave("Figs/num_pathway.pdf",width = 6,height = 4)

####差异基因profile相似性
sample_meta <- readRDS("data/all_sample_meta.rds")
rank_res <- vector("list",nrow(sample_meta))
for (i in seq_along(rank_res)){
  tmp <- readRDS(paste0("~/GEO_data/SRA/gse_out/",sample_meta$unid[i],"/deseq2_diff.rds"))
  tmp <- tmp %>% 
    select(gene_id, stat)
    #mutate(stat = ifelse(is.na(stat),min(stat,na.rm = T),stat))
  colnames(tmp)[2] <- paste0("stat_",i)
  rank_res[[i]] <- tmp
}
rank_res <- Reduce(inner_join,rank_res)
rownames(rank_res) <- rank_res$gene_id
rank_res$gene_id <- NULL
saveRDS(rank_res, "data/deseq2_stat.rds")

###计算相关性
library(doParallel)
library(foreach)
cal_ken_par <- function(var_dt, ncores){
  all_com <- combn(colnames(var_dt), 2, simplify = F)
  all_com <- sapply(all_com, paste, collapse = "-")
  ###对列进行并行
  my.cluster <- parallel::makeCluster(
    ncores, 
    type = "PSOCK"
  )
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  res <- foreach(
    i = all_com,
    .export = c("var_dt"),
    .packages = c("dplyr")
  ) %dopar% {
    com_tmp <- strsplit(i,"-")[[1]]
    cor_tmp <- kendallknight::kendall_cor(var_dt[[com_tmp[1]]],var_dt[[com_tmp[2]]])
    return(data.frame(coms = i, cor = cor_tmp))
  }
  parallel::stopCluster(cl = my.cluster)
  res <- bind_rows(res)
  return(res)
}

cor_res <- cal_ken_par(rank_res,ncores = 60)
saveRDS(cor_res,"~/DRIVE/data/deseq2_cor.rds")

####SMILE表示相似性
library(rcdk)
library(fingerprint)
drug_meta <- readxl::read_xlsx("data/drug_meta.xlsx")
drug_meta <- drug_meta %>% filter(!is.na(SMILE)) %>% select(-DrugType) %>% distinct_all()
drug_meta <- drug_meta %>% mutate(id = paste0("Drug_",1:nrow(.)))
mols <- parse.smiles(drug_meta$SMILE)
names(mols) <- drug_meta$id
mols <- mols[lengths(mols) != 0]
mols_can <- lapply(mols, get.smiles, flavor = smiles.flavors("Canonical")) %>% unlist()
mols <- sapply(mols_can, parse.smiles)
names(mols) <- names(mols_can)
fp <- lapply(mols, get.fingerprint, type='standard')
fp.sim <- fingerprint::fp.sim.matrix(fp, method='tanimoto')
rownames(fp.sim) <- names(mols)
colnames(fp.sim) <- names(mols)

fp.sim <- as.data.frame(fp.sim)
fp.sim$id1 <- rownames(fp.sim)
fp.sim <- fp.sim %>% tidyr::pivot_longer(cols = 1:248, names_to = "id2", values_to = "sim")
fp.sim <- left_join(fp.sim, drug_meta %>% select(Drug, id) %>% rename(id1 = id))
saveRDS(fp.sim, "data/drug_smiles_sim.rds")
###基于molformer对药物进行embedding，计算embedding 向量之间的相似性
dt <- mols_can %>% as.data.frame() %>% rename(smiles = 1)
dt$id <- rownames(dt)
write.csv(dt,"scripts/M2UMol/drug_smiles.csv",quote = F,row.names = F)

molformer_sim <- read.csv("data/drug_molformer_sim.csv")
colnames(molformer_sim)[1] <- "id1"
molformer_sim <- molformer_sim %>% 
  tidyr::pivot_longer(cols = 2:249, names_to = "id2",values_to = "sim")
molformer_sim <- left_join(molformer_sim, 
                           drug_meta %>% select(Drug, id) %>% rename(id1 = id))
saveRDS(molformer_sim, "data/drug_molformer_sim.rds")

#######
sample_meta <- readRDS("data/all_sample_meta.rds")
cor_res <- readRDS("data/deseq2_cor.rds")
drug_smiles_sim <- readRDS("~/DRIVE/data/drug_smiles_sim.rds")
drug_id <- drug_smiles_sim %>% select(id1,Drug) %>% distinct_all()
mapping <- sample_meta %>% 
  select(Drug, DrugClass,Tissue_Source2) %>% 
  mutate(ids = paste0("stat_",1:nrow(sample_meta))) %>% 
  rename(Tissue = 3) %>% 
  left_join(.,drug_id %>% rename(drug_id = id1)) %>% select(-Drug)

cor_res <- cor_res %>% 
  tidyr::separate_wider_delim(cols = coms, delim = "-", names = c("id1","id2"))
cor_res <- left_join(
  cor_res, mapping %>% 
    rename(id1 = ids, Drug_class1 = DrugClass, Tissue1 = Tissue, Drug1 = drug_id)
) %>% left_join(., mapping %>% 
                  rename(id2 = ids, Drug_class2 = DrugClass, Tissue2 = Tissue, Drug2 = drug_id))
###NA 是未知药物类型
dt <- cor_res %>% 
  filter(!is.na(Drug_class1)) %>% filter(!is.na(Drug_class2)) %>% 
  mutate(is_same_drug = ifelse(Drug_class1 == Drug_class2,"Y","N")) %>% 
  mutate(is_same_tissue = ifelse(Tissue1 == Tissue2,"Y","N"))
###相同类型的药物在同一组织来源的细胞系和不同组织来源的细胞系之间的相似性
dt1 <- dt %>% filter(is_same_drug == "Y")
dt1 <- dt %>% filter(Drug1 != Drug2)
ggviolin(dt1, x="is_same_tissue",y="cor",facet.by = "Drug_class1", 
         add = "boxplot",color = "is_same_tissue",
         palette = c("#80A9C9", "#89C7C1"),
         add.params = list(outliers = FALSE),
         xlab = "Same Type of Tissue",ylab = "Similarity",legend = "none")+
  stat_compare_means(label = "p.format")+
  theme(
    strip.text = element_text(size = 6, face = "bold"),
    strip.background = element_rect(fill = "lightblue")
  )
ggsave("Figs/tissue_exp_sim.pdf",width = 7,height = 8)

######
dt <- cor_res %>% 
  filter(Drug_class1 == Drug_class2)
dt1 <- dt %>% filter(Tissue1 == Tissue2) %>% filter(Drug1 != Drug2)
dt2 <- dt %>% filter(Tissue1 != Tissue2)

dt1_summ <- dt1 %>% group_by(Drug_class1) %>% 
  summarise(median_cor = median(cor)) %>% ungroup() %>% 
  arrange(median_cor)
dt1$Drug_class1 <- factor(dt1$Drug_class1, levels = dt1_summ$Drug_class1)
library(ggrain)
library(ggsci)
ggplot(dt1, aes(Drug_class1, cor, fill = Drug_class1)) +
  geom_rain(alpha = .5, point.args = list(alpha = 0),
            boxplot.args.pos = list(
              width = 0.1, position = position_nudge(x = 0.13)),
            violin.args.pos = list(
              side = "r",
              width = 0.7, position = position_nudge(x = 0.2))) +
  theme_pubr() +
  scale_fill_npg() +
  guides(fill = 'none', color = 'none') +
  coord_flip()+
  labs(title = 'Drug Treated on Same Type of Tissue',
       y = "Similarity",x="Drug Class")
ggsave("Figs/drug_class_sim_compare.pdf",width = 10,height = 8)
# 
# p1 <- ggplot(dt1, aes(x = cor, y = Drug_class1, fill = stat(x))) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#   scale_fill_viridis_c(name = "Transcriptome response similarity", option = "C") +
#   labs(title = 'Drug Treated on Same Type of Tissue',
#        x = "Similarity",y="Drug Class")+
#   theme_ridges(center_axis_labels = T)+
#   theme(legend.position = "none")
# p2 <- ggplot(dt2, aes(x = cor, y = Drug_class1, fill = stat(x))) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#   scale_fill_viridis_c(name = "Transcriptome response similarity", option = "C") +
#   labs(title = 'Drug Treated on Different Type of Tissue',
#        x = "Similarity",y="Drug Class")+
#   theme_ridges(center_axis_labels = T)+
#   theme(legend.position = "none")
# p1/p2
# ggsave("Figs/drug_class_sim_compare.pdf",width = 10,height = 10)
#####
####同一组织来源的细胞系中，不同类型的药物之间的相似性
####
dt <- cor_res %>% filter(Tissue1 == Tissue2)
dt <- dt %>% 
  filter(!is.na(Drug_class1)) %>% filter(!is.na(Drug_class2)) %>% 
  mutate(is_same_drug = ifelse(Drug_class1 == Drug_class2,"Same Drug Class",
                               "Not Same Drug Class")) %>% 
  mutate(cor_type = case_when(
    cor >= 0.5 ~ "High \nTranscriptome Response Similarity",
    cor < 0.5 ~ "Low \nTranscriptome Response Similarity"
  ))
###fisher检验
dt_summ <- dt %>%
  group_by(cor_type,is_same_drug) %>% summarise(counts=n()) %>% ungroup()

ggplot(data=dt_summ,aes(x=cor_type,y=counts,fill=is_same_drug))+
  geom_bar(stat = "identity",position="fill")+
  theme_pubr()+
  labs(y="Percent of cases (%)",title = "Chi-squared test, P-value < 2.2e-16")+
  scale_fill_manual(values=c("#3B4992FF","#EE0000FF"))+
  theme(axis.title.x = element_blank())+
  theme(legend.title = element_blank())
ggsave("Figs/same_tissue_diff_drug_sim.pdf",width = 8,height = 8)

####药物相似性
cor_res_smile <- left_join(cor_res,
                           drug_smiles_sim %>% 
                             rename(Drug1 = id1, Drug2 = id2) %>% 
                             select(Drug1, Drug2, sim))
cor_res_smile <- cor_res_smile %>% filter(!is.na(sim))

drug_molformer_sim <- readRDS("~/DRIVE/data/drug_molformer_sim.rds")
cor_res_molformer <- left_join(cor_res,
                               drug_molformer_sim %>%
                                 rename(Drug1 = id1, Drug2 = id2) %>%
                                 select(Drug1, Drug2, sim))
cor_res_molformer <- cor_res_molformer %>% filter(!is.na(sim))

dt1 <- cor_res_smile %>% 
  filter(Tissue1 == Tissue2) %>% 
  mutate(cor_type = case_when(
    cor >= 0.5 ~ "High",
    cor < 0.5 ~ "Low"
    ))
p1 <- ggboxplot(dt1, x="cor_type",y="sim",
                xlab = "Transcriptome Response Similarity",
                ylab = "SMILEs Similarity (Tanimoto)",
                fill = "cor_type",palette = c("#D9AE2C","#D88C27"))+
  theme(legend.position = "none")+
  stat_compare_means()

dt2 <- cor_res_molformer %>% 
  filter(Tissue1 == Tissue2) %>% 
  mutate(cor_type = case_when(
    cor >= 0.5 ~ "High",
    cor < 0.5 ~ "Low"
  ))
p2 <- ggboxplot(dt2, x="cor_type",y="sim",
                xlab = "Transcriptome Response Similarity",
                ylab = "SMILEs Embedding Similarity (MoLFormer)",
                fill = "cor_type",palette = c("#D9AE2C","#D88C27"))+
  theme(legend.position = "none")+
  stat_compare_means()

library(patchwork)
p1+p2
ggsave("Figs/trans_sim_smile_sim.pdf",width = 8,height = 6)

###用药物相似性能不能区分转录组相似性，ROC
dt1 <- dt1 %>% mutate(truth = factor(cor_type, levels = c("High","Low")))
dt2 <- dt2 %>% mutate(truth = factor(cor_type, levels = c("High","Low")))

p1 <- yardstick::roc_curve(dt1, truth, sim) %>% 
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_pubr()+
  annotate("text", x = 0.25, y = 0.9, 
           label = paste0("AUROC: ",round(yardstick::roc_auc(dt1, truth, sim)$.estimate,3)), 
           size = 5, color = "red")
p2 <- yardstick::roc_curve(dt2, truth, sim) %>% 
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_pubr()+
  annotate("text", x = 0.25, y = 0.9, 
           label = paste0("AUROC: ",round(yardstick::roc_auc(dt2, truth, sim)$.estimate,3)), 
           size = 5, color = "red")
p1 + p2
ggsave("Figs/drug_sim_cor_auc.pdf",width = 10,height = 6)
#####Figure3，主要做免疫
#################hallmarks 基因集
immune_res <- readRDS("data/immune_path_gsea.rds")
# met_gsea <- readRDS("data/met_gsea.rds")
sample_meta <- readRDS("data/all_sample_meta.rds")

plot_nes <- function(pathway_gsea, sample_meta_dt, pathway_name, need_drugs){
  pathway_gsea <- left_join(sample_meta_dt, pathway_gsea %>% rename(unid = id))
  dt_sig <- pathway_gsea %>%
    group_by(Drug) %>% mutate(treat_counts = length(unique(cell_line))) %>% 
    ungroup() %>% filter(treat_counts > 3) %>% 
    filter(Drug %in% need_drugs)
  #####
  all_drugs <- unique(dt_sig$Drug)
  drug_nes <- sapply(all_drugs,
                     function(x){
                       tmp <- dt_sig %>% filter(Drug == x)
                       tmp_dt <- data.frame(Drug = x, NES = tmp$NES) %>% arrange(desc(NES))
                       pc <- sum(tmp_dt$NES > 0)
                       nc <- sum(tmp_dt$NES < 0)
                       tmp_dt$y <- c(seq(from=0.5,by=1,length.out=pc),
                                     seq(from=-0.5,by=-1,length.out=nc))
                       return(tmp_dt)
                     },simplify = F)
  drug_nes <- bind_rows(drug_nes)
  
  ggplot(drug_nes, aes(x = Drug, y = y, fill = NES)) + 
    geom_tile(width = 0.5, height = 1, color = "black", size = 0.5)+
    scale_fill_gradient(low = "blue", high = "red") +  # 颜色渐变
    labs(x = "Drug",
         y = "Cells",
         fill = "NES",
         title = pathway_name)+
    theme(
      axis.text.y = element_blank(),    # 移除Y轴刻度标签
      axis.ticks.y = element_blank(),   # 移除Y轴刻度线
      axis.line.y = element_blank()     # 移除Y轴线
    )+
    geom_hline(yintercept = 0,  # Y=0的位置
               color = "black",    # 线条颜色
               linewidth = 1,     # 线条粗细
               linetype = "solid")+
    theme(
      panel.background = element_blank(),   # 移除主背景
      plot.background = element_blank(),    # 移除绘图区背景
      panel.grid = element_blank()         # 移除所有网格线
    )+
    theme(
      axis.text.x = element_text(
        size = 14,           # 字号
        angle = 90,          # 旋转角度
        hjust = 1,           # 水平对齐（旋转时通常用1）
        #vjust = 1            # 垂直对齐
      )
    )
}

dt <- left_join(sample_meta, immune_res %>% rename(unid = id))
dt_summ <- dt %>% 
  group_by(Drug,pathway) %>% 
  summarise(pc = sum(NES > 0), nc = sum(NES < 0)) %>% ungroup() %>% 
  rowwise() %>% 
  mutate(tc = pc+nc) %>% ungroup() %>% filter(tc > 3) %>% 
  rowwise() %>% 
  mutate(is_need = case_when(
    (pathway == "DOWN_Immune_Escape") & (pc > nc) ~ "yes",
    (pathway == "UP_Immune_Escape") & (pc < nc) ~ "yes",
    (pathway == "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION") & (pc > nc) ~ "yes",
    TRUE ~ "no"
  )) %>% ungroup() %>% 
  group_by(Drug) %>% summarise(yes_c = sum(is_need == "yes")) %>% ungroup() %>% 
  filter(yes_c == 3) %>% 
  filter(!(Drug %in% c("A","T")))

dt <- immune_res %>% 
  filter(pathway == "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION")
p1 <- plot_nes(dt, sample_meta, "APM", dt_summ$Drug)
dt <- immune_res %>% 
  filter(pathway == "DOWN_Immune_Escape")
p2 <- plot_nes(dt, sample_meta, "DOWN_Immune_Escape", dt_summ$Drug)
dt <- immune_res %>% 
  filter(pathway == "UP_Immune_Escape")
p3 <- plot_nes(dt, sample_meta, "UP_Immune_Escape", dt_summ$Drug)
library(patchwork)
p1 / p2 / p3
ggsave("Figs/nes.pdf",width = 8,height = 12)

#####可变剪切


