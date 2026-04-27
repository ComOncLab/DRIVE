library(dplyr)
##all_res
all_dir <- list.files("~/GEO_data/SRA/gse_out/",recursive = T,pattern = "A3SS.MATS.JCEC.txt",
                      full.names = T)
all_dir <- gsub("/rmats_out/.+","",all_dir)
as_type <- c("A3SS","A5SS","MXE","RI","SE")

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
  i = 1:length(all_dir),
  .export = c("all_dir","as_type"),
  .packages = c("dplyr")
) %dopar% {
  res_tmp <- vector("list")
  for (j in 1:length(as_type)){
    file_name <- paste0(all_dir[i],"/rmats_out/",as_type[j],".MATS.JC.txt")
    dt <- data.table::fread(file_name,data.table = F,check.names = T)
    if (is.character(dt$FDR)){
      dt <- data.table::fread(file_name,data.table = F,check.names = T,
                              colClasses=list(character=c("IJC_SAMPLE_1","IJC_SAMPLE_2",
                                                          "SJC_SAMPLE_1","SJC_SAMPLE_2",
                                                          "IncLevel1","IncLevel2"),
                                              numeric=c("PValue","FDR","IncLevelDifference")))
      dt$PValue <- as.numeric(dt$PValue)
      dt$FDR <- as.numeric(dt$FDR)
      dt$IncLevelDifference <- as.numeric(dt$IncLevelDifference)
    }
    dt$sample_id <- gsub(".+//","",all_dir[i])
    dt$as_type <- as_type[j]
    dt <- dt %>% 
      select(sample_id, as_type, ID, GeneID, geneSymbol, IJC_SAMPLE_1, IJC_SAMPLE_2,
             SJC_SAMPLE_1, SJC_SAMPLE_2, PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference)
    dt <- dt %>% 
      rowwise() %>% 
      mutate(reads_the = any((paste(IJC_SAMPLE_1,IJC_SAMPLE_2,
                                    SJC_SAMPLE_1,SJC_SAMPLE_2,sep = ",") %>% 
                                strsplit(.,",") %>% `[[`(1) %>% as.numeric()) > 10)) %>% 
      ungroup() %>% filter(reads_the) %>% filter(FDR < 0.05)
    dt <- dt %>% select(-reads_the)
    res_tmp <- c(res_tmp,list(dt))
  }
  res_tmp <- bind_rows(res_tmp[which(lengths(res_tmp) > 0)])
  return(res_tmp)
}
parallel::stopCluster(cl = my.cluster)
res <- bind_rows(res)
saveRDS(res,"data/rmats_res_fdr005.rds")

####save shiny results
rmats_res_fdr005 <- readRDS("~/Drug_splicing/data/rmats_res_fdr005.rds")
all_samples <- unique(rmats_res_fdr005$sample_id)
for (i in 1:length(all_samples)){
  dt <- rmats_res_fdr005 %>% filter(sample_id == all_samples[i])
  saveRDS(dt,paste0("scripts/Shiny/data/",all_samples[i],"_rmats.rds"))
}

######
res <- readRDS("data/rmats_res_fdr005.rds")
res <- res %>% filter(abs(IncLevelDifference) > 0.1)
sample_meta <- readRDS("data/all_sample_meta.rds")
res <- left_join(res, sample_meta %>% rename(sample_id = unid))
res_summ <- res %>% group_by(Tissue_Source2,cell_line,Drug) %>% 
  summarise(counts = n()) %>% ungroup() %>% 
  mutate(log2counts = log2(counts))
ord <- res_summ %>% group_by(Tissue_Source2) %>% summarise(med = median(log2counts)) %>% 
  ungroup() %>% arrange(desc(med))
p1 <- ggboxplot(res_summ,x="Tissue_Source2",y="log2counts",add = "jitter",
          color = "Tissue_Source2",order = ord$Tissue_Source2,xlab = F)+
  theme(legend.position = "none")+
  rotate_x_text(45)

res_summ <- res %>% group_by(DrugClass,cell_line,Drug) %>% 
  summarise(counts = n()) %>% ungroup() %>% 
  mutate(log2counts = log2(counts)) %>% 
  filter(!is.na(DrugClass))
ord <- res_summ %>% group_by(DrugClass) %>% summarise(med = median(log2counts)) %>% 
  ungroup() %>% arrange(desc(med))
p2 <- ggboxplot(res_summ,x="DrugClass",y="log2counts",add = "jitter",
          color = "DrugClass",order = ord$DrugClass,xlab = F)+
  theme(legend.position = "none")+
  rotate_x_text(90)
library(patchwork)
p1 / p2
ggsave("Figs/cell_drug_ascounts.pdf",width = 7,height = 8)

###比例
library(ggsci)
res_summ <- res %>% group_by(Tissue_Source2,as_type) %>% 
  summarise(counts = n()) %>% ungroup()
p1 <- ggplot(res_summ, aes(fill=as_type, y=counts, x=Tissue_Source2)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_npg()+
  theme_pubr()+
  labs(y="Proportion of splicing events",x=NULL)+
  guides(fill=guide_legend(title="Alternative splicing patterns"))+
  rotate_x_text(45)

res_summ <- res %>% group_by(DrugClass,as_type) %>% 
  summarise(counts = n()) %>% ungroup() %>% filter(!is.na(DrugClass))
p2 <- ggplot(res_summ, aes(fill=as_type, y=counts, x=DrugClass)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_npg()+
  theme_pubr()+
  labs(y="Proportion of splicing events",x=NULL)+
  guides(fill=guide_legend(title="Alternative splicing patterns"))+
  rotate_x_text(45)
p1 + p2
ggsave("Figs/cell_drug_as_type.pdf",width = 12,height = 6)

######
res <- readRDS("data/rmats_res_fdr005.rds")
sample_meta <- readRDS("data/all_sample_meta.rds")
res <- left_join(res, sample_meta %>% rename(sample_id = unid))
res <- res %>% filter(abs(IncLevelDifference) > 0.1)
###大于3个细胞系的drug
drugs <- res %>% group_by(Drug) %>% summarise(cell_counts = length(unique(cell_line))) %>% 
  ungroup()
drugs <- drugs %>% filter(cell_counts >= 3)
dt <- res %>% 
  filter(Drug %in% drugs$Drug)
dt_summ <- dt %>% group_by(sample_id, as_type) %>% 
  summarise(pos_c = mean(IncLevelDifference > 0),
            neg_c = mean(IncLevelDifference < 0)) %>% 
  ungroup() %>% left_join(., sample_meta %>% rename(sample_id = unid))
dt_stat <- dt_summ %>% 
  group_by(Drug, as_type) %>% 
  summarise(pvalue = wilcox.test(pos_c,neg_c)$p.value,
            diff = median(pos_c) - median(neg_c)) %>% ungroup()
dt_stat$padj <- p.adjust(dt_stat$pvalue,"fdr")
dt_stat_summ <- dt_stat %>% group_by(Drug) %>% 
  summarise(sig_c = sum(pvalue < 0.05),
            max_diff = max(abs(diff))) %>% ungroup() %>% 
  arrange(desc(sig_c),desc(max_diff))

####展示top5的药物
library(ggtext)
top5_drugs <- dt_stat_summ$Drug[1:5]
col_dt <- data.frame(
  as_type = unique(dt_stat$as_type),
  colors = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF") 
)
dt_sub <- dt_stat %>%
  filter(Drug %in% top5_drugs) %>%
  left_join(., col_dt) %>% 
  mutate(Drug = factor(Drug, levels = top5_drugs)) %>% 
  mutate(
    Y_label = sprintf("%s   <span style='color:%s;'>%s</span>", Drug, colors, as_type),
    Neglog10_pval = -log10(pvalue)
  ) %>% arrange(Drug) %>% 
  mutate(Y_label =  factor(Y_label, levels = rev(unique(Y_label))))
p_left <- ggplot(dt_sub, aes(x = diff, y = Y_label, fill = diff > 0)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_fill_manual(values = c("FALSE" = "#E41A1C", "TRUE" = "#377EB8"), guide = "none") +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  labs(x = "Median (Percentage of ΔPSI > 0)-Median (Percentage of ΔPSI < 0)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_markdown(size = 10),
    plot.margin = margin(5, 0, 5, 5)
  )
p_right <- ggplot(dt_sub, aes(x = Neglog10_pval, y = Y_label, fill = Neglog10_pval)) +
  geom_col(width = 0.7) +
  scale_fill_gradient(low = "#FEE08B", high = "#D73027", guide = "none") +
  labs(x = expression(-log[10](P~value)), y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),  
    plot.margin = margin(5, 5, 5, 0) 
  )
p <- p_left + p_right + plot_layout(widths = c(2, 1))
cairo_pdf("Figs/sp_top5_drugs.pdf",width = 7,height = 6)
print(p)
dev.off()

###展示所有的药物作为附图
dt_sig <- dt_stat_summ %>% filter(sig_c >= 1)
dt_part1 <- dt_sig[31:41,]
dt_sub <- dt_stat %>%
  filter(Drug %in% dt_part1$Drug) %>% 
  left_join(., col_dt) %>% 
  mutate(Drug = factor(Drug, levels = dt_part1$Drug)) %>% 
  mutate(
    Y_label = sprintf("%s   <span style='color:%s;'>%s</span>", Drug, colors, as_type),
    Neglog10_pval = -log10(pvalue)
  ) %>% arrange(Drug) %>% 
  mutate(Y_label =  factor(Y_label, levels = rev(unique(Y_label))))
p_left <- ggplot(dt_sub, aes(x = diff, y = Y_label, fill = diff > 0)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_fill_manual(values = c("FALSE" = "#E41A1C", "TRUE" = "#377EB8"), guide = "none") +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  labs(x = "Median (Percentage of ΔPSI > 0)-Median (Percentage of ΔPSI < 0)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_markdown(size = 10),
    plot.margin = margin(5, 0, 5, 5)
  )
p_right <- ggplot(dt_sub, aes(x = Neglog10_pval, y = Y_label, fill = Neglog10_pval)) +
  geom_col(width = 0.7) +
  scale_fill_gradient(low = "#FEE08B", high = "#D73027", guide = "none") +
  labs(x = expression(-log[10](P~value)), y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),  
    plot.margin = margin(5, 5, 5, 0) 
  )
p <- p_left + p_right + plot_layout(widths = c(2, 1))
cairo_pdf("Figs/sp_sig4_drugs.pdf",width = 7,height = 8)
print(p)
dev.off()


####Indisulam是已知的，作为control
dt1 <- dt %>% filter(Drug == "Indisulam") %>% 
  filter(as_type %in% c("RI","SE")) %>% 
  mutate(Diff_type = ifelse(IncLevelDifference > 0,
                            "ΔPSI > 0","ΔPSI < 0")) %>% 
  group_by(cell_line, as_type, Diff_type) %>% 
  summarise(counts = n()) %>% ungroup()
dt1$Diff_type <- factor(dt1$Diff_type,levels = c("ΔPSI > 0","ΔPSI < 0"))
p1 <- ggplot(dt1, aes(fill=Diff_type, y=counts, x=as_type)) + 
  geom_bar(position="fill", stat="identity")+
  facet_wrap(vars(cell_line), nrow  = 1)+
  scale_fill_npg()+
  theme_pubr()+
  labs(y="Proportion",x=NULL,title = "Indisulam")+
  guides(fill=guide_legend(title="ΔPSI"))
dt1 <- dt %>% filter(Drug == "Indisulam") %>% 
  filter(as_type %in% c("RI","SE")) %>% 
  group_by(sample_id,as_type) %>% 
  summarise(`ΔPSI > 0` = mean(IncLevelDifference > 0),
            `ΔPSI < 0` = mean(IncLevelDifference < 0)) %>% 
  ungroup() %>% 
  tidyr::pivot_longer(cols = c("ΔPSI > 0","ΔPSI < 0"), 
                      names_to = "type", values_to = "per")
p2 <- ggboxplot(dt1,x="type",y="per",facet.by = "as_type",add = "jitter",
          xlab = F, ylab = "Proportion", color = "type")+
  stat_compare_means()+
  scale_fill_npg()+
  theme(legend.position = "none")
library(patchwork)
# p <- p1+p2+plot_layout(widths = c(4,2))
cairo_pdf("Figs/Indisulam_as.pdf",width = 8,height = 5)
print(p1)
dev.off()

#####其他药物
dt1 <- dt %>% 
  filter(as_type %in% c("RI","SE")) %>% 
  filter(Drug == "KB-0742") %>% 
  group_by(sample_id,as_type) %>% 
  summarise(`ΔPSI > 0` = mean(IncLevelDifference > 0),
            `ΔPSI < 0` = mean(IncLevelDifference < 0)) %>% 
  ungroup() %>% 
  tidyr::pivot_longer(cols = c(`ΔPSI > 0`, `ΔPSI < 0`), names_to = "type", values_to = "per")
p <- ggboxplot(dt1,x="type",y="per",facet.by = "as_type",add = "jitter",nrow =1,
          color = "type",xlab = F,ylab = "Proportion",title = "KB-0742")+
  stat_compare_means()+
  scale_fill_npg()+
  theme(legend.position = "none")
cairo_pdf("Figs/KB_0742_as.pdf",width = 6,height = 6)
print(p)
dev.off()

#####剪切体GSEA
mapping <- readRDS("data/human_id_mapping.rds")
mapping <- mapping %>% select(gene_id_version, symbol) %>% distinct_all()
sp_path <- fgsea::gmtPathways("data/KEGG_SPLICEOSOME.v2026.1.Hs.gmt")

sample_meta <- readRDS("data/all_sample_meta.rds")
dt <- sample_meta %>% filter(Drug == "KB-0742")
res <- vector("list",nrow(dt))
for (i in 1:nrow(dt)){
  tmp <- readRDS(paste0("~/GEO_data/SRA/gse_out/",dt$unid[i],"/deseq2_diff.rds"))
  tmp <- tmp %>%
    filter(!is.na(stat)) %>% 
    left_join(.,mapping %>% rename(gene_id = gene_id_version) %>% 
                select(gene_id, symbol)) %>% 
    group_by(symbol) %>% 
    slice_max(abs(stat), with_ties = F) %>% ungroup() %>% 
    filter(nchar(symbol) > 0)
  ranks <- tmp$stat
  names(ranks) <- tmp$symbol
  ranks <- sort(ranks,decreasing = T)
  sp_path_dt <- data.frame(gs_name="KEGG_SPLICEOSOME",
                           genes=sp_path$KEGG_SPLICEOSOME)
  res_tmp <- clusterProfiler::GSEA(ranks, TERM2GENE = sp_path_dt)
  res[[i]] <- res_tmp
  names(res)[i] <- paste0(dt$cell_line[i],"-",dt$treat_time[i])
}
GseaVis::GSEAmultiGP(gsea_list = res,
                     geneSetID = "KEGG_SPLICEOSOME",
                     exp_name = gsub("GSE263153_","",dt$unid),
                     addPval = T,legend.position = "right",pvalY=0.98)
ggsave("Figs/splicesome_gsea.pdf",width = 12,height = 8)

####
res <- readRDS("data/rmats_res_fdr005.rds")
sample_meta <- readRDS("data/all_sample_meta.rds")
res <- left_join(res, sample_meta %>% rename(sample_id = unid))
res <- res %>% filter(abs(IncLevelDifference) > 0.1)
drugs <- res %>% group_by(Drug) %>% 
  summarise(cell_counts = length(unique(cell_line))) %>% ungroup()
drugs <- drugs %>% filter(cell_counts >= 3)
dt <- res %>% 
  filter(Drug %in% drugs$Drug)
dt <- dt %>% filter(Drug == "KB-0742")
dt_summ <- dt %>% group_by(sample_id, as_type) %>% summarise(counts =n()) %>% ungroup()
ggboxplot(dt_summ,x="as_type",y="counts",xlab = F,
          ylab = "Counts",add = "jitter",color = "as_type")+
  theme(legend.position = "none")+
  scale_color_npg()
ggsave("Figs/KB_0742_as_counts.pdf",width = 7,height = 7)

pos_gene <- dt %>% filter(as_type == "RI") %>% 
  filter(IncLevelDifference > 0)

pos_entrez <- clusterProfiler::bitr(pos_gene$geneSymbol,
                                    fromType = "SYMBOL",#现有的ID类型
                                    toType = "ENTREZID",#需转换的ID类型
                                    OrgDb = "org.Hs.eg.db")
KEGG_res <- clusterProfiler::enrichKEGG(gene = pos_entrez$ENTREZID,
                                        organism = "hsa", #物种Homo sapiens (智人)
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.05,
                                        pAdjustMethod = "BH",
                                        minGSSize = 10,
                                        maxGSSize = 500)
KEGG_res <- KEGG_res@result %>% filter(p.adjust < 0.05)
KEGG_res <- KEGG_res %>% filter(!grepl("disease",subcategory))

ggplot(data = KEGG_res, 
       aes(x = -log10(pvalue), y = Description, fill = Count)) +
  geom_bar(stat = "identity",width = 0.8) + 
  labs(x = "-log10(Pvalue)",
       y = "Pathway")+
  scale_fill_bs5("blue")+
  theme_pubr()
ggsave("Figs/KB_0742_res_kegg.pdf",width = 8,height = 8)




