library(dplyr)
library(Seurat)
sample_meta <- readxl::read_xlsx("data/CELL_metadata.xlsx",sheet = 2)
tn <- sample_meta %>% filter(`Treatement Timepoint` == "TN")
osi <- sample_meta %>% filter(Treatment == "osimertinib") %>% 
  filter(`Treatement Timepoint` != "TN")

load("data/IM01_Immune_Seurat_object_nodups-006.RData")
new_obj <- UpdateSeuratObject(tiss_immune)
sub_obj <- subset(new_obj, subset = sample_name %in% c(tn$`Sample Name`,osi$`Sample Name`))
sub_obj$sample_type <- case_when(
  sub_obj$sample_name %in% tn$`Sample Name` ~ "Naive",
  sub_obj$sample_name %in% osi$`Sample Name` ~ "Osimertinib Treatment"
)
sub_obj <- RunUMAP(sub_obj, dims = 1:50, reduction = "pca", reduction.name = "umap")
DimPlot(sub_obj, reduction = "umap", 
        group.by = c("sample_type","immune_subtype_annotation"),ncol=1)
ggsave("Figs/cell_umap.pdf",width = 7,height = 10)

obj_meta <- sub_obj@meta.data
obj_meta <- obj_meta %>% 
  select(sample_name, cell_id, immune_subtype_annotation, sample_type)
obj_summ <- obj_meta %>% 
  group_by(sample_type,sample_name) %>% 
  summarise(t_p = mean(immune_subtype_annotation == "T-cells")) %>% ungroup()
ggboxplot(obj_summ,x="sample_type",y="t_p",ylab = "Percentage of T cells",
          add = "jitter", fill = "sample_type",palette = c("#707DA6","#CCAD9D"),
          xlab = F) +
  stat_compare_means()+theme(legend.position = "none")
ggsave("Figs/Osimertinib_treat_T.pdf",width = 4,height = 5)

###所有药物
obj_meta <- new_obj@meta.data
obj_meta <- obj_meta %>% select(cell_id, sample_name,immune_subtype_annotation) %>% 
  left_join(sample_meta %>% select(`Sample Name`, `Treatement Timepoint`, Treatment) %>% 
              rename(sample_name = 1))
obj_meta <- obj_meta %>% 
  mutate(sample_type = case_when(
    `Treatement Timepoint` == "TN" ~ "Naive",
    TRUE ~ Treatment
  ))
obj_summ <- obj_meta %>% 
  group_by(sample_type, sample_name) %>% 
  summarise(t_p = mean(immune_subtype_annotation == "T-cells")) %>% ungroup()
ggboxplot(obj_summ,x="sample_type",y="t_p",ylab = "Percentage of T cells",
          add = "jitter", fill = "sample_type",
          xlab = F)

####cell rep med
exp <- readxl::read_xlsx("data/cell_reports_med_exp.xlsx") %>% as.data.frame()
meta <- readRDS("data/cell_rep_med_meta.rds")
meta <- meta %>% 
  tidyr::pivot_longer(cols = c(pre,post), names_to = "type", values_to = "sample") %>% 
  mutate(sample_id = paste0(sample,"R"))
meta <- meta %>% distinct_all()

rownames(exp) <- exp$Gene
exp$Gene <- NULL
exp_TPM <- exp %>% mutate(across(everything(), ~(./sum(.,na.rm = T))*10**6))
exp_TPM <- log2(as.matrix(exp_TPM) + 1) %>% as.data.frame()

library(GSVA)
cyt <- fgsea::gmtPathways("data/GOBP_T_CELL_MEDIATED_CYTOTOXICITY.v2026.1.Hs.gmt")
para <- gsvaParam(exp_TPM %>% as.matrix(), cyt,
                  minSize=5, maxSize=500, kcdf="Gaussian")
res <- gsva(para)
res <- t(res) %>% as.data.frame()
res$sample <- rownames(res)
colnames(res)[1] <- "score"

res <- left_join(res %>% rename(sample_id = sample), meta)
dt <- res %>% dplyr::filter(!is.na(sample))
dt_summ <- dt %>% group_by(patient) %>% summarise(counts = length(unique(type))) %>% 
  ungroup() %>% filter(counts == 2)
dt <- dt %>% filter(patient %in% dt_summ$patient)

dt <- dt %>%
  select(patient, type, score) %>% 
  tidyr::pivot_wider(names_from = type, values_from = score)
p1 <- ggpaired(dt, cond1 = "pre", cond2 = "post",
         fill = "condition", palette = "jco", 
         ylab = "GSVA score of T_CELL_MEDIATED_CYTOTOXICITY")+
  stat_compare_means(paired = T)+
  scale_fill_npg()+
  guides(fill = 'none', color = 'none')

###TIDE res
tide <- data.table::fread("~/tmp/cell_rep_med_res.txt",data.table = F)
res <- left_join(tide %>% rename(sample_id = V1), meta)
res <- res %>% dplyr::filter(!is.na(sample))
res_summ <- res %>% group_by(patient) %>% summarise(counts = length(unique(type))) %>% 
  ungroup() %>% filter(counts == 2)
res <- res %>% filter(patient %in% res_summ$patient)

dt <- res %>%
  select(patient, type, TIDE) %>% 
  rename(score = 3) %>% 
  tidyr::pivot_wider(names_from = type, values_from = score)
p2 <- ggpaired(dt, cond1 = "pre", cond2 = "post",
         fill = "condition", ylab = "TIDE score")+
  stat_compare_means(paired = T)+
  scale_fill_npg()+
  guides(fill = 'none', color = 'none')
p1 + p2 
ggsave("Figs/cell_rep_med_res.pdf",width = 8,height = 5)


dt <- res %>%
  select(patient, type, Dysfunction) %>% 
  rename(score = 3) %>% 
  tidyr::pivot_wider(names_from = type, values_from = score)
p3 <- ggpaired(dt, cond1 = "pre", cond2 = "post",
               fill = "condition", ylab = "T cell Dysfunction score")+
  stat_compare_means(paired = T)+
  scale_fill_npg()+
  guides(fill = 'none', color = 'none')

dt <- res %>%
  select(patient, type, Exclusion) %>% 
  rename(score = 3) %>% 
  tidyr::pivot_wider(names_from = type, values_from = score)
p4 <- ggpaired(dt, cond1 = "pre", cond2 = "post",
               fill = "condition", ylab = "T cell Exclusion score")+
  stat_compare_means(paired = T)+
  scale_fill_npg()+
  guides(fill = 'none', color = 'none')

library(patchwork)
p1 + p2 
ggsave("Figs/cell_rep_med_res.pdf",width = 8,height = 5)
