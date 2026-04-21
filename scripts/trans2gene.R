library(optparse)
option_list <- list(
  make_option(c("-i", "--input"), type = "character", 
              help = "salmon quant.sf 路径"),
  make_option(c("-o", "--output"), type = "character", 
              help = "output filename"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "显示详细信息")
)
# 解析参数
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)
salmon_exp <- args$input
outfile <- args$output

library(dplyr)
#edb <- EnsDb.Hsapiens.v115::EnsDb.Hsapiens.v115
# dt <- ensembldb::transcripts(edb,return.type="DataFrame",
#                              columns = c("tx_id_version","gene_id_version","gene_biotype","symbol")) %>% 
#   as.data.frame()
# saveRDS(dt,"data/human_id_mapping.rds")
# salmon_exp <- data.table::fread("~/GEO_data/SRA/out/GSM4120565/salmon_res/quant.sf",data.table = F)
# edb <- EnsDb.Mmusculus.v115::EnsDb.Mmusculus.v115
# dt <- ensembldb::transcripts(edb,return.type="DataFrame",
#                              columns = c("tx_id_version","gene_id_version","gene_biotype","symbol")) %>%
#   as.data.frame()
# saveRDS(dt, "data/mouse_id_mapping.rds")

library(tximport)
dt <- readRDS("~/Drug_splicing/data/human_id_mapping.rds")
txi <- tximport(salmon_exp, type = "salmon", tx2gene = dt %>% select(1,2))
salmon_gene <- txi$counts %>% as.data.frame() %>% 
  rename(counts = 1)
salmon_gene$gene_id <- rownames(salmon_gene)
salmon_gene <- salmon_gene %>% select(gene_id, counts)
saveRDS(salmon_gene, outfile)






