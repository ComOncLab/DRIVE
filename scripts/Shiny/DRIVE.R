library(dplyr)
library(shiny)
library(shinydashboard)
library(DT)     # 用于展示数据表格
library(ggplot2) # 用于画图
library(bslib)
library(shinyjs)
library(shinyalert)
library(ggsci)
library(ggpubr)
library(markdown)
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
cal_WTCS <- function(gene_sig, up_gene_name, down_gene_name){
  gene_sig <- gene_sig[order(gene_sig$FC,decreasing = TRUE,na.last = NA), ]
  sig_ranks <- gene_sig$FC
  names(sig_ranks) <- gene_sig$gene
  query_pathway <- list()
  query_pathway[["query_up_gene"]] <- up_gene_name
  query_pathway[["query_down_gene"]] <- down_gene_name
  fgseaRes <- fgsea::fgsea(query_pathway, sig_ranks)
  ES_up <- fgseaRes$ES[which(fgseaRes$pathway == "query_up_gene")]
  ES_up_padj <- fgseaRes$padj[which(fgseaRes$pathway == "query_up_gene")]
  ES_down <- fgseaRes$ES[which(fgseaRes$pathway == "query_down_gene")]
  ES_down_padj <- fgseaRes$padj[which(fgseaRes$pathway == "query_down_gene")]
  if (sum(sign(fgseaRes$ES)) == 0) {
    WTCS <- (ES_up - ES_down) / 2
  } else{
    WTCS <- 0
  }
  return(
    data.frame(
      WTCS = WTCS,
      ES_up = ES_up,
      ES_down = ES_down,
      ES_up_padj = ES_up_padj,
      ES_down_padj = ES_down_padj
    )
  )
}
cal_cor <- function(gene_sig, query_sig){
  ##gene_sig and query_sig is dataframe, first column is gene name, second is FC (log2FC)
  gene_sig <- gene_sig[order(gene_sig$FC,decreasing = TRUE,na.last = NA), ]
  query_sig <- query_sig[order(query_sig$FC,decreasing = TRUE,na.last = NA), ]
  query_sig <- query_sig %>% rename(query_FC = FC)
  merge_ref_query <- merge(gene_sig,
                           query_sig,
                           by = "gene",
                           all.X = TRUE)
  Spearman_whole <- cor(merge_ref_query$FC,
                        merge_ref_query$query_FC,
                        method = "spearman")
  Pearson_whole <- cor(merge_ref_query$FC,
                       merge_ref_query$query_FC,
                       method = "pearson")
  Cosine_whole <- as.vector(lsa::cosine(merge_ref_query$FC, merge_ref_query$query_FC))
  return(
    data.frame(
      Spearman_correlation = Spearman_whole,
      Pearson_correlation = Pearson_whole,
      Cosine_similarity = Cosine_whole
    )
  )
}
cal_sim_scores <- function(query_sig, K=100, ncores=20){
  ###query sig 用户输入的 signature 数据框，两列，一列是基因symbol，另一列是FC（log2FC）
  deseq_fc <- readRDS("~/Drug_splicing/data/deseq_fc.rds")
  xsum_per <- readRDS("~/Drug_splicing/data/Xsum_permutation.rds")
  css_per <- readRDS("~/Drug_splicing/data/Zhang_permutation.rds")
  query_sig <- query_sig[order(query_sig$FC,decreasing = TRUE,na.last = NA), ]
  query_up_gene <- head(query_sig$gene, K)
  query_down_gene <- tail(query_sig$gene, K)
  all_samples <- colnames(deseq_fc)
  library(doParallel)
  library(foreach)
  #create the cluster
  my.cluster <- parallel::makeCluster(
    ncores, 
    type = "PSOCK"
  )
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  res <- foreach(
    i = 1:length(all_samples),
    .export = c("cal_Xsum","cal_zhang","cal_WTCS","cal_cor","deseq_fc",
                "query_sig","query_up_gene","query_down_gene","xsum_per","css_per"),
    .packages = c("dplyr")
  ) %dopar% {
    ref_sig <- deseq_fc %>% select(all_samples[i]) %>% rename(FC = 1)
    ref_sig$gene <- rownames(ref_sig)
    wtcs_res <- cal_WTCS(ref_sig, up_gene_name = query_up_gene, down_gene_name = query_down_gene)
    xsum_res <- cal_Xsum(ref_sig,up_gene_name = query_up_gene, down_gene_name = query_down_gene)
    xsum_per_tmp <- xsum_per %>% select(all_samples[i])
    xsum_res_p <- sum(abs(xsum_per_tmp[,1]) >= abs(xsum_res)) / 10000
    css_res <- cal_zhang(ref_sig, up_gene_name = query_up_gene, down_gene_name = query_down_gene)
    css_per_tmp <- css_per %>% select(all_samples[i])
    css_p <- sum(abs(css_per_tmp[,1]) >= abs(css_res)) / 10000
    cor_res <- cal_cor(ref_sig, query_sig = query_sig)
    data.frame(
      signature_index = all_samples[i],
      WTCS = wtcs_res$WTCS,
      ES_up = wtcs_res$ES_up,
      ES_down = wtcs_res$ES_down,
      ES_up_padj = wtcs_res$ES_up_padj,
      ES_down_padj = wtcs_res$ES_down_padj,
      XSum = xsum_res,
      XSum_pvalue = xsum_res_p,
      CSS = css_res,
      CSS_pvalue = css_p,
      Spearman_correlation = cor_res$Spearman_correlation,
      Pearson_correlation = cor_res$Pearson_correlation,
      Cosine_similarity = cor_res$Cosine_similarity
    )
  }
  parallel::stopCluster(cl = my.cluster)
  all_res <- bind_rows(res)
  all_res <- all_res %>% rowwise() %>% 
    mutate(pos_counts = sum(c(WTCS, XSum, CSS, Spearman_correlation,
                              Pearson_correlation, Cosine_similarity) > 0)) %>% 
    mutate(type = ifelse(pos_counts > 3, "Pos","Neg")) %>% ungroup()
  pos <- all_res %>% filter(type == "Pos")
  neg <- all_res %>% filter(type == "Neg")
  ##前1% 赋值1
  Meta_score_pos <- as.numeric(rank(pos$WTCS,na.last = FALSE) > (nrow(pos) - round(nrow(pos) * 0.01))) +
    as.numeric(rank(pos$XSum,na.last = FALSE) > (nrow(pos) - round(nrow(pos) * 0.01))) +
    as.numeric(rank(pos$CSS,na.last = FALSE) > (nrow(pos) - round(nrow(pos) * 0.01))) +
    as.numeric(rank(pos$Spearman_correlation,na.last = FALSE) > (nrow(pos) - round(nrow(pos) *0.01))) +
    as.numeric(rank(pos$Pearson_correlation,na.last = FALSE) > (nrow(pos) - round(nrow(pos) * 0.01))) +
    as.numeric(rank(pos$Cosine_similarity,na.last = FALSE) > (nrow(pos) - round(nrow(pos) * 0.01)))
  Meta_score_neg <-  as.numeric(rank(neg$WTCS) <= round(nrow(neg) * 0.01)) +
    as.numeric(rank(neg$XSum) <= round(nrow(neg) * 0.01)) +
    as.numeric(rank(neg$CSS) <= round(nrow(neg) * 0.01)) +
    as.numeric(rank(neg$Spearman_correlation) <= round(nrow(neg) * 0.01)) +
    as.numeric(rank(neg$Pearson_correlation) <= round(nrow(neg) * 0.01)) +
    as.numeric(rank(neg$Cosine_similarity) <= round(nrow(neg) * 0.01))
  all_res <- bind_rows(pos, neg) %>% 
    bind_cols(., data.frame(Meta_score = c(Meta_score_pos, Meta_score_neg)))
  all_res <- all_res %>% select(signature_index, Meta_score, everything()) %>% arrange(desc(Meta_score))
  return(all_res)
}
radar_plot <- function(data) {
  assertthat::assert_that(nrow(data)==1)
  dat = data.frame(
    WTCS = data$WTCS,
    CSS = data$CSS,
    Spearman_correlation = data$Spearman_correlation,
    Pearson_correlation = data$Pearson_correlation,
    Cosine_similarity = data$Cosine_similarity
  )
  
  rownames(dat) <- data$signature_index
  
  # 添加最大和最小行（需根据实际范围调整）
  max_min <- data.frame(
    WTCS = c(1,-1),
    CSS = c(1,-1),
    Spearman_correlation = c(1,-1),
    Pearson_correlation = c(1,-1),
    Cosine_similarity = c(1,-1)
  )
  plot_data <- rbind(max_min, dat)
  
  # 绘制雷达图
  p = fmsb::radarchart(
    plot_data,
    title = paste0("Signature Name: ", data$signature_index),
    pcol = "#2D336B", # 数据线颜色
    plwd = 2, # 数据线宽度
    pfcol = rgb(0, 0, 1, 0.3), # 填充颜色
    cglcol = "grey", # 网格线颜色
    cglty = 1, # 网格线类型
    axislabcol = "grey", # 轴标签颜色
    vlcex = 1.2          # 变量标签字号
  )
  return(p)
}

# exp_query <- readRDS("~/Drug_splicing/scripts/Shiny/data/example_query_sig.rds")
# tt <- cal_sim_scores(exp_query, ncores = 20)
# radar_plot(tt[7,])

sample_meta <- readRDS("~/Drug_splicing/data/all_sample_meta.rds")
sample_hla <- readRDS("~/Drug_splicing/data/samples_with_HLA.rds")
all_rmats_samples <- readRDS("~/Drug_splicing/data/all_rmats_samples.rds")
all_deseq_samples <- readRDS("~/Drug_splicing/data/all_deseq_samples.rds")
escape_up_genes <- c("TGFB1", "IL10", "VEGFA", "CD274", "PVR","PDCD1LG2", 
                     "NECTIN2", "NECTIN3", "CD80", "CD86", 
                     "CD47", "TNFSF4", "TNFSF9", "CCL2","TNFRSF14",
                     "CSF1", "CEACAM1", "LGALS9",
                     ##LIN_TUMOR_ESCAPE_FROM_IMMUNE_ATTACK
                     "ACKR3","C14orf39", "CHD7", "DYNAP", "GJB4", "GPR149", 
                     "HTRA1", "IL2RG", "JAG1", "KHDRBS3", "MMD", "NPPB", 
                     "NRP2", "PLA2G7", "RGS16", "ST3GAL6", "SYCP3", 
                     "TNNT2", "VCAM1") %>% paste(.,collapse = ",") 
all_files <- list.files("~/Drug_splicing/scripts/Shiny/data/",pattern = "_pep.rds")
neo_files <- sample_meta %>% 
  filter(unid %in% gsub("_binding_pep.rds","",all_files))

get_drug_enrich <- function(sample_list, pathway_list, deseq_dir, ncores){
  library(doParallel)
  library(foreach)
  #create the cluster
  my.cluster <- parallel::makeCluster(
    ncores, 
    type = "PSOCK"
  )
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  res <- foreach(
    i = 1:length(sample_list),
    .export = c("sample_list","pathway_list"),
    .packages = c("dplyr","fgsea")
  ) %dopar% {
    tmp <- readRDS(paste0(deseq_dir,sample_list[i]))
    tmp <- tmp %>%
      filter(!is.na(stat)) %>% 
      group_by(symbol) %>% 
      slice_max(stat, with_ties = F) %>% ungroup() %>% 
      filter(nchar(symbol) > 0)
    ranks <- tmp$stat
    names(ranks) <- tmp$symbol
    gsea_res <- fgsea::fgsea(pathways = pathway_list, stats = ranks, nperm = 1000)
    gsea_res <- gsea_res %>% as.data.frame() %>% mutate(id = gsub("_deseq.rds","",sample_list[i]))
    return(gsea_res)
  }
  parallel::stopCluster(cl = my.cluster)
  res <- bind_rows(res) 
  res$padj <- p.adjust(res$pval,"fdr")
  return(res)
}
plot_nes <- function(pathway_gsea, sample_meta_dt, pathway_name, need_drugs){
  pathway_gsea <- left_join(sample_meta_dt, pathway_gsea %>% rename(unid = id))
  dt_sig <- pathway_gsea %>%
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
  drug_nes$Drug <- factor(drug_nes$Drug, levels = need_drugs)
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
# ============================
# 1. UI 部分 (前端界面)
# ============================
ui <- navbarPage(
  title = "DRIVE", # 顶部系统名称
  # 主题颜色（可选，让界面更好看）
  theme = bslib::bs_theme(version = 5, bootswatch = "journal"),
  # ---------------------------------
  # 主页面 1：数据概览
  # ---------------------------------
  tabPanel(
    title = "Transcriptome response",
    sidebarLayout(
      # 左侧：数据选择区
      sidebarPanel(
        width = 3, # 占比3/12,
        selectInput("sele_drug", "Select Drug:", 
                    choices = unique(sample_meta$Drug)),
        uiOutput("sele_cell_ui"),
        uiOutput("sele_sample_ui"),
        hr(),
        helpText("Please select drug and cell line combination")
      ),
      # 右侧：数据展示区
      mainPanel(
        width = 9, # 占比9/12
        tabsetPanel(
          tabPanel("Metadata",
                   hr(),
                   uiOutput("metainfo")
                   ),
          tabPanel("Gene Expression",
                   shinycssloaders::withSpinner(DTOutput("deseq_tab"),type = 3,color.background = "white"),
                   shinycssloaders::withSpinner(plotOutput("Volcano",height="1000px"),type = 3,color.background = "white")
                   ),
          tabPanel("Alternative Splicing",
                   uiOutput("as_ui")
                   ) 
        )
      )
    )
  ),
  tabPanel(
    title = "Drug Enrichment",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        textAreaInput("gene_list","Input Gene List", value = escape_up_genes,height="200px"),
        hr(),
        fileInput("gene_list_file", "Or choose Gene List File (.gmt or .txt)", accept = c(".txt",".gmt")),
        p("The uploaded txt file should contain one gene per line and no colnames. If If no file is uploaded, we will use the gene list in the text box."),
        br(),
        actionButton("drug_enrich_btn", "Enrich", class = "btn-primary")
      ),
      mainPanel(
        width = 9,
        shinycssloaders::withSpinner(DTOutput("drug_enrich_dt"),type = 3,color.background = "white"),
        hr(),
        shinycssloaders::withSpinner(plotOutput("drug_enrich_plot",height = "600px"),type = 3,
                                     color.background = "white")
      )
    )
  ),
  # ---------------------------------
  # ---------------------------------
  tabPanel(
    title = "Connectivity analysis",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput("gene_signature_file", "Upload a gene signature file (.txt)", accept = c(".txt")),
        p("The file should have two columns (space as separator), which the frist column is",
          strong(" gene (gene ID)")," and second column is ",strong(" FC (Log2FoldChange)")),
        actionButton("conn_analysis", "Run", class = "btn-primary"),
        downloadButton("down_example_sig", "Download Example File", class = "btn-info")
      ),
      mainPanel(
        width = 9,
        shinycssloaders::withSpinner(DTOutput("con_tab"),type = 3,color.background = "white"),
        layout_columns(
          card(card_header("Top 1 Positive Corrlated Signatures"),
               shinycssloaders::withSpinner(plotOutput("con_figs_pos",height = "500px"),
                                            type = 3,color.background = "white")),
          card(card_header("Top 1 Negative Corrlated Signatures"),
               shinycssloaders::withSpinner(plotOutput("con_figs_neg",height = "500px"),
                                            type = 3,color.background = "white")),
        )
      )
    )
  ),
  # ---------------------------------
  tabPanel(
    title = "AS-derived Neoantigen",
    sidebarLayout(
      sidebarPanel(
        width = 3, # 占比3/12,
        selectInput("sele_drug_neo", "Select Drug:", 
                    choices = unique(neo_files$Drug)),
        uiOutput("sele_cell_ui_neo"),
        uiOutput("sele_sample_ui_neo"),
        hr(),
        helpText("Please select drug and cell line combination")
      ),
      mainPanel(
        width = 9,
        card(card_header("Predicted %Rank for Peptides and HLA"),
             shinycssloaders::withSpinner(DTOutput("pep_rank"),
                                          type = 3,color.background = "white")),
        card(card_header("Information of Binding Peptides"),
             shinycssloaders::withSpinner(DTOutput("pep_meta"),
                                          type = 3,color.background = "white")),
        fluidRow(
          column(
            width = 8,
            shinycssloaders::withSpinner(plotOutput("neo_plot",height = "500px",width = "100%"),
                                         type = 3,color.background = "white")
          ),
          column(
            width = 4,
            wellPanel(
              style = "font-family: 'Microsoft YaHei', sans-serif; height: 400px; overflow-y: auto;",
              h5("This plot displays the differential alternative splicing events from the table above."),
              tags$ul(
                tags$li("The x-axis represents the differential expression log2FoldChange of the genes where the differential alternative splicing events are located"),
                tags$li("The y-axis shows the delta PSI of these events."),
                tags$li("The color indicates the reciprocal of the %Rank of the peptide with the minimum %Rank value among the peptides produced by the differential alternative splicing event."),
                tags$li("We have labeled the top 5 ranked events. The ID is Gene-rMATS Type-rMATS ID")
              )
            )
          )
        )
      )
    )
  ),
  tabPanel(
    title = "About",
    tags$iframe(src = "DRIVE_about.html", 
                width = "100%", 
                style = "height: 80vh;", # vh units ensure vertical fill
                frameBorder = "0")
  )
)
# ============================
# 2. Server 部分 (后端逻辑)
# ============================
server <- function(input, output, session) {
  shinyalert(
    title = "Welcome to DRIVE!",
    text = "Database of <b>D</b>rug induced <b>R</b>NA profile, <b>I</b>soform <b>V</b>ariability and <b>E</b>pitopes",
    imageUrl = "wel.png", # Image from www folder
    imageWidth = 700,
    imageHeight = 525,
    animation = T,size = "l",html = T
  )
  output$down_example_sig <- downloadHandler(
    filename = function() {
      "query_example.txt"
    },
    content = function(file) {
      file.copy("~/Drug_splicing/scripts/Shiny/data/example_query.txt", file)
    }
  )
  # --- 页面1：meta信息+转录组+可变剪切 ---
  output$sele_cell_ui <- renderUI({
    req(input$sele_drug)
    sele_drug <- input$sele_drug
    dt <- sample_meta %>% dplyr::filter(Drug == sele_drug)
    selectInput("sele_cell","Select Cell Line:",
                choices = unique(dt$cell_line))
  })
  output$sele_sample_ui <- renderUI({
    req(input$sele_cell)
    dt <- sample_meta %>% 
      dplyr::filter(Drug == input$sele_drug) %>% 
      dplyr::filter(cell_line == input$sele_cell)
    selectInput("sele_sample","Select Drug-Cell combination:",
                choices = unique(dt$unid))
  })
  
  output$metainfo <- renderUI({
    req(input$sele_sample)
    dt <- sample_meta %>% 
      dplyr::filter(unid == input$sele_sample) %>% 
      dplyr::mutate_all(as.character)
    dt[,which(is.na(dt[1,]))] <- "Not Available"
    dt[,which(grepl("Notspecified",dt[1,]))] <- "Not Available"
    dt_hla <- sample_hla %>% 
      dplyr::filter(unid == input$sele_sample)
    dt_hla <- dt_hla$HLA %>% strsplit(.,",") %>% `[[`(1)
    hlaa <- dt_hla[grep("HLA-A",dt_hla)]
    hlab <- dt_hla[grep("HLA-B",dt_hla)]
    hlac <- dt_hla[grep("HLA-C",dt_hla)]
    html_con <- readLines("~/Drug_splicing/scripts/Shiny/test.html")
    html_con <- paste(html_con, collapse = "\n")
    html_con <- gsub("GSEID",dt$study_alias,html_con,fixed = T)
    html_con <- gsub("CELLID",dt$cell_line,html_con,fixed = T)
    html_con <- gsub("TISSUEID",dt$Tissue_Source2,html_con,fixed = T)
    html_con <- gsub("DRUGID",dt$Drug,html_con,fixed = T)
    html_con <- gsub("DRUGCLASSID",dt$DrugClass,html_con,fixed = T)
    html_con <- gsub("DRUGSUBCLASSID",dt$DrugType,html_con,fixed = T)
    html_con <- gsub("TRECID",dt$treat_con,html_con,fixed = T)
    html_con <- gsub("TRETID",dt$treat_time,html_con,fixed = T)
    html_con <- gsub("SMILEID",dt$SMILE,html_con,fixed = T)
    
    hlaa <- sapply(hlaa,
                   function(x){gsub("gse",x,
                                    '<span class="sample-badge treatment-badge">gse</span>',
                                    fixed = T)}) %>% paste(.,collapse = "\n")
    hlab <- sapply(hlab,
                   function(x){gsub("gse",x,
                                    '<span class="sample-badge treatment-badge">gse</span>',
                                    fixed = T)}) %>% paste(.,collapse = "\n")
    hlac <- sapply(hlac,
                   function(x){gsub("gse",x,
                                    '<span class="sample-badge treatment-badge">gse</span>',
                                    fixed = T)}) %>% paste(.,collapse = "\n")
    
    html_con <- gsub("HLAA",hlaa,html_con,fixed = T)
    html_con <- gsub("HLAB",hlab,html_con,fixed = T)
    html_con <- gsub("HLAC",hlac,html_con,fixed = T)
    treatg <- sapply(strsplit(dt$treat_sample,",")[[1]],
                     function(x){gsub("gse",x,
                                      '<span class="sample-badge treatment-badge">gse</span>',
                                      fixed = T)}) %>% paste(.,collapse = "\n")
    html_con <- gsub("TREATGROUP",treatg,html_con,fixed = T)
    contrlg <- sapply(strsplit(dt$control_sample,",")[[1]],
                     function(x){gsub("gse",x,
                                      '<span class="sample-badge treatment-badge">gse</span>',
                                      fixed = T)}) %>% paste(.,collapse = "\n")
    html_con <- gsub("CONTROLGROUP",contrlg,html_con,fixed = T)
    writeLines(html_con, "tmp.html")
    includeHTML("tmp.html")
  })
  deseq_dt <- reactive({
    req(input$sele_sample)
    dt <- readRDS(paste0("~/Drug_splicing/scripts/Shiny/data/",input$sele_sample,"_deseq.rds"))
    dt
  })
  output$deseq_tab <- renderDT({
    datatable(deseq_dt()) %>% 
      formatRound(
        columns = c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj"),  # 指定要格式化的列
        digits = 3  # 保留2位小数
      )
  })
  output$Volcano <- renderPlot({
    EnhancedVolcano::EnhancedVolcano(
      deseq_dt(),
      x = "log2FoldChange",
      y = "pvalue",
      lab = deseq_dt()$symbol,
      #pCutoff = 0.05,
      FCcutoff = 1,
      title = input$sele_sample,
      subtitle = ""
    )
  })
  rmats_dt <- reactive({
    req(input$sele_sample)
    if (input$sele_sample %in% gsub("_rmats.rds","",all_rmats_samples)){
      dt <- readRDS(paste0("~/Drug_splicing/scripts/Shiny/data/",input$sele_sample,"_rmats.rds"))
      dt <- dt %>% mutate(geneSymbol = ifelse(is.na(geneSymbol),GeneID,geneSymbol))
      dt
    }else{
      NULL
    }
  })
  output$as_ui <- renderUI({
    if (is.null(rmats_dt())){
      textOutput("no_rmats")
    }else{
      genes <- rmats_dt()$geneSymbol %>% unique()
      tagList(
        fluidRow(DTOutput("rmats_tab")),
        hr(),
        fluidRow(
          column(width = 6, shinycssloaders::withSpinner(plotOutput("as_type",height="600px"),,type = 3,color.background = "white")),
          column(width = 6, plotOutput("as_gene",height="600px"))
        ),
        hr(),
        fluidRow(
          column(width = 4, selectInput("sele_as_gene","Which Gene to show ?",choices = genes)),
          column(width = 8, plotOutput("as_gene_select"))
        )
      )
    }
  })
  output$no_rmats <- renderText({
    "This Drug-Cell combination has no significant splicing events (FDR < 0.05)."
  })
  output$rmats_tab <- renderDT({
    datatable(rmats_dt() %>% select(-sample_id) %>% dplyr::rename(AS_type = as_type),
              options = list(
                autoWidth = TRUE,
                scrollX = TRUE
              ))
  })
  output$as_type <- renderPlot({
    dt_summ <- rmats_dt() %>% group_by(as_type) %>% 
      summarise(`ΔPSI > 0` = sum(IncLevelDifference > 0), `ΔPSI < 0` = sum(IncLevelDifference < 0)) %>% 
      ungroup() %>% 
      tidyr::pivot_longer(cols = c(`ΔPSI > 0`, `ΔPSI < 0`), names_to = "dpsi_type", values_to = "counts")
    ggbarplot(dt_summ, x="as_type",y="counts",fill = "dpsi_type",
              label = TRUE,lab.pos = "in",xlab = "Alternative Splicing Type",ylab = "Counts",
              title = "The number of alternative splicing types with significant changes.")+
      scale_fill_npg()+
      labs(fill = "ΔPSI = PSI(Treatment Group) - PSI(Control Group)")
  })
  output$as_gene <- renderPlot({
    dt_summ <- rmats_dt() %>% 
      group_by(geneSymbol,as_type) %>% summarise(counts = n()) %>% ungroup()
    dt_summ2 <- dt_summ %>% group_by(geneSymbol) %>% 
      summarise(t_c = sum(counts)) %>% ungroup() %>% arrange(desc(t_c)) %>% 
      slice_head(n=10)
    dt_summ <- dt_summ %>% filter(geneSymbol %in% dt_summ2$geneSymbol)
    ggbarplot(dt_summ, x="geneSymbol",y="counts",fill = "as_type",
              order = dt_summ2$geneSymbol,xlab = "Genes",ylab = "Counts",
              title = "Top 10 genes with the most number of significantly changed alternative splicing events.")+
      scale_fill_npg()+
      labs(fill = "AS Type")
  })
  output$as_gene_select <- renderPlot({
    req(input$sele_as_gene)
    gene_dt <- rmats_dt() %>% 
      filter(geneSymbol == input$sele_as_gene) %>% 
      select(as_type, IncLevel1, IncLevel2) %>% 
      tidyr::separate_wider_delim(cols = c(IncLevel1, IncLevel2),delim = ",",names_sep = "-") %>% 
      tidyr::pivot_longer(cols = 2:ncol(.), names_to = "type", values_to = "PSI") %>% 
      mutate(groups = gsub("-.+","",type)) %>% 
      mutate(groups = ifelse(groups == "IncLevel1", "PSI-Treat","PSI-Ctrl"))
    gene_dt$PSI <- as.numeric(gene_dt$PSI)
    ggboxplot(gene_dt,x="groups",y="PSI",facet.by = "as_type", fill = "groups")+
      stat_compare_means(label.y=1.01)+
      scale_fill_npg()+
      theme(legend.position = "none")
  })
  ##check upload file
  observeEvent(input$drug_enrich_btn, {
    file <- input$gene_list_file
    ext <- tools::file_ext(file$datapath)
    if (is.null(file)){
      gene_list <- input$gene_list
      if (nchar(gene_list) < 1){
        shinyalert("Oops!", "Must use text box or upload.", type = "error")
      }
      if (!grepl(",",gene_list)){
        shinyalert("Oops!", "Please use comma as separator in text box.", type = "error")
      }
    }
    req(file)
    if (ext == "txt"){
      tryCatch({
        dt <- read.table(file$datapath,sep = " ",header=F)
      }, error = function(e) {
        shinyalert("Oops!", "The uploaded txt file has some error.", type = "error")
      })
      if (ncol(dt) > 1){
        shinyalert("Oops!", "The uploaded txt file has more than one columns.", type = "error")
      }
    }else if (ext == "gmt"){
      tryCatch({
        dt <- fgsea::gmtPathways(file$datapath)
      }, error = function(e) {
        shinyalert("Oops!", "The uploaded gmt file has some error.", type = "error")
      })
    }else {
      shinyalert("Oops!", "Please upload a gmt or txt file.", type = "error")
    }
  })
  ###get gene list 
  gene_list_data <- eventReactive(input$drug_enrich_btn, {
    file <- input$gene_list_file
    if (is.null(file)){
      gene_list <- strsplit(input$gene_list,",")[[1]]
      gene_list <- list(USE_GENE = gene_list)
      return(gene_list)
    }else {
      ext <- tools::file_ext(file$datapath)
      if (ext == "txt"){
        dt <- read.table(file$datapath,sep = " ",header=F)
        gene_list <- list(USE_GENE = dt[,1])
        return(gene_list)
      }else{
        dt <- fgsea::gmtPathways(file$datapath)
        return(dt)
      }
    }
  })
  ####drug enrichment 
  drug_enrich_res <- reactive({
    req(gene_list_data())
    res <- get_drug_enrich(all_deseq_samples, gene_list_data(), "~/Drug_splicing/scripts/Shiny/data/", ncores = 30)
    return(res)
  })
  output$drug_enrich_dt <- renderDT({
    datatable(drug_enrich_res() %>% select(-1),
              options = list(
                autoWidth = TRUE,
                scrollX = TRUE
              )) %>% 
      formatRound(
        columns = c("pval", "padj","ES","NES"),  # 指定要格式化的列
        digits = 5
      )
  })
  ####画前15个正负比例差异最大的
  output$drug_enrich_plot <- renderPlot({
    dt <- left_join(sample_meta, drug_enrich_res() %>% rename(unid = id))
    dt_summ <- dt %>% 
      group_by(Drug, pathway) %>% 
      summarise(pc = sum(NES > 0), nc = sum(NES < 0)) %>% ungroup() %>% 
      rowwise() %>% 
      mutate(tc = pc+nc) %>% ungroup() %>% filter(tc > 3) %>% 
      rowwise() %>% mutate(diff = abs(pc - nc) / tc) %>% 
      ungroup() %>% 
      mutate(diff_type = ifelse((pc - nc) > 0,"pos","neg")) %>% 
      group_by(diff_type) %>% 
      slice_max(order_by = diff, n = 10, with_ties = F) %>% ungroup()
    plot_nes(drug_enrich_res(), sample_meta, 
             "The top 20 drugs with the greatest differences in NES \n (The first 10 drugs tend to have more negative NES, while the next 10 drugs tend to have more positive NES.)", 
             dt_summ$Drug)
  })
  #####check upload signature
  observeEvent(input$conn_analysis,{
    file <- input$gene_signature_file
    req(file)
    tryCatch({
      dt <- read.table(file$datapath,sep = " ",header=F)
    }, error = function(e) {
      shinyalert("Oops!", "The uploaded txt file has some error.", type = "error")
    })
    if (all(c("gene","FC") %in% colnames(dt))){
      shinyalert("Oops!", "column name gene and FC not exist in the uploaded file", type = "error")
    }
  })
  user_query_sig <- eventReactive(input$conn_analysis, {
    file <- input$gene_signature_file
    dt <- read.table(file$datapath,sep = " ",header=T)
    dt <- dt %>% filter(!is.na(FC))
    return(dt)
  })
  sim_res <- eventReactive(input$conn_analysis,{
    req(user_query_sig())
    cal_sim_scores(query_sig = user_query_sig(), K=100, ncores = 10)
  })
  output$con_tab <- renderDT({
    datatable(sim_res(),
              options = list(
                autoWidth = TRUE,
                scrollX = TRUE
              )) %>% 
      formatRound(
        columns = c("WTCS", "ES_up","ES_down","ES_up_padj", "ES_down_padj",
                    "XSum", "XSum_pvalue", "CSS", "CSS_pvalue",
                    "Spearman_correlation","Pearson_correlation","Cosine_similarity"),
        digits = 3
      )
  })
  output$con_figs_pos <- renderPlot({
    pos_res <- sim_res() %>% filter(type == "Pos")
    if (nrow(pos_res) > 1){
      radar_plot(pos_res[1,]) 
    }
  })
  output$con_figs_neg <- renderPlot({
    neg_res <- sim_res() %>% filter(type == "Neg")
    if (nrow(neg_res) > 1){
      radar_plot(neg_res[1,])
    }
  })
  ####
  output$sele_cell_ui_neo <- renderUI({
    req(input$sele_drug_neo)
    sele_drug <- input$sele_drug_neo
    dt <- neo_files %>% dplyr::filter(Drug == sele_drug)
    selectInput("sele_cell_neo","Select Cell Line:",
                choices = unique(dt$cell_line))
  })
  output$sele_sample_ui_neo <- renderUI({
    req(input$sele_cell_neo)
    dt <- neo_files %>% 
      dplyr::filter(Drug == input$sele_drug_neo) %>% 
      dplyr::filter(cell_line == input$sele_cell_neo)
    selectInput("sele_sample_neo","Select Drug-Cell combination:",
                choices = unique(dt$unid))
  })
  ##show tables for peptides
  pep_dt <- reactive({
    req(input$sele_sample_neo)
    dt <- readRDS(paste0("~/Drug_splicing/scripts/Shiny/data/",
                         input$sele_sample_neo,"_binding_pep.rds"))
    dt_meta <- readRDS(paste0("~/Drug_splicing/scripts/Shiny/data/",
                              input$sele_sample_neo,"_binding_pep_meta.rds"))
    dt_deg <- readRDS(paste0("~/Drug_splicing/scripts/Shiny/data/",
                             input$sele_sample_neo,"_deseq.rds"))
    return(list(dt,dt_meta,dt_deg))
  })
  
  output$pep_rank <- renderDT({
    datatable(pep_dt()[[1]],
              options = list(
                autoWidth = TRUE,
                scrollX = TRUE
              ))
  })
  output$pep_meta <- renderDT({
    datatable(pep_dt()[[2]] %>% select(-sample),
              options = list(
                autoWidth = TRUE,
                scrollX = TRUE
              ))
  })
  output$neo_plot <- renderPlot({
    dt <- pep_dt()[[1]]
    dt_meta <- pep_dt()[[2]]
    dt_deg <- pep_dt()[[3]]
    
    dt_rank <- apply(dt[,3:ncol(dt)],1,min)
    dt_rank <- data.frame(pep = dt$peptide, min_rank = dt_rank)
    dt_meta <- left_join(dt_meta, dt_rank %>% rename(seq = pep))
    dt_meta_summ <- dt_meta %>% 
      group_by(rMATS_type2, rMTS_ID) %>% 
      summarise(min_bind = min(min_rank),
                gene = unique(geneSymbol),
                dpsi = unique(IncLevelDifference)) %>% ungroup()
    
    dt_deg <- dt_deg %>% 
      filter(symbol %in% dt_meta_summ$gene) %>% 
      select(symbol, log2FoldChange) %>% 
      group_by(symbol) %>% 
      slice_max(abs(log2FoldChange), with_ties = F, na_rm = T) %>% ungroup() %>% 
      filter(nchar(symbol) > 0) %>% as.data.frame()
    
    dt_meta_summ <- left_join(dt_meta_summ, dt_deg %>% rename(gene = symbol))
    dt_meta_summ <- dt_meta_summ %>% mutate(abs_dpsi = abs(dpsi)) %>% 
      mutate(binding = 1/min_bind)
    ###label top ten events
    dt_sort <- dt_meta_summ %>% 
      arrange(desc(log2FoldChange), desc(abs_dpsi), desc(binding)) %>% 
      rowwise() %>% 
      mutate(ids = paste(gene, rMATS_type2, rMTS_ID, sep = "-")) %>% ungroup()
    if (nrow(dt_sort) > 5){
      dt_sort$ids[6:nrow(dt_sort)] <- ""
    }
    ggscatter(dt_sort,color ="binding",y="dpsi",x = "log2FoldChange",
              ylab = "Delta PSI", label = "ids", repel = TRUE,
              font.label = c(12, "plain", "#4c91c1"),size=4)+
      scale_color_gradient2(
        low = '#FEE5C1',
        high = '#B70503',
        mid = '#FC8C59',
        limits = c(range(dt_meta_summ$binding)[1], range(dt_meta_summ$binding)[2]),
        midpoint = (range(dt_meta_summ$binding)[2] - range(dt_meta_summ$binding)[1]) /2
      )+labs(color = "Binding")
  })
} 
# ============================
# 3. 运行应用
# ============================
shinyApp(ui = ui, server = server)