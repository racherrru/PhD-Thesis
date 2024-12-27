setwd("~/Desktop/R/04_Gene_Enrichment/Results")
pacman::p_load(tidyverse, readxl, janitor, ggrepel, ggfortify, patchwork, 
               clusterProfiler, DESeq2)
rm(list=ls())

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("dotplot", "clusterProfiler")

## Hallmark Analysis ----
### Overview ----

load("kegg.rda", verbose = TRUE)
load("viz.rda", verbose = TRUE)

raw.list <- lapply(raw, function(df){
  df@result %>% select(ID, NES, p.adjust) 
})

names(raw.list) <- gsub("^\\+ ", "", names(raw.list))

res.all <- bind_rows(raw.list, .id = "analysis") %>% 
  remove_rownames() %>% 
  mutate( 
    analysis = factor(analysis, levels = names(raw.list)),
    mylabel  = ifelse(p.adjust < 0.1, ID, "") )

## Volcano plot for GSEA
ggplot(res.all, aes(x = NES, y = -log10(p.adjust))) + 
  geom_point() + 
  geom_text_repel( aes(label = mylabel, col = factor(sign(NES))), max.overlaps = 5, size = 2) +
  facet_wrap(~ analysis, ncol = 1) +
  scale_color_manual( values = c("-1"="blue", "1"="red")) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 30),
        axis.text.x = element_text(size =30),
        strip.text = element_text(size = 30))+
  theme_bw()
ggsave("volcano_like_plot_kegg.png", height = 15, width = 6)


## Custom dotplot
top <- res.all %>% filter(p.adjust < 0.1, ) %>% distinct(ID) %>% pull()

res.sub <- res.all %>%
  filter(ID %in% top, p.adjust < 0.1)

levs <- res.sub %>% group_by(ID) %>% summarize( total = sum(NES) ) %>% arrange(total) %>% pull(ID)

res.sub <- res.sub %>% mutate(ID = factor(ID, levels = levs))

ggplot(res.sub, aes(x = analysis, y = ID, fill = NES, size = -log(p.adjust))) +
  geom_point(shape = 21, col = "grey80") + 
  scale_fill_gradient2(low = "blue", midpoint = 0, mid = "white", high = "red") +
  labs(x = NULL, y = NULL, title = "GSEA with KEGG", size = bquote(~-log[10] ~ italic(FDR))) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold.italic", hjust = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 0.5))
ggsave("GSEA_with_KEGG.png", height = 15, width = 10)

### Make pheatmap ----

for(analysis in names(raw)){
  
  print(analysis)
  
  out_folder <- file.path("pheatmap_KEGG", make_clean_names(analysis))
  
  res <- raw[[analysis]]
  
  sig.pathways <- res@result %>% 
    filter(p.adjust < 0.1) %>% 
    filter(ID != "REFERENCE_BCR_PLCG_CALCINEURIN_SIGNALING_PATHWAY") %>%
    filter(ID != "REFERENCE_CA2_CAM_CN_SIGNALING_PATHWAY") %>%
    filter(ID != "PATHOGEN_HIV_GP120_TO_CXCR4_GNAQ_PLCB_G_CALCINEURIN") %>%
    pull(ID)
  
  for(path in sig.pathways){
    
    goi <- res@result %>% 
      filter(ID %in% path) %>% 
      separate_rows(core_enrichment, sep = "/") %>% 
      select(gene_name = core_enrichment) %>% 
      left_join(gene_ann, by = "gene_name") %>% 
      filter(gene_type == "protein_coding") %>% 
      distinct(gene_name, .keep_all = TRUE)
    
    cat(nrow(goi), "\t", path, "\n")
    
    hm_mat <- goi %>%
      left_join(assay(vst_expr) %>% data.frame %>% rownames_to_column("ensid"), by = "ensid") %>% 
      column_to_rownames("gene_name") %>% 
      select(-ensid, -gene_type) %>% 
      data.matrix()
    
    phm <- pheatmap::pheatmap(hm_mat,
                              scale = "row", 
                              color = colorRampPalette(c("navy", "white", "red"))(20),
                              # border_color = "white",
                              
                              main = paste0("Core enrichment: ", path, " (", analysis, ")"),
                              
                              annotation_col    = ann_col,
                              annotation_colors = ann_colors,
                              
                              cluster_cols = FALSE,
                              # clustering_distance_cols = "correlation"                          
                              show_colnames = FALSE,
                              gaps_col = 9,
                              
                              cluster_rows = TRUE,
                              clustering_distance_rows = "correlation",                          
                              fontsize_row = min(500/nrow(goi), 12), 
                              
                              filename = paste0(out_folder, "/", path, ".png")
    )
    
    rm(path, goi, hm_mat, phm)
    
  }
  
}

sel <- res@result %>% filter(p.adjust < 0.05) %>% pull(ID)
clusterProfiler::dotplot(res)
dotplot(res, showCategory = sel)
