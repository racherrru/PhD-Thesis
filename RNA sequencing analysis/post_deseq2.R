# Setup ----
setwd("~/Desktop/R/03_Post_DESeq")
pacman::p_load(tidyverse, readxl, janitor, ggrepel, ggfortify,
               genefilter, car, edgeR, limma, DESeq2,
               pheatmap, EnhancedVolcano, patchwork, ggalt)
rm(list=ls())

load("data_w_uv.rda", verbose = TRUE)

## Number of significant genes affected by pollution
sapply(res, function(x) table(x$sig)) %>% t()
#            -1     0   1
# PM1PAH     56 17353  33
# PAH        26 16125  70
# UV|PM1PAH  28 19833  22
# UV|PAH    212 18861 200
# UV         34 16153  34

## Venn diagrams ----
hits <- lapply(res, function(x) x %>% filter(sig != 0) %>% pull(SYMBOL))
names(hits)

venn_lab <- c("+ PM1PAH", "+ PAH", "+ UV + PM1PAH", "+ UV + PAH")
venn_col <- c("orange2","magenta2","brown","red2")

ggvenn::ggvenn(hits,
               fill_color = venn_col,
               fill_alpha = 0.5,
               text_size = 3,
               set_name_size = 3)
ggsave("venn.png")

ggvenn::ggvenn(hits[c(2,4)],
               fill_alpha = 0.5,
               text_size = 3,
               set_name_size = 3)


rm(venn_col, venn_lab, hits)

## Merge and write out
res2 <- res
UV   <- res2[["UV"]]

for(i in names(res2)){
  colnames(res2[[i]])[4:7] <- paste0(i, "_", colnames(res2[[i]])[4:7])
}

res2 <- base::Reduce("full_join", res2)

openxlsx::write.xlsx(res2, file = "DEG_results_tmp.xlsx")
openxlsx::write.xlsx(res, file = "DEG_separated_tmp.xlsx")

rm(i)

## Impact of PAH
comb <- res2 %>% 
  dplyr::select(SYMBOL, 
                LFC.exposed   = "UV|PAH_LFC",
                LFC.unexposed = "PAH_LFC") %>% 
  mutate(diff = LFC.exposed - LFC.unexposed)

sel <- comb %>% 
  filter(abs(diff) > 1 | abs(LFC.exposed) > 2.5 | abs(LFC.unexposed) > 2.5) %>% 
  filter(str_starts(SYMBOL, "ENSG", negate = TRUE))

b1 <- ggplot(comb, aes(x=LFC.unexposed, y=LFC.exposed)) +
  geom_point(size = 0.5, col = "magenta3") +
  geom_text_repel(data = sel, aes(label = SYMBOL), col="black") +
  geom_abline(intercept = c(-1, 0, 1), linetype = c(3, 1, 3), col = "grey50") +
  theme_bw() +
  labs(title = "PAH vs. untreated (log2 Fold Change)", 
       x = "No UV Exposure", 
       y = "UV Exposure")

b1

# Impact of PM1PAH
comb2 <- res2 %>% 
  dplyr::select(SYMBOL, 
                LFC.exposed   = "UV|PM1PAH_LFC",
                LFC.unexposed = "PM1PAH_LFC") %>% 
  mutate(diff = LFC.exposed - LFC.unexposed)

sel2 <- comb2 %>% 
  filter(abs(diff) > 1 | abs(LFC.exposed) > 2.5 | abs(LFC.unexposed) > 2.5) %>% 
  filter(str_starts(SYMBOL, "ENSG", negate = TRUE))

b2 <- ggplot(comb2, aes(x=LFC.unexposed, y=LFC.exposed)) +
  geom_point(size = 0.5, col = "orange2") +
  geom_text_repel(data = sel, aes(label = SYMBOL), col="black") +
  geom_abline(intercept = c(-1, 0, 1), linetype = c(3, 1, 3), col = "grey50") +
  theme_bw() +
  labs(title = "PM1-PAH vs. untreated (log2 Fold Change)", 
       x = "No UV Exposure", 
       y = "UV Exposure")

b1/b2
ggsave("comparison.png", width = 7, height = 12)

rm(b1, b2, comb, comb2, comparison, sel, sel2)

## Investigate one gene ----
sapply( res, function(x) x %>% filter(SYMBOL == "EDNRA") ) %>% t()
#               ensid                SYMBOL    LFC      p            FDR          sig
# + PM1PAH      "ENSG00000101096.20" "NFATC2"  2.919675 5.018118e-16 1.75052e-12  1  
# + PAH         "ENSG00000101096.20" "NFATC2"  2.854621 2.159215e-15 3.502462e-12 1  
# + UV + PM1PAH "ENSG00000101096.20" "NFATC2"  1.909295 1.175467e-07 0.000116859  1  
# + UV + PAH    "ENSG00000101096.20" "NFATC2"  2.062533 1.031176e-08 1.195177e-06 1  

plot_genes <- function(genename){
  
  meta <- colData(vst_expr) %>% 
    data.frame() 
  id <- gene_ann %>% filter(gene_name == genename) %>% pull(ensid)
  
  meta$expr <- assay(vst_expr)[id, ]
  
  g <- ggplot(meta, aes(x=Treatment, y=expr, fill=UV_exposed)) +
    geom_dotplot(binaxis = "y", stackdir = "center") +
    facet_wrap( ~ UV_exposed, labeller = labeller(UV_exposed =  c("No" = "No UV Exposure",
                                                                  "Yes" = "UV Exposure"))) +
    labs(x = NULL , y = NULL, title = genename) + 
    theme_bw() +
    theme(legend.position = "none", 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.title.x = element_text(size = 18, colour = "black"),
          axis.title.y = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          plot.title = element_text(size = 18, colour = "black"))
  
  return(g)
}

# Function that scales all text within a plot up or down based on the scaling factor
safe_scale_text <- function(element, scale_factor) {
  if (!is.null(element$size)) {
    element$size <- element$size * scale_factor
  }
  
  return(element)
}

scale_plot_text <- function(plots, scale_factor=1) {
  scaled_plots <- lapply(plots, function(plot) { 
    current_theme <- plot$theme
    plot + theme(
      plot.title   = safe_scale_text(current_theme$plot.title, scale_factor),
      axis.title.x = safe_scale_text(current_theme$axis.title.x, scale_factor),
      axis.title.y = safe_scale_text(current_theme$axis.title.y, scale_factor),
      axis.text.x  = safe_scale_text(current_theme$axis.text.x, scale_factor),
      axis.text.y  = safe_scale_text(current_theme$axis.text.y, scale_factor)
    )
  })
  return(scaled_plots)
}

gene_list_down <- c("PALMD", "S1PR1", "STC1", "TNC", "DELEC1", "OR10AB1P", 
                    "WDR72", "OVCH2", "ROS1", "LINC01537", "CXCL14", 
                    "ZNF812P", "GBP2")
gene_list_up <- c( "CYP1B1", "ALDH1A3", "LINC01531", 
                   "CYP1A1", "NFATC2", "LGI3", "GDA", "WNT9A", "BIVM-ERCC5", 
                   "DCT", "EREG", "ADGRF1")


gene_plot_list <- list()

for (i in seq_along(gene_list_down)) {
  gene_plot_list[[i]] <- plot_genes(gene_list_down[[i]])
}

scaled_plot_list <- scale_plot_text(gene_plot_list, 1)
plots_col <- 2
tmp <- wrap_plots(scaled_plot_list, 
                             ncol = plots_col, 
                             nrow = ceiling(length(scaled_plot_list) / plots_col))
tmp + plot_annotation(title = "Significantly downregulated genes (VST Expression)", 
                      theme = theme(plot.title = element_text(size = 20, 
                                                              face = "bold", 
                                                              hjust = 0.5)))

ggsave("downreg.png", width =14, height = 20)

rm(scaled_plot_list)

gene_plot_list <- list()

for (i in seq_along(gene_list_up)) {
  gene_plot_list[[i]] <- plot_genes(gene_list_up[[i]])
}

scaled_plot_list <- scale_plot_text(gene_plot_list, 1)
plots_col <- 2
tmp <- patchwork::wrap_plots(scaled_plot_list, 
                             ncol = plots_col, 
                             nrow = ceiling(length(scaled_plot_list) / plots_col))
tmp + plot_annotation(title = "Significantly upregulated genes (VST Expression)", 
                      theme = theme(plot.title = element_text(size = 20, 
                                                              face = "bold", 
                                                              hjust = 0.5)))

ggsave("upreg.png", width =14, height = 20)

# To plot individual genes
plot_genes("GBP2")
plot_genes("ABCB1")
plot_genes("MYOF") #KEEP THIS
plot_genes("NFE2L2") #Transcription Factor

rm(scaled_plot_list, plots_col, tmp, gene_list_down, gene_list_up, plot_genes,
   safe_scale_text, safe_plot_text, scale_plot_text, gene_plot_list)


# Volcano Plot ----
gene_list <- c("PALMD", "S1PR1", "STC1", "TNC", "DELEC1", "OR10AB1P", 
               "WDR72", "OVCH2", "ROS1", "LINC01537", "CXCL14", 
               "ZNF812P", "GBP2", "CYP1B1", "ALDH1A3", "LINC01531", 
               "CYP1A1", "NFATC2", "LGI3", "GDA", "WNT9A", "BIVM-ERCC5", 
               "DCT", "EREG", "ADGRF1")

g_vol <- list()
for (i in names(res)) {
  print(i)
  g_vol[[i]] <- EnhancedVolcano(res[[i]],
                                lab = res[[i]]$SYMBOL, x = "LFC", y = "FDR",
                                pointSize = 1.0,
                                labSize = 4.0,
                                selectLab = paste0(gene_list),
                                #drawConnectors = TRUE,
                                FCcutoff = log2(2), pCutoff = 0.05, xlim = c(-5.5, 5.5),
                                title = i, subtitle = NULL, caption = NULL, 
                                axisLabSize = 12,
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE,
                                ylab = bquote(~-Log[10] ~ italic(FDR)), legendPosition = "none")
}

wrap_plots(g_vol, ncol = 2)
ggsave("volcano.png")

rm(g_vol)

# PHeatmap ----
ann_col <- colData(dds) %>%  
  data.frame() %>%  
  dplyr::select(-sizeFactor, -Batch, -group, -Condition)

gene_list <- ahr_list

glimpse(ann_col)
goi <- gene_ann %>%
  filter(gene_name %in% gene_list) %>% 
  filter((str_starts(gene_type, "protein_coding")))

hm_mat <- goi %>%
  left_join(assay(vst_expr) %>% data.frame %>% rownames_to_column("ensid")) %>% 
  column_to_rownames("gene_name") %>% 
  dplyr::select(-ensid, -gene_type) %>% 
  data.matrix()

ann_colors = list(
  Treatment =   c("None" = "grey80", "PM1PAH" = "grey30", "PAH" = "grey0"),
  UV_exposed =  c("No" = "grey0", "Yes" = "red2")
)

phm <- pheatmap(hm_mat, scale="row",
                annotation_col = ann_col,
                annotation_colors = ann_colors,
                #gaps_col = 9,
                #color = colorRampPalette(c("navy", "white", "red"))(20),
                border_color = "black",
                fontsize = 14, 
                show_colnames = FALSE,
                cluster_cols = FALSE,
                cluster_rows = TRUE
                # filename = "pheatmap.png"
                #clustering_distance_rows = "correlation",
                # clustering_distance_cols = "correlation"
)
ggsave("heatmap.png", phm, width = 20, height = 30)

rm (hm_mat, goi)
