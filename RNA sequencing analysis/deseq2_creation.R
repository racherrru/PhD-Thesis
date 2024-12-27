# Setup ----
setwd("~/Desktop/R/02_Differential_Gene_Analysis/")
pacman::p_load(tidyverse, janitor, edgeR, limma, DESeq2)
rm(list=ls())
load("processed.rda", verbose = TRUE)

# Make the DESeq object ----
identical( rownames(counts), rownames(geneinfo) )   #TRUE
identical( colnames(counts), rownames(sampleinfo) ) #TRUE

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleinfo, 
                              rowData = geneinfo,
                              design  = ~ -1)

rm(counts, geneinfo, sampleinfo, rs)

# PCA ----
vst_expr <- vst(dds)
design0  <- model.matrix( ~ Condition, data = colData(vst_expr))

assay(vst_expr) <- limma::removeBatchEffect(assay(vst_expr), vst_expr$Batch, design = design0)

# pc0 <- plotPCA(vst_expr, intgroup = c("Treatment","UV_exposed"), ntop = 500)
# pc  <- plotPCA(vst_expr, intgroup = c("Treatment", "UV_exposed"), returnData = TRUE, ntop = 500)

pc0 <- plotPCA(vst_expr, intgroup = "Condition", ntop = 500)
pc  <- plotPCA(vst_expr, intgroup = "Condition", returnData = TRUE, ntop = 500)

col_condition <- c("Unexposed_Untreated" = "black",
                   "Unexposed_PM1PAH"    = "black",
                   "Unexposed_PAH"       = "black", 
                   "Exposed_Untreated"   = "orange",
                   "Exposed_PM1PAH"      = "orange",
                   "Exposed_PAH"         = "orange")

col_label <- c("Untreated", "PM1PAH", "PAH", "UV", "UV|PM1PAH", "UV|PAH")

g <- ggplot(pc, aes(x = PC1, y = PC2, col = Condition, shape = Condition, label = name)) +
  geom_point(size = 5) +
  scale_color_manual(name = "Condition", values = col_condition, aesthetics = c("color", "fill"), labels = col_label) +
  scale_shape_manual(values = c(1, 10, 19, 0, 12, 15)) +
  labs(x = pc0$labels$x, y = pc0$labels$y, title = "Principal Component Analysis") +
  theme(text = element_text(size = 20)) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5))
g <- g + geom_encircle(aes(fill = Condition), alpha = 0.3, show.legend = FALSE, color = NA)
ggsave("pca.png", g)

rm(pc0, pc, design0, g, col_condition, col_label)

# Run the DE analysis ----
gene_ann <- rowData(dds) %>% 
  data.frame() %>% 
  select(gene_name, gene_type) %>% 
  # distinct(gene_name, .keep_all = TRUE) %>% 
  rownames_to_column("ensid")

table(dds$Condition, dds$Batch)
#                     B1 B2 B3
# Unexposed_Untreated  1  1  1
# Unexposed_PM1PAH     1  1  1
# Unexposed_PAH        1  1  1
# Exposed_Untreated    1  1  1
# Exposed_PM1PAH       1  1  1
# Exposed_PAH          1  1  1

dds@design <- ~ -1 + Condition + Batch

dds <- DESeq(dds)


# Extract the coefficients ----
resultsNames(dds)

get_results <- function(obj, name, num, denom, prefix = NULL, lfc_thr = log2(2), fdr_thr = 0.05){
  
  res <- results(dds, contrast = c(name, num, denom)) %>% 
    data.frame() %>% 
    rownames_to_column("ensid") %>%
    left_join(gene_ann) %>% 
    dplyr::select(ensid, SYMBOL=gene_name, baseMean, LFC=log2FoldChange, p=pvalue, FDR=padj) %>%
    filter(!is.na(LFC), !is.na(FDR)) %>%
    arrange(p) %>% 
    mutate (sig = sign(LFC) * (abs(LFC) > lfc_thr & FDR < fdr_thr)) 
  
  if( !is.null(prefix) ){
    colnames(res)[3:7] <- paste0(prefix, "_", colnames(res)[3:7])
  }
  
  return(res)  
}

## Genes affected by pollution and UV
res <- list()
res[["PM1PAH"]]      <- get_results(dds, "Condition", "Unexposed_PM1PAH", "Unexposed_Untreated")
res[["PAH"]]         <- get_results(dds, "Condition", "Unexposed_PAH",    "Unexposed_Untreated")
res[["UV|PM1PAH"]]   <- get_results(dds, "Condition", "Exposed_PM1PAH",   "Exposed_Untreated")
res[["UV|PAH"]]      <- get_results(dds, "Condition", "Exposed_PAH",      "Exposed_Untreated")
res[["UV"]]          <- get_results(dds, "Condition", "Exposed_Untreated","Unexposed_Untreated")

save(dds, res, gene_ann, vst_expr, file = "~/Desktop/R/03_Post_DESeq/data_w_uv.rda")


