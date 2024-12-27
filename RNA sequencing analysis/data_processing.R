# Setup ----
setwd("~/Desktop/R/01_Processing")
pacman::p_load(tidyverse, janitor)
rm(list=ls())

# Load counts file ----
counts <- read.csv("countsfile.csv", row.names=1, as.is=T) %>% data.matrix()
dim(counts) # 31479 rows 18 columns

rs <- rowSums(counts > 0)
hist(rs)
table(rs)

sampleinfo <- read.csv("metadata.csv", row.names = 1, as.is=T) %>% 
  mutate( Condition  = factor(Condition, 
                              levels = c("Unexposed_Untreated", "Unexposed_PM1PAH", "Unexposed_PAH",
                                         "Exposed_Untreated", "Exposed_PM1PAH", "Exposed_PAH")),
          Treatment  = factor(Treatment, levels = c("None", "PM1PAH", "PAH")),
          UV_exposed = factor(UV_exposed, levels = c("No", "Yes")),
          Batch      = factor(paste0("B", Batch), levels = c("B1", "B2", "B3")))

table(sampleinfo$UV_exposed, sampleinfo$Treatment)
#        None PM1PAH PAH
# No     3      3   3
# Yes    3      3   3

table(sampleinfo$Condition)
# Unexposed_Untreated    Unexposed_PM1PAH       Unexposed_PAH   Exposed_Untreated      Exposed_PM1PAH 
# 3                   3                   3                   3                   3 
# Exposed_PAH 
# 3 

geneinfo <- read.csv("annotation.csv", row.names = 6, as.is=T) %>%
  select(-gene_id, -ont)

save(counts, geneinfo, sampleinfo, rs, file = "~/Desktop/R/02_Differential_Gene_Analysis/processed.rda")
