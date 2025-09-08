load("/Users/mariagranados/TFM/LUNG/dndsout.RData")
load("/Users/mariagranados/TFM/LUNG/mats.RData")
load("/Users/mariagranados/TFM/LUNG/hipathia.RData")
load("/Users/mariagranados/TFM/LUNG/ko_egfr.RData")
load("/Users/mariagranados/TFM/LUNG/genes_info.RData")
load("/Users/mariagranados/TFM/dndsout.RData")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GOSemSim", "AnnotationHub", "org.Hs.eg.db",
                       "clusterProfiler", "DOSE", "meshes", "ReactomePA",
                       "pathview", "enrichplot","Rcpp","topGO"))
library(tidyverse)
library(dplyr)

#check logged data correctly
dndsout
genes_info
colnames(sig_change)

#------------------------------------------------------------------------------------
# This script is to track back samples with missense mutations in EGFR and assign
# the type of mutations present on each significant pathway
#------------------------------------------------------------------------------------

#selecting wmiss samples significant from dNdS analysis and relevant for egfr
# remember sig_genes <- sel_cv[sel_cv$qglobal_cv < 0.1, ] get a subset similarly but with sample names
sel_sig <- subset(sel_cv, qglobal_cv < 0.1)

#select egfr relevant
egfr_sel <- subset(sel_sig, gene_name == "EGFR")
egfr_sel
#get annotated mutations for those samples
egfr_muts <- subset(dndsout$annotmuts,
                    gene %in% egfr_sel$gene_name &
                      impact == "Missense")

# extract short IDs (first 12 characters)
egfr_samples <- unique(substr(egfr_muts$sampleID, 1, 12))
egfr_samples

# match expression/pathways samples
sig_change_short <- sig_change
colnames(sig_change_short) <- substr(colnames(sig_change), 1, 12)
common_samples <- intersect(colnames(sig_change_short), egfr_samples)

# Check how many overlap
common_samples
length(common_samples) #37 samples with egfr missense mutations and sig pathways
head(common_samples)

# Table showing those 37 samples specific mutations
egfr_muts$Sample <- substr(egfr_muts$sampleID, 1, 12)

# Samples with EGFR missense mutations and significant FC pathway changes
egfr_muts_sel <- egfr_muts %>%
  filter(Sample %in% common_samples) %>%
  dplyr::select(Sample, chr, pos, ref3_cod, mut3_cod,
         aachange, ntchange, codonsub)

View(egfr_muts_sel)
dim(egfr_muts_sel)
unique(egfr_muts_sel$codonsub)
write.csv(egfr_muts_sel, file = "egfr_missense_mutations_samples.csv", row.names = FALSE)

#------------------------------------------------------------------------------------
# 
colnames(sig_change_short)
# Reshape pathway results into long format (pathway x sample)
sig_change_unique <- sig_change_short
colnames(sig_change_unique) <- substr(colnames(sig_change), 1, 12)
colnames(sig_change_unique) <- make.unique(colnames(sig_change_unique)) 
View(sig_long)

sig_long <- sig_change_unique %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Pathway") %>%
  pivot_longer(-Pathway, names_to = "Sample", values_to = "FC") %>%
  filter(Sample %in% common_samples) %>%   # only 37 EGFR samples
  filter(FC > 2 | FC < -2)                 # only significant pathways

# detailed tables with samples, pathways and mutations
mut_pathways <- inner_join(
  egfr_muts_sel,
  sig_long,
  by = "Sample",
  relationship = "many-to-many"
)
View(mut_pathways)
write.csv(mut_pathways, file = "egfr_missense_mutations_pathways/sample.csv", row.names = FALSE)


# Final table with mutations per pathway indicating #samples per each 
mut_pathways_summary <- mut_pathways %>%
  group_by(Pathway, chr, pos,ref3_cod, mut3_cod,
           aachange, ntchange, codonsub) %>%
  summarise(n_samples = n_distinct(Sample), .groups = "drop")
View(mut_pathways_summary)
write.csv(mut_pathways_summary, file = "egfr_missense_mutations_pathways_summary.csv", row.names = FALSE)


