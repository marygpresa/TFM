load("/Users/mariagranados/TFM/LUNG/dndsout.RData")
load("/Users/mariagranados/TFM/LUNG/mats.RData")
load("/Users/mariagranados/TFM/LUNG/hipathia.RData")
load("/Users/mariagranados/TFM/LUNG/ko_egfr.RData")
load("/Users/mariagranados/TFM/LUNG/genes_info.RData")
load("/Users/mariagranados/TFM/dndsout.RData")
load("/Users/mariagranados/TFM/fc.RData")
load("/Users/mariagranados/TFM/sel_cv.RData")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GOSemSim", "AnnotationHub", "org.Hs.eg.db",
                       "clusterProfiler", "DOSE", "meshes", "ReactomePA",
                       "pathview", "enrichplot","Rcpp","topGO"))
library(tidyverse)
library(dplyr)
library(ggplot2)

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
# Genes under positive selection (likely cancer drivers)
sel_cv <- dndsout$sel_cv
sel_sig <- subset(sel_cv, qglobal_cv < 0.1)
View(sel_sig)

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
View(sig_change)
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


sig_long <- sig_change_unique %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Pathway") %>%
  pivot_longer(-Pathway, names_to = "Sample", values_to = "FC") %>%
  filter(Sample %in% common_samples) %>%   # only 37 EGFR samples
  filter(FC > 2 | FC < -2)                 # only significant pathways
View(sig_long)

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

renv::init()


#--------------------------------------------------------------------------------------
#adding metadata to compare per gender, tobacco and age
#--------------------------------------------------------------------------------------

BiocManager :: install(c("org.Hs.eg.db","edgeR","ComplexHeatmap","SummarizedExperiment","TCGAbiolinks","hpAnnot","sva"))
library(SummarizedExperiment)
library(TCGAbiolinks)

gdcprojects <- getGDCprojects()
clinical_data <- GDCquery_clinic("TCGA-LUAD", type = "clinical")
View(clinical_data)

#table with metadata, pathways and somatic mutations
mut_pathways_meta <- inner_join(
  mut_pathways,
  clinical_data %>% select(submitter_id, age_at_diagnosis, gender, exposure_type),
  by = c("Sample" = "submitter_id"),
  relationship = "many-to-many"
)
# will change age_at_diagnosis to years instead of days
mut_pathways_meta <- mut_pathways_meta %>%
  mutate(age_at_diagnosis = round(age_at_diagnosis / 365.25, 1))

View(mut_pathways_meta)
colnames(mut_pathways_meta)

#FC plot
library(reshape2)
cam_melt <- melt(sig_change)
bypass_plot <- ggplot(cam_melt, aes(x = Var1, y = value, fill = Var1)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(trans = "pseudo_log") +
  theme_minimal() +
  xlab("Pathway") +
  ylab("log2(Fold Change)") +
  labs(fill = "Pathways")  # rename legend title


colnames(cam_melt)
View(cam_melt)
colnames(sig_change)
View(sig_change)
View(clinical_data)
View(clinical)


# comparing FC paths with metadata info
grouping_FC <- cam_melt %>%
  left_join(clinical %>% select(barcode, gender, age_at_diagnosis, exposure_type),
            by = c("Var2" = "barcode"))
View(grouping_FC)

library(ggplot2)

bypass_plot <- ggplot(grouping_FC, aes(x = Var1, y = value, fill = gender)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(trans = "pseudo_log") +
  theme_minimal() +
  xlab("Pathway") +
  ylab("log2(Fold Change)") +
  labs(fill = "Gender")
print(bypass_plot)

#will filter for only significant values (><2 FC) for better visualisation
group_FC_filtered <- grouping_FC %>%
  filter(abs(value) > 2) 

# age in years not days
group_FC_filtered <- group_FC_filtered %>%
  mutate(age_years = round(age_at_diagnosis / 365.25, 1))

# Replace +/-Inf with pltos scale max and mins, otherwhise they wont appear in the plots
max_fc <- max(group_FC_filtered$value[is.finite(group_FC_filtered$value)], na.rm = TRUE)

group_FC_plot <- group_FC_filtered %>%
  mutate(value_plot = case_when(
    value == Inf ~ max_fc * 1.05,      # slightly above max
    value == -Inf ~ -max_fc * 1.05,   # slightly below min
    TRUE ~ value
  ))
View(group_FC_plot)

#gender and pathways
bypass_plot1 <- ggplot(group_FC_plot, aes(x = Var1, y = value_plot, fill = gender)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(trans = "pseudo_log") +
  theme_minimal() +
  xlab("Pathway") +
  ylab("log2(Fold Change)") +
  labs(fill = "Gender")
print(bypass_plot1) # wow mostly male


exposure_plot <- ggplot(group_FC_plot, aes(x = Var1, y = value_plot, fill = Var1)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(trans = "pseudo_log") +
  theme_minimal() +
  xlab("Pathway") +
  ylab("log2(Fold Change)") +
  facet_wrap(~ exposure_type)

print(exposure_plot) #all of them smoke/have been exposed to tobacco


#plot considering age + gender + pathways
ggplot(group_FC_plot, aes(x = age_years, y = value_plot, color = gender)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  xlab("Age at diagnosis (years)") +
  ylab("log2(Fold Change)") +
  labs(color = "Gender") +
  scale_y_continuous(trans = "pseudo_log") +
  facet_wrap(~ Var1, scales = "free_y") 


#boxplots highlight dispersion
ggplot(group_FC_plot, aes(x = gender, y = value_plot, fill = gender)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~ Var1, scales = "free_y") +
  theme_minimal() +
  ylab("log2(Fold Change)") +
  xlab("Gender")

save.image(file = "EGFRmutation-pathway.RData")
load("EGFRmutation-pathway.RData")

# oxitocin looks interesting
oxytocin_pathway <- subset(group_FC_plot, Var1 == "Oxytocin signaling pathway: CDKN1A")
View(oxytocin_pathway)

ggplot(oxytocin_pathway, aes(x = age_years, y = value_plot, color = gender)) +
  geom_jitter(width = 0, height = 0.1, size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  xlab("Age at diagnosis (years)") +
  ylab("log2(Fold Change)") +
  labs(color = "Gender") +
  ggtitle("Oxytocin Signaling Pathway: FC vs Age by Gender") +
  theme(plot.title = element_text(hjust = 0.5))
# plot shows no relation FC-age

ggplot(oxytocin_pathway, aes(x = gender, y = value_plot, fill = gender)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  theme_minimal() +
  ylab("log2(Fold Change)") +
  xlab("Gender") +
  ggtitle("Oxytocin Signaling Pathway FC Distribution") +
  theme(plot.title = element_text(hjust = 0.5))
# no differences per gender either


#comparing cohorts before and after KO
library(hipathia)

View(ko_egfr_vals)
ko_egfr_vals.names <- get_path_names(pathways, rownames(ko_egfr_vals))
ko_egfr_vals_cbind <- cbind(ko_egfr_vals.names, ko_egfr_vals)
ko_egfr_vals_cbind <- as.data.frame(ko_egfr_vals_cbind, stringsAsFactors = FALSE)
rownames(ko_egfr_vals_cbind) <- ko_egfr_vals_cbind$ko_egfr_vals.names
View(ko_egfr_vals_cbind)

View(path_vals)
path_vals.names <- get_path_names(pathways, rownames(path_vals))
path_vals_comp <- cbind(path_vals.names, path_vals)
path_vals_comp <- as.data.frame(path_vals_comp, stringsAsFactors = FALSE)
rownames(path_vals_comp) <- path_vals_comp$path_vals.names


oxytocin_pre <- path_vals_comp["Oxytocin signaling pathway: CDKN1A", , drop = FALSE]  # PRE KO
oxytocin_post <- ko_egfr_vals_cbind["Oxytocin signaling pathway: CDKN1A", , drop = FALSE]  # POST KO


oxytocin_long <- bind_rows(
  data.frame(Sample = colnames(oxytocin_pre),
             value = as.numeric(oxytocin_pre),
             Condition = "Pre-KO"),
  data.frame(Sample = colnames(oxytocin_post),
             value = as.numeric(oxytocin_post),
             Condition = "Post-KO") 
                  )
View(oxytocin_long)
oxytocin_long <- oxytocin_long[-1, ] #removing first row with column title

#add clinical data
oxytocin_long <- oxytocin_long %>%
  left_join(clinical %>% select(barcode, gender, age_at_diagnosis),
            by = c("Sample" = "barcode"))%>%
  mutate(age_years = age_at_diagnosis / 365.25)
View(oxytocin_long)


#plot pre-pos KO oxitocin
oxytocin_long$Condition <- factor(oxytocin_long$Condition, levels = c("Pre-KO", "Post-KO"))
ggplot(oxytocin_long, aes(x = Condition, y = value, color = gender)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  theme_minimal() +
  ylab("Pathway activation score") +
  xlab("") +
  labs(color = "Gender") +
  ggtitle("Oxytocin-CDKN1A Pathway activation Pre- and Post-EGFR KO") +
  theme(plot.title = element_text(hjust = 0.5))

View(oxytocin_post)



