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
  clinical_data %>% select(submitter_id, age_at_diagnosis, gender, exposure_type, ajcc_pathologic_stage),
  by = c("Sample" = "submitter_id"),
  relationship = "many-to-many"
)
# will change age_at_diagnosis to years instead of days
mut_pathways_meta <- mut_pathways_meta %>%
  mutate(age_at_diagnosis = round(age_at_diagnosis / 365.25, 1))

View(mut_pathways_meta)
colnames(mut_pathways_meta)

#FC plot
View(sig_change)
# Convert wide table to long format
cam_melt <- melt(sig_change)
bypass_plot <- ggplot(cam_melt, aes(x = Var1, y = value, fill = Var1)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(trans = "pseudo_log") +
  theme_minimal() +
  xlab("Pathway") +
  ylab("log2(Fold Change)") +
  labs(fill = "Pathways")  # rename legend title
print(bypass_plot)

# are we actually only seing 3-8 relevant patients for FoxO?
# Threshold for “strong change”
threshold <- 2

# Compute number of patients per pathway exceeding threshold
sig_df <- as.data.frame(sig_change)
sig_df$Pathway <- rownames(sig_change)
patient_counts <- sig_df %>%
  mutate(Pathway = rownames(.)) %>%
  rowwise() %>%
  mutate(NumPatients = sum(c_across(-Pathway) > threshold | c_across(-Pathway) < -threshold, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(NumPatients))

ggplot(patient_counts, aes(x = reorder(Pathway, NumPatients), y = NumPatients, fill = NumPatients)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Number of Patients with Strong log2(FC) per Pathway",
       x = "Pathway",
       y = "Number of Patients",
       fill = "Patients") +
  theme_minimal(base_size = 12)

foxo_counts <- subset(patient_counts, grepl("FoxO signaling pathway", Pathway))
ggplot(foxo_counts, aes(x = reorder(Pathway, NumPatients), y = NumPatients, fill = NumPatients)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Number of Patients with Strong log2(FC) per Pathway",
       x = "Pathway",
       y = "Number of Patients",
       fill = "Patients") +
  theme_minimal(base_size = 12)
# yep just three for FC>2 and 8 FC>1 but for oxytocin if about 400 patients affected

# comparing FC paths with metadata info
grouping_FC <- cam_melt %>%
  left_join(clinical %>% select(barcode, gender, age_at_diagnosis, exposure_type),
            by = c("Var2" = "barcode"))
View(grouping_FC)


bypass_plot <- ggplot(grouping_FC, aes(x = Var1, y = value, fill = gender)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(trans = "pseudo_log") +
  theme_minimal() +
  xlab("Pathway") +
  ylab("log2(Fold Change)") +
  labs(fill = "Gender")
print(bypass_plot)

#will filter for only significant values (><2 FC) for better visualization
group_FC_filtered <- grouping_FC %>%
  filter(abs(value) > 1) 

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
#almost all are foxO so will do one as representative

foxoFC <- subset(group_FC_plot, Var1 == "FoxO signaling pathway: FBXO32")
#all foxO results subset
foxoFC_temp <- subset(group_FC_plot, grepl("FoxO signaling pathway", Var1))
#avoiding patient duplicates, keep only unique sample IDs (Var2)
foxoFC <- foxoFC_temp[!duplicated(foxoFC_temp$Var2), ]
View(foxoFC)
table(clinical$ajcc_pathologic_stage[clinical$barcode %in% foxoFC$Var2])

# agrego stage
foxoFC_stage <- foxoFC %>%
  left_join(
    clinical %>% select(barcode, ajcc_pathologic_stage),
    by = c("Var2" = "barcode")
  )
View(foxoFC_stage)
ggplot(foxoFC, aes(x = gender, y = age_years, fill = gender)) +
  geom_boxplot() +
  geom_jitter(aes(color = gender),width = 0.2, alpha = 0.6) +
  theme_minimal() +
  ylab("Age (years)") +
  xlab("Sex") 

ggplot(foxoFC, aes(x = gender, y = value_plot, fill = gender)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_minimal() +
  ylab("log2(Fold Change)") +
  xlab("Sex") 


#boxplots highlight dispersion
ggplot(group_FC_plot, aes(x = gender, y = value_plot, fill = gender)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~ Var1, scales = "free_y") +
  theme_minimal() +
  ylab("log2(Fold Change)") +
  xlab("Sex")

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
  ggtitle("Oxytocin Signaling Pathway: FC vs Age by Sex") +
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
  left_join(clinical %>% select(barcode, gender, age_at_diagnosis,ajcc_pathologic_stage),
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



# adding expression after KO to oxytocin table
oxytocin_post_long <- data.frame(
  Sample_full = colnames(oxytocin_post),
  Expr_post_KO_oxy = as.numeric(as.character(oxytocin_post[1, ]))
)

#trim sample name to match
oxytocin_post_long <- oxytocin_post_long %>%
  mutate(Sample = sapply(strsplit(Sample_full, "-"), function(x) paste(x[1:3], collapse = "-")))

# avoid duplicates
oxytocin_post_long <- oxytocin_post_long %>%
  distinct(Sample, .keep_all = TRUE)
View(oxytocin_post_long)

# merge with mut_pathways_meta
mut_pathways_meta <- mut_pathways_meta %>%
  left_join(oxytocin_post_long %>% select(Sample, Expr_post_KO_oxy),
            by = "Sample")

View(mut_pathways_meta)
write.csv(mut_pathways_meta, file = "egfr_missense_mutations_pathways_metadata.csv", row.names = FALSE)

# pie charts
# gender
gender_counts <- mut_pathways_meta %>%
  filter(aachange == "L858R", Pathway == "Oxytocin signaling pathway: CDKN1A") %>%
  count(gender) %>%
  mutate(percent = 100 * n / sum(n))

ggplot(gender_counts, aes(x = "", y = n, fill = gender)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(percent, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 5) +
  labs(
    title = "Gender distribution for L858R (Oxytocin signaling pathway: CDKN1A)",
    x = NULL, y = NULL, fill = "Gender"
  ) +
  theme_void()

# tobacco pie chart
tobacco_counts <- mut_pathways_meta %>%
  filter(aachange == "L858R",
         Pathway == "Oxytocin signaling pathway: CDKN1A") %>%
  count(exposure_type) %>%
  mutate(percent = 100 * n / sum(n))

ggplot(tobacco_counts, aes(x = "", y = n, fill = exposure_type)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(percent, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 5) +
  labs(
    title = "Tobacco Use in L858R (Oxytocin signaling pathway: CDKN1A)",
    x = NULL, y = NULL, fill = "Smoker"
  ) +
  theme_void()

# stage
stage_counts <- mut_pathways_meta %>%
  filter(aachange == "L858R",
         Pathway == "Oxytocin signaling pathway: CDKN1A") %>%
  count(ajcc_pathologic_stage) %>%
  mutate(percent = 100 * n / sum(n))

ggplot(stage_counts, aes(x = "", y = n, fill = ajcc_pathologic_stage)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(percent, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 5) +
  labs(
    title = "Pathologic Stage in L858R (Oxytocin signaling pathway: CDKN1A)",
    x = NULL, y = NULL, fill = "Stage"
  ) +
  theme_void()


# age - gender
age_gender <- mut_pathways_meta %>%
  filter(aachange == "L858R",
         Pathway == "Oxytocin signaling pathway: CDKN1A")

ggplot(age_gender, aes(x = gender, y = age_at_diagnosis, fill = gender)) +
  geom_boxplot(alpha = 0.8) +
  labs(
    title = "Age distribution by gender (L858R, Oxytocin signaling pathway: CDKN1A)",
    x = "Gender",
    y = "Age (years)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

View(ko_egfr_comp) # has the p-calue and FDRp.value for pre-post KO comparison with wilcoxon





