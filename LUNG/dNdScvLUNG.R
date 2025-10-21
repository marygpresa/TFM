# source https://github.com/im3sanger/dndscv
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager :: install(c("SummarizedExperiment","TCGAbiolinks","GenomicRanges", "Biostrings", "IRanges", "BSgenome.Hsapiens.UCSC.hg38"))

install.packages("remotes")
remotes::install_github("im3sanger/dndscv")

library(TCGAbiolinks)
library(SummarizedExperiment)
library("dndscv")
library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")
library(dplyr)
library(ggplot2)
library(ggrepel) #truncate axis
library(grid) # for arrow in plot
library(scales)
data("dataset_simbreast", package="dndscv")

load("/Users/mariagranados/TFM/LUNG/dndsout.RData")
load("/Users/mariagranados/TFM/LUNG/maf.RData")

gdcprojects <- getGDCprojects()
getProjectSummary("TCGA-LUAD")
#I've got simple nucleotide variation data for lung adenocarcinoma

#to check the type of data
query_preview <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation"
)

preview <- getResults(query_preview)
str(preview)
unique(preview[, c("data_type", "analysis_workflow_type")])


# important to remember the TCGA-LUAD Genome of reference: hg38
# and the dndscv package requires hg19 unless I specify " refdb="hg38" as an argument in dndscv to use the default GRCh38/hg38 "

table(preview$data_type, preview$workflow_type)


query_somatic <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Simple Nucleotide Variation",
                  data.type = "Masked Somatic Mutation",
                  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query_somatic)
maf <- GDCprepare(query_somatic)
View(maf)
colnames(maf)
dim(maf)
length(unique(maf$Tumor_Sample_Barcode))
table(maf$Variant_Classification)

# need to filter conttiguous mutations within a sample (not single nucleotide)
unique(maf$Mutation_Status)
#returns "Somtic" meaning I don't need to filter out germline mutations

# create a data frame with only single nucleotide variants (SNVs) 
unique(maf$Variant_Type)
# still has "SNP" "DEL" "INS" "ONP" "TNP"
# need to filter out ONP: Oligo-nucleotide polymorphism  (4+ substitutions) and 
# TNP: Tri-nucleotide polymorphism (3substitutions)
maf <- maf[maf$Variant_Type %in% c("SNP", "DEL", "INS"), ]
save(query_somatic,maf, file = "maf.RData")

# in coding regions, which is ideal for dndscv
mutations <- data.frame(
  sampleID = maf$Tumor_Sample_Barcode,
  chromosome = maf$Chromosome,
  position = maf$Start_Position,
  ref = maf$Reference_Allele,
  mut = maf$Tumor_Seq_Allele2
)
View(mutations)
# need to change chromosome names to match the dndscv package removing "chr" before each
mutations$chromosome <- gsub("^chr", "", mutations$chromosome)

dndsout <- dndscv(mutations, refdb = "hg38")

# check outputs related to Driver discovery (positive selection) in cancer exomes/genome
globaldnds <- dndsout$globaldnds           # Global dN/dS values
print(dndsout$globaldnds)    # table with the global MLEs for the dN/dS ratios across all genes
# other outputs annotated table of coding mutations (*annotmuts*), MLEs of mutation rate parameters (*mle_submodel*), lists of samples and mutations excluded from the analysis and a table with the observed and expected number of mutations per gene (*genemuts*)
head(dndsout$annotmuts)
print(dndsout$nbreg$theta)

# Genes under positive selection plot
## Extract the data for plotting
#names_to_plot <- dndsout$globaldnds$name
#mle_values <- dndsout$globaldnds$mle

#plot MLE 
plotdf <- dndsout$globaldnds %>%
  filter(name != "wall") %>%  # remove overall summary
  mutate(name = recode(name,
                       "wmis" = "Missense",
                       "wnon" = "Nonsense",
                       "wspl" = "Splice-site",
                       "wtru" = "Truncating"))

# Plot ω (dN/dS MLE) with 95% CIs
ggplot(plotdf, aes(x = name, y = mle, fill = name)) +
  geom_col() +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh), width = 0.15) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Mutation Type",
    y = expression(omega~"(dN/dS MLE)"),
    title = "Selection Pressure by Mutation Type (dNdScv)"
  ) +
  theme_minimal(base_size = 14)

# Genes under positive selection (likely cancer drivers)
sel_cv <- dndsout$sel_cv
sig_genes <- sel_cv[sel_cv$qglobal_cv < 0.1, ]   # FDR < 10% (false discovery rate cutoff)
head(sig_genes)
write.csv(sel_cv, file = "sel_cv_dNdScv.csv", row.names = FALSE)
View(sel_cv)
# volcano plot significance vs effect
# Identify the extreme outlier(s)
extreme_gene <- sel_cv[which.max(sel_cv$wmis_cv), ]

# Volcano plot
sig_thresh <- 0.1 # Significance threshold
sel_cv$Extreme <- sel_cv$wmis_cv > 100 | -log10(sel_cv$qglobal_cv) > 15 # Add a flag for extreme points

ggplot(sel_cv, aes(x = wmis_cv, y = -log10(qglobal_cv))) +
  # all genes
  geom_point(alpha = 0.6, color = "gray40") + 
  # significant genes
  geom_point(data = subset(sel_cv, qglobal_cv < sig_thresh & !Extreme),
             aes(x = wmis_cv, y = -log10(qglobal_cv)),
             color = "red", size = 2) +
  # extreme significant points separately
  geom_point(data = subset(sel_cv, Extreme & qglobal_cv < sig_thresh),
             aes(x = wmis_cv, y = -log10(qglobal_cv)),
             color = "blue", size = 3, shape = 17) + # triangle for emphasis
  
  # significance line
  geom_hline(yintercept = -log10(sig_thresh), linetype = "dashed", color = "darkgray") +
  
  # labels for non-extreme significant points
  geom_text_repel(data = subset(sel_cv, qglobal_cv < sig_thresh & !Extreme),
                  aes(label = gene_name),
                  size = 3,
                  max.overlaps = 30) +
  
  # labels for extreme points in blue
  geom_text_repel(data = subset(sel_cv, Extreme & qglobal_cv < sig_thresh),
                  aes(label = gene_name),
                  color = "blue",
                  size = 3,
                  max.overlaps = 30) +
  
  # compress axes but keep all points
  scale_x_continuous(trans = "pseudo_log") +
  scale_y_continuous(trans = "pseudo_log") +
  
  theme_minimal(base_size = 14) +
  labs(title = "Genes under Selection in LUAD (dNdScv)",
       x = "wmis (Missense effect size)",
       y = "-log10(q-value)") +
  theme(legend.position = "none")

# table of selected genes
knitr::kable(sig_genes[, c("gene_name", "wmis_cv", "qglobal_cv")],
             caption = "Genes under positive selection in TCGA-LUAD (dNdScv)",
             digits = 3)

# Save only the columns of interest
write.csv(sig_genes[, c("gene_name", "wmis_cv", "qglobal_cv")],
          file = "significant_genes_dNdScv.csv",
          row.names = FALSE)

# types of mutations
dim(dndsout$annotmuts)
head(dndsout$annotmuts)
View(dndsout$annotmuts)
mut_counts <- table(dndsout$annotmuts$impact)
barplot(mut_counts, las = 2, col = "steelblue",
        main = "Distribution of mutation consequences",
        ylab = "Count")

#plot
mut_counts_df <- as.data.frame(mut_counts)
View(mut_counts_df)
write_csv(mut_counts_df, file = "mutation_spectrum_dNdScv.csv")
colnames(mut_counts_df) <- c("Impact", "Count")
ggplot(mut_counts_df, aes(x = Impact, y = Count, fill = Impact)) +
  geom_col() +
  theme_minimal(base_size = 14) +
  labs(title = "Mutation Spectrum in LUAD (dNdScv)")

save(mutations,dndsout, file = "dndsout.RData")
