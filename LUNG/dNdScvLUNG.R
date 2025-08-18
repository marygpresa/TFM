# source https://github.com/im3sanger/dndscv
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager :: install(c("SummarizedExperiment","TCGAbiolinks"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("remotes")
remotes::install_github("im3sanger/dndscv")
BiocManager::install(c("GenomicRanges", "Biostrings", "IRanges", "BSgenome.Hsapiens.UCSC.hg38"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library("dndscv")
library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")
library(dplyr)
data("dataset_simbreast", package="dndscv")

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
dndsout$globaldnds           # Global dN/dS values
print(dndsout$globaldnds)    # table with the global MLEs for the dN/dS ratios across all genes
# other outputs annotated table of coding mutations (*annotmuts*), MLEs of mutation rate parameters (*mle_submodel*), lists of samples and mutations excluded from the analysis and a table with the observed and expected number of mutations per gene (*genemuts*)
head(dndsout$annotmuts)
print(dndsout$nbreg$theta)

sel_cv_results <- dndsout$sel_cv
print(head(sel_cv_results))

# Genes under positive selection plot
## Extract the data for plotting
names_to_plot <- dndsout$globaldnds$name
mle_values <- dndsout$globaldnds$mle

## bar plot
dndsout_plot <- barplot(height = mle_values, names.arg = names_to_plot,
                        xlab = "dN/dS type", ylab = "MLE",
                        main = "Global dN/dS values",
                        ylim = c(0, 1.2), # Adjust y-axis limits as needed
                        col = "skyblue") # Choose a color



save(mutations,dndsout, file = "dndsout.RData")
