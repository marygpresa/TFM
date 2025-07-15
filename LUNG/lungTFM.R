# Descarga de datos
install.packages("installr")
library(installr)
updateR()

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager :: install(c("org.Hs.eg.db","edgeR","ComplexHeatmap","SummarizedExperiment","TCGAbiolinks"))

library(limma)
library(dplyr)
library(tidyr)
library(edgeR)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(TCGAbiolinks)
library(mafR)
library(DT)
library(biomaRt)


gdcprojects <- getGDCprojects()
getProjectSummary("TCGA-LUAD")

## creating a query
query_luad <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Transcriptome Profiling",
                      experimental.strategy = "RNA-Seq",
                      workflow.type = "STAR - Counts",
                      access = "open"
                      )
res <- getResults(query_luad)
View(res)

# download all quey data
GDCdownload(query_luad)

# obtain count matrix (gene expresion counts)
luad_data <- GDCprepare(query_luad, summarizedExperiment = TRUE)
## Check available assays
assayNames(luad_data)
mat <- assay(luad_data, "unstranded") #counts matrix
View(mat)
dim(mat) # 60660 genes 600 samples

#obtain and visualize clinical data detailed
##option 1
clinical_data <- GDCquery_clinic("TCGA-LUAD", type = "clinical")
View(clinical_data)
dim(clinical_data)

##option 2
# obtain metadata from genes (rowData) and clinical data (colData) (optional but will help get samples IDs)
gene_metadata <- as.data.frame(rowData(luad_data))
clinical <- as.data.frame(colData(luad_data))
View(clinical)
dim(clinical) # will use this df as clinical info, it has a bit more columns and info than the clinical_data generated.

# Find duplicates
table(duplicated(rownames(mat))) #  no duplicated genes
table(duplicated(colnames(mat))) # no duplicated samples
all(rownames(clinical) == colnames(mat))  #TRUE

#Set samples as rownames for clinical dataset
rownames(clinical) <- clinical$sample
dim(clinical)
length(colnames(mat)%in%clinical$sample) #600, all samples match
all(rownames(clinical)==colnames(mat)) # all samples are in same order

# Gene Annotation
# Conexión a Ensembl
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Limpiar las versiones del Ensembl ID en tu matriz (en mi caso los nombres de genes en la matriz tenian codgigoEnsembl.version)
# al tener la version puede no ser compatible con otras rutas en la anotacion, por lo tanto las eliminare de los ids
rownames(mat) <- gsub("\\..*", "", rownames(mat))
rownames(mat)

# Sacar metadata (anotacion ensembl, entrez y hgnc)
genes_info <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(mat),
  mart = mart
)

genes_info

#verify class
class(mat) #matrix
class(clinical) #data.frame

set.seed(234)

# Create DGEList
de.exp <- DGEList(counts = as.matrix(mat)) 

# Keep genes with at least 2 samples with CPM > 1
keep.exprs <- rowSums(cpm(de.exp) > 1) >= 2
de.exprs <- de.exp[keep.exprs, , keep.lib.sizes = FALSE] 

# TMM normalization
de.exprs_norm <- calcNormFactors(de.exprs, method = "TMM")
de.exprs_norm

# Log2 CPM-transformed normalized expression
tmmexp <- cpm(de.exprs_norm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 3)
hist(tmmexp, main = "TMM normalized log-CPM")

# Log2 CPM-transformed unnormalized expression
no_tmmexp <- cpm(de.exp, log = TRUE)
hist(no_tmmexp, main = "Raw log-CPM")

# Boxplots to compare before and after normalization
par(mfrow = c(1,2)) 
boxplot(no_tmmexp, las = 2, main = "A. Unnormalised data", ylab = "Log-cpm")
boxplot(tmmexp, las = 2, main = "B. Normalised data", ylab = "Log-cpm")
par(mfrow = c(1,1))


#design<-model.matrix(~tissue_type)
#design

#check pca
pca<-prcomp(t(tmmexp))
plot(pca$x[,1],pca$x[,2],col=(sinfo$Project.ID=="TCGA-LUAD")+1)
legend("topright",legend=c("Tumor","Normal"),col=c(1,2),pch=1)
title("before combat")
plot(pca$x[,1],pca$x[,2],col=(sinfo$Sample.Type=="Primary Tumor")+1)
legend("topright",legend=c("Normal","Tumor"),col=c(1,2),pch=1)
title("before combat")