# TFM 

#configuro todo para que se automatice y suba a github (de esta manera no pierdo
# archivos y la gente puede ver, recomendar y replicar mi trabajo)
system("git config --global user.name 'marygpresa'")
system("git config --global user.email 'mariagranadospresa@gmail.com'")

setwd("/Users/mariagranados/TFM")  # Asegura que estás en el directorio correcto
system("git status")  # Comprueba si Git funciona aquí

commit_message <- paste("Actualización:", Sys.time())
system("git add .")  
system(paste("git commit -m '", commit_message, "'", sep=""))  
system("git push origin main")

# Descarga de datos
install.packages("installr")
library(installr)

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager :: install(c("org.Hs.eg.db","edgeR","ComplexHeatmap","SummarizedExperiment","TCGAbiolinks"))

library(limma)
library(edgeR)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(TCGAbiolinks)


# Download data project: TARGET-AML from GDC (Genome Data Common)
## visualizing project
gdcprojects <- getGDCprojects()
getProjectSummary("TARGET-AML")

## creating a query
query_aml <- GDCquery(project = "TARGET-AML",
                      data.category = "Transcriptome Profiling",
                      experimental.strategy = "RNA-Seq",
                      workflow.type = "STAR - Counts",
                      access = "open"
                      )
res <- getResults(query_aml)
View(res)

# download all quey data
GDCdownload(query_aml)
# since data might change daily, I keep record here that the data used for this 
# analysis was downloaded from the GDC 27 March 2025 at 15:23.

# obtain count matrix (gene expresion counts)
aml_data <- GDCprepare(query_aml, summarizedExperiment = TRUE)
mat <- assay(aml_data, "unstranded") #counts matrix

#obtain and visualize clinical data
clinical_data <- GDCquery_clinic("TAGET-AML", type = "clinical")

# obtain metadata from genes (rowData) and clinical data (colData)
gene_metadata <- as.data.frame(rowData(aml_data))
clinical <- as.data.frame(colData(aml_data))

# save all data to avoid having to download it everytime I continue the analysis
# this is important because GDC webpage could crash or not function one day and I would
# not be able to continue the analysis, it also saves time
save(clinical, clinical_data, mat, gene_metadata, aml_data, file = "aml_data.RData")

# load the data when coming back to editing this code
load("aml_data.RData")

