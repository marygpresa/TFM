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

BiocManager :: install(c("org.Hs.eg.db","edgeR","ComplexHeatmap"))

library(limma)
library(edgeR)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(TCGAbiolinks)
