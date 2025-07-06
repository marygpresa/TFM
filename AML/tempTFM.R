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


# Paquetes
install.packages("installr")
library(installr)
installr::updateR() # actualizar R 4.5

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")

BiocManager :: install(c("org.Hs.eg.db","edgeR","ComplexHeatmap","SummarizedExperiment","TCGAbiolinks","biomaRt"))


# librerias
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

# load the data when coming back to editing this code
#load("aml_data.RData")

### STEP 1 DOWNLOAD DATA 
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
dim(res)

# download all query data
GDCdownload(query_aml)
# since data might change daily, I keep record here that the data TCGA-LAML used for this 
# analysis was downloaded from the GDC webpage on 12 May 2025 17:53
# ***esto era TARGET AML 27 March 2025 at 15:23.

# obtain count matrix (gene expresion counts)
#aml_data <- GDCprepare(query_aml, summarizedExperiment = TRUE)
#mat <- assay(aml_data, "unstranded") #counts matrix

#obtain and visualize clinical data
clinical_data <- GDCquery_clinic("TARGET-AML")
View(clinical_data)
colnames(clinical_data)

View(mat)
dim(mat)
# por ser sangre no hay muestreo pareado de tejido normal como en los tumores solidos
# no se puede hacer tumor vs normal
# alternativas:
#progression of recurense
#age at diagnosis
# vital status
# ctrl normal externo y lo normalizo todo junto?:GTEx
# sin ctrl y luego solo hago estadistica considerando alguna variable clinica

# Crear query para ese archivo
# Extraer los file_id disponibles
files <- query_aml$results[[1]]$file_id

# Revisar cuántos hay disponibles
length(files)

# Seleccionar aleatoriamente 1000 sin reemplazo
set.seed(123)  # para reproducibilidad
chunk_files <- sample(files, 1000, replace = FALSE)

# Crear una copia del query original
query_chunk <- query_aml

# Filtrar los resultados para que solo contenga los file_id seleccionados
query_chunk$results[[1]] <- query_chunk$results[[1]][query_chunk$results[[1]]$file_id %in% chunk_files, ]

# Verificar que ahora tenga solo 1000
nrow(query_chunk$results[[1]])  # debe dar 1000

#Preparar los datos descargados
aml_chunk <- GDCprepare(query_chunk, summarizedExperiment = TRUE)
mat <- assay(aml_chunk, "unstranded")
View(mat)








# obtain metadata from genes (rowData) and clinical data (colData)
gene_metadata <- as.data.frame(rowData(aml_data))
clinical <- as.data.frame(colData(aml_data))

# save all data to avoid having to download it everytime I continue the analysis
# this is important because GDC webpage could crash or not function one day and I would
# not be able to continue the analysis, it also saves time
save(clinical, clinical_data, mat, gene_metadata, aml_data, file = "aml_data.RData")

#agregar cuando logre correrlo  mat, gene_metadata, clinical

# load the data when coming back to editing this code
load("aml_data.RData")


# obtain clinical data
clinical_data <- GDC
clinical <- as.data.frame(colData(aml_data))


# loop para descargar data de TARGET-AML 3227 files
# Cambio la funcion de TCGA biolinks porque mi data no tiene una columna que esta intentando descargar
assignInNamespace(
  "colDataPrepare",
  function(cases) {
    library(dplyr)
    
    # Si no es data.frame, devolver tal cual
    if(!is.data.frame(cases)) {
      warning("cases no es un data.frame, se regresa tal cual.")
      return(cases)
    }
    
    # Normalizar nombres de columnas
    colnames(cases) <- tolower(colnames(cases))
    
    # Si no existe disease_response, agregarla con NA
    if (!"disease_response" %in% colnames(cases)) {
      cases$disease_response <- NA
    }
    
    # Preparar colData limpio y sin errores
    cases <- cases %>%
      dplyr::rename_at(vars(contains("last_follow_up")), ~"days_to_follow_up") %>%
      dplyr::group_by(submitter_id) %>%
      dplyr::filter(!is.na(submitter_id), !is.na(days_to_follow_up)) %>%
      dplyr::filter(dplyr::row_number() == which.max(days_to_follow_up)) %>%
      dplyr::ungroup() %>%
      dplyr::select(any_of(c("submitter_id", "days_to_follow_up", "disease_response")))
    
    return(cases)
  },
  ns = "TCGAbiolinks"
)

# Lista de todos los IDs
all_ids <- res$id

# Inicializar matriz de conteos y la de clinical
mat <- NULL
clinical <- data.frame()

# Definir tamaño de sub-chunk
chunk_size <- 100

# Calcular número total de chunks
n_chunks <- ceiling(length(all_ids) / chunk_size)

# Loop por cada file_id
for (file_id in res$id) {
  
  # Crear query para ese archivo
  query_chunk <- query_aml
  query_chunk$results[[1]] <- query_chunk$results[[1]][query_chunk$results[[1]]$file_id %in% file_id, ]
  
  # Descargar datos
  GDCdownload(query_chunk)
  
  # Preparar datos
  aml_chunk <- GDCprepare(query_chunk, summarizedExperiment = TRUE)
  
  # Extraer counts
  counts <- assay(aml_chunk, "unstranded")
  
  # Asignar nombres de muestra
  sample_id <- colData(aml_chunk)$submitter_id
  colnames(counts) <- sample_id
  
  # Combinar con la matriz general
  if (is.null(mat)) {
    mat <- counts
  } else {
    mat <- cbind(mat, counts)
  }
  
  # Extraer clinical de ese chunk
  chunk_clinical <- as.data.frame(colData(aml_chunk))
  clinical <- rbind(clinical, chunk_clinical)
}


# Revisar las dimensiones de la matriz final
dim(mat)
View(mat)
colnames(mat)

dim(clinical)
View(clinical)
save(mat,clinical, file = "mat.RData")
load("mat.RData")

# obtain Gene metadata 
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
