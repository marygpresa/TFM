# Descarga de datos
install.packages("installr")
library(installr)
updateR() # I am using R 4.5.1

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager :: install(c("org.Hs.eg.db","edgeR","ComplexHeatmap","SummarizedExperiment","TCGAbiolinks","hpAnnot"))


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
library(sva) #for combat, correcting batch effect
library(hipathia)
library(reshape2) #plots but shouldn't be necessary since it is masked in tidyr
library(hpAnnot)

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
unique(clinical$primary_site)
table(clinical$sample_type) 

# Find duplicates
table(duplicated(rownames(mat))) #  no duplicated genes (FALSE)
table(duplicated(colnames(mat))) # no duplicated samples (FALSE)
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


# Log2 CPM-transformed normalized expression
tmmexp <- cpm(de.exprs_norm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 3)
hist(tmmexp, main = "TMM normalized log-CPM")

# Log2 CPM-transformed unnormalized expression
no_tmmexp <- cpm(de.exp, log = TRUE)
hist(no_tmmexp, main = "Raw log-CPM")

save(luad_data,mat,clinical,de.exp,de.exprs_norm,tmmexp,no_tmmexp,file = "mats.RData")
load("mats.RData")


#PCA to check batch effects
pca<-prcomp(t(tmmexp))
plot(pca$x[,1], pca$x[,2], col = (clinical$sample_type == "Primary Tumor") + 1)
legend("topright", legend = c("Normal", "Tumor"), col = c(1, 2), pch = 1)
title("Before ComBat")

colnames(clinical)

# ComBat for batch effect correction?
table(clinical$is_ffpe) 
table(clinical$preservation_method)
#combat_data <- ComBat(dat = tmmexp, batch = clinical$ , par.prior = TRUE)


#----------------Hipathia calculations------------------------------

ls("package:hipathia") # view available functions in hipathia package
pathways <- load_pathways("hsa")
hipathia <- translate_data(tmmexp,"hsa") ##get entrez ids from matrix
hipathia <- normalize_data(hipathia) #normalize with hipathia
results <- hipathia(hipathia, pathways, verbose=TRUE) #perform hipathia analysis

# Extract the matrix of pathway activity scores
#path_vals <- get_paths_matrix(results) function get_pathways matrix doesn't exist anymore
path_vals <- assays(results)$paths
path_vals
path_vals <- normalize_paths(path_vals, pathways) 

#Group by sample type (normal vs tumor)
sample.group <-clinical$sample_type 
table(sample.group) # 539 Primary Tumor, 2 Recurrent Tumor, 59 Solid Tissue Normal)

# Perform Wilcoxon test for each pathway
comp <- do_wilcoxon(path_vals, sample.group, g1 = "Primary Tumor", g2 = "Solid Tissue Normal")
path.names <- get_path_names(pathways, rownames(comp))
comp <- cbind(path.names, comp)
table(comp$p.value<0.05) #1516 TRUE means 1516 pathways show significant differential activity (p < 0.05) between tumor and normal

save(pca, pathways,hipathia,path_vals,comp, file = "hipathia.RData")

#---------------KNOCK OUTS------------------------------
#Approach 1: EGFR (Entrez: 1956)
#Approach 2: group of genes?

#EGFR KO (1956)
hist(hipathia["1956",],100)
ko_egfr<-hipathia
ko_egfr["1956",sample.group=="Primary Tumor"]<-0  #KO
ko_egfr_res <- hipathia(ko_egfr, pathways, verbose=TRUE)
ko_egfr_vals <- assays(ko_egfr_res)$paths
ko_egfr_vals<-normalize_paths(ko_egfr_vals, pathways)

# wilcoxon with KO
ko_egfr_comp <- do_wilcoxon(ko_egfr_vals, sample.group, g1 = "Primary Tumor", g2 = "Solid Tissue Normal")
ko_egfr.names <- get_path_names(pathways, rownames(ko_egfr_comp))
ko_egfr_comp <- cbind(ko_egfr.names, ko_egfr_comp)
table(ko_egfr_comp$p.value<0.05) # 1553 TRUE means 1553 pathways show significant differential activity (p < 0.05) between tumor EGFR KO and normal

save(ko_egfr, ko_egfr_res, ko_egfr_vals, ko_egfr_comp, file = "ko_egfr.RData")


#---------------------PAIR-FOLD CHANGE ANALYSIS------------------------------
# fold-change base matrix
fc <- path_vals[,sample.group=="Primary Tumor"]
fc[,]<-0 #set to 0

#calculate fold change on tumor samples KOtumor/tumor
for(samp in colnames(fc)){
  fc[,samp]<-ko_egfr_vals[,samp]/path_vals[,samp]
}

dim(fc) #1876 539

#Change path names
rownames(fc)<-get_path_names(pathways, rownames(fc))
View(fc)

#Remove paths that don't change significally
no_change <- rownames(fc[apply(fc,1,function(x) all(x==1)),])
length(no_change) # 1140 out of the 1876 paths don't change (x==1)

#remove them
fc<-fc[!rownames(fc)%in%no_change,]
dim(fc) #941 539
#We have now 941 pathways with fc != 1
#Nevertheless 1876-1140 = 736 so the are extra pathways

# lets check again:
no_change2 <- rownames(fc[apply(fc,1,function(x) all(x==1)),])
length(no_change2) # 205 paths still there

#let's clean again
fc<-fc[!rownames(fc)%in%no_change2,]
dim(fc) #941 539, not working, they are probably NAs but will continue

# Get paths with max fc (the ones that should've completely turn off with the KO, x==0)
ko_egfr_total <- rownames(fc[apply(fc,1,function(x) all(x==0)),])
length(ko_egfr_total) # 23 pathways
# This means those paths depend completely on EGFR activity so there will be no bypass through other genes that might
# influence them if EGFR is inhibited, so we are not interested in them and will remove them.
fc<-fc[!rownames(fc)%in%ko_egfr_total,]
fc<-log2(fc)

#Find any FC>2 on any sample (person) (significant changes)
sig_change <- fc[apply(fc,1,function(x) any(x>2|x< -2)),,drop=F]
rownames(sig_change)
dim(sig_change) #371 539

#plot
cam_melt <- melt(sig_change)

violinplot_fc <- ggplot(cam_melt, aes(x=Var1, y=value)) + geom_violin(width=0.1) + 
                        theme(axis.title.x=element_blank(), 
                        axis.text.x=element_text(angle = 60,vjust = 1,hjust = 1),
                        plot.margin = unit(c(0.2,0.2,0.2,1.2),"cm"),
                        plot.title = element_text(hjust = 0.5)) +
                        ylab("log(FC)") +
                        ggtitle("Fold-change distribution of pathways affected in Lung Cancer")

print(violinplot_fc)

boxplot_fc <- ggplot(cam_melt, aes(x=Var1, y=value)) + geom_boxplot(width=0.1) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_text(angle = 60,vjust = 1,hjust = 1),
        plot.margin = unit(c(0.2,0.2,0.2,1.2),"cm"),
        plot.title = element_text(hjust = 0.5)) +
  ylab("log(FC)") +
  ggtitle("Fold-change distribution of pathways affected in Lung Cancer")

print(boxplot_fc)


#----------------------PATHWAY ANNOTATION------------------------------

ls("package:hpAnnot")
data("annofuns_GO_hsa", package = "hpAnnot")

