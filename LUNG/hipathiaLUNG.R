# Download data
install.packages("installr")
library(installr)
updateR() # I am using R 4.5.1

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager :: install(c("org.Hs.eg.db","edgeR","ComplexHeatmap","SummarizedExperiment","TCGAbiolinks","hpAnnot","sva"))


library(limma)
library(dplyr)
library(tidyr) # melt and others
library(edgeR)
library(ComplexHeatmap) #heatmaps
library(circlize) #heatmaps
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
library(reshape2) #for melt and plots but shouldn't be necessary since it is masked in tidyr
library(hpAnnot)
library(patchwork) # Arrange side by side plots
library(grid)  # for unit()

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

# Clean the Ensembl ID versions in the matrix (in my case, the gene names in the matrix had EnsemblID.version)
# Having the version may make them incompatible with other annotation paths, so I will remove the versions from the IDs
rownames(mat) <- gsub("\\..*", "", rownames(mat))
rownames(mat)

# Metadata (anotacion ensembl, entrez y hgnc)
genes_info <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(mat),
  mart = mart
)

genes_info
save(genes_info, file = "genes_info.RData")

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
ko_egfr_total <- ko_egfr_total[ko_egfr_total %in% rownames(fc)]
length(ko_egfr_total) # 22 pathways
# This means those paths depend completely on EGFR activity so there will be no bypass through other genes that might
# influence them if EGFR is inhibited, so we are not interested in them and will remove them.
fcoff<-fc[!rownames(fc)%in%ko_egfr_total,]
fcoff<-log2(fcoff) #log2 fold change

#Find any FC>2 on any sample (person) (significant changes)
sig_change <- fcoff[apply(fcoff,1,function(x) any(x>2|x< -2)),,drop=F]
rownames(sig_change)
dim(sig_change) #371 539
View(sig_change) #lots of NAs means not relevant for tumor, they come from ko_vals / tumor_vals where tumor_vals = 0
sig_change <- sig_change[!is.na(sig_change[,1]),] #remove NAs
dim(sig_change) #35 539 means 35 potential bypass pathways

save(fc, ko_egfr_total, fcoff, sig_change, file = "fc.RData")

#plots
# 22 EGFR-dependent pathways in Lung Cancer have FC=0 so can't graph that numerically
df <- data.frame(
  Pathway = ko_egfr_total,
  Status = rep("EGFR-dependent (turned off)", length(ko_egfr_total))
)

# Lollipop-style plot (not really nice since all values are 0)
ggplot(df, aes(x = Status, y = reorder(Pathway, Pathway))) +
  geom_point(color = "tomato", size = 4) +
  geom_segment(aes(x = 0, xend = 1, y = Pathway, yend = Pathway), color = "grey") +
  xlab("") + ylab("") +
  ggtitle("EGFR-dependent Pathways Completely Turned Off") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


#Fold-change distribution of pathways with significant change after inhibiting EGFR in Lung Cancer"
cam_melt <- melt(sig_change)
View(cam_melt)
boxplot_fc <- ggplot(cam_melt, aes(x=Var1, y=value)) + geom_boxplot(width=0.1) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_text(angle = 60,vjust = 1,hjust = 1),
        plot.margin = unit(c(0.2,0.2,0.2,1.2),"cm"),
        plot.title = element_text(hjust = 0.5)) +
  ylab("log(FC)") +
  ggtitle("Fold-change distribution of pathways with significant change after inhibiting EGFR in Lung Cancer")
print(boxplot_fc)


#cap extreme values for visualization only (there are logFCs really high and squeeze the graphs to make them fit)
cam_melt$logFC_capped <- pmin(pmax(cam_melt$value, -10), 10)
# Boxplot using the capped values
ggplot(cam_melt, aes(x = Var1, y = logFC_capped, fill = Var1)) +
  geom_boxplot() +
  coord_flip() +
  theme_minimal() +
  ylab("log2(Fold Change, capped at ±10)") +
  xlab("Pathway")


# to avoid caping values, graph changing scale gradually
bypass_plot <- ggplot(cam_melt, aes(x = Var1, y = value, fill = Var1)) +
                geom_boxplot() +
                coord_flip() +
                scale_y_continuous(trans = "pseudo_log") +
                theme_minimal() +
                xlab("Pathway") +
                ylab("log2(Fold Change)") +
                labs(fill = "Pathways")  # rename legend title
bypass_plot
save(bypass_plot, file = bypass_plot.png)
#---------------------------
# comparison of all groups plots and tables
# Tumor vs Normal wilcoxon
tumor_vs_normal_sig <- comp[comp$p.value < 0.05, ]
tumor_vs_normal_sig <- data.frame(
  Pathway = tumor_vs_normal_sig$path.names,
  p.value = tumor_vs_normal_sig$p.value
)
tumor_vs_normal_sig <- tumor_vs_normal_sig[order(tumor_vs_normal_sig$p.value), ]
View(tumor_vs_normal_sig) #table of significant pathways in tumor vs normal
write.csv(tumor_vs_normal_sig, file = "tumor_vs_normal_sig.csv", row.names = FALSE)


# EGFR KO vs Normal wilcoxon
ko_egfr_sig <- ko_egfr_comp[ko_egfr_comp$p.value < 0.05, ]
ko_egfr_sig <- data.frame(
  Pathway = ko_egfr_sig$ko_egfr.names,
  p.value = ko_egfr_sig$p.value
)
ko_egfr_sig <- ko_egfr_sig[order(ko_egfr_sig$p.value), ]
View(ko_egfr_sig) #table of significant pathways after EGFR KO
write.csv(ko_egfr_sig, file = "ko_egfr_sig.csv", row.names = FALSE)

#table for FC
# Completely turned off pathways
ko_egfr_total_table <- data.frame(
  Pathway = ko_egfr_total,
  Max_abs_FC = apply(fc[ko_egfr_total, , drop=FALSE], 1, function(x) max(abs(x), na.rm=TRUE)),
  Description = "Completely turned off"
)

# Potential bypass pathways
sig_change_table <- data.frame(
  Pathway = rownames(sig_change),
  Max_FC = apply(sig_change, 1, function(x) max(abs(x), na.rm=TRUE)),
  Description = "Potential bypass"
)

# Combine both tables
all_fc_table <- rbind(ko_egfr_total_table, sig_change_table)
View(sig_change_table)
# Sort by Max_FC descending
all_fc_table <- all_fc_table[order(-all_fc_table$Max_FC), ]
write.csv(all_fc_table, file = "all_fc_table.csv", row.names = FALSE)

# View table
View(all_fc_table) #table with relevant pathways
dim(all_fc_table) #57 paths
dim(tumor_vs_normal_sig) #1516 paths
table(all_fc_table%in% tumor_vs_normal_sig$Pathway) #54 TRUE 3 FALSE
# this means all but 3 (3 false) significant pathways according to FC after KO, were also significant pathways vs normal
setdiff(all_fc_table$Pathway, tumor_vs_normal_sig$Pathway)
#FoxO signaling pathway
#ErbB signaling pathway
#Pancreatic cancer

#--------------------------------------------------------------------------
#PLOTS of expression matrix before vs after KO
#--------------------------------------------------------------------------
#Pathways for each group
turned_off_paths <- ko_egfr_total_table$Pathway
changed_paths   <- sig_change_table$Pathway

# Make sure they exist in both expression matrices
turned_off_paths <- intersect(turned_off_paths, rownames(path_vals))
turned_off_paths <- intersect(turned_off_paths, rownames(ko_egfr_vals))

changed_paths <- intersect(changed_paths, rownames(path_vals))
changed_paths <- intersect(changed_paths, rownames(ko_egfr_vals))

# Extract expression matrices for tumor samples
expr_before_off <- path_vals[turned_off_paths, sample.group == "Primary Tumor", drop=FALSE]
expr_after_off  <- ko_egfr_vals[turned_off_paths, sample.group == "Primary Tumor", drop=FALSE]

expr_before_changed <- path_vals[changed_paths, sample.group == "Primary Tumor", drop=FALSE]
expr_after_changed  <- ko_egfr_vals[changed_paths, sample.group == "Primary Tumor", drop=FALSE]

# Remove pathways where expression is identical before and after
same_expr_off <- apply(expr_before_off == expr_after_off, 1, all, na.rm=TRUE)
expr_before_off <- expr_before_off[!same_expr_off, , drop=FALSE]
expr_after_off  <- expr_after_off[!same_expr_off, , drop=FALSE]

same_expr_changed <- apply(expr_before_changed == expr_after_changed, 1, all, na.rm=TRUE)
expr_before_changed <- expr_before_changed[!same_expr_changed, , drop=FALSE]
expr_after_changed  <- expr_after_changed[!same_expr_changed, , drop=FALSE]

# ---- Turned-off pathways ----
expr_before_long_off <- melt(expr_before_off, varnames = c("Pathway","Sample"), value.name = "Expression")
expr_before_long_off$Condition <- "Before"

expr_after_long_off <- melt(expr_after_off, varnames = c("Pathway","Sample"), value.name = "Expression")
expr_after_long_off$Condition <- "After"

expr_long_off <- rbind(expr_before_long_off, expr_after_long_off)

expr_long_off$Condition <- factor(expr_long_off$Condition, levels = c("Before", "After"))

egfr_dependant_plot <- ggplot(expr_long_off, aes(x = Pathway, y = Expression, fill = Condition)) +
                      geom_boxplot(position = position_dodge(width = 0.8)) +
                      theme_bw() +
                      theme(
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        plot.margin = margin(10, 15, 50, 40, "pt")  # top, right, bottom, left
                        ) +
                      labs(
                        title = "EGFR-dependent Pathways",
                        y = "Pathway Activity",
                        x = "Pathway"
                        )


egfr_dependant_plot
save(egfr_dependant_plot, file = "egfr_dependant_plot.png")


# ---- Significantly changed pathways ----
##pending,not completely necessary to visualize expression range is too wide so conclusion can be taken from fc plot

install.packages("renv")
renv::init()

