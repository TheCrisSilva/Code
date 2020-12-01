### Loading libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(tibble)
###loading data from the count matrix and making the seurat object

matrix.data<- GSE138852_counts_csv
df<- data.frame(GSE138852_counts_csv)
rownames(df) <- df[,1] ### first colum was assign as row names of the dataset

names(df)

df[,1] <- NULL # first column was removed from dataset because it has 13215 columns but it should only have 13214, 
#in order to compatible with feauture data
names(df)
###features data importing
feature.data <- GSE138852_covariates
head(feature.data)

####  setting up metadata and seurat object
memory.limit(size= 56000)
Grubman<- CreateSeuratObject(counts = df, project = "Grubman", min.cells = 3, min.features = 200) 
###min.cell and min.feature ???
#this is the setting from the tutorial


Grubman <- CreateSeuratObject(counts = df, project = "Grubman", meta.data = feature.data)
metadata <- feature.data$oupSample.batchCond
names(metadata) <- colnames(Grubman)
Grubman <- AddMetaData(Grubman, metadata, col.name = 'Condition')
metadata2 <- feature.data$oupSample.cellType
names(metadata2) <- colnames(Grubman)
Grubman <- AddMetaData(Grubman, metadata2, col.name = 'Cell_type')





###standard preprocess workflow

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Grubman[["percent.mt"]] <- PercentageFeatureSet(Grubman, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Grubman, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(Grubman, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Grubman, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

###this data set was downloaded after the QC therefore I'm not going to do anything

###Identification of highly variable features (feature selection)

Grubman <- FindVariableFeatures(Grubman, selection.method = "vst", nfeatures = 2500)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Grubman), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Grubman)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

top10
write.csv(top10, file = "top10markers.csv")

###Scaling the data
#paper says scaling was done for the dataset however without scaling we can't run PCA

all.genes <- rownames(Grubman)
Grubman <- ScaleData(Grubman, features = all.genes)



###Perform linear dimension reduction
Grubman <- RunPCA(Grubman, features = VariableFeatures(object = Grubman))
DimHeatmap(Grubman, dims = 1:10, cells = 13214, balanced = TRUE)

###Determine dimensionality of dataset

ElbowPlot(Grubman)

###Cluster the cells

Grubman <- FindNeighbors(Grubman, dims = 1:10)
Grubman <- FindClusters(Grubman, resolution = 0.5)
head(Idents(Grubman), 5)

#UMAP
Grubman <- RunUMAP(object = Grubman, dims = 1:10)
# Plot results
DimPlot(object = Grubman, reduction = "umap", label = TRUE)

DimPlot(Grubman, reduction = "umap", group.by = "Condition", label = TRUE)
DimPlot(Grubman, reduction = "umap", group.by = "Cell_type", label = TRUE)

saveRDS(Grubman, file = "E://All Research projects//Alzeihmer research//Progress//11282020//Grubman_10pc.rds")

Grubman <- readRDS(file="E://All Research projects//Alzeihmer research//Progress//11302020//Grubman_10pc.rds")
################################################################################################################
###Subsetting Oligocytes and DEG


Oligocytes <- subset(Grubman, idents = c("0","1","2","3","6","8"))

#Identification of highly variable genes
Oligocytes <- FindVariableFeatures(Oligocytes, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes in Astrocytes
top10Oligocytes <- head(VariableFeatures(Oligocytes), 10)
top10Oligocytes
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Oligocytes)
plot2 <- LabelPoints(plot = plot1, points = top10Oligocytes, repel = TRUE)
plot1 + plot2

write.csv(top10Oligocytes, file = "top10Oligocytevariablegenes.csv")

#Scale data
all.genes <- rownames(Oligocytes)
Oligocytes <- ScaleData(Oligocytes, features = all.genes)

#Linear dimension reduction
Oligocytes <- RunPCA(Oligocytes, features = VariableFeatures(object = Oligocytes))
DimHeatmap(Oligocytes, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(Oligocytes)

#clustering
Oligocytes <- FindNeighbors(Oligocytes, dims = 1:4)
Oligocytes <- FindClusters(Oligocytes, resolution = 0.5)

Oligocytes <- RunUMAP(Oligocytes, dims = 1:4)

DimPlot(object = Oligocytes, reduction = "umap", label = TRUE)

DimPlot(Oligocytes, reduction = "umap", group.by = "Condition", label = TRUE)
DimPlot(Oligocytes, reduction = "umap", group.by = "Cell_type", label = TRUE) ###some doublets and unId s are here



#finding all DEG between disease and healthy brain Oligocytes

Oligocytes_AD_Ct_markers <- FindMarkers(Oligocytes, ident.1 = c(4,5, 7,3,1), ident.2 = c(0,2,6) , min.pct = 0.25, loffc.threshold= log(2))
head(Oligocytes_AD_Ct_markers, n = 5)
write.csv (Oligocytes_AD_Ct_markers, file = "all_DEG_in_Oligocytes.csv")

top20upregulated_AD_Oligo_markers <- Oligocytes_AD_Ct_markers%>% top_n(n = 20, wt = avg_logFC)
write.csv(top20upregulated_AD_Oligo_markers, file = "top20upregulatedgenes_AD_Oligocytes.csv")
#############################################################################################################
###Subsetting Astrocystes and DEG


Astrocytes <- subset(Grubman, idents = c("10","7","4"))

#Identification of highly variable genes
Astrocytes <- FindVariableFeatures(Astrocytes, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes in Astrocytes
top10Astrocytes <- head(VariableFeatures(Astrocytes), 10)
top10Astrocytes
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Astrocytes)
plot2 <- LabelPoints(plot = plot1, points = top10Astrocytes, repel = TRUE)
plot1 + plot2

write.csv(top10Astrocytes, file = "top10Astrocytevariablegenes.csv")

#Scale data
all.genes <- rownames(Astrocytes)
Astrocytes <- ScaleData(Astrocytes, features = all.genes)

#Linear dimension reduction
Astrocytes <- RunPCA(Astrocytes, features = VariableFeatures(object = Astrocytes))
DimHeatmap(Astrocytes, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(Astrocytes)


Astrocytes <- FindNeighbors(Astrocytes, dims = 1:9)
Astrocytes <- FindClusters(Astrocytes, resolution = 0.5)

Astrocytes <- RunUMAP(Astrocytes, dims = 1:9)

DimPlot(object = Astrocytes, reduction = "umap", label = TRUE)

DimPlot(Astrocytes, reduction = "umap", group.by = "Condition", label = TRUE)
DimPlot(Astrocytes, reduction = "umap", group.by = "Cell_type", label = TRUE) ###some doublets and unId s are here



#finding all DEG between disease and healthy brain astrocytes

AD_Astrocytes_upregulatedgenes <- FindMarkers(Astrocytes, ident.1 = c(3,7, 5,8), ident.2 = c(0,1,2,4,6) , min.pct = 0.25, loffc.threshold= log(2))
head(AD_Astrocytes_upregulatedgenes, n = 5)
write.csv (Astrocytes_AD_Ct_markers, file = "all_DEG_in_astrocytes.csv")

top20upregulated_AD_Astro_markers <- Astrocytes_AD_Ct_markers%>% top_n(n = 20, wt = avg_logFC)
write.csv(top20upregulated_AD_Astro_markers, file = "top20upregulatedgenes_AD_Astrocytes.csv")
##############################################################################################################
###Subsetting neurones

Neuron <- subset(Grubman, idents = c("9","17","12"))

#Identification of highly variable genes
Neuron <- FindVariableFeatures(Neuron, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes in Astrocytes
top10Neuron <- head(VariableFeatures(Neuron), 10)
top10Neuron
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Neuron)
plot2 <- LabelPoints(plot = plot1, points = top10Neuron, repel = TRUE)
plot1 + plot2

write.csv(top10Neuron, file = "top10Neuronvariablegenes.csv")

#Scale data
all.genes <- rownames(Neuron)
Neuron <- ScaleData(Neuron, features = all.genes)

#Linear dimension reduction
Neuron <- RunPCA(Neuron, features = VariableFeatures(object = Neuron))
DimHeatmap(Neuron, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(Neuron)


Neuron <- FindNeighbors(Neuron, dims = 1:7)
Neuron <- FindClusters(Neuron, resolution = 0.5)

Neuron <- RunUMAP(Neuron, dims = 1:7)

DimPlot(object = Neuron, reduction = "umap", label = TRUE)

DimPlot(Neuron, reduction = "umap", group.by = "Condition", label = TRUE)
DimPlot(Neuron, reduction = "umap", group.by = "Cell_type", label = TRUE) ###some doublets and unId s are here



#finding all DEG between disease and healthy brain Neuron

Neuron_AD_Ct_markers <- FindMarkers(Neuron, ident.1 = c(0,3,7,2), ident.2 = c(5,8,9,4,6,1) , min.pct = 0.25, loffc.threshold= log(2))
head(Neuron_AD_Ct_markers, n = 5)
write.csv (Neuron_AD_Ct_markers, file = "all_DEG_in_Neuron.csv")

top20upregulated_AD_Neuron_markers <- Neuron_AD_Ct_markers%>% top_n(n = 20, wt = avg_logFC)
write.csv(top20upregulated_AD_Neuron_markers, file = "top20upregulatedgenes_AD_Neuron.csv")

######################################################################################################

###Subsetting OPC

OPC <- subset(Grubman, idents = c("5","14"))

#Identification of highly variable genes
OPC <- FindVariableFeatures(OPC, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes in Astrocytes
top10OPC <- head(VariableFeatures(OPC), 10)
top10OPC
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(OPC)
plot2 <- LabelPoints(plot = plot1, points = top10OPC, repel = TRUE)
plot1 + plot2

write.csv(top10OPC, file = "top10OPCvariablegenes.csv")

#Scale data
all.genes <- rownames(OPC)
OPC <- ScaleData(OPC, features = all.genes)

#Linear dimension reduction
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC))
DimHeatmap(OPC, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(OPC)


OPC <- FindNeighbors(OPC, dims = 1:6)
OPC <- FindClusters(OPC, resolution = 0.5)

OPC <- RunUMAP(OPC, dims = 1:6)

DimPlot(object = OPC, reduction = "umap", label = TRUE)

DimPlot(OPC, reduction = "umap", group.by = "Condition", label = TRUE)
DimPlot(OPC, reduction = "umap", group.by = "Cell_type", label = TRUE) ###some doublets and unId s are here



#finding all DEG between disease and healthy brain Neuron

OPC_AD_Ct_markers <- FindMarkers(OPC, ident.1 = c(5,3,4), ident.2 = c(0,2,1,6) , min.pct = 0.25, loffc.threshold= log(2))
head(OPC_AD_Ct_markers, n = 5)
write.csv (OPC_AD_Ct_markers, file = "all_DEG_in_OPC.csv")

top20upregulated_AD_OPC_markers <- OPC_AD_Ct_markers%>% top_n(n = 20, wt = avg_logFC)
write.csv(top20upregulated_AD_OPC_markers, file = "top20upregulatedgenes_AD_OPC.csv")
###################################################################################################################

#Microglia?
