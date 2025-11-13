install.packages("Seurat")
install.packages("viridis")
install.packages("scCustomize")

library(Seurat)
library(viridis)
library(scCustomize)

# Path to data
data_dir <- "~/Documents/Comp Bio Project/Data/Drosophila brain /"

# Load dataset --> need to have exact seurat names!
L1_brain_drosophila <- Read10X(data.dir = data_dir)

# Initialize the Seurat object with the raw (non-normalized data)
L1_brain_drosophila <- CreateSeuratObject(counts = L1_brain_drosophila, project = "brain drosophila", min.cells = 3, min.features = 200)

# QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

# Mitochondrial contamination
L1_brain_drosophila[["percent.mt"]] <- PercentageFeatureSet(L1_brain_drosophila, pattern = "^mt:")

# Q/C filter 
VlnPlot(L1_brain_drosophila, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
L1_brain_drosophila <- subset(L1_brain_drosophila, subset = nFeature_RNA > 1500 & nFeature_RNA < 4000 & percent.mt < 10)
VlnPlot(L1_brain_drosophila, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalizing data 
L1_brain_drosophila <- NormalizeData(L1_brain_drosophila, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features 
L1_brain_drosophila <- FindVariableFeatures(L1_brain_drosophila, selection.method = "vst", nfeatures = 2000)

# Scaling the data 
all.genes <- rownames(L1_brain_drosophila)
L1_brain_drosophila<- ScaleData(L1_brain_drosophila, features = all.genes)

# First principal component --> list of genes with the most positive and negative loadings (= correlation or anti-correlation) 
L1_brain_drosophila <- RunPCA(L1_brain_drosophila, features = VariableFeatures(object = L1_brain_drosophila))

# Find min PC possible to include variability over sample
stdv <- L1_brain_drosophila[["pca"]]@stdev
percent_stdv <- (stdv/sum(stdv)) * 100
cumulative <- cumsum(percent_stdv)
co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                     percent_stdv[2:length(percent_stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min_pc <- min(co1, co2)
print(min_pc) ß
# --> 11

# Determine dimensionality:
#‘Elbow plot’ = ranking of principle components based on the percentage of variance explained by each one
ElbowPlot(L1_brain_drosophila, ndims = 100)
# --> Here I'm going to choose 31 (I would've probably just chose 30 but they used 31 in the paper)

# Cluster the cells:
# FindNeighbors() function --> construct a KNN graph based on the euclidean distance in PCA space & refine the edge weights between any two cells based on the shared overlap in their local neighborhoods 
L1_brain_drosophila <- FindNeighbors(L1_brain_drosophila, dims = 1:31)
# FindClusters() function  --> Louvain algorithm (default) to group cells together, resolution parameter of 0.4-1.2 typically returns good results for single-cell datasets (increased values leading to a greater number of clusters)
L1_brain_drosophila <- FindClusters(L1_brain_drosophila, resolution = 2)
# First UMAP
L1_brain_drosophila<- RunUMAP(L1_brain_drosophila, dims = 1:31)

# Determining clusters of low quality --> nFeature /cluster 
L1_brain_drosophila.data <- L1_brain_drosophila@meta.data
mean_nFeature <- mean(L1_brain_drosophila.data$nFeature_RNA, na.rm = TRUE)
sd_nFeature <- sd(L1_brain_drosophila.data$nFeature_RNA, na.rm = TRUE)
clusterfilter <- mean_nFeature - (1.5*(sd_nFeature))
print(clusterfilter)
mean_nFeature_clusters <- L1_brain_drosophila.data%>%
  group_by(seurat_clusters) %>%
  summarize(mean_nFeature_cluster = mean(nFeature_RNA, na.rm = TRUE))
print(mean_nFeature_clusters, n = 100)
# clusterfilter = 1659.27
# No clusters of significant low quality! 

# Look at generated clusters
DimPlot(L1_brain_drosophila, reduction = "umap")

# Highlight individual cluster
cluster_to_highlight <- c("18", "19")

DimPlot(
  L1_brain_drosophila,
  reduction = "umap",
  group.by = "seurat_clusters",
  cols = ifelse(
    levels(L1_brain_drosophila$seurat_clusters) == cluster_to_highlight,
    "red",       # color for highlighted cluster
    "lightgrey"  # color for other clusters
  )
)

cluster_to_highlight <- "17"

DimPlot(
  L1_brain_drosophila,
  reduction = "umap",
  group.by = "seurat_clusters",
  cols = ifelse(
    levels(L1_brain_drosophila$seurat_clusters) == cluster_to_highlight,
    "red",       # color for highlighted cluster
    "lightgrey"  # color for other clusters
  )
)

# Determine identity of clusters:
pal <- viridis(10, option = "turbo") # I like this color way we can change this later though!

# Big picture --> neurons, neural progenitors, glial cells, undifferentiated neurons 
# Neurons = characterized by elav, & ChAT expression 
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "elav", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "ChAT", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
cluster16.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 16) # ChAT !
cluster9.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 9) # other neural markers
cluster4.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 4) # elav & GABA neuron markers!!!
cluster13.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 13) #other neural markers
cluster17.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 17) #elav!
# Clusters 0, 1, 2, 3, 4, 5, 6, 7, 9, 12, 13, 16, & 17 

# Neural progenitors = characterized by N, & dpn expression 
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "N", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "dpn", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
# CLusters 8, 10, 14, 15

# Glia = characterized by repo expression (+ loco and vir-1, the paper used repo but i wantred to confirm bc it's kind of lowly expressed)
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "repo", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "loco", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "vir-1", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
cluster18.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 18) #repo included and other markers included
cluster19.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 19) #repo included and other markers included
cluster11.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 11) #repo and other markers included
# Clusters 11, 18 ,& 19 

Idents(L1_brain_drosophila) <- "seurat_clusters" 
L1_brain_drosophila<-RenameIdents(L1_brain_drosophila, '11' = "Glial cells",'18' = "Glial cells", '19' = "Glial cells")
L1_brain_drosophila<-RenameIdents(L1_brain_drosophila, '0' = "Neurons", '1' = "Neurons", '2' = "Neurons", '3' = "Neurons", '4' = "Neurons", '5' = "Neurons", '6' = "Neurons", '7' = "Neurons", '9' = "Neurons", '12' = "Neurons", '13' = "Neurons", '16' = "Neurons", '17' = "Neurons")
L1_brain_drosophila<-RenameIdents(L1_brain_drosophila, '8' = "Neural Progenitors",'10' = "Neural Progenitors",'14' = "Neural Progenitors", '15' = "Neural Progenitors")

L1_brain_drosophila[["assigned.ident"]] <- Idents(object = L1_brain_drosophila)
DimPlot(L1_brain_drosophila, reduction = "umap", group.by = c("assigned.ident"))

#Our genes: 
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Nrx-1", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "rg", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "kis", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "vnd", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "SK", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "14-3-3ζ", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")

