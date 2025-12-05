install.packages("Seurat")
install.packages("viridis")
install.packages("scCustomize")

library(Seurat)
library(viridis)
library(scCustomize)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

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
print(min_pc) 
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
mean_nFeature_clusters <- L1_brain_drosophila.data %>%
  group_by(seurat_clusters) %>%
  summarize(mean_nFeature_clusters = mean(nFeature_RNA, na.rm = TRUE))
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

cluster_to_highlight <- "16"

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

# Labels to make --> neurons (+ subtypes), neural progenitors, glial cells, undifferentiated neurons 
# Neurons = characterized by elav, & ChAT expression 
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "elav", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "ChAT", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
# Clusters 0, 1, 2, 3, 4, 5, 6, 7, 9, 12, 13, 16, & 17 

# Monoaminergic neurons = characterized by Vmat expression 
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Vmat", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
cluster17.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 17) 
# Cluster 17 

# GABAergic neurons = characterized by Gad1 expression 
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Gad1", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
cluster4.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 4) # GABA neuron markers!!!
# Cluster 1, and 4

# Glutamatergic neurons = characterized by VGlut expression 
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "VGlut", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
cluster9.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 9)
# Cluster 5, 9, and  0 

# Cholinergic neurons = characterized by VAChT expression 
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "VAChT", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
cluster16.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 16) # VAChT!
# Cluster 7, 3, 2, 6, 16, & 12 

# Peptidergic neurons = characterized by VAChT expression 
FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "twit", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 0.5, reduction = "umap")
cluster13.markers <- FindMarkers(L1_brain_drosophila, only.pos = TRUE, ident.1 = 13) 
# Cluster 13

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
L1_brain_drosophila<-RenameIdents(L1_brain_drosophila, '13' = "Peptidergic neurons")
L1_brain_drosophila<-RenameIdents(L1_brain_drosophila, '7' = "Cholinergic neurons", '3' = "Cholinergic neurons", '2' = "Cholinergic neurons", '6' = "Cholinergic neurons", '16' = "Cholinergic neurons", '12' = "Cholinergic neurons")
L1_brain_drosophila<-RenameIdents(L1_brain_drosophila, '5' = "Glutamatergic neurons", '9' = "Glutamatergic neurons", '0' = "Glutamatergic neurons")
L1_brain_drosophila<-RenameIdents(L1_brain_drosophila, '1' = "GABAergic neurons", '4' = "GABAergic neurons")
L1_brain_drosophila<-RenameIdents(L1_brain_drosophila, '17' = "Monoaminergic neurons")
L1_brain_drosophila<-RenameIdents(L1_brain_drosophila, '8' = "Neural progenitors",'10' = "Neural progenitors",'14' = "Neural progenitors", '15' = "Neural progenitors")
L1_brain_drosophila[["assigned.ident"]] <- Idents(object = L1_brain_drosophila)

# Save Object 
saveRDS(L1_brain_drosophila, file = "/Users/cmdb/Documents/Comp Bio Project/L1_brain_drosophila")

#UMAP for presentation 
my_colors <- brewer.pal(n = 8, name = "Set2")  
DimPlot(
  L1_brain_drosophila,
  reduction = "umap",
  group.by = "assigned.ident",
  pt.size = 2.5,
  cols = my_colors
)  + 
  theme(
    legend.text = element_text(size = 17),   
  )


#Related genes:

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "vnd", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "vnd")
#Uniformly lowly expressed

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "SK", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "SK")
#Uniformly expressed in neurons 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Abl", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Abl")
#Uniformly expressed

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "CG8665", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "CG8665")
#Uniformly not expressed 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "CG9650", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "CG9650")
#Highly expressed in cholinergic, glutamatergic, & GABAergic neurons 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Rep", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Rep")
#Uniformly lowly expressed

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Chc", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Chc")
#Highly expressed in cholinergic, glutamatergic, GABAergic neurons, and Glial cells 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Chc", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Chc")
#Highly expressed in cholinergic, glutamatergic, GABAergic neurons, and Glial cells 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "CRMP", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "CRMP")
#Most highly expressed in cholinergic neurons BUT also expressed in monomergic, glutamergic and GABAergic neurons

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Sra-1", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Sra-1")
#Medium/Low expression in most cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "shi", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "shi")
#Highly expressed in monoaminergic, glutamatergic, cholinergic, and GABAergic neurons 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Eph", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Eph")
#Medium/Low expression in all cell types, most highly expressed in monoaminergic, glutamatergic, cholinergic, and GABAergic neurons  

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "cher", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "cher")
#Very low expression in most cell types, highest in glial cells 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "CG8916", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "CG8916")
#Very low expression in most cell types

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Ufd4", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Ufd4")
#Medium expression in most cell types, highest in monoaminergic, cholinergic, and GABAergic neurons as well as in glial cells 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "if", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "if")
#Very low expression in most cell types

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "dop", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "dop")
#Low/Medium expression in most cell types, highest in peptidergic and monoaminergic neurons 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "CG15439", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "CG15439")
#Low expression in all cell types except peptidergic neurons 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Rap1", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Rap1")
#High expression in all cell types with the lowest expression in neural progenitors 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "rhea", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "rhea")
#High expression in peptidergic neurons 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "TER94", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "TER94")
#Medium/High expression in all cell types, with highest expression in cholinergic neurons and neural progenitors 

#Causal genes: 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Nrx-1", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Nrx-1")
#High expression in all neuronal cell types except peptidergic neurons where it is expressed a bit lower 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "rg", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "rg")
#Very low expression in all cell types, except in glial cells where it is not expressed

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Act5C", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Act5C")
#Very high expression in all cell types

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Ssadh", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Ssadh")
#Medium expression in all cell types with highest expression in glial cells and GABAergic cells

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Raf", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Raf")
#Medium/low expression in all cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "sti", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "sti")
#Highly expressed in neural progenitors and not other cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "bru3", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "bru3")
#Highly expressed in all neural cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "cic", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "cic")
#Medium expression in all cell types, highest being in glutamatergic, cholinergic, and GABAergic neurons 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Cul3", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Cul3")
#Low/medium expression in all cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "DIP2", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "DIP2")
#Low/medium expression in all cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "fne", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "fne")
#High expression in all neuronal cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Lcch3", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Lcch3")
#Medium expression in all neuronal cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Pits", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Pits")
#Medium/high expression in all cell types, highest in GABAergic, cholinergic, glutamatergic, and monoaminergic neurons 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Chi", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Chi")
#Low/medium expression in all cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "skd", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "skd")
#Medium expression in all cell types 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "hiw", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "hiw")
#Low expression in all cell types, highest in cholinergic neurons

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Nup154", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Nup154")
#Low expression in all cell types, highest in neural progenitors

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Rfx", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Rfx")
#Low expression in all cell types, except glutamatergic neurons where it is mediumly expressed 

FeaturePlot_scCustom(seurat_object = L1_brain_drosophila, features = "Sos", order = T, colors_use = pal, label = T, na_cutoff = 0.5, pt.size = 1, reduction = "umap")
VlnPlot(L1_brain_drosophila, features = "Sos")
#Low expression in all cell types

#Generating table of genes of interest and their average expression in each cluster:

#Define the gene lists
related_genes <- c("vnd", "SK", "Abl", "CG8665", "CG9650", "Rep", "Chc", "CRMP", "Sra-1", "shi", "Eph", "cher", "CG8916", "Ufd4", "if", "dop", "CG15439", "Rap1", "rhea", "TER94")
causal_genes <- c("Nrx-1", "rg", "Act5C", "Ssadh", "Raf", "stj", "bru3", "cic", "Cul3", "DIP2", "fne", "Lcch3", "Pits", "Chi", "skd", "hiw", "Nup154", "Rfx", "Sos")

#Compute average expression per cluster for each list
avg_rel <- AverageExpression(L1_brain_drosophila, features = related_genes)$RNA
avg_caus <- AverageExpression(L1_brain_drosophila, features = causal_genes)$RNA

#Average expression across genes within each list 
cluster_avg_rel <- colMeans(avg_rel, na.rm = TRUE)
cluster_avg_caus <- colMeans(avg_caus, na.rm = TRUE)

# Difference between lists
diff_relcaus <- cluster_avg_caus - cluster_avg_rel

# Plot expression differences 
desired_order <- c("Neural progenitors", "Monoaminergic neurons", "Glutamatergic neurons", "Cholinergic neurons", "Peptidergic neurons", "Glial cells", "GABAergic neurons")

results$cluster <- factor(results$cluster, levels = desired_order)

expr_clust_causrel <- ggplot(results, 
                             aes(x = cluster, y = observed_difference, fill = cluster)) +
  geom_col() +
  geom_text(aes(label = gsub(" ", "\n", cluster)),
            vjust = -0.5,
            size = 7,
            lineheight = 0.9) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) + 
  theme_bw(base_size = 25) +
  labs(
    y = "Δ Avg Expr (Causal − Related)",
    x = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none"
  )
ggsave("expression_casrel.png", plot = expr_clust_causrel, width = 15, height = 8, dpi = 600)

# Combine both gene lists
all_genes <- unique(c(related_genes, causal_genes))
n_rel <- length(related_genes)
n_caus <- length(causal_genes)

# Save observed difference
T_obs <- diff_relcaus 

# Permutation test
n_perm <- 5000
set.seed(123) #random number generator 

avg_all <- AverageExpression(L1_brain_drosophila, features = all_genes)$RNA

# Actual differences
obs_rel <- colMeans(avg_all[related_genes, , drop=FALSE], na.rm=TRUE)
obs_caus <- colMeans(avg_all[causal_genes, , drop=FALSE], na.rm=TRUE)
T_obs <- obs_caus - obs_rel

T_perm <- matrix(NA, nrow=n_perm, ncol=ncol(avg_all))
colnames(T_perm) <- colnames(avg_all)

for (i in 1:n_perm) {
  perm <- sample(all_genes)
  rel_perm <- perm[1:n_rel]
  caus_perm <- perm[(n_rel+1):(n_rel+n_caus)]
  rel_vals <- colMeans(avg_all[rel_perm, , drop=FALSE])
  caus_vals <- colMeans(avg_all[caus_perm, , drop=FALSE])
  T_perm[i, ] <- caus_vals - rel_vals
}

# Compute p-values per cluster
p_vals <- sapply(seq_along(T_obs), function(j) {
  mean(T_perm[, j] >= T_obs[j])
})
names(p_vals) <- names(T_obs)

# Combine results
results <- data.frame(
  cluster = names(T_obs),
  observed_difference = T_obs,
  p_value = p_vals
)

# Shorten very small p-values for display
results$label <- ifelse(results$p_value < 0.001, "<0.001", sprintf("p=%.3g", results$p_value))

# Plot with p-values on top of bars
causrel_pvalue <- ggplot(results, aes(x = cluster, y = observed_difference, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = my_colors) +
  theme_bw(base_size = 25) +
  labs(
    y = "Δ Avg Expr (Causal − Related)",
    x = NULL
  ) +
geom_text(
  aes(label = label),
  vjust = -0.5,
  size = 6,
  fontface = "bold"
) +
scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
theme(
  legend.position = "none",
  axis.text.x = element_text(size = 20, vjust = 1)
)
ggsave("p-value_expression_casrel.png", plot = causrel_pvalue, width = 15, height = 8, dpi = 600)
