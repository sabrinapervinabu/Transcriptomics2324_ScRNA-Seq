# http://159.149.160.56/GT_2023_sc/Seurat.html
# https://panglaodb.se/view_data.php?sra=SRA713577&srs=SRS3363004
# SRA713577	SRS3363004	Peripheral blood mononuclear cells #Homo sapiens

install.packages("tidyverse")
install.packages('Seurat')
library(Seurat)
    ###term “features” to denote genes

library(patchwork)
library(ggplot2)

load("~/Downloads/trascrittomica/scRNA-seq/SRA713577_SRS3363004.sparse.RData")

#### Rename the rows of the table with just the gene symbol ####
#split the name in two part, one with the correct one and one with the _ENSG part that i want delete
splitnames = strsplit(rownames(sm), "_ENSG")
# la funzione sapply() simplifies the split into a vector, in this case I create a vector with the first part of the
# all the splits using '[',1
splitnames = sapply(splitnames, `[`, 1)
# I need to replace the _ with - because Seurat give me a warning that say:
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Error in validObject(.Object) :
# invalid class “LogMap” object: Duplicate rownames not allowed
splitnames = gsub("_", "-", splitnames)
# Here, I want that all the gene symbols in the splitnames vector are unique because having duplicate
# row names in your matrix can cause issues when creating the Seurat object
unique_gene_symbols = make.unique(splitnames)
# I add the correct row names and I create the Seurat object
rownames(sm) = unique_gene_symbols #Renaming the rows of the table with just the gene symbo

#### Initialize the Seurat object with the raw (non-normalized) data ####
PB_data <- CreateSeuratObject(counts = sm,
                              project = "scRNA-Seq", min.cells = 3, min.features = 200)
PB_data
    #An object of class Seurat 
    #21070 features across 3568 samples within 1 assay 
    #Active assay: RNA (21070 features, 0 variable features)
    #1 layer present: counts

# Lets examine the counts for a few genes in the first thirty cells. They are in sm
sm[c("CD3D", "TCL1A", "MS4A1"), 1:30]
    #3 x 30 sparse Matrix of class "dgCMatrix"
    #[[ suppressing 30 column names ‘AAACCTGAGAACAACT’, ‘AAACCTGAGCGTGTCC’, ‘AAACCTGAGTACACCT’ ... ]]

    #CD3D  . 1 1 . 1 . 1 . . 1 2 . . . . 1 1 . . . 5 . . 3 . . . . 5 2
    #TCL1A . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #MS4A1 . . . . . . . . . . 1 . . . . . . . . . . 3 . . . . . . . .

    #The . values in the matrix represent zeroes (no reads for that gene in that cell).
    #most values in an scRNA-seq matrix are 0, #Seurat uses a sparse-matrix representation whenever possible
    #The data type is “dgCMatrix", resulting in significant memory and speed savings 

#Column names are the original barcodes of each cell
head(colnames(PB_data))

#### Cell quality control #### 
# main quality control parameters
# - the number of unique genes detected in each cell: low-quality cells or empty droplets will often have very few genes 
# - cell doublets or multiplets exhibit an aberrantly high gene count 
# - the total number of molecules detected within a cell (correlates strongly with unique genes) 
# - the percentage of reads that map to the mitochondrial genome: low-quality / dying cells often exhibit extensive mitochondrial contamination

# We can calculate mitochondrial QC metrics with the PercentageFeatureSet() function, 
# which calculates the percentage of counts originating from a selected set of features (genes)
grep("^MT-",rownames(PB_data),value = TRUE)
    # [1] "MT-ATP6" "MT-ATP8" "MT-CO1"  "MT-CO2"  "MT-CO3"  "MT-CYB"  "MT-ND1"  "MT-ND2"  "MT-ND3"  "MT-ND4" 
    # [11] "MT-ND4L" "MT-ND5"  "MT-ND6"  "MT-RNR1" "MT-RNR2" "MT-TP"  

# In scRNA-Seq count tables MT gene names start with “MT-” 
# (do not forget the “-” symbol, there are nuclear genes with name starting with just MT without the dash!).
# Also remember that in mouse gene names are usually in lowercase (Mt- or even mt-).

# The [[ operator can add columns to object metadata. This is a great place to store additional info/data
# “Pattern” defines the criterion with which we select the genes by name. 
# “^MT-” means “starting with MT-” (starting is the ^ symbol)
PB_data[["percent.mt"]] <- PercentageFeatureSet(PB_data, pattern = "^MT-")

#ribosomal protein genes “eat up” a lot of reads because highly expressed. 
#Their gene symbol usually starts by RPL or RPS
grep("^RP[LS]",rownames(PB_data),value = TRUE)
PB_data[["percent.rbp"]] <- PercentageFeatureSet(PB_data, pattern = "^RP[LS]")

# The number of unique genes (called here features) and total molecules (reads after UMI filtering) 
# are automatically calculated during CreateSeuratObject(). 
# You can find them stored in the object meta data, together with the values just computed
head(PB_data@meta.data, 5)

#### Quality control violin plots ####
# Visualize QC metrics as violin plots - also adding the RPL genes
VlnPlot(PB_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4)
#without dots:
VlnPlot(PB_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4, pt.size=0)

# We can check if the different parameters are correlated with one another

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(PB_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PB_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot3 <- FeatureScatter(PB_data, feature1 = "nCount_RNA", feature2 = "percent.rbp")
plot3
plot1 + plot2 + plot3

# All in all the only visible correlation is between the number of reads and the number of genes detected. 
# When they are too low, the droplet was empty. Too high, probably a doublet. 
# On the basis of these plot, we have to decide thresholds for cell quality control.

PB_data <- subset(PB_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#The thresholds are on the number of genes (and not “features”!!!!!) detected (between 200 and 2500), 
#and on the MT DNA/RNA (5%), chosen by looking at the distributions above. 
#Cells removed are “outilers” in the distributions. 
#Notice also that the % of ribosomal protein genes is not usually employed at this stage to filter out cells, 
#it is a “natural” parameter that we can use later on.

# Let us see how many samples (cells) remain:
PB_data
    #An object of class Seurat 
    # 21070 features across #3095 samples within 1 assay VS 21070 features across #3568 samples within 1 assay 
    #Active assay: RNA (21070 features, 0 variable features)
    #1 layer present: counts

#Almost 500 cells have been removed. 

#### Normalizing the data ####
#10x data are usually just transformed into counts per 10,000 reads. 
#But, the final “expression estimate” it’s given by the log of the normalized counts.
PB_data <- NormalizeData(PB_data, normalization.method = "LogNormalize", scale.factor = 10000)

# The original and normalized counts are buried inside the Seurat object PB_data. Let us find them:
PB_data@assays
    #$RNA
    #Assay (v5) data with 21070 features for 3095 cells
    #First 10 features:
    #  A1BG-AS1, A1BG, A2M-AS1, A2M, A2MP1, AAAS, AACS, AAED1, AAGAB, AAK1 
    #Layers:
    #  counts, data 
PB_data@assays$RNA
    #Assay (v5) data with 21070 features for 3095 cells
    #First 10 features:
    #  A1BG-AS1, A1BG, A2M-AS1, A2M, A2MP1, AAAS, AACS, AAED1, AAGAB, AAK1 
    #Layers:
    #  counts, data 

#raw counts are here, uncomment to take a look - two possible notations
#LayerData(PB_data, assay = "RNA", layer = "counts")
#PB_data[["RNA"]]$counts
#normalized counts are here - two possible notations
#LayerData(PB_data, assay = "RNA", layer = "data")
#PB_data[["RNA"]]$data
#normalized counts after scaling are here
#LayerData(PB_data, assay = "RNA", layer = "scale.data")
#PB_data[["RNA"]]$scale.data

#Now that we have normalized the expression values, we can take a look to the genes 
# that have the highest mean expression across our cells:
apply(PB_data[["RNA"]]$data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
head(gene.expression, n=50)
#   MALAT1    RPS27      B2M   TMSB4X    RPL10    RPL13    RPL21    RPL34   RPL13A    RPS18    RPS14    RPL32 
# 6.145262 4.888452 4.884248 4.842857 4.635024 4.519836 4.489505 4.456312 4.451034 4.432674 4.259496 4.234354 

#And to the distribution of the expression of MALAT1 and of another typical housekeeping gene (GAPDH yok):
VlnPlot(PB_data, features = c("MALAT1","GAPDH"))
VlnPlot(PB_data, features = c("MALAT1","RPL13A"))  #[alternative housekeeping genes: ACTB", "RPL13A", "HPRT1"]
VlnPlot(PB_data, features = c("MALAT1","ACTB")) #prettiest - codifica per la proteina beta-actina
VlnPlot(PB_data, features = c("MALAT1","HPRT1")) #yok :\

cc.genes.updated.2019
CellCycleScoring(PB_data, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> PB_data
PB_data[[]]
PB_data <- FindVariableFeatures(PB_data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(PB_data), 10)
plot1 <- VariableFeaturePlot(PB_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#The idea is to shift the expression of each gene, so that, across cells, the mean expression is 0 and the variance is 1.
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate. 
# In practice, values are sort of “binarized”, or rather “ternarized” -> 
      # 0 “high expression”, 0 “average expression” <0 “under expression, or no expression at all”. 
#Notice that this is done for all the genes - not only the most variable ones.
all.genes <- rownames(PB_data)
PB_data <- ScaleData(PB_data, features = all.genes)
#Let us see the scaled values for a gene:
PB_data@assays$RNA
    #Assay (v5) data with 21070 features for 3095 cells
    #First 10 features:
     # A1BG-AS1, A1BG, A2M-AS1, A2M, A2MP1, AAAS, AACS, AAED1, AAGAB, AAK1 
    #Layers:
     # counts, data, scale.data 

#The results of this are stored in 
#PB_data[["RNA"]]$scale.data
#example
#PB_data[["RNA"]]$scale.data["A1BG",]

#PERO' ONLINE DICE COSI' - USA @ E NON $
#The results of this are stored in 
#pbmc[["RNA"]]@scale.data
#example
#pbmc[["RNA"]]@scale.data["MS4A1",]


#### DIMENSION REDUCTION #### WHY DO WE DO IT
#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#The recommended method for 10x data is PCA. 
#Notice that it is performed on the “variable features” only (the 2000 most variable genes), not on all the genes.
# Also notice: the variable genes are identified considering the log(x+1) of counts, 
# while the PCA is performed after data have been scaled as explained above.
PB_data <- FindVariableFeatures(PB_data, selection.method = "vst", nfeatures = 2000)
variable.features <- VariableFeatures(PB_data)
print(length(variable.features)) 

PB_data <- RunPCA(PB_data, features = VariableFeatures(object = PB_data))
# Examine and visualize PCA results a few different ways
print(PB_data[["pca"]], dims = 1:5, nfeatures = 5)
    # PC_ 1 
    #Positive:  CSTA, CST3, CLEC7A, FCN1, MNDA
    #Negative:  RPS18, RPS4X, AC073861.1, EEF1A1, RPL34 
    # PC_ 2 
    #Positive:  RPL39, RPL34, RPS18, RPLP1, RPS8 
    #Negative:  ACRBP, CMTM5, TMEM40, CAVIN2, GP9 
    #PC_ 3 
    #Positive:  LINC00926, CD79A, MS4A1, TNFRSF13C, IGHM 
    #Negative:  S100A4, S100A6, TYROBP, NKG7, SRGN 
    #PC_ 4 
    #Positive:  NKG7, GZMB, KLRD1, CST7, FGFBP2 
    #Negative:  IL7R, LTB, NOSIP, TRAC, PRKCQ-AS1 
    #PC_ 5 
    #Positive:  G0S2, CMTM2, FCGR3B, NAMPT, CXCR2   
    #Negative:  PPP1R14B, SERPINF1, NRP1, RPS4X, CLEC4C

VizDimLoadings(PB_data, dims = 1:2, reduction = "pca")
#NO Cell cycle phase grouping.
#Notice that the cells were colored according to the CC phase
DimPlot(PB_data, reduction = "pca")
#The top principal components therefore represent a robust compression 
# of the dataset. However, how many components should we choose to include?

# In the Seurat tutorial the “Jackstraw” procedure is explained. 
# Another more intuitive way is to represent the % of variance explained by each PC
#with ndims we can choose how many PC to plot
ElbowPlot(PB_data, ndims=50)

# Esegui il PCA sul tuo dataset
PB_data <- RunPCA(PB_data, features = VariableFeatures(object = PB_data))
# Estrai la varianza spiegata per ciascun componente principale
variance_explained <- PB_data[["pca"]]@stdev^2
variance_explained <- variance_explained / sum(variance_explained)
# Calcola la somma cumulativa della varianza spiegata
cumulative_variance_explained <- cumsum(variance_explained)
# Trova il numero di componenti principali che spiegano almeno il 75% della varianza
num_pcs_75_variance <- which(cumulative_variance_explained >= 0.75)[1]
# Stampa il risultato
num_pcs_75_variance
ElbowPlot(PB_data, ndims = 50) + 
  geom_vline(xintercept = num_pcs_75_variance, color = "blue", linetype = "dashed") + 
  geom_vline(xintercept = 10, color = "red", linetype = "dashed") +
  ggtitle("Elbow Plot with Variance Explained and Chosen PCs")

ElbowPlot(PB_data, ndims = 50) + 
  geom_vline(xintercept = num_pcs_75_variance, color = "blue", linetype = "dashed") + 
  ggtitle("Elbow Plot with Chosen PCs")

# In this example, all three approaches yielded similar results, 
# but we might have been justified in choosing anything between PC 7-12 as a cutoff.

# Performing downstream analyses with only 5 PCs does significantly and adversely affect results. 
# On the other hand, too many PCs can introduce noise and split cells by technical variation only.

# Another very effective rule is to keep all the PCs until 70-75% of the variance is explained
pc.touse <- (PB_data$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.75))
pc.touse
#[1] 20 il num di PC da usare 
#So in this way we would use many more dimensions than in the vignette…! >>> e quindi???

#### CLUSTERING ####
# Seurat first constructs a kNN graph based on the euclidean distance in PCA space, and refines the edge 
# weights between any two cells based on the shared overlap in their local neighborhoods

#This step is performed using the FindNeighbors() function, and takes as input the previously defined 
# dimensionality of the dataset (first 10 PCs).By default, the “k” of kNN is set to 20.

PB_data <- FindNeighbors(PB_data, dims = 1:20)

# the Louvain algorithm (default) are applied to iteratively group cells together, 
                    # with the goal of optimizing the standard modularity function.
#The FindClusters() function implements this procedure, and contains a resolution parameter 
# that sets the ‘granularity’ of the downstream clustering, 
# with increased values leading to a greater number of clusters

# Seurat authors find that setting this parameter between 0.4-1.2 typically returns good results 
# for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.

# In practice, you have to try different combinations of PCA dimensions and resolution, 
              # until you are satisfied by the results.
#The clusters can be found using the Idents() function, or in the field “seurat_clusters” of the pbmc object.

PB_data <- FindClusters(PB_data, resolution = 0.5)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

  ## Number of nodes: 3095
  ## Number of edges: 129454

  ## Running Louvain algorithm...
  ## Maximum modularity in 10 random starts: 0.8821
  ## Number of communities: 11
  ## Elapsed time: 0 seconds

# Look at cluster IDs of the first 5 cells
head(Idents(PB_data), 5)

## AAACCTGAGTACACCT AAACCTGAGTGTACGG AAACCTGCATCACAAC AAACCTGGTAAATGAC AAACCTGTCAGCTCGG 
## 1                2                0                1                5 
## Levels: 0 1 2 3 4 5 6 7 8 9 10

head(PB_data[[]],5)
##                  orig.ident nCount_RNA nFeature_RNA percent.mt percent.rbp      S.Score   G2M.Score Phase old.ident RNA_snn_res.0.5 seurat_clusters
## AAACCTGAGTACACCT  scRNA-Seq       2662         1039   2.329076    38.65515  0.009254571 -0.05628027     S scRNA-Seq               1               1
## AAACCTGAGTGTACGG  scRNA-Seq       2535         1048   4.733728    29.78304  0.048466029  0.05695519   G2M scRNA-Seq               2               2
## AAACCTGCATCACAAC  scRNA-Seq       4348         1331   2.782889    46.11316  0.022314370 -0.01959941     S scRNA-Seq               0               0
## AAACCTGGTAAATGAC  scRNA-Seq       4140         1467   2.584541    40.84541  0.012167490  0.02551380   G2M scRNA-Seq               1               1
## AAACCTGTCAGCTCGG  scRNA-Seq       1027          508   1.752678    43.62220 -0.016401227  0.11994602   G2M scRNA-Seq               5               5

# In all 11 clusters were found, numbered from 0 to 10. 
# We can plot them in the space of the first two PCA components

DimPlot(PB_data, reduction = "pca")

# We can see how the points are indeed spread by the first two PCs, 
# but there is still some overlap among clusters. Remember that clustering is done on 20 PCs, not only 2.

# Or we select to look at the projection along any two of the 20 PC as we want:

DimPlot(PB_data,reduction="pca", dims=c(4,9))
DimPlot(PB_data,reduction="pca", dims=c(4,5))

# But we know that for visualization and 2D plotting there are better strategies. 
# t_SNE, always on the PC dimensions chosen for clustering:
PB_data <- RunTSNE(PB_data, dims=1:20)
p <- DimPlot(PB_data, reduction = "tsne")

# Estrazione delle informazioni sui cluster
cluster_counts <- table(Idents(PB_data))
cluster_counts_df <- as.data.frame(cluster_counts)
colnames(cluster_counts_df) <- c("cluster", "count")

# Estrazione delle coordinate centrali dei cluster
tsne_data <- as.data.frame(PB_data@reductions$tsne@cell.embeddings)
tsne_data$cluster <- as.factor(Idents(PB_data))
centroids <- aggregate(cbind(tSNE_1, tSNE_2) ~ cluster, tsne_data, mean)

# Unione dei dati dei cluster con i conteggi
centroids$count <- cluster_counts_df$count[match(centroids$cluster, cluster_counts_df$cluster)]

# Aggiunta delle annotazioni al plot
p + geom_text(data=centroids, aes(x=tSNE_1, y=tSNE_2, label=paste(count)),
              color="black", size=4, vjust=-1, fontface="bold")

# Or UMAP, that as of today is the preferred method:
library(ggplot2)
library(reticulate)
reticulate::py_install(packages ='umap-learn')
PB_data <- RunUMAP(PB_data, dims = 1:20)
umap_plot <- DimPlot(PB_data, reduction = "umap")

# Estrazione delle informazioni sui cluster
cluster_counts <- table(Idents(PB_data))
cluster_counts_df <- as.data.frame(cluster_counts)
colnames(cluster_counts_df) <- c("cluster", "count")

# Estrazione delle coordinate centrali dei cluster
umap_data <- as.data.frame(PB_data@reductions$umap@cell.embeddings)
umap_data$cluster <- as.factor(Idents(PB_data))
centroids <- aggregate(cbind(umap_1, umap_2) ~ cluster, umap_data, mean)

# Unione dei dati dei cluster con i conteggi
centroids$count <- cluster_counts_df$count[match(centroids$cluster, cluster_counts_df$cluster)]

# Aggiunta delle annotazioni al plot
umap_plot <- umap_plot + geom_text(data=centroids, aes(x=umap_1, y=umap_2, label=paste(count)),
                                   color="black", size=4, vjust=-1, fontface="bold")

# Visualizzazione del plot
print(umap_plot)

# We can also check whether some of the critical quality parameters influenced the clustering we got:
VlnPlot(PB_data,features="nCount_RNA")
VlnPlot(PB_data,features="nFeature_RNA")
VlnPlot(PB_data,features="percent.mt")
VlnPlot(PB_data,features="percent.rbp")
# As we can see, there seem to be two clusters (7 and 8) in which we can notice relevant differences in 
# library sizes and number of expressed genes. 

# On the other hand (see later on) we will be anyway able to assign a cell identity to both. 
# So, its is not a “techical problem”. They indeed are cells that contain less RNA than the others. 
      # anche x il mio dataset ????

#maybe the problem can come from a specific cell-cycle phase  
install.packages("dplyr")
library(dplyr)
library(ggplot2)

PB_data@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")
# A sizable number of cells of cluster 8 is in the G2M phase, 
# with respect to the other clusters. Once again, if we are able to assign a cell type to the cluster 
# this is not a problem.

#### Finding “marker” genes and assigning cell types to clusters ####
# Seurat includes a function that can be used to find genes 
        # a) over expressed between two clusters  
        # b) overexpressed in one cluster with respect to all the others. 
# The function permits to employ different tests, including those used for bulk RNA-Seq. 
# For 10x data, the choice is to employ a non parametric test (once again, the Wilcoxon test!) 
  # which is the default. 
# Notice also another parameter (min.pct): it means that a gene has to be found expressed (counts > 0)
# in at least 25% of the cells of the cluster.

# find all markers of cluster 2 versus all the others
cluster2.markers <- FindMarkers(PB_data, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
head(cluster2.markers, n = 5)
#              p_val avg_log2FC pct.1 pct.2     p_val_adj
#CCL5  4.371877e-168   2.203083 0.990 0.465 9.211544e-164
#TRGC2 4.198333e-123   3.322982 0.505 0.084 8.845888e-119
#CD8A  5.450401e-116   3.159085 0.480 0.081 1.148400e-111
#LYAR  2.381836e-114   2.272279 0.683 0.191 5.018528e-110
#CST7  2.326201e-107   1.491196 0.782 0.233 4.901306e-103

# The above are the markers of cluster 2, with fold change, %of cells expressing it, pvalue and FDR. 
# Sorted by increasing FDR. 
# Could you guess which is the corresponding cell type?

cluster0.markers <- FindMarkers(PB_data, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
head(cluster0.markers, n = 5)
#              p_val avg_log2FC pct.1 pct.2     p_val_adj
#RPL21 2.487040e-139  0.4830454     1 0.984 5.240193e-135
#RPL34 2.889045e-126  0.4844558     1 0.983 6.087217e-122
#RPS3A 5.439967e-119  0.4909968     1 0.975 1.146201e-114
#RPL32 1.697302e-111  0.4544616     1 0.984 3.576215e-107
#RPS27 8.000789e-102  0.3839355     1 0.989  1.685766e-97

#non va bene continuare così perchè ne escono i RibosomalProteins 


# The comparison can be made against one or more of the other clusters. 
# Cluster 2 versus 0 and 1:
cluster2_01.markers <- FindMarkers(PB_data, ident.1 = 2, ident.2 = c(0, 1), min.pct = 0.25)
head(cluster2_01.markers, n = 5)

#               p_val avg_log2FC pct.1 pct.2     p_val_adj
#CCL5  1.176365e-199   3.882564 0.990 0.401 2.478600e-195
#CST7  1.625789e-168   4.028134 0.782 0.123 3.425537e-164
#NKG7  1.233227e-167   4.289905 0.911 0.343 2.598410e-163
#KLRD1 1.342835e-141   5.042117 0.587 0.049 2.829353e-137
#GZMA  2.207277e-132   2.850336 0.814 0.223 4.650732e-128

# Cluster-cluster comparisons might be useful in case the “one versus all” analysis returns 
# some common genes for more than one cluster, or if at the end of the annotation 
# we have two or more clusters attributed to the same cell type

##### confronto dei geni tra tutti i cluster, tutti vs tutti ####
#The one vs. all analysis can be iterated automatically:
# we return only genes "over expressed", found in at least 25% of the cells, 
# and with a logFC threshold of at least 0.25

PB_data.markers <- FindAllMarkers(PB_data, only.pos = TRUE, min.pct = 0.25, 
                                  logfc.threshold = 0.25)

#And we can output the top n (in this case 5) genes for each cluster. 
#Notice that here they are sorted by logFC - more informative than “p_val_adj”, 
#since a lot of genes will have a FDR close to zero with smallest changes
PB_data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) 
## A tibble: 55 × 7
## Groups:   cluster [11]
#.     p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene      
#.     <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>     
# 1 2.15e-80       3.04 0.259 0.033  4.53e-76 0       NOG       
# 2 6.71e-60       2.38 0.251 0.047  1.41e-55 0       FHIT      
# 3 2.63e-79       2.13 0.419 0.112  5.54e-75 0       IL6ST     
# 4 2.34e-99       1.91 0.595 0.198  4.93e-95 0       CCR7      
# 5 3.50e-62       1.87 0.41  0.132  7.38e-58 0       MAL       
# 6 4.65e-84       2.92 0.281 0.036  9.80e-80 1       TTC39C-AS1
# 7 3.63e-54       1.84 0.333 0.089  7.65e-50 1       AQP3      
# 8 5.70e-38       1.50 0.309 0.102  1.20e-33 1       CORO1B    
# 9 2.74e-31       1.49 0.255 0.084  5.78e-27 1       TNFRSF25  
#10 2.61e-42       1.43 0.353 0.117  5.50e-38 1       TRAT1     
## ℹ 45 more rows
## ℹ Use `print(n = ...)` to see more rows

PB_data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print (n=55) #avrai tt i geni
#    p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene  
# 1 2.15e- 80       3.04 0.259 0.033 4.53e- 76 0       NOG      0 e 6 Extravillous trophoblasts 
# 2 6.71e- 60       2.38 0.251 0.047 1.41e- 55 0       FHIT !!! ok, xo dropouts Oligodendrocytes/Oligodendrocyte precursor cells (brain)
# 3 2.63e- 79       2.13 0.419 0.112 5.54e- 75 0       IL6ST    Glandular and luminal cells	
# 4 2.34e- 99       1.91 0.595 0.198 4.93e- 95 0       CCR7     0 e 6 Langerhans cells(skin)/Syncytiotrophoblasts(placenta)
# 5 3.50e- 62       1.87 0.41  0.132 7.38e- 58 0       MAL       

VlnPlot(PB_data, features = c("FHIT", "S100A4")) #CLUSTER 0: FHIT Oligodendrocytes

# 6 4.65e- 84       2.92 0.281 0.036 9.80e- 80 1       TTC39C-AS1     ok
# 7 3.63e- 54       1.84 0.333 0.089 7.65e- 50 1       AQP3           ok Squamous epithelial cells
# 8 5.70e- 38       1.50 0.309 0.102 1.20e- 33 1       CORO1B         ok not marker Distal enterocytes
# 9 2.74e- 31       1.49 0.255 0.084 5.78e- 27 1       TNFRSF25   !!! ok T-cells
#10 2.61e- 42       1.43 0.353 0.117 5.50e- 38 1       TRAT1          0 e 1 T-cells

VlnPlot(PB_data, features = c("TNFRSF25", "S100A4"))

#11 4.20e-123       3.32 0.505 0.084 8.85e-119 2       TRGC  !!! ok T-cells/Nk-cells
#12 5.45e-116       3.16 0.48  0.081 1.15e-111 2       CD8A      2 e 6  T-cells/Nk-cells
#13 1.49e- 77       3.10 0.302 0.041 3.13e- 73 2       GZMK  !!! ok T-cells/Nk-cells
#14 4.31e- 80       2.49 0.515 0.144 9.08e- 76 2       GZMH      2 e 3  Nk-cells/T-cells
#15 1.33e- 76       2.48 0.426 0.096 2.80e- 72 2       KLRG1 !!! ok Nk-cells/Granulocytes/T-cells

VlnPlot(PB_data, features = c("TRGC2", "S100A4"))

#16 5.90e-186       6.25 0.374 0.011 1.24e-181 3       IGFBP7     ok Peritubular cells nah
#17 0               5.33 0.805 0.055 0         3       KLRF1  !!! very ok Nk-cells NETTAMENTE up regulated 
#18 1.11e-111       5.19 0.264 0.013 2.34e-107 3       SH2D1B     ok Nk-cells NETTAMENTE up regulated 
#19 1.44e-233       5.12 0.629 0.057 3.03e-229 3       SPON2      very ok Nk-cells NETTAMENTE up regulated 
#20 7.95e-170       4.89 0.442 0.031 1.67e-165 3       KLRC1      ok Nk-cells/T-cells

VlnPlot(PB_data, features = c("KLRF1", "S100A4"))

#21 1.49e-232       6.69 0.42  0.006 3.15e-228 4  !!!  IGHG2     ok Plasma cells - B-cells
#22 1.22e-128       6.57 0.286 0.01  2.57e-124 4       IGHGP     meh
#23 1.28e-112       6.34 0.272 0.012 2.69e-108 4       IGHA2     meh
#24 9.14e-235       6.16 0.457 0.01  1.93e-230 4       IGHG4     ok Plasma cells - B-cells
#25 1.30e-259       6.02 0.558 0.018 2.74e-255 4       IGHG3     best Plasma cells - B-cells

VlnPlot(PB_data, features = c("IGHG2", "S100A4"))

#26 1.70e-  5       5.72 0.192 0.403 3.59e-  1 5   -   HBB       Erythroid cells (liver and bon marrow)
#27 4.24e- 20       2.31 0.444 0.24  8.93e- 16 5   ?   B4GALT1    Monocytes/Langerhans cells
#28 1.00e- 41       2.19 0.697 0.423 2.11e- 37 5   ?   HNRNPH1    Langerhans cells
#29 3.40e-  8       1.91 0.308 0.196 7.16e-  4 5   ?   APOBEC3C   Nk-cells
#30 2.91e- 13       1.86 0.423 0.272 6.14e-  9 5   ?   C16orf54   Nk-cells

VlnPlot(PB_data, features = c("HNRNPH1", "S100A4"))

#31 3.99e-219       4.22 0.718 0.046 8.41e-215 6       LINC02446     ok. T cells
#32 6.05e- 66       3.03 0.381 0.049 1.28e- 61 6       CA6           ok. Serous glandular cells
#33 1.60e-131       2.68 0.867 0.152 3.37e-127 6       CD8B          2 e 6  
#34 2.65e- 44       2.66 0.254 0.032 5.58e- 40 6       LRRN3         meh 
#35 4.72e- 35       2.10 0.37  0.083 9.94e- 31 6       NELL2         ok  Excitatory neurons

VlnPlot(PB_data, features = c("NELL2", "S100A4"))

#36 0               8.32 0.837 0.017 0         7       LYZ     !!! Langerhans cells
#37 3.39e-205       6.53 0.416 0.008 7.14e-201 7       MAFB        Monocytes 
#38 0               6.44 0.882 0.032 0         7       CSTA        7 e 10   Suprabasal keratinocytes
#39 1.93e-158       6.24 0.303 0.004 4.06e-154 7       HNMT        Hofbauer cell
#40 4.78e-198       6.11 0.416 0.009 1.01e-193 7       LILRA5      Monocytes

VlnPlot(PB_data, features = c("LILRA5", "S100A4"))

#41 0               6.26 0.893 0.026 0         8       TCL1A    !!! Dendritic cells
#42 4.11e-269       5.53 0.878 0.047 8.66e-265 8       IGHD         8 e 4  Plasma cells
#43 5.83e-128       5.23 0.351 0.011 1.23e-123 8       COL19A1      Skeletal myocytes
#44 1.39e-174       5.00 0.595 0.031 2.93e-170 8       FCER2        B-cells
#45 3.11e-103       4.73 0.351 0.017 6.54e- 99 8       PCDH9        Oligodendrocytes

VlnPlot(PB_data, features = c("TCL1A", "S100A4"))

#46 0              10.1  0.765 0.007 0         9       TMEM40       Syncytiotrophoblasts placenta
#47 3.45e-210       9.81 0.471 0.004 7.28e-206 9       AC090409.1
#48 4.17e-170       9.77 0.353 0.002 8.78e-166 9       CLDN5  !!!   Adipocytes
#49 4.71e-237       9.68 0.588 0.007 9.93e-233 9       CTTN         not marker 
#50 2.09e-264       9.59 0.882 0.019 4.41e-260 9       ACRBP        spermatids

VlnPlot(PB_data, features = c("ACRBP", "S100A4"))

#51 6.88e- 86      10.6  0.278 0.002 1.45e- 81 10      CXCL1   !!!  Basal respiratory cells
#52 1.54e-120       9.93 0.444 0.004 3.25e-116 10      LINC01506 
#53 5.69e-208       9.78 0.833 0.008 1.20e-203 10      FCGR3B       Serous glandular cells
#54 1.18e-301       9.76 0.889 0.005 2.48e-297 10      CMTM2        spermatids
#55 5.72e-201       9.62 0.667 0.005 1.20e-196 10      CXCR2        Nk-cells


#CLUSTER 0: FHIT Oligodendrocytes
#CLUSTER 1: TNFRSF25 T-cells
#CLUSTER 2: KLRG1 Granulocytes (T-cells/Nk-cells)
#CLUSTER 3: KLRF1 Nk-cells (NETTAMENTE)
#CLUSTER 4: IGHG2 B-cells
#CLUSTER 5: B4GALT1 Monocytes
#CLUSTER 6: NELL2 Excitatory neurons
#CLUSTER 7: LYZ Langerhans cells
#CLUSTER 8: TCL1A Dendritic cells
#CLUSTER 9: CLDN5 Adipocytes
#CLUSTER 1O: CXCL1 Basal respiratory cells

#the “real” markers (sorted by logFC, not by FDR) we can plot their expression with a heatmap:
FeaturePlot(PB_data, features = c("FHIT", "TNFRSF25", "KLRG1", "KLRF1", "IGHG2", "B4GALT1", 
                                    "NELL2", "LYZ", "TCL1A", "CLDN5", "CXCL1"))

#Or in single cells grouped by cluster:
PB_data.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(PB_data, features = top10$gene) + NoLegend()
#Here we can observe that some clusters seem to share part of their over-expressed genes, 
      # e.g.  0 and 1, 4 and 8

# Let us see the top markers for cluster 0, sorted by fold change and NOT by FDR:
cluster0.markers <- FindMarkers(PB_data, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster0.markers <- cluster0.markers[order(-cluster0.markers$avg_log2FC),]
head(cluster0.markers, n = 10)
 
 cluster2.markers <- FindMarkers(PB_data, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
 cluster2.markers <- cluster2.markers[order(-cluster2.markers$avg_log2FC),]
 head(cluster2.markers, n = 10)
 
 cluster0AND2.markers <- FindMarkers(PB_data, ident.1 = c(0,2), min.pct = 0.25, test.use = "wilcox")
 cluster0AND2.markers <- cluster0AND2.markers[order(-cluster0AND2.markers$avg_log2FC),]
 head(cluster0AND2.markers, n = 10)
 
 cluster20.markers <- FindMarkers(PB_data, ident.1 = 2, ident.2 = 0, min.pct = 0.25, test.use = "wilcox")
 cluster20.markers <- cluster20.markers[order(-cluster20.markers$avg_log2FC),]
 head(cluster20.markers, n = 10)
 
 cluster51.markers <- FindMarkers(PB_data, ident.1 = 1, ident.2 = 5, min.pct = 0.25, test.use = "wilcox")
 cluster51.markers <- cluster51.markers[order(-cluster51.markers$avg_log2FC),]
 head(cluster51.markers, n = 10)
 
 cluster51.markers <-
   cluster51.markers[order(cluster51.markers$avg_log2FC),]
 head(cluster51.markers, n = 10)

 VlnPlot(PB_data, features = c("FHIT", "S100A4")) #CLUSTER0  
 VlnPlot(PB_data, features = c("TNFRSF25", "S100A4")) #CLUSTER1
 VlnPlot(PB_data, features = c("KLRG1", "S100A4")) #CLUSTER2 
 VlnPlot(PB_data, features = c("KLRF1", "S100A4")) #CLUSTER3 
 VlnPlot(PB_data, features = c("IGHG2", "S100A4")) #CLUSTER4 
 VlnPlot(PB_data, features = c("B4GALT1", "S100A4")) #CLUSTER5 
 VlnPlot(PB_data, features = c("NELL2", "S100A4")) #CLUSTER6 
 VlnPlot(PB_data, features = c("LYZ", "S100A4")) #CLUSTER7 
 VlnPlot(PB_data, features = c("TCL1A", "S100A4")) #CLUSTER8 
 VlnPlot(PB_data, features = c("CLDN5", "S100A4")) #CLUSTER9
 VlnPlot(PB_data, features = c("CXCL1", "S100A4")) #CLUSTER10
 
 #CLUSTER 0: FHIT Oligodendrocytes
 #CLUSTER 1: TNFRSF25 T-cells
 #CLUSTER 2: KLRG1 Granulocytes 
 #CLUSTER 3: KLRF1 Nk-cells (NETTAMENTE)
 #CLUSTER 4: IGHG2 B-cells
 #CLUSTER 5: B4GALT1 Monocytes
 #CLUSTER 6: NELL2 Excitatory neurons
 #CLUSTER 7: LYZ Langerhans cells
 #CLUSTER 8: TCL1A Dendritic cells
 #CLUSTER 9: CLDN5 Adipocytes
 #CLUSTER 1O: CXCL1 Basal respiratory cells
 
 DotPlot(PB_data, features = c("FHIT", "TNFRSF25", "KLRG1", "KLRF1", "IGHG2", "B4GALT1",
                               "NELL2", "LYZ", "TCL1A", "CLDN5", "CXCL1"), main='DotPlot: expression of the markers')
 
 new.cluster.ids <- c("Oligodendrocytes", "T-cells", "Granulocytes", "Nk-cells", 
                      "B-cells", "Monocytes", "Excitatory neurons", "Langerhans cells", 
                      "Dendritic cells", "Adipocytes", "Basal respiratory cells")
 names(new.cluster.ids) <- levels(PB_data)
 PB_data <- RenameIdents(PB_data, new.cluster.ids)
 DimPlot(PB_data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
 DimPlot(PB_data, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
 
 DimPlot(PB_data, reduction = "umap", label = TRUE, pt.size = 0.5) 
 DimPlot(PB_data, reduction = "tsne", label = TRUE, pt.size = 0.5) 