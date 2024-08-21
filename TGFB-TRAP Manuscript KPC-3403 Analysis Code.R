######################

# Manuscript: Targeting heterogeneous tumor microenvironments in pancreatic cancer mouse models of metastasis by TGFB depletion
# Author: Heng-Chung Kung
# Code for KPC-3403 Lung Metastasis Model Single Nuclear RNA Sequencing Analysis
# For All Cell Types, CAFs, Lymphoid Cells, Myeloid Cells

######################


# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(Matrix)
library(ggpubr)

# Function to read in .h5 file after being processed by CellBender
ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (hdf5r::existsGroup(infile, 'matrix')) {
    # cellranger version 3
    message('CellRanger version 3+ format H5')
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    message('CellRanger version 2 format H5')
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}


setwd("~/Desktop/LZ01JHU521_IncludeIntrons/LZ01JHU521_CellBender_Output/")

# Sample Shortened Names
# Control: sn1_4
# TGFB-TRAP: sn2_4
# anti-IL-1B: sn3_4
# anti-PD-1: sn4_4
# anti-PD-1 + TGFB-TRAP: sn5_4
# anti-PD-1 + anti-Il-1B: sn6_4

# Read in filtered file from CellBender for Control (sn1_4)
sn1_4 <- hdf5r::H5File$new('./7_1_out_filtered.h5', mode = 'r')

sn1_4 <- Matrix::sparseMatrix(
  i = sn1_4[['matrix/indices']][],
  p = sn1_4[['matrix/indptr']][],
  x = sn1_4[['matrix/data']][],
  dimnames = list(
    sn1_4[['matrix/features/name']][],
    sn1_4[['matrix/barcodes']][]
  ),
  dims = sn1_4[['matrix/shape']][],
  index1 = FALSE
)


# Read in filtered file from CellBender for TGFB-TRAP (sn2_4)
sn2_4 <- hdf5r::H5File$new('./8_2_out_filtered.h5', mode = 'r')

sn2_4 <- Matrix::sparseMatrix(
  i = sn2_4[['matrix/indices']][],
  p = sn2_4[['matrix/indptr']][],
  x = sn2_4[['matrix/data']][],
  dimnames = list(
    sn2_4[['matrix/features/name']][],
    sn2_4[['matrix/barcodes']][]
  ),
  dims = sn2_4[['matrix/shape']][],
  index1 = FALSE
)


# Read in filtered file from CellBender for anti-IL-1B (sn3_4)
sn3_4 <- hdf5r::H5File$new('./9_1_out_filtered.h5', mode = 'r')

sn3_4 <- Matrix::sparseMatrix(
  i = sn3_4[['matrix/indices']][],
  p = sn3_4[['matrix/indptr']][],
  x = sn3_4[['matrix/data']][],
  dimnames = list(
    sn3_4[['matrix/features/name']][],
    sn3_4[['matrix/barcodes']][]
  ),
  dims = sn3_4[['matrix/shape']][],
  index1 = FALSE
)


# Read in filtered file from CellBender for antiPD-1 (sn4_4)
sn4_4 <- hdf5r::H5File$new('./10_1_out_filtered.h5', mode = 'r')

sn4_4 <- Matrix::sparseMatrix(
  i = sn4_4[['matrix/indices']][],
  p = sn4_4[['matrix/indptr']][],
  x = sn4_4[['matrix/data']][],
  dimnames = list(
    sn4_4[['matrix/features/name']][],
    sn4_4[['matrix/barcodes']][]
  ),
  dims = sn4_4[['matrix/shape']][],
  index1 = FALSE
)


# Read in filtered file from CellBender for anti-PD-1 + TGFB-TRAP (sn5_4)
sn5_4 <- hdf5r::H5File$new('./11_1_out_filtered.h5', mode = 'r')

sn5_4 <- Matrix::sparseMatrix(
  i = sn5_4[['matrix/indices']][],
  p = sn5_4[['matrix/indptr']][],
  x = sn5_4[['matrix/data']][],
  dimnames = list(
    sn5_4[['matrix/features/name']][],
    sn5_4[['matrix/barcodes']][]
  ),
  dims = sn5_4[['matrix/shape']][],
  index1 = FALSE
)


# Read in filtered file from CellBender for anti-PD-1 + anti-IL1B (sn6_4)
sn6_4 <- hdf5r::H5File$new('./12_4_out_filtered.h5', mode = 'r')

sn6_4 <- Matrix::sparseMatrix(
  i = sn6_4[['matrix/indices']][],
  p = sn6_4[['matrix/indptr']][],
  x = sn6_4[['matrix/data']][],
  dimnames = list(
    sn6_4[['matrix/features/name']][],
    sn6_4[['matrix/barcodes']][]
  ),
  dims = sn6_4[['matrix/shape']][],
  index1 = FALSE
)



# Initialize the Seurat object with the raw (non-normalized) data.
sn1_4 <- CreateSeuratObject(counts = sn1_4, project = "sn1_4", min.cells = 3, min.features = 200)
sn2_4 <- CreateSeuratObject(counts = sn2_4, project = "sn2_4", min.cells = 3, min.features = 200)
sn3_4 <- CreateSeuratObject(counts = sn3_4, project = "sn3_4", min.cells = 3, min.features = 200)
sn4_4 <- CreateSeuratObject(counts = sn4_4, project = "sn4_4", min.cells = 3, min.features = 200)
sn5_4 <- CreateSeuratObject(counts = sn5_4, project = "sn5_4", min.cells = 3, min.features = 200)
sn6_4 <- CreateSeuratObject(counts = sn6_4, project = "sn6_4", min.cells = 3, min.features = 200)


# set working directory to drive
setwd("~/Desktop/LZ01JHU521_IncludeIntrons/")

# Merge data (merge all the samples within one batch together into one Seurat object called PDAC_combined)
PDAC_combined <- merge(sn1_4, y = c(sn2_4, sn3_4, sn4_4, sn5_4, sn6_4), add.cell.ids = c("1_4", "2_4", "3_4", "4_4", "5_4", "6_4"), project = "520_Combined")
#Display
PDAC_combined

# Remove counts for subsequent analysis given overlap with rRNA overlap based on https://www.nature.com/articles/s41467-021-27035-8
counts <- GetAssayData(PDAC_combined, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('Gm42418','AY036118'))),]
PDAC_combined <- subset(PDAC_combined, features = rownames(counts))

# QC-ing
PDAC_combined[["percent.mt"]] <- PercentageFeatureSet(PDAC_combined, pattern = "^MT-")
head(PDAC_combined@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(PDAC_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize QC metrics as scatter plots: RNA counts vs mitochondrial% and RNA counts vs features
plot1 <- FeatureScatter(PDAC_combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PDAC_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 

# Filter out cells with cells with unique feature counts <200 and >2,500 AND >10% mitochondrial counts
filtered_PDAC_combined <- subset(PDAC_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)


# Normalizing the Data
norm_PDAC_combined <- NormalizeData(filtered_PDAC_combined)
# Identify highly variable features, feature selection
varfeat_PDAC_combined <- FindVariableFeatures(norm_PDAC_combined, selection.method = "vst", nfeatures = 2000)

# Scale the data
# Shifts expression such that the mean expression across all cells is 0
# Scales data so the variance across cells is 1
all.genes <- rownames(varfeat_PDAC_combined)
scaled_PDAC_combined<- ScaleData(varfeat_PDAC_combined, features = all.genes)

# Linear Dimensional Reduction of Scaled Data
scaled_PDAC_combined <- RunPCA(scaled_PDAC_combined, features = VariableFeatures(object = scaled_PDAC_combined))
# View the key differential features of first five principal components (5 pos enriched, 5 neg), of course you can change this
print(scaled_PDAC_combined[["pca"]], dims = 1:5, nfeatures = 5)
# Visualize top genes associated with reduction components
VizDimLoadings(scaled_PDAC_combined, dims = 1:5, reduction = "pca")
# View 2D scatter plot of desired PCs
DimPlot(scaled_PDAC_combined, reduction = "pca")
# Explore correlated feature sets via heatmap. This code will display heatmaps for the first 15 PCs but you can adjust to just 1 (dim = x) if desired
DimHeatmap(y, dims = 1:20, cells = 500, balanced = TRUE)

# Determine the dimensionality of data. What PCs are most statistically significant in their variation? 
PDAC_combined <- JackStraw(scaled_PDAC_combined, dims = 20)
PDAC_combined <- ScoreJackStraw(PDAC_combined, dims = 1:20)

#Visualize JackStraw plot. Distribution of p-values for each PC. Dashed line is a uniform distribubtion.
JackStrawPlot(PDAC_combined, dims = 1:20)
# Compare JackStraw Plot to Elbow plot to see where most variation arises. As plots level off, there is less variation coming from these populations. 
ElbowPlot(PDAC_combined, n = 20)

# Cluster the cells based on dimensions of interest from JackStraw and ElbowPlot above
PDAC_combined <- FindNeighbors(PDAC_combined, dims = 1:9)
PDAC_combined <- FindClusters(PDAC_combined, resolution = 1.2)


# Run non-linear dimensional reduction. Create UMAP/tSNE. Specifiy dimensionality based on significant PCs. 
PDAC_combined <- RunUMAP(PDAC_combined, dims = 1:9)
# Visualize UMAP
DimPlot(PDAC_combined, reduction = "umap", label = TRUE)
DimPlot(PDAC_combined, label = TRUE, ncol = 3, split.by = "orig.ident")

# Save the data here as an RDS file. This can be loaded back into R or shared with others for analysis at this step without having to redo the computationally-heavy previous steps. 
saveRDS(PDAC_combined, file = "~/Desktop/LZ01JHU521_IncludeIntrons/LZ01JHU521_Combined_CellBender_rRNA.rds")

# Identify the differentially expressed features for each cluster. 
Idents(PDAC_combined) <- "seurat_clusters"
PDAC_combined.markers <- FindAllMarkers(PDAC_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PDAC_combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
# Save all markers as a spreadsheet for easier analysis
write.csv (PDAC_combined.markers %>% group_by(cluster), file ="~/Desktop/LZ01JHU521_IncludeIntrons/LZ01JHU521_Combined_Markers_nFeature_Cellbender_rRNA.csv", quote =TRUE)



# Visualizing expression of different markers for specific cell types with UMAPs and Violin plots
# Acinar
FeaturePlot(PDAC_combined, features = c("Prss2", "Ctrb1", "Try5", "Pnliprp1"))
VlnPlot(PDAC_combined, features = c("Prss2", "Ctrb1", "Try5", "Pnliprp1"))


# CAF
FeaturePlot(PDAC_combined, features = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Vim", "Fap", "Pdpn"))
VlnPlot(PDAC_combined, features = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Vim", "Fap", "Pdpn"))

# myCAF
FeaturePlot(PDAC_combined, features = c("Mmp11", "Myl9", "Hopx", "Postn", "Tpm1", "Thy1", "Col12a1", "Acta2" ))
VlnPlot(PDAC_combined, features = c("Mmp11", "Myl9", "Hopx", "Postn", "Tpm1", "Thy1", "Col12a1", "Acta2"))

# iCAF
FeaturePlot(PDAC_combined, features = c("Dpt", "Lmna", "Agtr1a", "Has1", "Cxcl1", "Cxcl2", "Ccl2", "Clec3b", "Col14a1", "Il6", "Il6ra"))
VlnPlot(PDAC_combined, features = c("Dpt", "Lmna", "Agtr1a", "Has1", "Il6ra"))
VlnPlot(PDAC_combined, features = c("Cxcl1", "Cxcl2", "Ccl2", "Clec3b", "Col14a1", "Il6"))

# apCAF
FeaturePlot(PDAC_combined, features = c("H2-Ab1", "Cd74", "Saa3", "Slpi"))
VlnPlot(PDAC_combined, features = c("H2-Ab1", "Cd74", "Saa3", "Slpi"))

# Schwann
FeaturePlot(PDAC_combined, features = c("Sox10", "S100b"))
VlnPlot(PDAC_combined, features = c("Sox10", "S100b"))

# Ductal Malignant
FeaturePlot(PDAC_combined, features = c("Krt19", "Krt7", "Krt17", "Slpi", "Epcam", "Tspan8"))
FeaturePlot(PDAC_combined, features = c("Krt18", "Sox9", "Krt8", "Krt19", "Krt7", "Epcam"))
VlnPlot(PDAC_combined, features = c("Krt18", "Sox9", "Krt8", "Krt19", "Krt7", "Epcam"))

#Epithelial Normal
FeaturePlot(PDAC_combined, features = c("Cftr"))
VlnPlot(PDAC_combined, features = c('Cftr'))

# Endothelial
FeaturePlot(PDAC_combined, features = c("Cdh5", "Plvap", "Vwf", "Cldn5", "Pecam1"))
VlnPlot(PDAC_combined, features = c("Cdh5", "Plvap", "Vwf", "Cldn5", "Pecam1"))

# Adipocytes
FeaturePlot(PDAC_combined, features = c("Plin1", "Lpl", "Pparg", "Fabp4", "Adipoq"))

#Hepatoctye 
FeaturePlot(PDAC_combined, features = c("Alb", "Apoa1", "Pck1"))
VlnPlot(PDAC_combined, features = c("Alb", "Apoa1", "Pck1"))

#Acinar
FeaturePlot(renamed_PDAC_combined, features = c("Prss2", "Ctrb1", "Try5", "Pnliprp1"))
VlnPlot(renamed_PDAC_combined, features = c("Prss2", "Ctrb1", "Try5", "Pnliprp1"))

# Immune 
FeaturePlot(PDAC_combined, features = c("Ptprc"))
VlnPlot(PDAC_combined, features = c("Ptprc"))

# B cell
FeaturePlot(PDAC_combined, features = c("Ms4a1", "Cd79a", "Cd79b", "Cd19", "Bank1"))
VlnPlot(PDAC_combined, features = c("Ms4a1", "Cd79a", "Cd79b", "Cd19", "Bank1"))

# NK Cell
FeaturePlot(PDAC_combined, features = c("Plin1", "Lpl"))
VlnPlot(PDAC_combined, features = c("Plin1", "Lpl"))

#Lymphoid T Cell
FeaturePlot(PDAC_combined, features = c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Gzma", "Foxp3", "Cd96", "Themis"))
VlnPlot(PDAC_combined, features = c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Gzma", "Foxp3", "Cd96", "Themis"))
FeaturePlot(PDAC_combined, features = c("Nr4a3", "Prdm1", "C1qbp"))
VlnPlot(PDAC_combined, features = c("Nr4a3", "Prdm1", "C1qbp"))

#Myeloid cell
FeaturePlot(PDAC_combined, features = c("Cd14", "Fcgr2b", "Cd68", "Cd52", "Cd74", "Itgam"))
VlnPlot(PDAC_combined, features = c("Cd14", "Fcgr2b", "Cd68", "Cd52", "Cd74", "Itgam"))
FeaturePlot(PDAC_combined, features = c("Trem1", "Trem2", "Marco", "C1qa", "C1qb"))
VlnPlot(PDAC_combined, features = c("Trem1", "Trem2", "Marco", "C1qa", "C1qb"))




# Label clusters
new.cluster.ids <- c("Ductal", "Ductal", "Myeloid", "3", "Ductal", "Fibroblast", "Ductal", "Ductal",
                     "Ductal", "Fibroblast", "Lymphoid", "Endothelial", "12","Lymphoid","Myeloid", "Ductal", "Fibroblast", "Fibroblast", "Adipocyte")
names(new.cluster.ids) <- levels(PDAC_combined)
renamed_PDAC_combined <- RenameIdents(PDAC_combined, new.cluster.ids)
renamed_PDAC_combined$cell_type <- Idents(renamed_PDAC_combined)
# Figure S15A
DimPlot(renamed_PDAC_combined, reduction = "umap", cols = c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "darkblue", "orange", "purple"),order = c("12","3","Adipocyte", "Endothelial", "Myeloid", "Lymphoid", "Fibroblast", "Ductal")) + ggtitle("KPC-3403 (All Samples)") + theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Desktop/LZ01JHU521_IncludeIntrons/Figures for Publication/KPC3404 (All Treatments).png", width = 25, height = 20, units = "cm")




#################################

# CAF Analysis

LZ01JHU521_Fibroblast <- subset(renamed_PDAC_combined, idents = c("Fibroblast"))

# Display object
LZ01JHU521_Fibroblast

# Recluster CAFs
varfeat_CAF_combined <- FindVariableFeatures(LZ01JHU521_Fibroblast, selection.method = "vst", nfeatures = 2000)

# Scale the data
# Shifts expression such that the mean expression across all cells is 0
# Scales data so the variance across cells is 1
CAF.all.genes <- rownames(varfeat_CAF_combined)
scaled_CAF_combined<- ScaleData(varfeat_CAF_combined, features = CAF.all.genes)

# Linear Dimensional Reduction of Scaled Data
scaled_CAF_combined <- RunPCA(scaled_CAF_combined, features = VariableFeatures(object = scaled_CAF_combined))

# Determine the dimensionality of data. What PCs are most statistically significant in their variation? 
CAF_combined <- JackStraw(scaled_CAF_combined, dims = 20)
CAF_combined <- ScoreJackStraw(CAF_combined, dims = 1:20)

#Visualize JackStraw plot. Distribution of p-values for each PC. Dashed line is a uniform distribubtion.
JackStrawPlot(CAF_combined, dims = 1:20)
# Compare JackStraw Plot to Elbow plot to see where most variation arises. As plots level off, there is less variation coming from these populations. 
ElbowPlot(CAF_combined, ndims = 30)

# Cluster the cells
# Determine the dimensions of interest from JackStraw and ElbowPlot above, and add to code below.
CAF_combined <- FindNeighbors(CAF_combined, dims = 1:8)
CAF_combined <- FindClusters(CAF_combined, resolution = 0.5)
# Cluster IDs for the first 5 cells
head(Idents(CAF_combined), 5)

# Run non-linear dimensional reduction. Create UMAP/tSNE. Specifiy dimensionality based on significant PCs. 
CAF_combined <- RunUMAP(CAF_combined, dims = 1:8)
# Visualize UMAP
DimPlot(CAF_combined, reduction = "umap", label = TRUE)


# Find all markers for every cluster compared to all remaining cells. This will report only the positive values. 
CAF_combined.markers <- FindAllMarkers(CAF_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CAF_combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


# Fibroblast Markers
#Fibroblast
FeaturePlot(CAF_combined, features = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Vim", "Fap", "Pdpn"))
VlnPlot(CAF_combined, features = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Vim", "Fap", "Pdpn"))

#myCAF from Tuveson et al. and Zheng et al. 
FeaturePlot(CAF_combined, features = c("Acta2", "Myl9", "Tpm1", "Postn", "Tagln", "Thy1", "Col12a1", "Mmp11", "Hopx", "Igfbp3","Thbs2", "Tpm2", "Ctgf" ))
VlnPlot(CAF_combined, features = c("Mmp11", "Myl9", "Hopx", "Postn", "Tpm1", "Thy1", "Col12a1", "Acta2", "Tagln", "Igfbp3","Thbs2", "Tpm2", "Ctgf"))

FeaturePlot(CAF_combined, features = c("Gjb2", "Col2a1", "Pkp1", "Cdh10", "Crabp1", "Gjb3", "Dock8", "Acsbg1", "Itih2", "Shank2", "Car6", "Klhdc8a", "Bcat1", "Col8a1" ))
VlnPlot(CAF_combined, features = c("Gjb2", "Col2a1", "Pkp1", "Cdh10", "Crabp1", "Gjb3", "Dock8", "Acsbg1", "Itih2", "Shank2", "Car6", "Klhdc8a", "Bcat1", "Col8a1"))

FeaturePlot(CAF_combined, features = c("Cdh3", "Lipg", "Serpine2", "Col15a1", "Tnc", "Sparcl1", "Spp1", "Tgfb1", "Crlf1",  "Col12a1", "Sdc1", "Cthrc1"))
VlnPlot(CAF_combined, features = c("Cdh3", "Lipg", "Serpine2", "Col15a1", "Tnc", "Sparcl1", "Spp1", "Tgfb1", "Crlf1",  "Col12a1", "Sdc1", "Cthrc1"))

#iCAF from Tuveson et al. and Zheng et al. 
FeaturePlot(CAF_combined, features = c("Dpt", "Lmna", "Agtr1a", "Has1", "Has2", "Cxcl1", "Cxcl2", "Cxcl12", "Ccl2", "Clec3b", "Col14a1", "Il6", "Ly6c1", "Cfd", "Csf3", "Lif"))
VlnPlot(CAF_combined, features = c("Dpt", "Lmna", "Agtr1a", "Has1", "Has2", "Cxcl1", "Cxcl2", "Cxcl12", "Ccl2", "Clec3b", "Col14a1", "Il6", "Ly6c1", "Cfd", "Csf3", "Lif"))

FeaturePlot(CAF_combined, features = c("Il21", "Il1a", "Cxcl3", "H2-Q10", "Il11", "Prdm1", "Pnpla3", "Ism1", "Ifi205", "Gsn", "Svep1", "Plpp3", "Scara3", "Pcolce2", "Efemp1", "Ly6a","Spink5", "Adamts5"))
VlnPlot(CAF_combined, features = c("Il21", "Il1a", "Cxcl3", "H2-Q10", "Il11", "Prdm1", "Pnpla3", "Ism1", "Ifi205", "Gsn", "Svep1", "Plpp3", "Scara3", "Pcolce2", "Efemp1", "Ly6a", "Spink5", "Adamts5"))

FeaturePlot(CAF_combined, features = c("S100a8", "Saa1", "Il13ra2", "Hp", "Il1a", "Gpr84", "Saa4", "Ndufa4l2", "Nos2", "Ism1", "Ifi27l2a", "Tnxb"))
VlnPlot(CAF_combined, features = c("S100a8", "Saa1", "Il13ra2", "Hp", "Il1a", "Gpr84", "Saa4", "Ndufa4l2", "Nos2", "Ism1", "Ifi27l2a", "Tnxb"))

FeaturePlot(CAF_combined, features = c("Scara5", "Sfrp2", "S100a10", "S100a4", "Fbln2", "Gpx3", "Pnrc1", "Fbn1", "Sgk1", "Adh1", "Ccdc80", "Abl2"))
VlnPlot(CAF_combined, features = c("Scara5", "Sfrp2", "S100a10", "S100a4", "Fbln2", "Gpx3", "Pnrc1", "Fbn1", "Sgk1", "Adh1", "Ccdc80", "Abl2"))

FeaturePlot(CAF_combined, features = c("C4b", "C3"))
VlnPlot(CAF_combined, features = c("C4b", "C3"))


# Isolate TGFB-TRAP samples 
Idents(CAF_combined) <- "orig.ident"
CAF_combined_TGFB <- subset(CAF_combined,  idents = c("sn1_4", "sn2_4", "sn4_4", "sn5_4"))
Idents(CAF_combined_TGFB) <- "seurat_clusters"
DimPlot(CAF_combined_TGFB) 

# Organize into 3 main clusters: CAF_c1, CAF_c2, CAF_c3
CAF.new.cluster.ids <- c("CAF_c2", "CAF_c1", "CAF_c1", "CAF_c3", "CAF_c3", "CAF_c1", "CAF_c2", "CAF_c2", "CAF_c2")
names(CAF.new.cluster.ids) <- levels(CAF_combined)
renamed_CAF_combined <- RenameIdents(CAF_combined, CAF.new.cluster.ids)
# Display labeled UMAP
# Figure S15B
DimPlot(renamed_CAF_combined, reduction = "umap", label = FALSE, pt.size = 0.5, order = c("CAF_c3", "CAF_c2", "CAF_c1"))



# CAF Phenotype myCAF and iCAF Violin Plots (Figure S15C-Q)

library(tidyr)
library(forcats)

Idents(CAF_combined_TGFB) <- "seurat_clusters"
DimPlot(CAF_combined_TGFB, label = TRUE)
DimPlot(CAF_combined_TGFB, split.by = "orig.ident",ncol = 2, label = TRUE)

# All TGFB-TRAP Samples
Idents(CAF_combined_TGFB) <- "seurat_clusters"

# All TGFB-TRAP Samples: CAF_c1
CAF_all_TGFB_c1 <- subset(renamed_CAF_combined, idents = c("CAF_c1"))
CAF_all_TGFB_c1_iCAF_genes = data.frame(t(data.frame(CAF_all_TGFB_c1[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_all_TGFB_c1_iCAF_genes <- CAF_all_TGFB_c1_iCAF_genes %>% gather(Genes, Expression) 
CAF_all_TGFB_c1_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_all_TGFB_c1_myCAF_genes = data.frame(t(data.frame(CAF_all_TGFB_c1[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_all_TGFB_c1_myCAF_genes <- CAF_all_TGFB_c1_myCAF_genes %>% gather(Genes, Expression) 
CAF_all_TGFB_c1_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_all_TGFB_c1_combined_df <- rbind(CAF_all_TGFB_c1_iCAF_genes, CAF_all_TGFB_c1_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_all_TGFB_c1_combined_df[, "Expression"])) / 100000
CAF_all_TGFB_c1_combined_df$Expression <- CAF_all_TGFB_c1_combined_df$Expression + abs(noise)
ggplot(CAF_all_TGFB_c1_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 23, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 All TGFB-TRAP Samples Cluster 1 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF all TGFB-TRAP c1 Genes.png", width = 30, height = 20, units = "cm")

# All TGFB-TRAP Samples: CAF_c2
Idents(CAF_combined_TGFB) <- "seurat_clusters"
CAF_all_TGFB_c2 <- subset(renamed_CAF_combined, idents = c("CAF_c2"))
CAF_all_TGFB_c2_iCAF_genes = data.frame(t(data.frame(CAF_all_TGFB_c2[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_all_TGFB_c2_iCAF_genes <- CAF_all_TGFB_c2_iCAF_genes %>% gather(Genes, Expression) 
CAF_all_TGFB_c2_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_all_TGFB_c2_myCAF_genes = data.frame(t(data.frame(CAF_all_TGFB_c2[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_all_TGFB_c2_myCAF_genes <- CAF_all_TGFB_c2_myCAF_genes %>% gather(Genes, Expression) 
CAF_all_TGFB_c2_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_all_TGFB_c2_combined_df <- rbind(CAF_all_TGFB_c2_iCAF_genes, CAF_all_TGFB_c2_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_all_TGFB_c2_combined_df[, "Expression"])) / 100000
CAF_all_TGFB_c2_combined_df$Expression <- CAF_all_TGFB_c2_combined_df$Expression + abs(noise)
ggplot(CAF_all_TGFB_c2_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 23, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 All TGFB-TRAP Samples Cluster 2 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF all TGFB-TRAP c2 Genes.png", width = 30, height = 20, units = "cm")

# All TGFB-TRAP Samples: CAF_c3
Idents(CAF_combined_TGFB) <- "seurat_clusters"
CAF_all_TGFB_c3 <- subset(renamed_CAF_combined_TGFB, idents = c("CAF_c3"))
CAF_all_TGFB_c3_iCAF_genes = data.frame(t(data.frame(CAF_all_TGFB_c3[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_all_TGFB_c3_iCAF_genes <- CAF_all_TGFB_c3_iCAF_genes %>% gather(Genes, Expression) 
CAF_all_TGFB_c3_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_all_TGFB_c3_myCAF_genes = data.frame(t(data.frame(CAF_all_TGFB_c3[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_all_TGFB_c3_myCAF_genes <- CAF_all_TGFB_c3_myCAF_genes %>% gather(Genes, Expression) 
CAF_all_TGFB_c3_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_all_TGFB_c3_combined_df <- rbind(CAF_all_TGFB_c3_iCAF_genes, CAF_all_TGFB_c3_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_all_TGFB_c3_combined_df[, "Expression"])) / 100000
CAF_all_TGFB_c3_combined_df$Expression <- CAF_all_TGFB_c3_combined_df$Expression + abs(noise)
ggplot(CAF_all_TGFB_c3_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 23, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 All TGFB-TRAP Samples Cluster 3 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF all TGFB-TRAP c3 Genes.png", width = 30, height = 20, units = "cm")


# Control (sn1_4)
# Control: CAF_c1
CAF_sn1_4 <- subset(CAF_combined, idents = "sn1_4")
Idents(CAF_sn1_4) <- "seurat_clusters"
DimPlot(CAF_sn1_4, label = T)
CAF_sn1_c1 <- subset(CAF_sn1_4, idents = c("1", "2", "5"))
CAF_sn1_c1_iCAF_genes = data.frame(t(data.frame(CAF_sn1_c1[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn1_c1_iCAF_genes <- CAF_sn1_c1_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c1_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn1_c1_myCAF_genes = data.frame(t(data.frame(CAF_sn1_c1[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn1_c1_myCAF_genes <- CAF_sn1_c1_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c1_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn1_c1_combined_df <- rbind(CAF_sn1_c1_iCAF_genes, CAF_sn1_c1_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn1_c1_combined_df[, "Expression"])) / 100000
CAF_sn1_c1_combined_df$Expression <- CAF_sn1_c1_combined_df$Expression + abs(noise)
ggplot(CAF_sn1_c1_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 Control Cluster 1 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn1_4 c1 Genes.png", width = 30, height = 20, units = "cm")

# Control: CAF_c2
Idents(CAF_sn1_4) <- "seurat_clusters"
CAF_sn1_c2 <- subset(CAF_sn1_4, idents = c( "0", "6"))
CAF_sn1_c2_iCAF_genes = data.frame(t(data.frame(CAF_sn1_c2[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn1_c2_iCAF_genes <- CAF_sn1_c2_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c2_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn1_c2_myCAF_genes = data.frame(t(data.frame(CAF_sn1_c2[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn1_c2_myCAF_genes <- CAF_sn1_c2_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c2_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn1_c2_combined_df <- rbind(CAF_sn1_c2_iCAF_genes, CAF_sn1_c2_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn1_c2_combined_df[, "Expression"])) / 100000
CAF_sn1_c2_combined_df$Expression <- CAF_sn1_c2_combined_df$Expression + abs(noise)
ggplot(CAF_sn1_c2_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 Control Cluster 2 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn1_4 c2 Genes.png", width = 30, height = 20, units = "cm")

# Control: CAF_c3
Idents(CAF_sn1_4) <- "seurat_clusters"
CAF_sn1_c3 <- subset(CAF_sn1_4, idents = c("4"))
CAF_sn1_c3_iCAF_genes = data.frame(t(data.frame(CAF_sn1_c3[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn1_c3_iCAF_genes <- CAF_sn1_c3_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c3_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn1_c3_myCAF_genes = data.frame(t(data.frame(CAF_sn1_c3[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn1_c3_myCAF_genes <- CAF_sn1_c3_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c3_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn1_c3_combined_df <- rbind(CAF_sn1_c3_iCAF_genes, CAF_sn1_c3_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn1_c3_combined_df[, "Expression"])) / 100000
CAF_sn1_c3_combined_df$Expression <- CAF_sn1_c3_combined_df$Expression + abs(noise)
ggplot(CAF_sn1_c3_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 Control Cluster 3 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn1_4 c3 Genes.png", width = 30, height = 20, units = "cm")



# TGFB-TRAP (sn2_4)
# TGFB-TRAP CAF_c1
CAF_sn2_4 <- subset(CAF_combined, idents = "sn2_4")
Idents(CAF_sn2_4) <- "seurat_clusters"
DimPlot(CAF_sn2_4, label = T)
CAF_sn2_c1 <- subset(CAF_sn2_4, idents = c("1", "2", "5"))
CAF_sn2_c1_iCAF_genes = data.frame(t(data.frame(CAF_sn2_c1[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn2_c1_iCAF_genes <- CAF_sn2_c1_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c1_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn2_c1_myCAF_genes = data.frame(t(data.frame(CAF_sn2_c1[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn2_c1_myCAF_genes <- CAF_sn2_c1_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c1_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn2_c1_combined_df <- rbind(CAF_sn2_c1_iCAF_genes, CAF_sn2_c1_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn2_c1_combined_df[, "Expression"])) / 100000
CAF_sn2_c1_combined_df$Expression <- CAF_sn2_c1_combined_df$Expression + abs(noise)
ggplot(CAF_sn2_c1_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 TGFB-TRAP Cluster 1 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn2_4 c1 Genes.png", width = 30, height = 20, units = "cm")

# TGFB-TRAP CAF_c2
Idents(CAF_sn2_4) <- "seurat_clusters"
CAF_sn2_c2 <- subset(CAF_sn2_4, idents = c("0", "6","7","8"))
CAF_sn2_c2_iCAF_genes = data.frame(t(data.frame(CAF_sn2_c2[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn2_c2_iCAF_genes <- CAF_sn2_c2_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c2_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn2_c2_myCAF_genes = data.frame(t(data.frame(CAF_sn2_c2[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn2_c2_myCAF_genes <- CAF_sn2_c2_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c2_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn2_c2_combined_df <- rbind(CAF_sn2_c2_iCAF_genes, CAF_sn2_c2_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn2_c2_combined_df[, "Expression"])) / 100000
CAF_sn2_c2_combined_df$Expression <- CAF_sn2_c2_combined_df$Expression + abs(noise)
ggplot(CAF_sn2_c2_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 TGFB-TRAP Cluster 2 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn2_4 c2 Genes.png", width = 30, height = 20, units = "cm")

# TGFB-TRAP CAF_c3
Idents(CAF_sn2_4) <- "seurat_clusters"
CAF_sn2_c3 <- subset(CAF_sn2_4, idents = c("3", "4"))
CAF_sn2_c3_iCAF_genes = data.frame(t(data.frame(CAF_sn2_c3[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn2_c3_iCAF_genes <- CAF_sn2_c3_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c3_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn2_c3_myCAF_genes = data.frame(t(data.frame(CAF_sn2_c3[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn2_c3_myCAF_genes <- CAF_sn2_c3_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c3_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn2_c3_combined_df <- rbind(CAF_sn2_c3_iCAF_genes, CAF_sn2_c3_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn2_c3_combined_df[, "Expression"])) / 100000
CAF_sn2_c3_combined_df$Expression <- CAF_sn2_c3_combined_df$Expression + abs(noise)
ggplot(CAF_sn2_c3_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 TGFB-TRAP Cluster 3 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn2_4 c3 Genes.png", width = 30, height = 20, units = "cm")



# anti-PD-1 (sn4_4)
CAF_sn4_4 <- subset(CAF_combined, idents = "sn4_4")
# anti-PD-1 CAF_c1
Idents(CAF_sn4_4) <- "seurat_clusters"
DimPlot(CAF_sn4_4, label = TRUE)
CAF_sn4_c1 <- subset(CAF_sn4_4, idents = c("1", "2", "5"))
CAF_sn4_c1_iCAF_genes = data.frame(t(data.frame(CAF_sn4_c1[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn4_c1_iCAF_genes <- CAF_sn4_c1_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c1_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn4_c1_myCAF_genes = data.frame(t(data.frame(CAF_sn4_c1[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn4_c1_myCAF_genes <- CAF_sn4_c1_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c1_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn4_c1_combined_df <- rbind(CAF_sn4_c1_iCAF_genes, CAF_sn4_c1_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn4_c1_combined_df[, "Expression"])) / 100000
CAF_sn4_c1_combined_df$Expression <- CAF_sn4_c1_combined_df$Expression + abs(noise)
ggplot(CAF_sn4_c1_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 anti-PD1 Cluster 1 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn4_4 c1 Genes.png", width = 30, height = 20, units = "cm")

# anti-PD-1 CAF_c2
Idents(CAF_sn4_4) <- "seurat_clusters"
CAF_sn4_c2 <- subset(CAF_sn4_4, idents = c("0", "6"))
CAF_sn4_c2_iCAF_genes = data.frame(t(data.frame(CAF_sn4_c2[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn4_c2_iCAF_genes <- CAF_sn4_c2_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c2_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn4_c2_myCAF_genes = data.frame(t(data.frame(CAF_sn4_c2[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn4_c2_myCAF_genes <- CAF_sn4_c2_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c2_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn4_c2_combined_df <- rbind(CAF_sn4_c2_iCAF_genes, CAF_sn4_c2_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn4_c2_combined_df[, "Expression"])) / 100000
CAF_sn4_c2_combined_df$Expression <- CAF_sn4_c2_combined_df$Expression + abs(noise)
ggplot(CAF_sn4_c2_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 anti-PD1 Cluster 2 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn4_4 c2 Genes.png", width = 30, height = 20, units = "cm")

# anti-PD-1 CAF_c3
Idents(CAF_sn4_4) <- "seurat_clusters"
CAF_sn4_c3 <- subset(CAF_sn4_4, idents = c("3", "4"))
CAF_sn4_c3_iCAF_genes = data.frame(t(data.frame(CAF_sn4_c3[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn4_c3_iCAF_genes <- CAF_sn4_c3_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c3_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn4_c3_myCAF_genes = data.frame(t(data.frame(CAF_sn4_c3[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn4_c3_myCAF_genes <- CAF_sn4_c3_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c3_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn4_c3_combined_df <- rbind(CAF_sn4_c3_iCAF_genes, CAF_sn4_c3_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn4_c3_combined_df[, "Expression"])) / 100000
CAF_sn4_c3_combined_df$Expression <- CAF_sn4_c3_combined_df$Expression + abs(noise)
ggplot(CAF_sn4_c3_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 anti-PD1 Cluster 3 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn4_4 c3 Genes.png", width = 30, height = 20, units = "cm")



# anti-PD-1 + TGFB-TRAP (sn5_4)
# anti-PD-1 + TGFB-TRAP CAF_c1
CAF_sn5_4 <- subset(CAF_combined, idents = c("sn5_4"))
Idents(CAF_sn5_4) <- "seurat_clusters"
DimPlot(CAF_sn5_4, label = T)
CAF_sn5_c1 <- subset(CAF_sn5_4, idents = c("1", "2", "5"))
CAF_sn5_c1_iCAF_genes = data.frame(t(data.frame(CAF_sn5_c1[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn5_c1_iCAF_genes <- CAF_sn5_c1_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c1_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn5_c1_myCAF_genes = data.frame(t(data.frame(CAF_sn5_c1[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn5_c1_myCAF_genes <- CAF_sn5_c1_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c1_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn5_c1_combined_df <- rbind(CAF_sn5_c1_iCAF_genes, CAF_sn5_c1_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn5_c1_combined_df[, "Expression"])) / 100000
CAF_sn5_c1_combined_df$Expression <- CAF_sn5_c1_combined_df$Expression + abs(noise)
ggplot(CAF_sn5_c1_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 anti-PD1 + TGFB-TRAP Cluster 1 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn5_4 c1 Genes.png", width = 30, height = 20, units = "cm")

# anti-PD-1 + TGFB-TRAP CAF_c2
Idents(CAF_sn5_4) <- "seurat_clusters"
CAF_sn5_c2 <- subset(CAF_sn5_4, idents = c("0"))
CAF_sn5_c2_iCAF_genes = data.frame(t(data.frame(CAF_sn5_c2[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn5_c2_iCAF_genes <- CAF_sn5_c2_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c2_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn5_c2_myCAF_genes = data.frame(t(data.frame(CAF_sn5_c2[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn5_c2_myCAF_genes <- CAF_sn5_c2_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c2_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn5_c2_combined_df <- rbind(CAF_sn5_c2_iCAF_genes, CAF_sn5_c2_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn5_c2_combined_df[, "Expression"])) / 100000
CAF_sn5_c2_combined_df$Expression <- CAF_sn5_c2_combined_df$Expression + abs(noise)
ggplot(CAF_sn5_c2_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 anti-PD1 + TGFB-TRAP Cluster 2 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn5_4 c2 Genes.png", width = 30, height = 20, units = "cm")

# anti-PD-1 + TGFB-TRAP CAF_c3
Idents(CAF_sn5_4) <- "seurat_clusters"
CAF_sn5_c3 <- subset(CAF_sn5_4, idents = c("3", "4"))
CAF_sn5_c3_iCAF_genes = data.frame(t(data.frame(CAF_sn5_c3[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn5_c3_iCAF_genes <- CAF_sn5_c3_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c3_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn5_c3_myCAF_genes = data.frame(t(data.frame(CAF_sn5_c3[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc"),])))
CAF_sn5_c3_myCAF_genes <- CAF_sn5_c3_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c3_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn5_c3_combined_df <- rbind(CAF_sn5_c3_iCAF_genes, CAF_sn5_c3_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3",  "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2",  "Col12a1", "Dock8", "Col15a1", "Col8a1", "Tnc")))
noise <- rnorm(n = length(x = CAF_sn5_c3_combined_df[, "Expression"])) / 100000
CAF_sn5_c3_combined_df$Expression <- CAF_sn5_c3_combined_df$Expression + abs(noise)
ggplot(CAF_sn5_c3_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-3403 anti-PD1 + TGFB-TRAP Cluster 3 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/KPC3403/KPC3403 CAF sn5_4 c3 Genes.png", width = 30, height = 20, units = "cm")


############################

# Lymphoid Analysis

# Isolate all lymphoid cells
LZ01JHU521_Lymphoid <- subset(renamed_PDAC_combined, idents = c("Lymphoid"))
DimPlot(LZ01JHU521_Lymphoid, label = TRUE)

# Recluster lymphoid cells
Lymphoid_combined <- FindVariableFeatures(LZ01JHU520_T_cells, selection.method = "vst", nfeatures = 2000)
# Scale the data
lymphoid.all.genes <- rownames(Lymphoid_combined)
Lymphoid_combined<- ScaleData(Lymphoid_combined, features = lymphoid.all.genes)
# Linear Dimensional Reduction of Scaled Data
Lymphoid_combined <- RunPCA(Lymphoid_combined, features = VariableFeatures(object = Lymphoid_combined))

# Determine the dimensionality of data. What PCs are most statistically significant in their variation? 
Lymphoid_combined <- JackStraw(Lymphoid_combined, dims = 20)
Lymphoid_combined <- ScoreJackStraw(Lymphoid_combined, dims = 1:20)

#Visualize JackStraw plot. Distribution of p-values for each PC. Dashed line is a uniform distribubtion.
JackStrawPlot(Lymphoid_combined, dims = 1:20)
# Compare JackStraw Plot to Elbow plot to see where most variation arises. As plots level off, there is less variation coming from these populations. 
ElbowPlot(Lymphoid_combined)

# Cluster the cells
# Determine the dimensions of interest from JackStraw and ElbowPlot above, and add to code below.
Lymphoid_combined <- FindNeighbors(Lymphoid_combined, dims = 1:11)
Lymphoid_combined <- FindClusters(Lymphoid_combined, resolution = 0.8)

# Run non-linear dimensional reduction. Create UMAP/tSNE. Specifiy dimensionality based on significant PCs. 
Lymphoid_combined <- RunUMAP(Lymphoid_combined, dims = 1:11)
# Visualize UMAP
DimPlot(Lymphoid_combined, reduction = "umap", label = TRUE)

# Find all differentially expressed genes for each cluster
Lymphoid_combined.markers <- FindAllMarkers(Lymphoid_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Lymphoid_combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# CD8 T cells
FeaturePlot(Lymphoid_combined, features = c("Cd8a", "Cd8b1", "Nkg7", "Gzmk"))
FeaturePlot(Lymphoid_combined, features = c("Gzma", "Gzmb", "Prf1", "Ccl5"))

# CD4 T cells and Tregs
FeaturePlot(Lymphoid_combined, features = c("Cd4", "Ctla4", "Il2ra", "Ltb", "Il7r", "Foxp3"))

# B cells
FeaturePlot(Lymphoid_combined, features = c("Ms4a1", "Cd79a", "Bank1", "Cd79b", "Cd19"))


#Labeling
DimPlot(Lymphoid_combined, label = TRUE)
Lymphoid.cluster.ids = c("CD8", "Treg", "B", "CD4", "CD4", "CD4", "6" )
names(Lymphoid.cluster.ids) <- levels(Lymphoid_combined)
renamed_Lymphoid_combined <- RenameIdents(Lymphoid_combined, Lymphoid.cluster.ids)
renamed_Lymphoid_combined$cell_type <- Idents(renamed_Lymphoid_combined)
DimPlot(renamed_Lymphoid_combined, label = FALSE, cols = c("coral2", "goldenrod2", "green3", "turquoise2", "purple"), order = rev(c("CD8", "CD4", "Treg", "B", "6"))) +xlim(-8,13)+ylim(-9,8) + ggtitle(" KPC-3403 Lymphoid (All Samples)") + theme(plot.title = element_text(hjust = 0.5))


# Isolate TGFB Samples
Idents(renamed_Lymphoid_combined) <- "orig.ident"
treatment_ids <- c("sn1_4 (Control)", "sn2_4 (TGFB-TRAP)", "sn3_4 (anti-IL1B)","sn4_4 (anti-PD1)", "sn5_4 (anti-PD1 + TGFB-TRAP)", "sn6_4 (anti-PD1 + anti-IL1B)")
names(treatment_ids) <- levels(renamed_Lymphoid_combined)
renamed_Lymphoid_combined <- RenameIdents(renamed_Lymphoid_combined, treatment_ids)
renamed_Lymphoid_combined$treatment <- Idents(object = renamed_Lymphoid_combined)
# Isolate TGFB Samples
Idents(renamed_Lymphoid_combined) <- "orig.ident"
renamed_Lymphoid_TGFB <- subset(renamed_Lymphoid_combined,  idents = c("sn1_4", "sn2_4", "sn4_4", "sn5_4"))
# TGFB Samples Only Figure S16B
Idents(renamed_Lymphoid_TGFB) <- "cell_type"
DimPlot(renamed_Lymphoid_TGFB, label = FALSE, cols = c("coral2", "goldenrod2", "green3", "turquoise2", "purple"), order = rev(c("CD8", "CD4", "Treg", "B", "6"))) +xlim(-8,13)+ylim(-9,8) + ggtitle(" KPC-3403 Lymphoid (TGFB-TRAO Samples)") + theme(plot.title = element_text(hjust = 0.5))



#########################

# Myeloid Analysis
LZ01JHU521_Myeloid <- subset(renamed_PDAC_combined, idents = c("Myeloid"))
DimPlot(LZ01JHU521_Myeloid, label = TRUE)

# Recluster Myeloid Cells
Myeloid_combined <- FindVariableFeatures(LZ01JHU520_Myeloid, selection.method = "vst", nfeatures = 2000)
# Scale the data
# Shifts expression such that the mean expression across all cells is 0
# Scales data so the variance across cells is 1
myeloid.all.genes <- rownames(Myeloid_combined)
Myeloid_combined<- ScaleData(Myeloid_combined, features = myeloid.all.genes)

# Linear Dimensional Reduction of Scaled Data
Myeloid_combined <- RunPCA(Myeloid_combined, features = VariableFeatures(object = Myeloid_combined))

# Determine the dimensionality of data. What PCs are most statistically significant in their variation? 
Myeloid_combined <- JackStraw(Myeloid_combined, dims = 20)
Myeloid_combined <- ScoreJackStraw(Myeloid_combined, dims = 1:20)

#Visualize JackStraw plot. Distribution of p-values for each PC. Dashed line is a uniform distribubtion.
JackStrawPlot(Myeloid_combined, dims = 1:20)
# Compare JackStraw Plot to Elbow plot to see where most variation arises. As plots level off, there is less variation coming from these populations. 
ElbowPlot(Myeloid_combined)

# Cluster the cells
# Determine the dimensions of interest from JackStraw and ElbowPlot above, and add to code below.
Myeloid_combined <- FindNeighbors(Myeloid_combined, dims = 1:19)
Myeloid_combined <- FindClusters(Myeloid_combined, resolution = 0.8)


# Run non-linear dimensional reduction. Create UMAP/tSNE. Specifiy dimensionality based on significant PCs. 
Myeloid_combined <- RunUMAP(Myeloid_combined, dims = 1:19)
# Visualize UMAP
DimPlot(Myeloid_combined, reduction = "umap", label = TRUE)

# Find differentially expressed genes for each cluster
Idents(Myeloid_combined) <- "seurat_clusters"
Myeloid_combined.markers <- FindAllMarkers(Myeloid_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Myeloid_combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Macrophage/Monocyte Markers
FeaturePlot(Myeloid_combined, features = c("Itgam", "Cd163", "Cd33", "Cd14", "Cd68", "Tlr2", "Mrc1", "Cd80", "Cd86", "Tgfb1", "Csf1"))
VlnPlot(Myeloid_combined, features = c("Itgam", "Cd163", "Cd33", "Cd14", "Cd68", "Tlr2", "Mrc1", "Cd80", "Cd86", "Tgfb1", "Csf1"))
FeaturePlot(Myeloid_combined, features = c("C1qa", "C1qb", "Trem2", "Apoe" ,"Mrc1", "C1qc", "H2-Eb1", "H2-Aa", "H2-Ab1", "AW112010", "Cxcl9", "Ccl8"))
FeaturePlot(Myeloid_combined, features = c("Nos2", "Tnf", "Ly6c2", "Plac8","S100a4", "Thbs1", "Chil3", "Ptgs2", "Cxcl3", "Cxcl2", "Ml2", "Cxcl1", "Serpine1", "Spp1"))

#Granulocyte Markers
FeaturePlot(Myeloid_combined, features = c("S100a8", "Fcgr4", "Fut4", "Enpp3", "Itgam", "Cd63", "Cd33", "Cd14", "Cd69"))
FeaturePlot(Myeloid_combined, features = c("Itgax", "Fcgr3", "Cxcr2", "Fcgr2b",  "Csf3r", "Itgal", "Itgb2", "Itgb1", "Itga1"))
FeaturePlot(Myeloid_combined, features = c("Cd44", "Icam1", "Ncam1", "Icam2", "Vcam1", "Mcam", "Alcam", "F11r", "Jam2"))
FeaturePlot(Myeloid_combined, features = c("Kit", "Ifngr1", "Il1r1", "Il3ra", "Il4ra", "Il6ra", "Il7r", "Pdgfra", "Pdgfrb"))


# Label Cell Types
Idents(Myeloid_combined) <- "seurat_clusters"
Myeloid.cluster.ids = c("Macrophage", "Granulocyte", "Macrophage", "Macrophage", "Granulocyte", "Macrophage")
names(Myeloid.cluster.ids) <- levels(Myeloid_combined)
renamed_Myeloid_combined <- RenameIdents(Myeloid_combined, Myeloid.cluster.ids)
DimPlot(renamed_Myeloid_combined) + ggtitle("Myeloid (All Samples)") + theme(plot.title = element_text(hjust = 0.5))
renamed_Myeloid_combined$cell_type <- Idents(object = renamed_Myeloid_combined)

# Isolate TGFB-TRAP samples 
Idents(renamed_Myeloid_combined) <- "orig.ident"
Myeloid_TGFB <- subset(renamed_Myeloid_combined,  idents = c("sn1_4", "sn2_4", "sn4_4", "sn5_4"))
Idents(Myeloid_TGFB) <- "cell_type"
# Figure S16E
DimPlot(Myeloid_TGFB)


