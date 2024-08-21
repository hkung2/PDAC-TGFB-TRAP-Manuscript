######################

# Manuscript: Targeting heterogeneous tumor microenvironments in pancreatic cancer mouse models of metastasis by TGFB depletion
# Author: Heng-Chung Kung
# Code for KPC-4545 Liver Metastasis Model Single Nuclear RNA Sequencing Analysis
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


# Sample Shortened Names
# Control: sn1_4
# TGFB-TRAP: sn2_4
# anti-IL-1B: sn3_4
# anti-PD-1: sn4_4
# anti-PD-1 + TGFB-TRAP: sn5_4
# anti-PD-1 + anti-Il-1B: sn6_4


# Read in filtered file from CellBender for Control (sn1_4)
sn1_4 <- hdf5r::H5File$new('/Volumes/ZhengData/Rex/Single Cell Data/LZ01JHU520/LZ01JHU520_000_analysis/cellranger/count_includeIntrons/1_4_LZ_5DGE/1_4_out_filtered.h5', mode = 'r')

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
sn2_4 <- hdf5r::H5File$new('/Volumes/ZhengData/Rex/Single Cell Data/LZ01JHU520/LZ01JHU520_000_analysis/cellranger/count_includeIntrons/2_4_LZ_5DGE/2_4_out_filtered.h5', mode = 'r')

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
sn3_4 <- hdf5r::H5File$new('/Volumes/ZhengData/Rex/Single Cell Data/LZ01JHU520/LZ01JHU520_000_analysis/cellranger/count_includeIntrons/3_4_LZ_5DGE/3_4_out_filtered.h5', mode = 'r')

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
sn4_4 <- hdf5r::H5File$new('/Volumes/ZhengData/Rex/Single Cell Data/LZ01JHU520/LZ01JHU520_000_analysis/cellranger/count_includeIntrons/4_4_LZ_5DGE/4_4_out_filtered.h5', mode = 'r')

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
sn5_4 <- hdf5r::H5File$new('/Volumes/ZhengData/Rex/Single Cell Data/LZ01JHU520/LZ01JHU520_000_analysis/cellranger/count_includeIntrons/5_4_LZ_5DGE/5_4_out_filtered.h5', mode = 'r')

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
sn6_4 <- hdf5r::H5File$new('/Volumes/ZhengData/Rex/Single Cell Data/LZ01JHU520/LZ01JHU520_000_analysis/cellranger/count_includeIntrons/6_4_LZ_5DGE/6_4_out_filtered.h5', mode = 'r')

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
setwd("/Volumes/ZhengData/Rex/Single Cell Analysis/LZ01JHU520_IncludeIntrons/")

# Merge data (merge all the samples within one batch together into one Seurat object called PDAC_combined)
PDAC_combined <- merge(sn1_4, y = c(sn2_4, sn3_4, sn4_4, sn5_4, sn6_4), add.cell.ids = c("1_4", "2_4", "3_4", "4_4", "5_4", "6_4"), project = "520_Combined")

# Display
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
PDAC_combined <- FindClusters(PDAC_combined, resolution = 0.8)

# Run non-linear dimensional reduction. Create UMAP
PDAC_combined <- RunUMAP(PDAC_combined, dims = 1:9)
# Visualize UMAP
DimPlot(PDAC_combined, reduction = "umap", label = TRUE)
# Save the data here as an RDS file. This can be loaded back into R or shared with others for analysis at this step without having to redo the computationally-heavy previous steps. 
saveRDS(PDAC_combined, file = "/Volumes/ZhengData/Rex/Single Cell Analysis/LZ01JHU520_IncludeIntrons/LZ01JHU520_Combined_CellBender_rRNA.rds")

# Identify the differentially expressed features for each cluster. 
PDAC_combined.markers <- FindAllMarkers(PDAC_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PDAC_combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


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
# In this example, X, Y, and Z are the cell types corresponding to clusters 0, 1, and 2, respectively. 
# These are only the basic cell types, can do more identification in the future
PDAC_combined <- readRDS("~/Desktop/LZ01JHU520_IncludeIntrons/PDAC_Combined.RDS")
Idents(PDAC_combined) <- "seurat_clusters"
new.cluster.ids <- c("Ductal", "Ductal", "Fibroblast", "Ductal", "Lymphoid", "Ductal", "Myeloid", "Fibroblast",
                     "Lymphoid", "Myeloid", "Myeloid", "Fibroblast", "Ductal", "Endothelial", "Acinar", "Hepatocyte")
names(new.cluster.ids) <- levels(PDAC_combined)
renamed_PDAC_combined <- RenameIdents(PDAC_combined, new.cluster.ids)


renamed_PDAC_combined$cell_type <- Idents(object = renamed_PDAC_combined)

# Display labeled UMAP
Idents(renamed_PDAC_combined) <- "cell_type"
# Figure 4A
DimPlot(renamed_PDAC_combined, label = FALSE) + ggtitle("KPC-4545 (All Samples)") + theme(plot.title = element_text(hjust = 0.5))

# Label samples by their treatment 
Idents(renamed_PDAC_combined) <- "orig.ident"
treatment_ids <- c("sn1_4 (Control)", "sn2_4 (TGFB-TRAP)", "sn3_4 (anti-IL1B)","sn4_4 (anti-PD1)", "sn5_4 (anti-PD1 + TGFB-TRAP)", "sn6_4 (anti-PD1 + anti-IL1B)")
names(treatment_ids) <- levels(renamed_PDAC_combined)
renamed_PDAC_combined <- RenameIdents(renamed_PDAC_combined, treatment_ids)
renamed_PDAC_combined$treatment <- Idents(object = renamed_PDAC_combined)

# Extract TGFB-TRAP samples for presentation and analysis (sn1_4, sn2_4, sn4_4, sn5_4)
Idents(renamed_PDAC_combined) <- "orig.ident"
renamed_PDAC_combined_TGFB <- subset(renamed_PDAC_combined,  idents = c("sn1_4", "sn2_4", "sn4_4", "sn5_4"))
Idents(renamed_PDAC_combined_TGFB) <- "cell_type"
# Figure 4B
DimPlot(renamed_PDAC_combined_TGFB, ncol = 2, split.by = "treatment")

# Create heatmap for notable markers of identified cell type for TGFB-TRAP samples
combined_averages_TGFB <- AverageExpression(renamed_PDAC_combined_TGFB, return.seurat = TRUE)
dotplot_markers_TGFB <- c("Krt7", "Krt18", "Epcam", "Col1a1", "Col3a1", "Dcn", "Ptprc", "Themis", "Cd96", "Cd3e", "Bank1", 
                          "Cd74", "Itgam", "Fcgr2b", "Cdh5", "Pecam1", "Vwf", "Prss2", "Ctrb1", "Pnliprp1", "Alb" , "Pck1")
marker_heatmap_matrix_TGFB <- as.data.frame((combined_averages_TGFB@assays$RNA@scale.data)[dotplot_markers_TGFB,])
marker_heatmap_matrix_TGFB <-  marker_heatmap_matrix_TGFB[,c("Ductal", "Fibroblast", "Lymphoid", "Myeloid", "Endothelial", "Acinar", "Hepatocyte")]
# Figure 4C
jpeg(paste0("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/TGFB Renamed_Cluster_Heatmap_top3.jpeg"),  width = 25, height = 9, units='cm', res = 300)
Heatmap(t(marker_heatmap_matrix_TGFB), name = "Normalized \nExpression", rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()


#####################################################

# Ductal Analysis 
PDAC_combined_Ductal_TGFB <- subset(renamed_PDAC_combined_TGFB, idents = c("Ductal"))
DimPlot(PDAC_combined_Ductal_TGFB, label = TRUE)

# Create Grouped Violin Plots with TGFB Downstream and Signaling Genes for Ductal cells
TGFB_Signaling_Violin_Plot_Ductal <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(PDAC_combined_Ductal_TGFB, features = signature,
            pt.size = 0.1, 
            assay = "RNA",
            group.by = "orig.ident", 
            y.max = y_max)    + stat_compare_means(comparisons = test_sign, label = "p.signif", method = "wilcox.test") + theme(legend.position = 'none', axis.title.x = element_blank())
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 2) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 14, height = 10)
}

# Comparisons to make
my_comparisons <- rev(list(c("sn1_4","sn2_4"),c("sn4_4","sn5_4")))

# TGFB Downstream Signaling Genes
TGFB_downstream_genes_of_interest_Ductal <- c("Runx2","Ccn2", "Serpine1", "Zeb1", "Fn1", "Vim", "Thbs1")
# Figure 4D
TGFB_Signaling_Violin_Plot_Ductal(gene_signature = TGFB_downstream_genes_of_interest_Ductal, file_name = "~/Desktop/Ductal_TGFB_downstream_genes_of_interest", test_sign = my_comparisons)

# TGFB Pathway Genes
TGFB_pathway_genes_of_interest <- c("Tgfbr1", "Tgfbr2", "Tgfbr3", "Smad1","Smad2", "Smad3", "Smad4", "Smad5", "Smad6", "Smad7")
# Supplemental Figure 9A
TGFB_Signaling_Violin_Plot_Ductal(gene_signature = TGFB_pathway_genes_of_interest, file_name = "~/Desktop/Ductal_1245_TGFB_Pathway", test_sign = my_comparisons)


#####################################################

# CAF Analysis 

# Isolate CAFs for Downstream analysis
LZ01JHU520_Fibroblast <- subset(renamed_PDAC_combined, idents = c("Fibroblast"))
#Display object
LZ01JHU520_Fibroblast

# Recluster CAFs
varfeat_CAF_combined <- FindVariableFeatures(LZ01JHU520_Fibroblast, selection.method = "vst", nfeatures = 2000)
all.genes.CAF <- rownames(varfeat_CAF_combined)
scaled_CAF_combined<- ScaleData(varfeat_CAF_combined, features = all.genes.CAF)

# Linear Dimensional Reduction of Scaled Data
scaled_CAF_combined <- RunPCA(scaled_CAF_combined, features = VariableFeatures(object = scaled_CAF_combined))

# Determine the dimensionality of data. What PCs are most statistically significant in their variation? 
CAF_combined <- JackStraw(scaled_CAF_combined, dims = 20)
CAF_combined <- ScoreJackStraw(CAF_combined, dims = 1:20)

#Visualize JackStraw plot. Distribution of p-values for each PC. Dashed line is a uniform distribubtion.
JackStrawPlot(CAF_combined, dims = 1:20)
# Compare JackStraw Plot to Elbow plot to see where most variation arises. As plots level off, there is less variation coming from these populations. 
ElbowPlot(CAF_combined)

# Cluster the cells
# Determine the dimensions of interest from JackStraw and ElbowPlot above, and add to code below.
CAF_combined<- FindNeighbors(CAF_combined, dims = 1:4)
CAF_combined <- FindClusters(CAF_combined, resolution = 0.5)

# Run non-linear dimensional reduction. Create UMAP/tSNE. Specifiy dimensionality based on significant PCs. 
CAF_combined <- RunUMAP(CAF_combined, dims = 1:4)
# Visualize UMAP
DimPlot(CAF_combined, reduction = "umap", label = TRUE)

Idents(object = CAF_combined) <- "seurat_clusters"
CAF.new.cluster.ids <- c("s0", "s1", "s2", "s3", "s4", "s5", "s6")
names(CAF.new.cluster.ids) <- levels(CAF_combined)
renamed_CAF_combined <- RenameIdents(CAF_combined, CAF.new.cluster.ids)
# Display labeled UMAP
# Figure 5A
DimPlot(renamed_CAF_combined, reduction = "umap", label = TRUE, pt.size = 0.5)

# Table for CAF Distribution
CAF_distribution_table <- table(Idents(renamed_CAF_combined), renamed_CAF_combined$orig.ident)
CAF_distribution_table

# Find all markers for every cluster compared to all remaining cells. This will report only the positive values. 
CAF_combined.markers <- FindAllMarkers(CAF_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CAF_combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
# Save all markers as a spreadsheet for easier analysis
write.csv (CAF_combined.markers %>% group_by(cluster), file ="~/Desktop/LZ01JHU520_IncludeIntrons/CAF/CAF_combined_markers.csv", quote =TRUE)

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

# Figure 5B
Idents(CAF_combined) <- "orig.ident"
# Isolate TGFB-TRAP samples
CAF_TGFB <- subset(CAF_combined,  idents = c("sn1_4", "sn2_4", "sn4_4", "sn5_4"))
Idents(CAF_TGFB) <- "seurat_clusters"
combined.cluster.ids <- c("s0", "s1", "s2", "s3", "s4", "s5", "s6")
names(combined.cluster.ids) <- levels(CAF_TGFB)
renamed_CAF_TGFB <- RenameIdents(CAF_TGFB, combined.cluster.ids)
DimPlot(renamed_CAF_TGFB, ncol = 2, split.by = "treatment", label = TRUE) 


# Lrrc15 Expression 
# Figure 5C
FeaturePlot(CAF_combined, features = "Lrrc15")
VlnPlot(CAF_combined, features = "Lrrc15")
FeaturePlot(CAF_sn1_4, features = "Lrrc15") + ggtitle("sn1_4 (Control) Lrrc15") + theme(plot.title = element_text(hjust = 0.5))
FeaturePlot(CAF_sn2_4, features = "Lrrc15")+ ggtitle("sn2_4 (TGFB-TRAP) Lrrc15") + theme(plot.title = element_text(hjust = 0.5))
FeaturePlot(CAF_sn4_4, features = "Lrrc15")+ ggtitle("sn4_4 (anti-PD1) Lrrc15") + theme(plot.title = element_text(hjust = 0.5))
FeaturePlot(CAF_sn5_4, features = "Lrrc15")+ ggtitle("sn5_4 (anti-PD1 + TGFB-TRAP) Lrrc15") + theme(plot.title = element_text(hjust = 0.5))


# Extract clusters 2 and 4 relative count
Idents(CAF_combined) <- "seurat_clusters"
# Isolate CAF subcluster 2 and subcluster 4
CAF_c2_c4 <- subset(CAF_combined, idents = c("2", "4"))
meta.data.cluster <- unique(x = CAF_c2_c4@meta.data$orig.ident)
Idents(CAF_c2_c4) <- "orig.ident"
for(group in meta.data.cluster) {
  group.cells <- WhichCells(object = CAF_c2_c4 , idents  = group)
  data_to_write_out <- as.data.frame(x = as.matrix(x = CAF_c2_c4@assays$RNA@scale.data[, group.cells]))
  data_to_write_out <- t(as.data.frame(data_to_write_out['Lrrc15',]))
  write.csv(x = data_to_write_out, row.names = TRUE, file = paste0("~/Desktop/LZ01JHU520_IncludeIntrons/CAF/Lrrc15", "/", group, "c2_c4_Lrrc15_scaled_count.csv"))
}

# Figure 5E based on Dominguez et al. Cancer Discovery Figure 3I Left
figure_3I_genes_left <- c("Tgfbi", "Igfbp3", "Tagln", "Acta2", "Cnn1",  "Sema7a", "Ctgf", "Tpm1", "Fstl3", "Sh3pxd2a", "Adam19", "Hsbp1", "Adam12", "Tns1", "Pxdc1", "Col4a1", "Actg2")
# Figure 5E 
Idents(CAF_c2_c4) <- "orig.ident"
# Isolate TGFB-TRAP samples 
CAF_c2_c4_TGFB <- subset(CAF_c2_c4, idents = c("sn1_4", "sn2_4", "sn4_4", "sn5_4"))
fig_3I_left_c2_c4_between_samples_TGFB <- AverageExpression(CAF_c2_c4_TGFB, features = figure_3I_genes_left, group.by = "orig.ident", slot = "data")
fig_3I_left_c2_c4_between_samples_TGFB <- as.data.frame(fig_3I_left_c2_c4_between_samples_TGFB$RNA)
fig_3I_left_c2_c4_between_samples_TGFB<- as.data.frame(t(apply(fig_3I_left_c2_c4_between_samples_TGFB, 1, cal_z_score)))
fig_3I_left_c2_c4_between_samples_TGFB <- as.matrix(t(fig_3I_left_c2_c4_between_samples_TGFB))
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_fun(seq(-3, 3))
Heatmap(fig_3I_left_c2_c4_between_samples_TGFB, name = "mat", col = col_fun, cluster_rows = FALSE)


# Thsd4 Expression 
# Figure S11
FeaturePlot(CAF_combined_TGFB, features = "Thsd4", split.by = "orig.ident") +
  patchwork::plot_layout(ncol = 2, nrow = 2)
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/TGFB-TRAP Manuscript V7 for JCI Insight/PhotoShop/CAF s5 Thsd4 Expression.png", width = 22, height = 20, units = "cm")


# Create Grouped Violin Plots with TGFB Downstream and Signaling Genes for Ductal cells
# Figure 4E
TGFB_Signaling_Violin_Plot_CAF <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(CAF_TGFB, features = signature,
            pt.size = 0.1, 
            assay = "RNA",
            group.by = "orig.ident", 
            y.max = y_max)    + stat_compare_means(comparisons = test_sign, label = "p.signif", method = "wilcox.test") + theme(legend.position = 'none', axis.title.x = element_blank())
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 2) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 14, height = 10)
}

# Comparisons to be made
my_comparisons <- rev(list(c("sn1_4","sn2_4"),c("sn4_4","sn5_4")))
# TGFB pathway genes
TGFB_pathway_genes_of_interest <- c("Tgfbr1", "Tgfbr2", "Tgfbr3", "Smad1","Smad2", "Smad3", "Smad4", "Smad5", "Smad6", "Smad7")
# TGFB downstream signaling genes for CAFs
TGFB_downstream_genes_of_interest_CAF <- c("Runx2","Ccn2", "Zeb1", "Fn1", "Col1a1", "Col3a1", "Vim", "Thbs1")
# Figure 4E
TGFB_Signaling_Violin_Plot_CAF(gene_signature = TGFB_downstream_genes_of_interest_CAF, file_name = "~/Desktop/CAF_TGFB_downstream_genes_of_interest", test_sign = my_comparisons)
# Figure S9B
TGFB_Signaling_Violin_Plot_CAF(gene_signature = TGFB_pathway_genes_of_interest, file_name = "~/Desktop/CAF_1245_TGFB_Pathway", test_sign = my_comparisons)


# Organizing CAF subclusters (s0-s6) into CAF cluster 1 (CAF_c1), CAF cluster 2 (CAF_c2), CAF cluster 3 (CAF_c3)
Idents(CAF_TGFB) <- "seurat_clusters"
CAF_ids <- c("CAF_c1", "CAF_c1", "CAF_c2",  "CAF_c3", "CAF_c2", "CAF_c1", "CAF_c3")
names(CAF_ids) <- levels(renamed_CAF_TGFB)
renamed_CAF_TGFB <- RenameIdents(renamed_CAF_TGFB, CAF_ids)
# Figure S10A
DimPlot(renamed_CAF_TGFB, ncol = 2, split.by = "treatment") + ggtitle("KPC-4545 CAFs (TGFB-TRAP Samples Separated)") + theme(plot.title = element_text(hjust = 0.5)) 
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/CAF Renamed TGFB-TRAP 4 samples separated.png", width = 25, height = 20, units = "cm")
# Figure S10B
DimPlot(renamed_CAF_TGFB, label = FALSE) + ggtitle("KPC-4545 CAFs (TGFB-TRAP Samples Combined)") + theme(plot.title = element_text(hjust = 0.5)) 
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/CAF Renamed TGFB-TRAP 4 Samples Only Combined.png", width = 25, height = 20, units = "cm")


# Differential Gene Expression for CAFs based on CAF_c1, CAF_c2, and CAF_c3 between different treatments and Volcano Plots
# Figure 6A
# Figure S13A-F

# CAF_c1
# CAF_c1 sn1_4 vs sn2_4
CAF_c1 <- subset(renamed_CAF_combined, idents = "CAF_c1")
DimPlot(CAF_c1, label = TRUE)
Idents(CAF_c1) <- "orig.ident"
CAF_c1_sn1_v_sn2 <- FindMarkers(CAF_c1, ident.1 = "sn2_4", ident.2 = "sn1_4", test.use = "MAST", logfc.threshold = 0 )
EnhancedVolcano(CAF_c1_sn1_v_sn2,
                lab = rownames(CAF_c1_sn1_v_sn2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CAF_c1: sn1_4 (Control) vs sn2_4 (TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/CAF_c1 sn1_4 vs sn2_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
renamed_CAF_c1_sn1_v_sn2_upregulated <- CAF_c1_sn1_v_sn2[CAF_c1_sn1_v_sn2$avg_log2FC > 0.5 & CAF_c1_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c1_sn1_v_sn2_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c1_sn1_v_sn2_upregulated_p_adj.csv")
renamed_CAF_c1_sn1_v_sn2_downregulated <- CAF_c1_sn1_v_sn2[CAF_c1_sn1_v_sn2$avg_log2FC < -0.5 & CAF_c1_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c1_sn1_v_sn2_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c1_sn1_v_sn2_downregulated_p_adj.csv")
renamed_CAF_c1_sn1_v_sn2_up <- CAF_c1_sn1_v_sn2[CAF_c1_sn1_v_sn2$avg_log2FC > 0, ]
renamed_CAF_c1_sn1_v_sn2_down <- CAF_c1_sn1_v_sn2[CAF_c1_sn1_v_sn2$avg_log2FC < 0, ]
write.csv(renamed_CAF_c1_sn1_v_sn2_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c1_sn1_v_sn2_up.csv")
write.csv(renamed_CAF_c1_sn1_v_sn2_down,"~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c1_sn1_v_sn2_down.csv")


# CAF_c1 sn4_4 vs sn5_4
Idents(CAF_c1) <- "orig.ident"
CAF_c1_sn4_v_sn5 <- FindMarkers(CAF_c1, ident.1 = "sn5_4", ident.2 = "sn4_4", test.use = "MAST", logfc.threshold = 0 )
EnhancedVolcano(CAF_c1_sn4_v_sn5,
                lab = rownames(CAF_c1_sn4_v_sn5),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CAF_c1: sn4_4 (anti-PD1) vs sn2_4 (anti-PD1 + TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/CAF_c1 sn4_4 vs sn5_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
renamed_CAF_c1_sn4_v_sn5_upregulated <- CAF_c1_sn4_v_sn5[CAF_c1_sn4_v_sn5$avg_log2FC > 0.5 & CAF_c1_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c1_sn4_v_sn5_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c1_sn4_v_sn5_upregulated_p_adj.csv")
renamed_CAF_c1_sn4_v_sn5_downregulated <- CAF_c1_sn4_v_sn5[CAF_c1_sn4_v_sn5$avg_log2FC < -0.5 & CAF_c1_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c1_sn4_v_sn5_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c1_sn4_v_sn5_downregulated_p_adj.csv")
renamed_CAF_c1_sn4_v_sn5_up <- CAF_c1_sn4_v_sn5[CAF_c1_sn4_v_sn5$avg_log2FC > 0, ]
renamed_CAF_c1_sn4_v_sn5_down <- CAF_c1_sn4_v_sn5[CAF_c1_sn4_v_sn5$avg_log2FC < 0, ]
write.csv(renamed_CAF_c1_sn4_v_sn5_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c1_sn4_v_sn5_up.csv")
write.csv(renamed_CAF_c1_sn4_v_sn5_down,"~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c1_sn4_v_sn5_down.csv")


# CAF_c2
# CAF_c2 sn1_4 vs sn2_4
CAF_c2 <- subset(renamed_CAF_combined, idents = "CAF_c2")
DimPlot(CAF_c2, label = TRUE)
Idents(CAF_c2) <- "orig.ident"
CAF_c2_sn1_v_sn2 <- FindMarkers(CAF_c2, ident.1 = "sn2_4", ident.2 = "sn1_4", test.use = "MAST", logfc.threshold = 0 )
EnhancedVolcano(CAF_c2_sn1_v_sn2,
                lab = rownames(CAF_c2_sn1_v_sn2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CAF_c2: sn1_4 (Control) vs sn2_4 (TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/CAF_c2 sn1_4 vs sn2_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
renamed_CAF_c2_sn1_v_sn2_upregulated <- CAF_c2_sn1_v_sn2[CAF_c2_sn1_v_sn2$avg_log2FC > 0.5 & CAF_c2_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c2_sn1_v_sn2_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c2_sn1_v_sn2_upregulated_p_adj.csv")
renamed_CAF_c2_sn1_v_sn2_downregulated <- CAF_c2_sn1_v_sn2[CAF_c2_sn1_v_sn2$avg_log2FC < -0.5 & CAF_c2_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c2_sn1_v_sn2_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c2_sn1_v_sn2_downregulated_p_adj.csv")
renamed_CAF_c2_sn1_v_sn2_up <- CAF_c2_sn1_v_sn2[CAF_c2_sn1_v_sn2$avg_log2FC > 0, ]
renamed_CAF_c2_sn1_v_sn2_down <- CAF_c2_sn1_v_sn2[CAF_c2_sn1_v_sn2$avg_log2FC < 0, ]
write.csv(renamed_CAF_c2_sn1_v_sn2_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c2_sn1_v_sn2_up.csv")
write.csv(renamed_CAF_c2_sn1_v_sn2_down,"~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c2_sn1_v_sn2_down.csv")



# CAF_c2 sn4_4 vs sn5_4
Idents(CAF_c2) <- "orig.ident"
CAF_c2_sn4_v_sn5 <- FindMarkers(CAF_c2, ident.1 = "sn5_4", ident.2 = "sn4_4", test.use = "MAST", logfc.threshold = 0 )
EnhancedVolcano(CAF_c2_sn4_v_sn5,
                lab = rownames(CAF_c2_sn4_v_sn5),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CAF_c2: sn4_4 (anti-PD1) vs sn2_4 (anti-PD1 + TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/CAF_c2 sn4_4 vs sn5_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
renamed_CAF_c2_sn4_v_sn5_upregulated <- CAF_c2_sn4_v_sn5[CAF_c2_sn4_v_sn5$avg_log2FC > 0.5 & CAF_c2_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c2_sn4_v_sn5_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c2_sn4_v_sn5_upregulated_p_adj.csv")
renamed_CAF_c2_sn4_v_sn5_downregulated <- CAF_c2_sn4_v_sn5[CAF_c2_sn4_v_sn5$avg_log2FC < -0.5 & CAF_c2_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c2_sn4_v_sn5_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c2_sn4_v_sn5_downregulated_p_adj.csv")
renamed_CAF_c2_sn4_v_sn5_up <- CAF_c2_sn4_v_sn5[CAF_c2_sn4_v_sn5$avg_log2FC > 0, ]
renamed_CAF_c2_sn4_v_sn5_down <- CAF_c2_sn4_v_sn5[CAF_c2_sn4_v_sn5$avg_log2FC < 0, ]
write.csv(renamed_CAF_c2_sn4_v_sn5_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c2_sn4_v_sn5_up.csv")
write.csv(renamed_CAF_c2_sn4_v_sn5_down,"~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c2_sn4_v_sn5_down.csv")



# CAF_c3
# CAF_c3 sn1_4 vs sn2_4
CAF_c3 <- subset(renamed_CAF_combined, idents = "CAF_c3")
DimPlot(CAF_c3, label = TRUE)
Idents(CAF_c3) <- "orig.ident"
CAF_c3_sn1_v_sn2 <- FindMarkers(CAF_c3, ident.1 = "sn2_4", ident.2 = "sn1_4", test.use = "MAST", logfc.threshold = 0 )
EnhancedVolcano(CAF_c3_sn1_v_sn2,
                lab = rownames(CAF_c3_sn1_v_sn2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CAF_c3: sn1_4 (Control) vs sn2_4 (TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/CAF_c3 sn1_4 vs sn2_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
renamed_CAF_c3_sn1_v_sn2_upregulated <- CAF_c3_sn1_v_sn2[CAF_c3_sn1_v_sn2$avg_log2FC > 0.5 & CAF_c3_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c3_sn1_v_sn2_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c3_sn1_v_sn2_upregulated_p_adj.csv")
renamed_CAF_c3_sn1_v_sn2_downregulated <- CAF_c3_sn1_v_sn2[CAF_c3_sn1_v_sn2$avg_log2FC < -0.5 & CAF_c3_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c3_sn1_v_sn2_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c3_sn1_v_sn2_downregulated_p_adj.csv")
renamed_CAF_c3_sn1_v_sn2_up <- CAF_c3_sn1_v_sn2[CAF_c3_sn1_v_sn2$avg_log2FC > 0, ]
renamed_CAF_c3_sn1_v_sn2_down <- CAF_c3_sn1_v_sn2[CAF_c3_sn1_v_sn2$avg_log2FC < 0, ]
write.csv(renamed_CAF_c3_sn1_v_sn2_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c3_sn1_v_sn2_up.csv")
write.csv(renamed_CAF_c3_sn1_v_sn2_down,"~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c3_sn1_v_sn2_down.csv")



# CAF_c3 sn4_4 vs sn5_4
Idents(CAF_c3) <- "orig.ident"
CAF_c3_sn4_v_sn5 <- FindMarkers(CAF_c3, ident.1 = "sn5_4", ident.2 = "sn4_4", test.use = "MAST", logfc.threshold = 0 )
EnhancedVolcano(CAF_c3_sn4_v_sn5,
                lab = rownames(CAF_c3_sn4_v_sn5),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CAF_c3: sn4_4 (anti-PD1) vs sn2_4 (anti-PD1 + TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/CAF_c3 sn4_4 vs sn5_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
renamed_CAF_c3_sn4_v_sn5_upregulated <- CAF_c3_sn4_v_sn5[CAF_c3_sn4_v_sn5$avg_log2FC > 0.5 & CAF_c3_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c3_sn4_v_sn5_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c3_sn4_v_sn5_upregulated_p_adj.csv")
renamed_CAF_c3_sn4_v_sn5_downregulated <- CAF_c3_sn4_v_sn5[CAF_c3_sn4_v_sn5$avg_log2FC < -0.5 & CAF_c3_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(renamed_CAF_c3_sn4_v_sn5_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/CAF/DGE/CAF_c3_sn4_v_sn5_downregulated_p_adj.csv")
renamed_CAF_c3_sn4_v_sn5_up <- CAF_c3_sn4_v_sn5[CAF_c3_sn4_v_sn5$avg_log2FC > 0, ]
renamed_CAF_c3_sn4_v_sn5_down <- CAF_c3_sn4_v_sn5[CAF_c3_sn4_v_sn5$avg_log2FC < 0, ]
write.csv(renamed_CAF_c3_sn4_v_sn5_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c3_sn4_v_sn5_up.csv")
write.csv(renamed_CAF_c3_sn4_v_sn5_down,"~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CAF_c3_sn4_v_sn5_down.csv")




# myCAF and iCAF gene panels in Figure S10C-N
library(dplyr)
library(tidyr)
library(forcats)
Idents(CAF_combined) <- "orig.ident"
CAF_sn1_4 <- subset(CAF_combined, idents = c("sn1_4"))
CAF_sn2_4 <- subset(CAF_combined, idents = c("sn2_4"))
CAF_sn3_4 <- subset(CAF_combined, idents = c("sn3_4"))
CAF_sn4_4 <- subset(CAF_combined, idents = c("sn4_4"))
CAF_sn5_4 <- subset(CAF_combined, idents = c("sn5_4"))
CAF_sn6_4 <- subset(CAF_combined, idents = c("sn6_4"))

# sn1_4 Cluster 1
Idents(CAF_sn1_4) <- "seurat_clusters"
CAF_sn1_c0 <- subset(CAF_sn1_4, idents = c("0", "1", "5"))
CAF_sn1_c0_iCAF_genes = data.frame(t(data.frame(CAF_sn1_c0[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn1_c0_iCAF_genes <- CAF_sn1_c0_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c0_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn1_c0_myCAF_genes = data.frame(t(data.frame(CAF_sn1_c0[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn1_c0_myCAF_genes <- CAF_sn1_c0_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c0_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn1_c0_combined_df <- rbind(CAF_sn1_c0_iCAF_genes, CAF_sn1_c0_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn1_c0_combined_df[, "Expression"])) / 100000
CAF_sn1_c0_combined_df$Expression <- CAF_sn1_c0_combined_df$Expression + abs(noise)
ggplot(CAF_sn1_c0_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 Control Cluster 1 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/IL1B/Figures for Publication/Supplemental Figure/CAF sn1_4 c0 Genes.png", width = 30, height = 20, units = "cm")

# sn1_4 Cluster 2
Idents(CAF_sn1_4) <- "seurat_clusters"
CAF_sn1_c2 <- subset(CAF_sn1_4, idents = c("2", "4"))
CAF_sn1_c2_iCAF_genes = data.frame(t(data.frame(CAF_sn1_c2[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn1_c2_iCAF_genes <- CAF_sn1_c2_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c2_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn1_c2_myCAF_genes = data.frame(t(data.frame(CAF_sn1_c2[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn1_c2_myCAF_genes <- CAF_sn1_c2_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c2_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn1_c2_combined_df <- rbind(CAF_sn1_c2_iCAF_genes, CAF_sn1_c2_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn1_c2_combined_df[, "Expression"])) / 100000
CAF_sn1_c2_combined_df$Expression <- CAF_sn1_c2_combined_df$Expression + abs(noise)
ggplot(CAF_sn1_c2_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 Control Cluster 2 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/IL1B/Figures for Publication/Supplemental Figure/CAF sn1_4 c2 Genes.png", width = 30, height = 20, units = "cm")

# sn1_4 Cluster 3
Idents(CAF_sn1_4) <- "seurat_clusters"
CAF_sn1_c3 <- subset(CAF_sn1_4, idents = c("3", "6"))
CAF_sn1_c3_iCAF_genes = data.frame(t(data.frame(CAF_sn1_c3[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn1_c3_iCAF_genes <- CAF_sn1_c3_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c3_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn1_c3_myCAF_genes = data.frame(t(data.frame(CAF_sn1_c3[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn1_c3_myCAF_genes <- CAF_sn1_c3_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn1_c3_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn1_c3_combined_df <- rbind(CAF_sn1_c3_iCAF_genes, CAF_sn1_c3_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn1_c3_combined_df[, "Expression"])) / 100000
CAF_sn1_c3_combined_df$Expression <- CAF_sn1_c3_combined_df$Expression + abs(noise)
ggplot(CAF_sn1_c3_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 Control Cluster 3 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/IL1B/Figures for Publication/Supplemental Figure/CAF sn1_4 c3 Genes.png", width = 30, height = 20, units = "cm")


# sn2_4 Cluster 1
Idents(CAF_sn2_4) <- "seurat_clusters"
CAF_sn2_c0 <- subset(CAF_sn2_4, idents = c("0", "1", "5"))
CAF_sn2_c0_iCAF_genes = data.frame(t(data.frame(CAF_sn2_c0[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn2_c0_iCAF_genes <- CAF_sn2_c0_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c0_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn2_c0_myCAF_genes = data.frame(t(data.frame(CAF_sn2_c0[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn2_c0_myCAF_genes <- CAF_sn2_c0_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c0_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn2_c0_combined_df <- rbind(CAF_sn2_c0_iCAF_genes, CAF_sn2_c0_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn2_c0_combined_df[, "Expression"])) / 100000
CAF_sn2_c0_combined_df$Expression <- CAF_sn2_c0_combined_df$Expression + abs(noise)
ggplot(CAF_sn2_c0_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 TGFB-TRAP Cluster 1 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/CAF sn2_4 c0 Genes.png", width = 30, height = 20, units = "cm")

# sn2_4 Cluster 2
Idents(CAF_sn2_4) <- "seurat_clusters"
CAF_sn2_c2 <- subset(CAF_sn2_4, idents = c("2", "4"))
CAF_sn2_c2_iCAF_genes = data.frame(t(data.frame(CAF_sn2_c2[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn2_c2_iCAF_genes <- CAF_sn2_c2_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c2_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn2_c2_myCAF_genes = data.frame(t(data.frame(CAF_sn2_c2[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn2_c2_myCAF_genes <- CAF_sn2_c2_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c2_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn2_c2_combined_df <- rbind(CAF_sn2_c2_iCAF_genes, CAF_sn2_c2_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn2_c2_combined_df[, "Expression"])) / 100000
CAF_sn2_c2_combined_df$Expression <- CAF_sn2_c2_combined_df$Expression + abs(noise)
ggplot(CAF_sn2_c2_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 TGFB-TRAP Cluster 2 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/CAF sn2_4 c2 Genes.png", width = 30, height = 20, units = "cm")

# sn2_4 Cluster 3
Idents(CAF_sn2_4) <- "seurat_clusters"
CAF_sn2_c3 <- subset(CAF_sn2_4, idents = c("3", "6"))
CAF_sn2_c3_iCAF_genes = data.frame(t(data.frame(CAF_sn2_c3[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn2_c3_iCAF_genes <- CAF_sn2_c3_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c3_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn2_c3_myCAF_genes = data.frame(t(data.frame(CAF_sn2_c3[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn2_c3_myCAF_genes <- CAF_sn2_c3_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn2_c3_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn2_c3_combined_df <- rbind(CAF_sn2_c3_iCAF_genes, CAF_sn2_c3_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn2_c3_combined_df[, "Expression"])) / 100000
CAF_sn2_c3_combined_df$Expression <- CAF_sn2_c3_combined_df$Expression + abs(noise)
ggplot(CAF_sn2_c3_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 TGFB-TRAP Cluster 3 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/CAF sn2_4 c3 Genes.png", width = 30, height = 20, units = "cm")



# sn4_4 Cluster 1
Idents(CAF_sn4_4) <- "seurat_clusters"
CAF_sn4_c0 <- subset(CAF_sn4_4, idents = c("0", "1", "5"))
CAF_sn4_c0_iCAF_genes = data.frame(t(data.frame(CAF_sn4_c0[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn4_c0_iCAF_genes <- CAF_sn4_c0_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c0_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn4_c0_myCAF_genes = data.frame(t(data.frame(CAF_sn4_c0[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn4_c0_myCAF_genes <- CAF_sn4_c0_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c0_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn4_c0_combined_df <- rbind(CAF_sn4_c0_iCAF_genes, CAF_sn4_c0_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn4_c0_combined_df[, "Expression"])) / 100000
CAF_sn4_c0_combined_df$Expression <- CAF_sn4_c0_combined_df$Expression + abs(noise)
ggplot(CAF_sn4_c0_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 anti-PD1 Cluster 1 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/IL1B/Figures for Publication/Supplemental Figure/CAF sn4_4 c0 Genes.png", width = 30, height = 20, units = "cm")

# sn4_4 Cluster 2
Idents(CAF_sn4_4) <- "seurat_clusters"
CAF_sn4_c2 <- subset(CAF_sn4_4, idents = c("2", "4"))
CAF_sn4_c2_iCAF_genes = data.frame(t(data.frame(CAF_sn4_c2[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn4_c2_iCAF_genes <- CAF_sn4_c2_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c2_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn4_c2_myCAF_genes = data.frame(t(data.frame(CAF_sn4_c2[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn4_c2_myCAF_genes <- CAF_sn4_c2_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c2_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn4_c2_combined_df <- rbind(CAF_sn4_c2_iCAF_genes, CAF_sn4_c2_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn4_c2_combined_df[, "Expression"])) / 100000
CAF_sn4_c2_combined_df$Expression <- CAF_sn4_c2_combined_df$Expression + abs(noise)
ggplot(CAF_sn4_c2_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 anti-PD1 Cluster 2 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/IL1B/Figures for Publication/Supplemental Figure/CAF sn4_4 c2 Genes.png", width = 30, height = 20, units = "cm")

# sn4_4 Cluster 3
Idents(CAF_sn4_4) <- "seurat_clusters"
CAF_sn4_c3 <- subset(CAF_sn4_4, idents = c("3", "6"))
CAF_sn4_c3_iCAF_genes = data.frame(t(data.frame(CAF_sn4_c3[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn4_c3_iCAF_genes <- CAF_sn4_c3_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c3_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn4_c3_myCAF_genes = data.frame(t(data.frame(CAF_sn4_c3[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn4_c3_myCAF_genes <- CAF_sn4_c3_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn4_c3_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn4_c3_combined_df <- rbind(CAF_sn4_c3_iCAF_genes, CAF_sn4_c3_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn4_c3_combined_df[, "Expression"])) / 100000
CAF_sn4_c3_combined_df$Expression <- CAF_sn4_c3_combined_df$Expression + abs(noise)
ggplot(CAF_sn4_c3_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 anti-PD1 Cluster 3 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/IL1B/Figures for Publication/Supplemental Figure/CAF sn4_4 c3 Genes.png", width = 30, height = 20, units = "cm")



# sn5_4 Cluster 1
Idents(CAF_sn5_4) <- "seurat_clusters"
CAF_sn5_c0 <- subset(CAF_sn5_4, idents = c("0", "1", "5"))
CAF_sn5_c0_iCAF_genes = data.frame(t(data.frame(CAF_sn5_c0[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn5_c0_iCAF_genes <- CAF_sn5_c0_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c0_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn5_c0_myCAF_genes = data.frame(t(data.frame(CAF_sn5_c0[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn5_c0_myCAF_genes <- CAF_sn5_c0_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c0_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn5_c0_combined_df <- rbind(CAF_sn5_c0_iCAF_genes, CAF_sn5_c0_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn5_c0_combined_df[, "Expression"])) / 100000
CAF_sn5_c0_combined_df$Expression <- CAF_sn5_c0_combined_df$Expression + abs(noise)
ggplot(CAF_sn5_c0_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 anti-PD1 + TGFB-TRAP Cluster 1 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/CAF sn5_4 c0 Genes.png", width = 30, height = 20, units = "cm")

# sn5_4 Cluster 2
Idents(CAF_sn5_4) <- "seurat_clusters"
CAF_sn5_c2 <- subset(CAF_sn5_4, idents = c("2", "4"))
CAF_sn5_c2_iCAF_genes = data.frame(t(data.frame(CAF_sn5_c2[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn5_c2_iCAF_genes <- CAF_sn5_c2_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c2_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn5_c2_myCAF_genes = data.frame(t(data.frame(CAF_sn5_c2[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn5_c2_myCAF_genes <- CAF_sn5_c2_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c2_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn5_c2_combined_df <- rbind(CAF_sn5_c2_iCAF_genes, CAF_sn5_c2_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn5_c2_combined_df[, "Expression"])) / 100000
CAF_sn5_c2_combined_df$Expression <- CAF_sn5_c2_combined_df$Expression + abs(noise)
ggplot(CAF_sn5_c2_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 anti-PD1 + TGFB-TRAP Cluster 2 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/CAF sn5_4 c2 Genes.png", width = 30, height = 20, units = "cm")

# sn5_4 Cluster 3
Idents(CAF_sn5_4) <- "seurat_clusters"
CAF_sn5_c3 <- subset(CAF_sn5_4, idents = c("3", "6"))
CAF_sn5_c3_iCAF_genes = data.frame(t(data.frame(CAF_sn5_c3[["RNA"]]@data[c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3", "Abl2", "Fbn2"),])))
CAF_sn5_c3_iCAF_genes <- CAF_sn5_c3_iCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c3_iCAF_genes['PhenotypeMarker'] <- "iCAF"
CAF_sn5_c3_myCAF_genes = data.frame(t(data.frame(CAF_sn5_c3[["RNA"]]@data[c("Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1"),])))
CAF_sn5_c3_myCAF_genes <- CAF_sn5_c3_myCAF_genes %>% gather(Genes, Expression) 
CAF_sn5_c3_myCAF_genes['PhenotypeMarker'] <- "myCAF"
CAF_sn5_c3_combined_df <- rbind(CAF_sn5_c3_iCAF_genes, CAF_sn5_c3_myCAF_genes) %>% 
  mutate(Genes = fct_relevel(Genes, c("Svep1", "Plpp3", "Fbn1", "Ccdc80", "C3",  "Abl2", "Fbn2", "Tpm1", "Igfbp3", "Thbs2", "Spp1", "Col12a1", "Dock8", "Col15a1")))
noise <- rnorm(n = length(x = CAF_sn5_c3_combined_df[, "Expression"])) / 100000
CAF_sn5_c3_combined_df$Expression <- CAF_sn5_c3_combined_df$Expression + abs(noise)
ggplot(CAF_sn5_c3_combined_df, aes(x = Genes, y = Expression)) +  geom_violin(trim=TRUE, scale = "width", aes(fill = PhenotypeMarker)) + 
  geom_jitter(width = 0.2) +theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold" ), panel.border = element_blank(), axis.text.y=element_text(size=20),axis.text.x=element_text(size=20, angle = 45, hjust = 1), axis.title = element_text(size = 20),legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) + ylim(0, 5) + ggtitle("KPC-4545 anti-PD1 + TGFB-TRAP Cluster 3 CAF Gene Panel") + ylab("Expression Level")
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Supplemental Figure/CAF sn5_4 c3 Genes.png", width = 30, height = 20, units = "cm")


#####################################################

# Lymphoid Analysis

# Isolate all lymphoid cells
LZ01JHU520_Lymphoid <- subset(renamed_PDAC_combined, idents = c("Lymphoid"))
DimPlot(LZ01JHU520_Lymphoid, label = TRUE)

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
Lymphoid_combined <- FindNeighbors(Lymphoid_combined, dims = 1:15)
Lymphoid_combined <- FindClusters(Lymphoid_combined, resolution = 0.8)


# Run non-linear dimensional reduction. Create UMAP/tSNE. Specifiy dimensionality based on significant PCs. 
Lymphoid_combined <- RunUMAP(Lymphoid_combined, dims = 1:15)
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
Idents(Lymphoid_combined) <- "seurat_clusters"
Lymphoid.cluster.ids = c("CD8", "CD4", "Treg", "CD8", "CD8", "CD4", "B", "7", "CD4", "CD4", "Plasma" )
names(Lymphoid.cluster.ids) <- levels(Lymphoid_combined)
renamed_Lymphoid_combined <- RenameIdents(Lymphoid_combined, Lymphoid.cluster.ids)
# Figure 5F
DimPlot(renamed_Lymphoid_combined, cols = c("coral2", "goldenrod2", "green3", "turquoise2", "royalblue1", "maroon2"))  + ggtitle("KPC-4545 Lymphoid (All Samples)") + theme(plot.title = element_text(hjust = 0.5))
renamed_Lymphoid_combined$cell_type <- Idents(object = renamed_Lymphoid_combined)

# Isolate TGFB Samples
Idents(renamed_Lymphoid_combined) <- "orig.ident"
treatment_ids <- c("sn1_4 (Control)", "sn2_4 (TGFB-TRAP)", "sn3_4 (anti-IL1B)","sn4_4 (anti-PD1)", "sn5_4 (anti-PD1 + TGFB-TRAP)", "sn6_4 (anti-PD1 + anti-IL1B)")
names(treatment_ids) <- levels(renamed_Lymphoid_combined)
renamed_Lymphoid_combined <- RenameIdents(renamed_Lymphoid_combined, treatment_ids)
renamed_Lymphoid_combined$treatment <- Idents(object = renamed_Lymphoid_combined)
# Isolate TGFB Samples
Idents(renamed_Lymphoid_combined) <- "orig.ident"
renamed_Lymphoid_TGFB <- subset(renamed_Lymphoid_combined,  idents = c("sn1_4", "sn2_4", "sn4_4", "sn5_4"))
# Figure 5G
DimPlot(renamed_Lymphoid_TGFB, ncol = 2, label = FALSE,split.by = "treatment", cols = c("coral2", "goldenrod2", "green3", "turquoise2", "royalblue1", "maroon2"))
# TGFB Samples Only Figure S16A
DimPlot(renamed_Lymphoid_TGFB, label = TRUE)


# Differential Gene Expression for CD8 and CD4 T cells

# CD8
# Figure S13G-H
DimPlot(CD8_combined, label = TRUE)

# CD8: sn1_4 vs sn2_4
Idents(CD8_combined) <- "orig.ident"
CD8_sn1_v_sn2 <- FindMarkers(CD8_combined, ident.1 = "sn2_4", ident.2 = "sn1_4", test.use = "MAST", logfc.threshold = 0 )
EnhancedVolcano(CD8_sn1_v_sn2,
                lab = rownames(CD8_sn1_v_sn2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CD8: sn1_4 (Control) vs sn2_4 (TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Lymphoid CD8 sn1_4 vs sn2_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
CD8_sn1_v_sn2_upregulated <- CD8_sn1_v_sn2[CD8_sn1_v_sn2$avg_log2FC > 0.5 & CD8_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(CD8_sn1_v_sn2_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Lymphoid/DGE/all/CD8_sn1_v_sn2_sn1_v_sn2_upregulated_p_adj.csv")
CD8_sn1_v_sn2_downregulated <- CD8_sn1_v_sn2[CD8_sn1_v_sn2$avg_log2FC < -0.5 & CD8_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(CD8_sn1_v_sn2_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Lymphoid/DGE/all/CD8_sn1_v_sn2_sn1_v_sn2_downregulated_p_adj.csv")
CD8_sn1_v_sn2_up <- CD8_sn1_v_sn2[CD8_sn1_v_sn2$avg_log2FC > 0, ]
CD8_sn1_v_sn2_down <- CD8_sn1_v_sn2[CD8_sn1_v_sn2$avg_log2FC < 0, ]
write.csv(CD8_sn1_v_sn2_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CD8_sn1_v_sn2_sn1_v_sn2_up.csv")
write.csv(CD8_sn1_v_sn2_down, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CD8_sn1_v_sn2_sn1_v_sn2_down.csv")


# CD8: sn4_4 vs sn5_4
CD8_combined_sn4_v_sn5 <- FindMarkers(CD8_combined, ident.1 = "sn5_4", ident.2 = "sn4_4", test.use = "MAST", logfc.threshold = 0)
EnhancedVolcano(CD8_combined_sn4_v_sn5,
                lab = rownames(CD8_combined_sn4_v_sn5),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CD8: sn4_4 (anti-PD1) vs sn5_4 (anti-PD1 + TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Lymphoid CD8 sn4_4 vs sn5_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
CD8_sn4_v_sn5_upregulated <- CD8_combined_sn4_v_sn5[CD8_combined_sn4_v_sn5$avg_log2FC > 0.5 & CD8_combined_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(CD8_sn4_v_sn5_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Lymphoid/DGE/all/CD8_sn4_v_sn5_upregulated_p_adj.csv")
CD8_sn4_v_sn5_downregulated <- CD8_combined_sn4_v_sn5[CD8_combined_sn4_v_sn5$avg_log2FC < -0.5 & CD8_combined_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(CD8_sn4_v_sn5_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Lymphoid/DGE/all/CD8_sn4_v_sn5_downregulated_p_adj.csv")
CD8_sn4_v_sn5_up <- CD8_combined_sn4_v_sn5[CD8_combined_sn4_v_sn5$avg_log2FC > 0, ]
CD8_sn4_v_sn5_down <- CD8_combined_sn4_v_sn5[CD8_combined_sn4_v_sn5$avg_log2FC < 0, ]
write.csv(CD8_sn4_v_sn5_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CD8_sn4_v_sn5_up.csv")
write.csv(CD8_sn4_v_sn5_down, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CD8_sn4_v_sn5_down.csv")




# CD4
# Figure 6B
# Figure S13I-J

# CD4: sn1_4 vs sn2_4
DimPlot(CD4_combined, label = TRUE)
Idents(CD4_combined) <- "orig.ident"
CD4_sn1_v_sn2 <- FindMarkers(CD4_combined, ident.1 = "sn2_4", ident.2 = "sn1_4", test.use = "MAST", logfc.threshold = 0 )
EnhancedVolcano(CD4_sn1_v_sn2,
                lab = rownames(CD4_sn1_v_sn2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CD4: sn1_4 (Control) vs sn2_4 (TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)
ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Lymphoid CD4 sn1_4 vs sn2_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
CD4_sn1_v_sn2_upregulated <- CD4_sn1_v_sn2[CD4_sn1_v_sn2$avg_log2FC > 0.5 & CD4_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(CD4_sn1_v_sn2_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Lymphoid/DGE/all/CD4_sn1_v_sn2_sn1_v_sn2_upregulated_p_adj.csv")
CD4_sn1_v_sn2_downregulated <- CD4_sn1_v_sn2[CD4_sn1_v_sn2$avg_log2FC < -0.5 & CD4_sn1_v_sn2$p_val_adj < 0.05, ]
write.csv(CD4_sn1_v_sn2_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Lymphoid/DGE/all/CD4_sn1_v_sn2_sn1_v_sn2_downregulated_p_adj.csv")
CD4_sn1_v_sn2_up <- CD4_sn1_v_sn2[CD4_sn1_v_sn2$avg_log2FC > 0, ]
CD4_sn1_v_sn2_down <- CD4_sn1_v_sn2[CD4_sn1_v_sn2$avg_log2FC < 0, ]
write.csv(CD4_sn1_v_sn2_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CD4_sn1_v_sn2_up.csv")
write.csv(CD4_sn1_v_sn2_down, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CD4_sn1_v_sn2_down.csv")


# CD4: sn4_4 vs sn5_4
CD4_combined_sn4_v_sn5 <- FindMarkers(CD4_combined, ident.1 = "sn5_4", ident.2 = "sn4_4", test.use = "MAST", logfc.threshold = 0)
EnhancedVolcano(CD4_combined_sn4_v_sn5,
                lab = rownames(CD4_combined_sn4_v_sn5),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CD4: sn4_4 (anti-PD1) vs sn5_4 (anti-PD1 + TGFB-TRAP)',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

ggsave("~/Desktop/LZ01JHU520_IncludeIntrons/Figures for Publication/Lymphoid CD4 sn4_4 vs sn5_4 Volcano Plot.png", width = 25, height = 20, units = "cm")
# Find differentially expressed genes to feed into pathway analysis
CD4_sn4_v_sn5_upregulated <- CD4_combined_sn4_v_sn5[CD4_combined_sn4_v_sn5$avg_log2FC > 0.5 & CD4_combined_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(CD4_sn4_v_sn5_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Lymphoid/DGE/all/CD4_sn4_v_sn5_upregulated_p_adj.csv")
CD4_sn4_v_sn5_downregulated <- CD4_combined_sn4_v_sn5[CD4_combined_sn4_v_sn5$avg_log2FC < -0.5 & CD4_combined_sn4_v_sn5$p_val_adj < 0.05, ]
write.csv(CD4_sn4_v_sn5_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Lymphoid/DGE/all/CD4_sn4_v_sn5_downregulated_p_adj.csv")
CD4_sn4_v_sn5_up <- CD4_combined_sn4_v_sn5[CD4_combined_sn4_v_sn5$avg_log2FC > 0, ]
CD4_sn4_v_sn5_down <- CD4_combined_sn4_v_sn5[CD4_combined_sn4_v_sn5$avg_log2FC < 0, ]
write.csv(CD4_sn4_v_sn5_up, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CD4_sn4_v_sn5_up.csv")
write.csv(CD4_sn4_v_sn5_down, "~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/DGE Tables/CD4_sn4_v_sn5_down.csv")


## T cell exhaustion gene list for CD8 T cells
# Figure S12
Idents(renamed_Lymphoid_TGFB) <- "cell_type"
renamed_Lymphoid_TGFB_CD8 <- subset(renamed_Lymphoid_TGFB, idents = "CD8")
Idents(renamed_Lymphoid_TGFB_CD8) <- "orig.ident"

# https://www.nature.com/articles/s41467-022-35238-w
T_ex_gene_list <- c("Tigit","Havcr2", "Lag3", "Pdcd1", "Itgae", "Entpd1", "Tnfrsf9", "Tox", "Tbx21")
Idents(renamed_Lymphoid_TGFB_CD8) <- "orig.ident"
renamed_Lymphoid_TGFB_CD8_no_sn2 <- subset(renamed_Lymphoid_TGFB_CD8, idents = c("sn1_4", "sn4_4", "sn5_4"))

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
# Average expression
renamed_Lymphoid_TGFB_CD8_ex <- AverageExpression(renamed_Lymphoid_TGFB_CD8_no_sn2, features = T_ex_gene_list, group.by = "orig.ident", slot = "data")
renamed_Lymphoid_TGFB_CD8_ex <- as.data.frame(renamed_Lymphoid_TGFB_CD8_ex$RNA)
renamed_Lymphoid_TGFB_CD8_ex <- as.data.frame(t(apply(renamed_Lymphoid_TGFB_CD8_ex, 1, cal_z_score)))
renamed_Lymphoid_TGFB_CD8_ex <- as.matrix(t(renamed_Lymphoid_TGFB_CD8_ex))

library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_fun(seq(-3, 3))
png("~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/TGFB-TRAP Manuscript V7 for JCI Insight/PhotoShop/TGFB-TRAP Figure S10/Cd8 T Cell Exhaustion Panel.png",width=15,height=7,units="cm", res = 1200)
Heatmap(renamed_Lymphoid_TGFB_CD8_ex, name = "z-score", col = col_fun)
dev.off()



#####################################################

# Myeloid Analysis
LZ01JHU520_Myeloid <- subset(renamed_PDAC_combined, idents = c("Myeloid"))
DimPlot(LZ01JHU520_Myeloid, label = TRUE)

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
Myeloid_combined <- FindNeighbors(Myeloid_combined, dims = 1:15)
Myeloid_combined <- FindClusters(Myeloid_combined, resolution = 0.8)

# Run non-linear dimensional reduction. Create UMAP/tSNE. Specifiy dimensionality based on significant PCs. 
Myeloid_combined <- RunUMAP(Myeloid_combined, dims = 1:15)
# Visualize UMAP
DimPlot(Myeloid_combined, reduction = "umap", label = TRUE)
DimPlot(Myeloid_combined, label = TRUE, ncol = 3, split.by = "orig.ident")

# Find differentially expressed genes for each cluster
Idents(Myeloid_combined) <- "seurat_clusters"
Myeloid_combined.markers <- FindAllMarkers(Myeloid_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Myeloid_combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
# Save all markers as a spreadsheet for easier analysis
write.csv (Myeloid_combined.markers %>% group_by(cluster), file ="~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Myeloid_combined_markers_NEW_6_9_10.csv", quote =TRUE)

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

# Label Cell Type
Idents(Myeloid_combined) <- "seurat_clusters"
DimPlot(Myeloid_combined, label = TRUE)
Myeloid.cluster.ids = c("Macrophage", "Granulocyte", "Macrophage", "Macrophage", "Granulocyte", "Macrophage", "Macrophage", "Macrophage", "Macrophage", "Macrophage", "Macrophage", "Macrophage" )
names(Myeloid.cluster.ids) <- levels(Myeloid_combined)
renamed_Myeloid_combined <- RenameIdents(Myeloid_combined, Myeloid.cluster.ids)
# Figure 5H
DimPlot(renamed_Myeloid_combined, reduction = "umap",  pt.size = 0.5) + ggtitle("KPC-4545 Myeloid (All Samples)") + theme(plot.title = element_text(hjust = 0.5))
renamed_Myeloid_combined$cell_type <- Idents(object = renamed_Myeloid_combined)

# Isolate TGFB-TRAP Samples
Idents(renamed_Myeloid_combined) <- "orig.ident"
Myeloid_TGFB <- subset(renamed_Myeloid_combined,  idents = c("sn1_4", "sn2_4", "sn4_4", "sn5_4"))
Idents(Myeloid_TGFB) <- "cell_type"
# Figure 5I
DimPlot(Myeloid_TGFB, ncol = 2, split.by = "orig.ident")




# Differential Gene Expression for Granucoytes and Macrophage 

#  Granulocyte
DimPlot(renamed_Myeloid_TGFB, label = TRUE)
renamed_Myeloid_TGFB_Granulocyte <- subset(renamed_Myeloid_TGFB, idents = "Granulocyte")

# Granulocyte: sn1_4 vs sn2_4
Idents(renamed_Myeloid_TGFB_Granulocyte) <- "orig.ident"
renamed_Granulocyte_sn1_v_sn2 <- FindMarkers(renamed_Myeloid_TGFB_Granulocyte, ident.1 = "sn1_4", ident.2 = "sn2_4", test.use = "MAST")
EnhancedVolcano(renamed_Granulocyte_sn1_v_sn2,
                lab = rownames(renamed_Granulocyte_sn1_v_sn2),
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'Granulocyte: Control vs TGFB-TRAP',
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 6.0)
renamed_Granulocyte_sn1_v_sn2_upregulated <- renamed_Granulocyte_sn1_v_sn2[renamed_Granulocyte_sn1_v_sn2$avg_log2FC >= 0.58 & renamed_Granulocyte_sn1_v_sn2$p_val < 0.05, ]
write.csv(renamed_Granulocyte_sn1_v_sn2_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Myeloid/DGE/Granulocyte/renamed_Granulocyte_sn1_v_sn2_upregulated.csv")
renamed_Granulocyte_sn1_v_sn2_downregulated <- renamed_Granulocyte_sn1_v_sn2[renamed_Granulocyte_sn1_v_sn2$avg_log2FC <= -0.58 & renamed_Granulocyte_sn1_v_sn2$p_val < 0.05, ]
write.csv(renamed_Granulocyte_sn1_v_sn2_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Myeloid/DGE/Granulocyte/renamed_Granulocyte_sn1_v_sn2_downregulated.csv")


# Granulocyte: sn4_4 vs sn5_4
renamed_Granulocyte_sn4_v_sn5 <- FindMarkers(renamed_Myeloid_TGFB_Granulocyte, ident.1 = "sn4_4", ident.2 = "sn5_4", test.use = "MAST")
EnhancedVolcano(renamed_Granulocyte_sn4_v_sn5,
                lab = rownames(renamed_Granulocyte_sn4_v_sn5),
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'Granulocyte: anti-PD1 vs anti-PD1 + TGFB-TRAP',
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 6.0)
renamed_Granulocyte_sn4_v_sn5_upregulated <- renamed_Granulocyte_sn4_v_sn5[renamed_Granulocyte_sn4_v_sn5$avg_log2FC >= 0.58 & renamed_Granulocyte_sn4_v_sn5$p_val < 0.05, ]
write.csv(renamed_Granulocyte_sn4_v_sn5_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Myeloid/DGE/Granulocyte/renamed_Granulocyte_sn4_v_sn5_upregulated.csv")
renamed_Granulocyte_sn4_v_sn5_downregulated <- renamed_Granulocyte_sn4_v_sn5[renamed_Granulocyte_sn4_v_sn5$avg_log2FC <= -0.58 & renamed_Granulocyte_sn4_v_sn5$p_val < 0.05, ]
write.csv(renamed_Granulocyte_sn4_v_sn5_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Myeloid/DGE/Granulocyte/renamed_Granulocyte_sn4_v_sn5_downregulated.csv")



#  Macrophage
DimPlot(renamed_Myeloid_TGFB, label = TRUE)
renamed_Myeloid_TGFB_Macrophage <- subset(renamed_Myeloid_TGFB, idents = "Macrophage")

# Macrophage: sn1_4 vs sn2_4
Idents(renamed_Myeloid_TGFB_Macrophage) <- "orig.ident"
renamed_Macrophage_sn1_v_sn2 <- FindMarkers(renamed_Myeloid_TGFB_Macrophage, ident.1 = "sn1_4", ident.2 = "sn2_4", test.use = "MAST")
EnhancedVolcano(renamed_Macrophage_sn1_v_sn2,
                lab = rownames(renamed_Macrophage_sn1_v_sn2),
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'Macrophage: Control vs TGFB-TRAP',
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 6.0)
renamed_Macrophage_sn1_v_sn2_upregulated <- renamed_Macrophage_sn1_v_sn2[renamed_Macrophage_sn1_v_sn2$avg_log2FC >= 0.58 & renamed_Macrophage_sn1_v_sn2$p_val < 0.05, ]
write.csv(renamed_Macrophage_sn1_v_sn2_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Myeloid/DGE/Macrophage/renamed_Macrophage_sn1_v_sn2_upregulated.csv")
renamed_Macrophage_sn1_v_sn2_downregulated <- renamed_Macrophage_sn1_v_sn2[renamed_Macrophage_sn1_v_sn2$avg_log2FC <= -0.58 & renamed_Macrophage_sn1_v_sn2$p_val < 0.05, ]
write.csv(renamed_Macrophage_sn1_v_sn2_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Myeloid/DGE/Macrophage/renamed_Macrophage_sn1_v_sn2_downregulated.csv")


# Macrophage: sn4_4 vs sn5_4
renamed_Macrophage_sn4_v_sn5 <- FindMarkers(renamed_Myeloid_TGFB_Macrophage, ident.1 = "sn4_4", ident.2 = "sn5_4", test.use = "MAST")
EnhancedVolcano(renamed_Macrophage_sn4_v_sn5,
                lab = rownames(renamed_Macrophage_sn4_v_sn5),
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'Macrophage: anti-PD1 vs anti-PD1 + TGFB-TRAP',
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 6.0)
renamed_Macrophage_sn4_v_sn5_upregulated <- renamed_Macrophage_sn4_v_sn5[renamed_Macrophage_sn4_v_sn5$avg_log2FC >= 0.58 & renamed_Macrophage_sn4_v_sn5$p_val < 0.05, ]
write.csv(renamed_Macrophage_sn4_v_sn5_upregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Myeloid/DGE/Macrophage/renamed_Macrophage_sn4_v_sn5_upregulated.csv")
renamed_Macrophage_sn4_v_sn5_downregulated <- renamed_Macrophage_sn4_v_sn5[renamed_Macrophage_sn4_v_sn5$avg_log2FC <= -0.58 & renamed_Macrophage_sn4_v_sn5$p_val < 0.05, ]
write.csv(renamed_Macrophage_sn4_v_sn5_downregulated, "~/Desktop/LZ01JHU520_IncludeIntrons/Immune Analysis/Myeloid/DGE/Macrophage/renamed_Macrophage_sn4_v_sn5_downregulated.csv")





###################################
# Ligand Receptor Interaction with iTALK

library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)


PDAC_combined_raw_matrix <- as.matrix(GetAssayData(PDAC_combined, slot = "counts"))

# Create Annotations for Lymphoid Cells
cluster_letters <- Idents(object = renamed_Lymphoid_combined)
names(cluster_letters) <- colnames(x = renamed_Lymphoid_combined)
renamed_Lymphoid_combined <- AddMetaData(
  object = renamed_Lymphoid_combined,
  metadata = cluster_letters,
  col.name = 'labeled_ident'
)
Lymphoid_annotations <- renamed_Lymphoid_combined@meta.data %>% 
  rownames_to_column("barcodes")  %>% 
  dplyr::select("barcodes","labeled_ident")
Lymphoid_annotations <- Lymphoid_annotations %>% column_to_rownames(., var = "barcodes")


# Create Annotations for Myeloid Cells
DimPlot(renamed_Myeloid_combined, reduction = "umap", label = TRUE, pt.size = 0.5)

# Write granulocyte annotations onto matrix
renamed_Granulocyte_combined <- subset(renamed_Myeloid_combined, idents = c("Granulocyte"))
library(tibble)
cluster_letters <- Idents(object = renamed_Granulocyte_combined)
names(cluster_letters) <- colnames(x = renamed_Granulocyte_combined)
renamed_Granulocyte_combined <- AddMetaData(
  object = renamed_Granulocyte_combined,
  metadata = cluster_letters,
  col.name = 'labeled_ident'
)
Granulocyte_annotations <- renamed_Granulocyte_combined@meta.data %>% 
  rownames_to_column("barcodes")  %>% 
  dplyr::select("barcodes","labeled_ident")
Granulocyte_annotations <- Granulocyte_annotations %>% column_to_rownames(., var = "barcodes")

# Write Macrophage annotations into matrix
renamed_Macrophage_combined <- subset(renamed_Myeloid_combined, idents = c("Macrophage"))
DimPlot(renamed_Macrophage_combined, reduction = "umap", label = TRUE, pt.size = 0.5)
library(tibble)
cluster_letters <- Idents(object = renamed_Macrophage_combined)
names(cluster_letters) <- colnames(x = renamed_Macrophage_combined)
renamed_Macrophage_combined <- AddMetaData(
  object = renamed_Macrophage_combined,
  metadata = cluster_letters,
  col.name = 'labeled_ident'
)
Macrophage_annotations <- renamed_Macrophage_combined@meta.data %>% 
  rownames_to_column("barcodes")  %>% 
  dplyr::select("barcodes","labeled_ident")
Macrophage_annotations <- Macrophage_annotations %>% column_to_rownames(., var = "barcodes")


# Write CAF cluster annotations onto matrix (CAF_c1, CAF_c2, CAF_c3)
CAF.new.cluster.id <- c("CAF_c1", "CAF_c1", "CAF_c2", "CAF_c3", "CAF_c2", "CAF_c1", "CAF_c3")
names(CAF.new.cluster.id) <- levels(CAF_combined)
renamed_CAF_combined <- RenameIdents(CAF_combined, CAF.new.cluster.id)
DimPlot(renamed_CAF_combined, label = TRUE)
CAF_cluster_ident <- Idents(object = renamed_CAF_combined)
names(CAF_cluster_ident) <- colnames(x = renamed_CAF_combined)
renamed_CAF_combined <- AddMetaData(
  object = renamed_CAF_combined,
  metadata = CAF_cluster_ident,
  col.name = 'labeled_ident'
)
CAF_annotations <- renamed_CAF_combined@meta.data %>% 
  rownames_to_column("barcodes")  %>% 
  dplyr::select("barcodes","labeled_ident")
CAF_annotations <- CAF_annotations %>% column_to_rownames(., var = "barcodes")

library(tibble)  # for `rownames_to_column` and `column_to_rownames`

# CAF annotations
CAF_annotations_sn1 <- CAF_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(1:307) %>%
  column_to_rownames('barcode')

CAF_annotations_sn2 <- CAF_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(308:969) %>%
  column_to_rownames('barcode')

CAF_annotations_sn3 <- CAF_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(970:2088) %>%
  column_to_rownames('barcode')

CAF_annotations_sn4 <- CAF_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(2089:2349) %>%
  column_to_rownames('barcode')

CAF_annotations_sn5 <- CAF_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(2350:2585) %>%
  column_to_rownames('barcode')

CAF_annotations_sn6 <- CAF_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(2586:2787) %>%
  column_to_rownames('barcode')



# Granulocyte Anntoations
Granulocyte_annotations_sn1 <- Granulocyte_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(1:131) %>%
  column_to_rownames('barcode')

Granulocyte_annotations_sn2 <- Granulocyte_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(132:181) %>%
  column_to_rownames('barcode')

Granulocyte_annotations_sn3 <- Granulocyte_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(182:444) %>%
  column_to_rownames('barcode')

Granulocyte_annotations_sn4 <- Granulocyte_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(445:608) %>%
  column_to_rownames('barcode')

Granulocyte_annotations_sn5 <- Granulocyte_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(609:662) %>%
  column_to_rownames('barcode')

Granulocyte_annotations_sn6 <- Granulocyte_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(663:674) %>%
  column_to_rownames('barcode')


# Macrophage Annotations
Macrophage_annotations_sn1 <- Macrophage_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(1:500) %>%
  column_to_rownames('barcode')

Macrophage_annotations_sn2 <- Macrophage_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(501:689) %>%
  column_to_rownames('barcode')

Macrophage_annotations_sn3 <- Macrophage_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(690:1356) %>%
  column_to_rownames('barcode')

Macrophage_annotations_sn4 <- Macrophage_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(1357:1476) %>%
  column_to_rownames('barcode')

Macrophage_annotations_sn5 <- Macrophage_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(1477:1577) %>%
  column_to_rownames('barcode')

Macrophage_annotations_sn6 <- Macrophage_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(1578:1634) %>%
  column_to_rownames('barcode')


# Lymphoid Annotations
Lymphoid_annotations_sn1 <- Lymphoid_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(1:84) %>%
  column_to_rownames('barcode')

Lymphoid_annotations_sn2 <- Lymphoid_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(85:90) %>%
  column_to_rownames('barcode')

Lymphoid_annotations_sn3 <- Lymphoid_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(91:434) %>%
  column_to_rownames('barcode')

Lymphoid_annotations_sn4 <- Lymphoid_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(435:672) %>%
  column_to_rownames('barcode')

Lymphoid_annotations_sn5 <- Lymphoid_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(673:834) %>%
  column_to_rownames('barcode')

Lymphoid_annotations_sn6 <- Lymphoid_annotations %>%
  rownames_to_column('barcode') %>%
  dplyr::slice(835:836) %>%
  column_to_rownames('barcode')

Anno_CAF_Immune_sn1 <- rbind(CAF_annotations_sn1, Lymphoid_annotations_sn1, Granulocyte_annotations_sn1, Macrophage_annotations_sn1)
Anno_CAF_Immune_sn2 <- rbind(CAF_annotations_sn2, Lymphoid_annotations_sn2, Granulocyte_annotations_sn2, Macrophage_annotations_sn2)
Anno_CAF_Immune_sn3 <- rbind(CAF_annotations_sn3, Lymphoid_annotations_sn3, Granulocyte_annotations_sn3, Macrophage_annotations_sn3)
Anno_CAF_Immune_sn4 <- rbind(CAF_annotations_sn4, Lymphoid_annotations_sn4, Granulocyte_annotations_sn4, Macrophage_annotations_sn4)
Anno_CAF_Immune_sn5 <- rbind(CAF_annotations_sn5, Lymphoid_annotations_sn5, Granulocyte_annotations_sn5, Macrophage_annotations_sn5)
Anno_CAF_Immune_sn6 <- rbind(CAF_annotations_sn6, Lymphoid_annotations_sn6, Granulocyte_annotations_sn6, Macrophage_annotations_sn6)


# Function to convert mice genes to human homologs
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_mouse_to_human <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  return (output)
}



#sn1_4 (Control): Figure 6A, B; Figure S13A, B
Total.exp_CAF_Immune_sn1 <- PDAC_combined_raw_matrix[,rownames(Anno_CAF_Immune_sn1)]
iTalk_data_sn1 <- as.data.frame(t(Total.exp_CAF_Immune_sn1))
iTalk_data_sn1$cell_type <- Anno_CAF_Immune_sn1$labeled_ident
unique(iTalk_data_sn1$cell_type)

# find top highly expressed genes to conduct ligand-receptor analysis later
highly_exprs_genes_sn1<-rawParse(iTalk_data_sn1,top_genes=50,stats='mean')

# convert mice genes to human homologs
highly_exprs_genes_sn1$human_homolog <- NA
for (i in 1: nrow(highly_exprs_genes_sn1)) {
  homolog = convert_mouse_to_human(highly_exprs_genes_sn1[i, "gene"])
  if (is.null(homolog)) {
    highly_exprs_genes_sn1[i, "human_homolog"] = "none"
  } else {
    highly_exprs_genes_sn1[i, "human_homolog"] = homolog[1]
  }
}
writexl::write_xlsx(highly_exprs_genes_sn1, "./highly_exprs_genes_sn1.xlsx")
df_sn1 <- highly_exprs_genes_sn1
df_sn1 <- df_sn1[df_sn1$human_homolog != "none",] 
df_sn1['gene'] <- df_sn1['human_homolog']

library(ggsci)
# Color choice for circo plot
cell_col<-structure(c(pal_npg()(10), pal_igv()(9), pal_uchicago("light")(9), pal_futurama()(12), pal_aaas()(10)),names=levels(iTalk_data_sn1$cell_type))

#find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')

#Visualization for Manuscript 
res_sn1<-NULL
res_filtered_sn1 <- NULL
for(comm_type in comm_list){
  res_cat<-FindLR(df_sn1,datatype='mean count',comm_type = comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]

  # exclude interactions between different clusters of CAFs and from CD8 T cells to CAFs and granulocytes
  res_cat_filtered <- res_cat %>% filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c2")) %>% 
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CD8" & cell_to == "Granulocyte"))
  
  #top 30 ligand-receptor pairs
  jpeg(paste0("~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/TGFB-TRAP Manuscript V3 Edited 03282024/PhotoShop/TGFB-TRAP Figure 7/Figure 7 Individual/sn1_4_",comm_type,"_top_30.jpeg"),  width = 350, height = 350, units='mm', res = 300)
  par(cex = 2.95)
  LRPlot(res_cat_filtered[1:30,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat_filtered$cell_from_mean_exprs[1:30],link.arr.width=res_cat_filtered$cell_to_mean_exprs[1:30], track.height_1=circlize::uh(1,'mm'),track.height_2=circlize::uh(48,'mm'),text.vjust = '0.7cm')
  dev.off()
  res_sn1<-rbind(res_sn1,res_cat)
  res_filtered_sn1 <- rbind(res_filtered_sn1, res_cat_filtered)
}

# Unfiltered top ligand receptor pairs
writexl::write_xlsx(res_sn1, "./sn1_4_ligand_receptor_results.xlsx")

# Filtered top ligand receptor pairs
writexl::write_xlsx(res_filtered_sn1, "./sn1_4_ligand_receptor_results_filtered.xlsx")
#######







#sn2_4 (TGFB-TRAP): Figure 6C, d; Figure S13C, D
Total.exp_CAF_Immune_sn2 <- PDAC_combined_raw_matrix[,rownames(Anno_CAF_Immune_sn2)]
iTalk_data_sn2 <- as.data.frame(t(Total.exp_CAF_Immune_sn2))
iTalk_data_sn2$cell_type <- Anno_CAF_Immune_sn2$labeled_ident
unique(iTalk_data_sn2$cell_type)

# find top  percent highly expressed genes for ligand-receptor pair determination
highly_exprs_genes_sn2<-rawParse(iTalk_data_sn2,top_genes=50,stats='mean')

# convert mice genes to human homologs
highly_exprs_genes_sn2$human_homolog <- NA
for (i in 1: nrow(highly_exprs_genes_sn2)) {
  homolog = convert_mouse_to_human(highly_exprs_genes_sn2[i, "gene"])
  if (is.null(homolog)) {
    highly_exprs_genes_sn2[i, "human_homolog"] = "none"
  } else {
    highly_exprs_genes_sn2[i, "human_homolog"] = homolog[1]
  }
}
writexl::write_xlsx(highly_exprs_genes_sn2, "./highly_exprs_genes_sn2.xlsx")
df_sn2 <- highly_exprs_genes_sn2
df_sn2 <- df_sn2[df_sn2$human_homolog != "none",] 
df_sn2['gene'] <- df_sn2['human_homolog']

# find top highly expressed genes to conduct ligand-receptor analysis later
comm_list<-c('growth factor','other','cytokine','checkpoint')

# Color Palette for Circos plot
cell_col<-structure(c(pal_npg()(10), pal_igv()(9), pal_uchicago("light")(9), pal_futurama()(12), pal_aaas()(10)),names=levels(iTalk_data_sn2$cell_type))


#Visualization for Manuscript
res_sn2<-NULL
res_filtered_sn2 <- NULL
for(comm_type in comm_list){
  res_cat<-FindLR(df_sn2,datatype='mean count',comm_type = comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
 
  # exclude interactions between different clusters of CAFs and from CD8 T cells to CAFs and granulocytes
  res_cat_filtered <- res_cat %>% filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c2")) %>% 
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CD8" & cell_to == "Granulocyte"))
  #top 30 ligand-receptor pairs
  jpeg(paste0("~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/TGFB-TRAP Manuscript V3 Edited 03282024/PhotoShop/TGFB-TRAP Figure 7/Figure 7 Individual/sn2_4_",comm_type,"_top_30.jpeg"),  width = 350, height = 350, units='mm', res = 300)
  par(cex = 2.95)
  LRPlot(res_cat_filtered[1:30,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat_filtered$cell_from_mean_exprs[1:30],link.arr.width=res_cat_filtered$cell_to_mean_exprs[1:30], track.height_1=circlize::uh(1,'mm'),track.height_2=circlize::uh(48,'mm'),text.vjust = '0.7cm')
  dev.off()
  res_sn2<-rbind(res_sn2,res_cat)
  res_filtered_sn2 <- rbind(res_filtered_sn2, res_cat_filtered)
}

# Unfiltered top ligand receptor pairs
writexl::write_xlsx(res_sn2, "./sn2_4_ligand_receptor_results.xlsx")

# Filtered top ligand receptor pairs
writexl::write_xlsx(res_filtered_sn2, "./sn2_4_ligand_receptor_results_filtered.xlsx")

#######


#sn4_4 (anti-PD-1): Figure 6E, F; Figure S13E, F
Total.exp_CAF_Immune_sn4 <- PDAC_combined_raw_matrix[,rownames(Anno_CAF_Immune_sn4)]
iTalk_data_sn4 <- as.data.frame(t(Total.exp_CAF_Immune_sn4))
iTalk_data_sn4$cell_type <- Anno_CAF_Immune_sn4$labeled_ident
unique(iTalk_data_sn4$cell_type)

# find top highly expressed genes to conduct ligand-receptor analysis later
highly_exprs_genes_sn4<-rawParse(iTalk_data_sn4,top_genes=50,stats='mean')

# convert mice genes to human homologs
highly_exprs_genes_sn4$human_homolog <- NA
for (i in 1: nrow(highly_exprs_genes_sn4)) {
  homolog = convert_mouse_to_human(highly_exprs_genes_sn4[i, "gene"])
  if (is.null(homolog)) {
    highly_exprs_genes_sn4[i, "human_homolog"] = "none"
  } else {
    highly_exprs_genes_sn4[i, "human_homolog"] = homolog[1]
  }
}
writexl::write_xlsx(highly_exprs_genes_sn4, "./highly_exprs_genes_sn4.xlsx")
df_sn4 <- highly_exprs_genes_sn4
df_sn4 <- df_sn4[df_sn4$human_homolog != "none",] 
df_sn4['gene'] <- df_sn4['human_homolog']

#find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')

# color palette
cell_col<-structure(c(pal_npg()(10), pal_igv()(9), pal_uchicago("light")(9), pal_futurama()(12), pal_aaas()(10)),names=levels(iTalk_data_sn4$cell_type))

#Visualization for Manuscript
res_sn4<-NULL
res_filtered_sn4 <- NULL
for(comm_type in comm_list){
  res_cat<-FindLR(df_sn4,datatype='mean count',comm_type = comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
  
  # exclude interactions between different clusters of CAFs and from CD8 T cells to CAFs and granulocytes
  res_cat_filtered <- res_cat %>% filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c2")) %>% 
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CD8" & cell_to == "Granulocyte"))
  #top 30 ligand-receptor pairs
  jpeg(paste0("~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/TGFB-TRAP Manuscript V3 Edited 03282024/PhotoShop/TGFB-TRAP Figure 7/Figure 7 Individual/sn4_4_",comm_type,"_top_30.jpeg"),  width = 350, height = 350, units='mm', res = 300)
  par(cex = 2.95)
  LRPlot(res_cat_filtered[1:30,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat_filtered$cell_from_mean_exprs[1:30],link.arr.width=res_cat_filtered$cell_to_mean_exprs[1:30], track.height_1=circlize::uh(1,'mm'),track.height_2=circlize::uh(48,'mm'),text.vjust = '0.7cm')
  dev.off()
  res_sn4<-rbind(res_sn4,res_cat)
  res_filtered_sn4 <- rbind(res_filtered_sn4, res_cat_filtered)
}

# Unfiltered top ligand receptor pairs
writexl::write_xlsx(res_sn4, "./sn4_4_ligand_receptor_results.xlsx")

# Filtered top ligand receptor pairs
writexl::write_xlsx(res_filtered_sn4, "./sn4_4_ligand_receptor_results_filtered.xlsx")

#######


#sn5_4 (anti-PD-1 + TGFB-TRAP): Figure 6G, H; Figure S13G, H
Total.exp_CAF_Immune_sn5 <- PDAC_combined_raw_matrix[,rownames(Anno_CAF_Immune_sn5)]
iTalk_data_sn5 <- as.data.frame(t(Total.exp_CAF_Immune_sn5))
iTalk_data_sn5$cell_type <- Anno_CAF_Immune_sn5$labeled_ident
unique(iTalk_data_sn5$cell_type)

# find top highly expressed genes to conduct ligand-receptor analysis later
highly_exprs_genes_sn5<-rawParse(iTalk_data_sn5,top_genes=50,stats='mean')

# convert mice genes to human homologs
highly_exprs_genes_sn5$human_homolog <- NA
for (i in 1: nrow(highly_exprs_genes_sn5)) {
  homolog = convert_mouse_to_human(highly_exprs_genes_sn5[i, "gene"])
  if (is.null(homolog)) {
    highly_exprs_genes_sn5[i, "human_homolog"] = "none"
  } else {
    highly_exprs_genes_sn5[i, "human_homolog"] = homolog[1]
  }
}
writexl::write_xlsx(highly_exprs_genes_sn5, "./highly_exprs_genes_sn5.xlsx")
df_sn5 <- highly_exprs_genes_sn5
df_sn5 <- df_sn5[df_sn5$human_homolog != "none",] 
df_sn5['gene'] <- df_sn5['human_homolog']

#find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')

# Color palette
cell_col<-structure(c(pal_npg()(10), pal_igv()(9), pal_uchicago("light")(9), pal_futurama()(12), pal_aaas()(10)),names=levels(iTalk_data_sn5$cell_type))

#Visualization for Manuscript
res_sn5<-NULL
res_filtered_sn5 <- NULL
for(comm_type in comm_list){
  res_cat<-FindLR(df_sn5,datatype='mean count',comm_type = comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
  
  # exclude interactions between different clusters of CAFs and from CD8 T cells to CAFs and granulocytes
  res_cat_filtered <- res_cat %>% filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CAF_c1" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CAF_c2" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c2")) %>% 
    filter(!(cell_from == "CAF_c3" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c1")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c2")) %>%
    filter(!(cell_from == "CD8" & cell_to == "CAF_c3")) %>%
    filter(!(cell_from == "CD8" & cell_to == "Granulocyte"))
  #top 30 ligand-receptor pairs
  jpeg(paste0("~/Desktop/LZ01JHU520_IncludeIntrons/Manuscript Iterations/TGFB-TRAP Manuscript V3 Edited 03282024/PhotoShop/TGFB-TRAP Figure 7/Figure 7 Individual/sn5_4_",comm_type,"_top_30.jpeg"),  width = 350, height = 350, units='mm', res = 300)
  par(cex = 2.95)
  LRPlot(res_cat_filtered[1:30,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat_filtered$cell_from_mean_exprs[1:30],link.arr.width=res_cat_filtered$cell_to_mean_exprs[1:30], track.height_1=circlize::uh(1,'mm'),track.height_2=circlize::uh(48,'mm'),text.vjust = '0.7cm')
  dev.off()
  res_sn5<-rbind(res_sn5,res_cat)
  res_filtered_sn5 <- rbind(res_filtered_sn5, res_cat_filtered)
}

# Unfiltered top ligand receptor pairs
writexl::write_xlsx(res_sn5, "./sn5_4_ligand_receptor_results.xlsx")

# Filtered top ligand receptor pairs
writexl::write_xlsx(res_filtered_sn5, "./sn5_4_ligand_receptor_results_filtered.xlsx")
#######











