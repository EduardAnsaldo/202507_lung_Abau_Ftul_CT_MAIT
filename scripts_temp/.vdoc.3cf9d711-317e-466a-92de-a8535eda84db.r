#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# General R and plotting 
library(conflicted)
library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)
library(cowplot)
library(gridExtra)
library(ggrepel)
library(stringr)
library(VennDiagram)
library(pheatmap)
library(viridis)
library(here)
library(readxl)

# Single Cell Analysis Packages
library(Seurat)
library(scRepertoire)
library(circlize)
library(scCustomize)
library(SingleR)
library(celldex)
library(UCell)
library(presto)
library(scDblFinder)
library(immunarch)

# DEG, pathway enrichment and visualization packages
library(DESeq2)
library(gprofiler2)
library(clusterProfiler)
library(DOSE)
library(pathview)
library(org.Mm.eg.db)
library(scRepertoire)
library(enrichplot)
library(msigdbr)
library(rlang)

# Python Version 3.9.6!!! 

knitr::dep_prev()
# peconflicts_scout()
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(base::setdiff)
conflicts_prefer(dplyr::rename)

#
#
#
#| results: hide 
# Setting working directory and seed
# renv::update()
# renv::snapshot()
# renv::init()
# renv::status()
set.seed(3514)
i_am('scripts/Analysis.qmd')
here()
results_path <- here('results/')
dir.create(results_path)
results_path <- here('results/')
dir.create(results_path)
figures_path <- here('results/figures/')
dir.create(figures_path)
# data_path <- here('data/')
# dir.create(data_path)
tables_path <- here('results/tables/')
dir.create(tables_path)

# Loading custom functions
source(here('scripts/function_template.r'))
source(here('scripts/gProfiler2_functions.r'))
source(here('scripts/circlize_functions.r'))

# Setting up color palettes
diverging_palette <- hcl.colors(n = 20,'Purple-Green',rev = T)
sequential_palette_dotplot <- hcl.colors(n = 20,'YlGn',rev = T)
sequential_palette <- hcl.colors(n = 20,'Mako',rev = T)
#
#
#
#
#
#
#
#
#
#
#
#

#| results: hide
#| #Load the dataset from the cellranger outs
#| results: hide
#Load the dataset from the cellranger outs
scdata <- Read10X(data.dir = here("data/cluster_processing_2/aggr_all/outs/count/filtered_feature_bc_matrix"))

#Initialize the seurat object with the raw (non-normalized data)
seurat <- CreateSeuratObject(counts=scdata$'Gene Expression', min.cells = 3)

# #Add HTO data as a new assay independent from RNA
HTO <- CreateAssayObject(counts = scdata$'Antibody Capture')
seurat[["HTO"]] <- HTO

#
#
#
#
#
#
#
#
#

seurat$barcode <- Cells(seurat)
seurat@meta.data <- seurat@meta.data |>
    mutate(GEM = str_extract(barcode, '\\d')) |>
    mutate(GEM = case_when(
            GEM == '1' ~ 'CT_MAIT',
            GEM == '2' ~ 'Ftula_MAIT',
            GEM == '3' ~ 'Abau_MAIT'
    )) 
# Splitting the dataset into the three groups
seurat_list <- SplitObject(seurat, split.by = 'GEM')

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| results: hide

seurat <- seurat_list$Abau_MAIT

#
#
#
#
#
#
#
print('Antibody Capture total reads:')
seurat@assays$HTO$counts |> as.matrix() |>
        rowSums()  |> as.data.frame() 

HTO_data <- seurat@assays$HTO$counts

# Extracting HTO and biotin data
HTO_data <-  HTO_data[-c(3), ]
# colnames(scdata$'HTO') <- colnames(scdata$'Antibody Capture')
# rownames(scdata$'HTO') <- rownames(scdata$'Antibody Capture')[-c(2, 3, 4, 5, 6, 7, 11, 17, 18, 19)]

# #Add HTO data as a new assay independent from RNA
HTO <- CreateAssayObject(counts = HTO_data)
seurat[["HTO"]] <- HTO
#
#
#
#
#
#
#
#
#
#
#
#| layout-ncol: 2
#Quantifying percentage mitochondria
Idents(seurat) <- 'Tissue'
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, assay="RNA",pattern = "mt-")

VlnPlot(seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol=3, pt.size = 0.01)

# Setting thresholds:
nCount_RNA_threshold <- 1500
nFeature_RNA_threshold <- 500
percent_mt_threshold <- 6


p1 <- seurat@meta.data %>% 
  	ggplot(aes(x=nCount_RNA)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = nCount_RNA_threshold)

p2 <- seurat@meta.data %>% 
  	ggplot(aes(x=nFeature_RNA)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = nFeature_RNA_threshold)

p3 <- seurat@meta.data %>% 
  	ggplot(aes(x=percent.mt)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = percent_mt_threshold)

grid.arrange(p1, p2, p3)    

# # Visualize feature relationships
# # plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")

# # plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # plot1+plot2

cell_number <- list(nrow(seurat[[]]))
#
#
#
#
#
# Filter data

cell_number <- list(nrow(seurat[[]]))
seurat <- subset(seurat, subset = nFeature_RNA  > nFeature_RNA_threshold & percent.mt < percent_mt_threshold & nCount_RNA > nCount_RNA_threshold)
cell_number <- append(cell_number, nrow(seurat[[]]))

print("Cell number before and after filterings")
print(cell_number)
#
#
#
#
#
#
#
# Normalize RNA data with SCTransform
# seurat <- SCTransform(seurat, verbose = T)

# Normalize HTO data with CLR
seurat <- NormalizeData(seurat, assay = "HTO",normalization.method = "CLR")
#
#
#
#
#
#
#
#| fig-width: 14
#| fig-height: 10
plot1 <- seurat@meta.data %>% 
  	ggplot(aes(x=nCount_HTO)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") 

plot2 <- seurat@meta.data %>% 
  	ggplot(aes(x=nFeature_HTO)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") 

plot1+plot2

seurat <- HTODemux(seurat, assay = "HTO", positive.quantile = 0.99)

table(seurat$HTO_classification.global)

#tSNE-visualization
seurat_subset <- subset(seurat, idents = "Negative", invert =TRUE)
DefaultAssay(seurat_subset) <- "HTO"
seurat_subset <- ScaleData(seurat_subset, assay = "HTO",  features = rownames(seurat_subset), verbose = FALSE )
seurat_subset <- RunPCA(seurat_subset, assay = "HTO", rownames(seurat_subset), approx = FALSE)
seurat_subset <- RunTSNE(seurat_subset, assay = "HTO", dims = 1:24,  perplexity = 100, check_duplicates=FALSE)
DimPlot(seurat_subset)
#
#
#
#
print("Number of Singlets")
print(nrow(seurat@meta.data |> filter( HTO_classification.global  == 'Singlet')))

print("Number of Doublets")
print(nrow(seurat@meta.data |> filter( HTO_classification.global  == 'Doublet')))

# Filtering singlets only
Idents(seurat) <- "HTO_classification.global"
seurat <-subset(seurat, idents = "Singlet")
#
#
#
#
#
#
#

seurat@meta.data <- seurat@meta.data |>
    mutate(Groups = 'Abau') |>
    mutate(Samples = str_c(Groups, hash.ID, sep = '_')) |> 
    mutate(Samples = case_when(
            hash.ID == 'Mouse-4' ~ 'Abau_Mouse-3',
            hash.ID == 'Mouse-5' ~ 'Abau_Mouse-4',
            TRUE ~ Samples
    )) |>
    mutate(Samples = factor(Samples, levels = c(
            'Abau_Mouse-1',
            'Abau_Mouse-2',
            'Abau_Mouse-3',
            'Abau_Mouse-4')))

#
#
#
seurat_Abau <- seurat
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| results: hide

seurat <- seurat_list$Ftula_MAIT

#
#
#
#
#
#
#
print('Antibody Capture total reads:')
seurat@assays$HTO$counts |> as.matrix() |>
        rowSums()  |> as.data.frame() 

HTO_data <- seurat@assays$HTO$counts

# Extracting HTO and biotin data
# HTO_data <-  HTO_data[-c(3), ]
# colnames(scdata$'HTO') <- colnames(scdata$'Antibody Capture')
# rownames(scdata$'HTO') <- rownames(scdata$'Antibody Capture')[-c(2, 3, 4, 5, 6, 7, 11, 17, 18, 19)]

# #Add HTO data as a new assay independent from RNA
HTO <- CreateAssayObject(counts = HTO_data)
seurat[["HTO"]] <- HTO
#
#
#
#
#
#
#
#
#
#
#
#| layout-ncol: 2
#Quantifying percentage mitochondria
Idents(seurat) <- 'Tissue'
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, assay="RNA",pattern = "mt-")

VlnPlot(seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol=3, pt.size = 0.01)

# Setting thresholds:
nCount_RNA_threshold <- 1500
nFeature_RNA_threshold <- 500
percent_mt_threshold <- 6


p1 <- seurat@meta.data %>% 
  	ggplot(aes(x=nCount_RNA)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = nCount_RNA_threshold)

p2 <- seurat@meta.data %>% 
  	ggplot(aes(x=nFeature_RNA)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = nFeature_RNA_threshold)

p3 <- seurat@meta.data %>% 
  	ggplot(aes(x=percent.mt)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = percent_mt_threshold)

grid.arrange(p1, p2, p3)    

# # Visualize feature relationships
# # plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")

# # plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # plot1+plot2

cell_number <- list(nrow(seurat[[]]))
#
#
#
#
#
# Filter data

cell_number <- list(nrow(seurat[[]]))
seurat <- subset(seurat, subset = nFeature_RNA  > nFeature_RNA_threshold & percent.mt < percent_mt_threshold & nCount_RNA > nCount_RNA_threshold)
cell_number <- append(cell_number, nrow(seurat[[]]))

print("Cell number before and after filterings")
print(cell_number)
#
#
#
#
#
#
#
# Normalize RNA data with SCTransform
# seurat <- SCTransform(seurat, verbose = T)

# Normalize HTO data with CLR
seurat <- NormalizeData(seurat, assay = "HTO",normalization.method = "CLR")
#
#
#
#
#
#
#
#| fig-width: 14
#| fig-height: 10
plot1 <- seurat@meta.data %>% 
  	ggplot(aes(x=nCount_HTO)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") 

plot2 <- seurat@meta.data %>% 
  	ggplot(aes(x=nFeature_HTO)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") 

plot1+plot2

seurat <- HTODemux(seurat, assay = "HTO", positive.quantile = 0.97)

table(seurat$hash.ID)
table(seurat$HTO_classification.global)

#tSNE-visualization
seurat_subset <- subset(seurat, idents = "Negative", invert =TRUE)
DefaultAssay(seurat_subset) <- "HTO"
seurat_subset <- ScaleData(seurat_subset, assay = "HTO",  features = rownames(seurat_subset), verbose = FALSE )
seurat_subset <- RunPCA(seurat_subset, assay = "HTO", rownames(seurat_subset), approx = FALSE)
seurat_subset <- RunTSNE(seurat_subset, assay = "HTO", dims = 1:24,  perplexity = 100, check_duplicates=FALSE)
DimPlot(seurat_subset)
#
#
#
#
print("Number of Singlets")
print(nrow(seurat@meta.data |> filter( HTO_classification.global  == 'Singlet')))

print("Number of Doublets")
print(nrow(seurat@meta.data |> filter( HTO_classification.global  == 'Doublet')))

# Filtering singlets only
Idents(seurat) <- "HTO_classification.global"
seurat <-subset(seurat, idents = "Singlet")
#
#
#
#
#
#
#

seurat@meta.data <- seurat@meta.data |>
    mutate(Groups = 'Ftula') |>
    mutate(Samples = str_c(Groups, hash.ID, sep = '_')) |> 
    mutate(Samples = factor(Samples, levels = c(
            'Ftula_Mouse-1',
            'Ftula_Mouse-2',
            'Ftula_Mouse-3',
            'Ftula_Mouse-4',
            'Ftula_Mouse-5')))

#
#
#
seurat_Ftula <- seurat
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| results: hide

seurat <- seurat_list$CT_MAIT

#
#
#
#
#
#
#
print('Antibody Capture total reads:')
seurat@assays$HTO$counts |> as.matrix() |>
        rowSums()  |> as.data.frame() 

HTO_data <- seurat@assays$HTO$counts

# Extracting HTO and biotin data
HTO_data <-  HTO_data[-c(3), ]
# colnames(scdata$'HTO') <- colnames(scdata$'Antibody Capture')
# rownames(scdata$'HTO') <- rownames(scdata$'Antibody Capture')[-c(2, 3, 4, 5, 6, 7, 11, 17, 18, 19)]

# #Add HTO data as a new assay independent from RNA
HTO <- CreateAssayObject(counts = HTO_data)
seurat[["HTO"]] <- HTO
#
#
#
#
#
#
#
#
#
#
#
#| layout-ncol: 2
#Quantifying percentage mitochondria
Idents(seurat) <- 'Tissue'
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, assay="RNA",pattern = "mt-")

VlnPlot(seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol=3, pt.size = 0.01)

# Setting thresholds:
nCount_RNA_threshold <- 1500
nFeature_RNA_threshold <- 500
percent_mt_threshold <- 6


p1 <- seurat@meta.data %>% 
  	ggplot(aes(x=nCount_RNA)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = nCount_RNA_threshold)

p2 <- seurat@meta.data %>% 
  	ggplot(aes(x=nFeature_RNA)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = nFeature_RNA_threshold)

p3 <- seurat@meta.data %>% 
  	ggplot(aes(x=percent.mt)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = percent_mt_threshold)

grid.arrange(p1, p2, p3)    

# # Visualize feature relationships
# # plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")

# # plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # plot1+plot2

cell_number <- list(nrow(seurat[[]]))
#
#
#
#
#
# Filter data

cell_number <- list(nrow(seurat[[]]))
seurat <- subset(seurat, subset = nFeature_RNA  > nFeature_RNA_threshold & percent.mt < percent_mt_threshold & nCount_RNA > nCount_RNA_threshold)
cell_number <- append(cell_number, nrow(seurat[[]]))

print("Cell number before and after filterings")
print(cell_number)
#
#
#
#
#
#
#
# Normalize RNA data with SCTransform
# seurat <- SCTransform(seurat, verbose = T)

# Normalize HTO data with CLR
seurat <- NormalizeData(seurat, assay = "HTO",normalization.method = "CLR")
#
#
#
#
#
#
#
#| fig-width: 14
#| fig-height: 10
plot1 <- seurat@meta.data %>% 
  	ggplot(aes(x=nCount_HTO)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") 

plot2 <- seurat@meta.data %>% 
  	ggplot(aes(x=nFeature_HTO)) + 
  	geom_density() + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") 

plot1+plot2

seurat <- HTODemux(seurat, assay = "HTO", positive.quantile = 0.99)

table(seurat$hash.ID)
table(seurat$HTO_classification.global)

#tSNE-visualization
seurat_subset <- subset(seurat, idents = "Negative", invert =TRUE)
DefaultAssay(seurat_subset) <- "HTO"
seurat_subset <- ScaleData(seurat_subset, assay = "HTO",  features = rownames(seurat_subset), verbose = FALSE )
seurat_subset <- RunPCA(seurat_subset, assay = "HTO", rownames(seurat_subset), approx = FALSE)
seurat_subset <- RunTSNE(seurat_subset, assay = "HTO", dims = 1:24,  perplexity = 100, check_duplicates=FALSE)
DimPlot(seurat_subset)
#
#
#
#
print("Number of Singlets")
print(nrow(seurat@meta.data |> filter( HTO_classification.global  == 'Singlet')))

print("Number of Doublets")
print(nrow(seurat@meta.data |> filter( HTO_classification.global  == 'Doublet')))

# Filtering singlets only
Idents(seurat) <- "HTO_classification.global"
seurat <-subset(seurat, idents = "Singlet")
#
#
#
#
#
#
#

seurat@meta.data <- seurat@meta.data |>
    mutate(Groups = 'CT') |>
    mutate(Samples = str_c(Groups, hash.ID, sep = '_')) |> 
    mutate(Samples=case_when(
            hash.ID == 'Mouse-4' ~ 'CT_Mouse-3',
            hash.ID == 'Mouse-5' ~ 'CT_Mouse-4',
            TRUE ~ Samples
    )) |>  
    mutate(Samples = factor(Samples, levels = c(
            'CT_Mouse-1',
            'CT_Mouse-2',
            'CT_Mouse-3',
            'CT_Mouse-4')))

#
#
#
seurat_CT <- seurat
#
#
#
#
#
#
#
data_list <- list(CT_MAIT = seurat_CT, Abau_MAIT = seurat_Abau, Ftula_MAIT = seurat_Ftula)
# Merging the datasets
seurat <- Merge_Seurat_List(data_list, 
                       add.cell.ids = names(data_list),
                       project = "MAIT_Abau_Ftula_CT_Lungs")

seurat <- JoinLayers(seurat, assay = 'RNA')

#
#
#
#
#
#
#
#
#

seurat@meta.data <- seurat@meta.data |>
    mutate(Groups = factor(Groups, levels = c('CT', 'Abau', 'Ftula'))) |>
    mutate(Samples = factor(Samples, levels = c(
            'CT_Mouse-1',
            'CT_Mouse-2',
            'CT_Mouse-3',
            'CT_Mouse-4',
            'Abau_Mouse-1',
            'Abau_Mouse-2',
            'Abau_Mouse-3',
            'Abau_Mouse-4',
            'Ftula_Mouse-1',
            'Ftula_Mouse-2',
            'Ftula_Mouse-3',
            'Ftula_Mouse-4',
            'Ftula_Mouse-5')))

seurat <- RenameCells(seurat, new.names=seurat$barcode)
#
#
#
#
#
#
#
#
#
contigs <- read.csv(here('data/cluster_processing_2/aggr_all/outs/vdj_t/filtered_contig_annotations.csv'))

#Demultiplexing VDJ libraries
contig_list <- createHTOContigList(contigs, seurat, group.by = "Samples")
contig_list <- contig_list[c(levels(seurat$Samples))]
# Combining the Contigs
combined <- combineTCR(contig_list, removeNA = TRUE, samples=c(levels(seurat$Samples)), filterMulti = TRUE)

# # #Adding Samples as variables
combined <- addVariable(combined, variable.name='Samples', variables=c(unique(seurat$Samples)))                                    

## Integrating data with seurat object
cell_names <- Cells(seurat)
groups_merge <- pull(seurat@meta.data, Samples)

#Changing barcodes on the seurat object for the merge with VDJ data
new_cell_names <- paste(groups_merge, cell_names, sep='_')
seurat$barcode <- cell_names
seurat <- RenameCells(seurat, new.names=new_cell_names)

# #Combining VDJ and Seurat Object
seurat <- combineExpression(combined, seurat, proportion = TRUE, cloneCall='aa', group.by='Samples', chain = 'both')
# Renaming cells back to original barcodes
seurat <- RenameCells(seurat, new.names=seurat$barcode)
#
#
#
#
#
#
#
#
#
# Normalize RNA data with SCTransform
seurat <- SCTransform(seurat, verbose = F) 
# seurat@assays$SCT
# seurat <- PrepSCTFindMarkers(seurat)
#
#
#
#
#
#
#
#
#
#| fig-width: 4
#| fig-height: 4
#| fig-cap: 'Elbow plot showing the percentage of variability represented by by each PC to select number of dimensions '

#Dimensionality reduction
seurat <- RunPCA(seurat,npcs = 100)

#Determining dimensionality of the dataset
ElbowPlot(seurat, ndims = 100)
#
#
#
#
#
#
#
#
#
# Describe number of dimensions

dimensions  <- 25 
dimensions
#
#
#
#| fig-cap: UMAP plots showing the clustering results at different resolutions
#| fig-width: 10
#| fig-height: 10

seurat <- RunUMAP(seurat, dims=1:dimensions, verbose = F, seed.use = 3514L)
# resolutions <- c(0.25, 0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5, 1.75)
resolutions <- c(0.05, 0.1, 0.15, 0.25, 0.375, 0.5)
seurat <- FindNeighbors(seurat, dims = 1:dimensions, verbose = F)
seurat <- FindClusters(seurat, resolution = resolutions, verbose = F)

p <- list()
i <- 1
for (resolution in resolutions ) {
    Idents(seurat) <- paste0('SCT_snn_res.', resolution) 
    p2 <- DimPlot(seurat, reduction = "umap", label = TRUE) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")    
    p[[i]] <- p2
    i <- i+1
}

plot <- grid.arrange(grobs = p)
# ggsave('initial_clustering_results_by_resolution.pdf', path = here('result/figures'), plot = plot)
#
#
#
#
#
#
#
#
#
resolution <- 0.375
resolution
#
#
#
object_annotations <- 'full_object_merged_pre_filtering'
#
#
#
#
#
#| layout-ncol: 2
Idents(seurat) <- paste0('SCT_snn_res.', resolution)
seurat[['seurat_clusters']]<- Idents(seurat)
DimPlot(seurat, reduction = "umap", label = TRUE) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
ggsave(paste0('UMAP_clusters_R_', resolution, '_', object_annotations, '.pdf'), path = figures_path)
DimPlot(seurat, reduction = "umap", label = TRUE, split.by = 'Groups', ncol = 2) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
ggsave(paste0('UMAP_clusters_by_group_R_', resolution, '_', object_annotations, '.pdf'), path = figures_path, width = 5, height = 5)
# DimPlot(seurat, reduction = "umap", label = TRUE, split.by = 'Samples', ncol = 4) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
# ggsave(paste0('UMAP_clusters_by_sample_R_', resolution, '_', object_annotations, '.pdf'), path = figures_path)
DimPlot(seurat, reduction = "umap", label = FALSE, group.by = 'Groups') + ggtitle(paste0('R ', resolution)) + theme(legend.position = "right")
#
#
#
#
#
#| layout-ncol: 3
# seurat <- PrepSCTFindMarkers(seurat)
plot1 <- top_genes_per_cluster(seurat, n_genes_to_plot = 3, object_annotations = object_annotations, tables_path = tables_path, figures_path = figures_path, results_path = results_path, run_pathway_enrichment = F) 
#
#
#
#| fig-width: 8
#| fig-height: 12
print(plot1)
ggsave(paste0(figures_path, 'DotPlot_Top3_per_cluster', object_annotations, '.pdf'), width = 8, height = 12)
```
#
#
#
#
#
#
#
local_path <- paste0(figures_path, object_annotations,'_cell_type_annotations')
unlink(local_path,recursive = T)
dir.create(local_path)

# Normalize and scale data
# seurat <- JoinLayers(seurat, assay = 'RNA')
seurat <- NormalizeData(seurat, assay = "RNA", normalization.method = "LogNormalize", scale.fct = 10000)
seurat <- ScaleData(seurat, assay = 'RNA')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| fig-cap: UMAP plot showing the clustering results with SingleR fine annotations on a per cluster basis
#| fig-width: 9
#| fig-height: 6
seurat <- annotate_seurat_with_SingleR_Eduard(seurat, local_path, database = 'ImmGen', annotation_basis = 'cluster_coarse', split_by_groups = FALSE)
#
#
#
#
#
#
#
#
#
#| fig-cap: UMAP plot showing the clustering results with SingleR coarse annotations on a per cell basis

annotate_seurat_with_SingleR_Eduard(seurat, local_path, database = 'ImmGen', annotation_basis = 'cell_coarse', split_by_groups = FALSE)
#
#
#
#
#
#| fig-cap: Feature plots showing the expression of T cell markers
#| fig-width: 8
#| fig-height: 10
FeaturePlot_scCustom(seurat, features = c('Cd3e', 'Cd4', 'Cd8a', 'Cd8b1', 'Trdc', 'Trac'), pt.size = 0.01, colors_use = sequential_palette) +
    theme(legend.position = 'none')
#
#
#
#| fig-cap: Feature plots showing the expression of myeloid cell markers
#| fig-width: 8
#| fig-height: 10
FeaturePlot_scCustom(seurat, features = c('Sirpa', 'Xcr1', 'Itgae', 'Itgax', 'Csf1r', 'Csf2ra'), pt.size = 0.1, colors_use = sequential_palette) +
    theme(legend.position = 'none')
#
#
#
#| fig-cap: Feature plots showing the expression of neuron Markers
#| #| fig-width: 8
#| fig-height: 13
FeaturePlot_scCustom(seurat, features = c('Tubb3', 'Snap25', 'Gad1', 'Slc17a7'), pt.size = 0.1, colors_use = sequential_palette, num_columns = 2) 
#
#
#
#
#
#
#
Idents(seurat) <- 'seurat_clusters'
seurat <- subset(seurat, idents = c('10'), invert = TRUE)
DimPlot(seurat, reduction = "umap", label = TRUE) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
#
#
#
#
#
#
#
#
#
# Normalize RNA data with SCTransform
seurat <- SCTransform(seurat, verbose = F) 
# seurat@assays$SCT
# seurat <- PrepSCTFindMarkers(seurat)
#
#
#
#
#
#
#
#
#
#
#| fig-width: 4
#| #| fig-height: 4
#| fig-cap: 'Elbow plot showing the percentage of variability represented by by each PC to select number of dimensions '

seurat <- quietVDJgenes(seurat)

#Dimensionality reduction
seurat <- RunPCA(seurat,npcs = 100)

#Determining dimensionality of the dataset
ElbowPlot(seurat, ndims = 100)
#
#
#
#
#
#
#
#
#
# Describe number of dimensions

dimensions  <- 25 
dimensions
#
#
#
#| fig-cap: UMAP plots showing the clustering results at different resolutions
#| fig-width: 10
#| fig-height: 10

seurat <- RunUMAP(seurat, dims=1:dimensions, verbose = F, seed.use = 3514L)
# resolutions <- c(0.25, 0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5, 1.75)
resolutions <- c(0.05, 0.1, 0.15, 0.175, 0.2, 0.225, 0.25, 0.375, 0.5)
seurat <- FindNeighbors(seurat, dims = 1:dimensions, verbose = F)
seurat <- FindClusters(seurat, resolution = resolutions, verbose = F)

p <- list()
i <- 1
for (resolution in resolutions ) {
    Idents(seurat) <- paste0('SCT_snn_res.', resolution) 
    p2 <- DimPlot(seurat, reduction = "umap", label = TRUE) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")    
    p[[i]] <- p2
    i <- i+1
}

plot <- grid.arrange(grobs = p)
# ggsave('initial_clustering_results_by_resolution.pdf', path = here('result/figures'), plot = plot)
#
#
#
#
#
#
#
#
#
resolution <- 0.15
resolution
#
#
#
#
#
#
#
#
#
saveRDS(seurat, file = here('data/MAIT_Abau_Ftula_CT_Lungs.rds'))
#
#
#
#
#
seurat <- readRDS(here('data/MAIT_Abau_Ftula_CT_Lungs.rds'))
object_annotations <- 'full_object_filtered'
resolution <- 0.15
#
#
#
#
#
#| layout-ncol: 2
Idents(seurat) <- paste0('SCT_snn_res.', resolution)
seurat[['seurat_clusters']]<- Idents(seurat)
DimPlot(seurat, reduction = "umap", label = TRUE) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
ggsave(paste0('UMAP_clusters_R_', resolution, '_', object_annotations, '.pdf'), path = figures_path)
DimPlot(seurat, reduction = "umap", label = TRUE, split.by = 'Groups', ncol = 2) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
ggsave(paste0('UMAP_clusters_by_group_R_', resolution, '_', object_annotations, '.pdf'), path = figures_path, width = 5, height = 5)
DimPlot(seurat, reduction = "umap", label = TRUE, split.by = 'Samples', ncol = 4) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
ggsave(paste0('UMAP_clusters_by_sample_R_', resolution, '_', object_annotations, '.pdf'), path = figures_path)
DimPlot(seurat, reduction = "umap", label = FALSE, group.by = 'Groups') + ggtitle(paste0('R ', resolution)) + theme(legend.position = "right")
#
#
#
#
#
extract_cell_counts(seurat, seurat_clusters, figures_path, tables_path, object_annotations=object_annotations)
#
#
#
#
#
#| layout-ncol: 3
# seurat <- PrepSCTFindMarkers(seurat)
plot1 <- top_genes_per_cluster(seurat, n_genes_to_plot=5, object_annotations, tables_path = tables_path, figures_path = figures_path, results_path = results_path, run_pathway_enrichment = T) 
#
#
#
#| fig-width: 8
#| fig-height: 12
print(plot1)
ggsave(paste0(figures_path, 'DotPlot_Top3_per_cluster', object_annotations, '.pdf'), width = 8, height = 12)
#
#
#
#
#
#| layout-ncol: 2
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, assay="RNA",pattern = "^Rp[sl]")
FeaturePlot_scCustom(seurat, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo'), colors_use = sequential_palette)
FeaturePlot_scCustom(seurat, features = c('Zbtb16', '', 'Icos'), colors_use = sequential_palette)
VlnPlot(seurat, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo'), group.by = 'seurat_clusters', pt.size = 0.1) + NoLegend()
VlnPlot(seurat, features = c('Zbtb16', 'Ubc', 'Icos'), group.by = 'seurat_clusters', pt.size = 0.1) + NoLegend()
FeatureScatter(seurat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', split.by = 'seurat_clusters', ncol = 3)
#
#
#
#
#
#
#
#
#
#| fig-width: 4
#| #| fig-height: 4
#| fig-cap: 'Elbow plot showing the percentage of variability represented by by each PC to select number of dimensions '

# Regressing out ribosomal protein content
seurat <- SCTransform(seurat, vars.to.regress = c('percent.ribo', 'percent.mt'), verbose = T)
seurat <- quietVDJgenes(seurat)

#Dimensionality reduction
seurat <- RunPCA(seurat,npcs = 100)

#Determining dimensionality of the dataset
ElbowPlot(seurat, ndims = 100)
#
#
#
#
#
#
#
#
#
# Describe number of dimensions

dimensions  <- 25 
dimensions
#
#
#
#| fig-cap: UMAP plots showing the clustering results at different resolutions
#| fig-width: 10
#| fig-height: 10

seurat <- RunUMAP(seurat, dims=1:dimensions, verbose = F, seed.use = 3514L)
# resolutions <- c(0.25, 0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5, 1.75)
resolutions <- c(0.05, 0.1, 0.15, 0.175, 0.2, 0.225, 0.25, 0.375, 0.5)
seurat <- FindNeighbors(seurat, dims = 1:dimensions, verbose = F)
seurat <- FindClusters(seurat, resolution = resolutions, verbose = F)

p <- list()
i <- 1
for (resolution in resolutions ) {
    Idents(seurat) <- paste0('SCT_snn_res.', resolution) 
    p2 <- DimPlot(seurat, reduction = "umap", label = TRUE) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")    
    p[[i]] <- p2
    i <- i+1
}

plot <- grid.arrange(grobs = p)
# ggsave('initial_clustering_results_by_resolution.pdf', path = here('result/figures'), plot = plot)
#
#
#
#
#
#
#
#
#
resolution <- 0.15
resolution
#
#
#
#
#
#
#
#
#
saveRDS(seurat, file = here('data/MAIT_Abau_Ftula_CT_Lungs.rds'))
#
#
#
#
#
seurat <- readRDS(here('data/MAIT_Abau_Ftula_CT_Lungs.rds'))
object_annotations <- 'full_object_filtered'
resolution <- 0.15
#
#
#
#
#
#| layout-ncol: 2
Idents(seurat) <- paste0('SCT_snn_res.', resolution)
seurat[['seurat_clusters']]<- Idents(seurat)
DimPlot(seurat, reduction = "umap", label = TRUE) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
ggsave(paste0('UMAP_clusters_R_', resolution, '_', object_annotations, '.pdf'), path = figures_path)
DimPlot(seurat, reduction = "umap", label = TRUE, split.by = 'Groups', ncol = 2) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
ggsave(paste0('UMAP_clusters_by_group_R_', resolution, '_', object_annotations, '.pdf'), path = figures_path, width = 5, height = 5)
DimPlot(seurat, reduction = "umap", label = TRUE, split.by = 'Samples', ncol = 4) + ggtitle(paste0('R ', resolution)) + theme(legend.position = "none")
ggsave(paste0('UMAP_clusters_by_sample_R_', resolution, '_', object_annotations, '.pdf'), path = figures_path)
DimPlot(seurat, reduction = "umap", label = FALSE, group.by = 'Groups') + ggtitle(paste0('R ', resolution)) + theme(legend.position = "right")
#
#
#
#
#
extract_cell_counts(seurat, seurat_clusters, figures_path, tables_path, object_annotations=object_annotations)
#
#
#
#
#
#| layout-ncol: 3
#seurat <- PrepSCTFindMarkers(seurat)
plot1 <- top_genes_per_cluster(seurat, n_genes_to_plot=5, object_annotations, tables_path = tables_path, figures_path = figures_path, results_path = results_path, run_pathway_enrichment = T) 
#
#
#
#| fig-width: 8
#| fig-height: 12
print(plot1)
ggsave(paste0(figures_path, 'DotPlot_Top3_per_cluster', object_annotations, '.pdf'), width = 8, height = 12)
#
#
#
#
#
#
#

DefaultAssay(seurat) <- 'SCT'
signatures_table <- read_excel(here('data/Hagoulet_signatures.xlsx'), skip = 1)
signatures_table <- signatures_table |>
    dplyr::filter(!(row_number() == 1)) |>
    dplyr::select(-Signature) 

vector <- signatures_table$ISG |> purrr::discard(is.na)

signatures <- map(.x = signatures_table, .f = ~ purrr::discard(.x, is.na))

# Removing genes that are not detected
blacklist_elements <- c('H2afz','Smsc2','Hist1h2ap','Flg','Flg2','Bmp2','Bmp4','Ngf','Fgf7','Fgf10','Fgf22','Ctgf','Mmp3','Mmp13','Defb1','Defb6','Epgn','Mmp2','Angpt1','Hgf','Cyr61','Inhbb','Btc','Ereg','Lep','Vegfc','Vegfd','Pdgfc','Cxcl12','Shh','Ihh','Dll3','Dll4','Wnt1','Wnt2','Wnt3a','Wnt5a','Wnt6','Wnt7a','Wnt7b','Wnt8a','Wnt8b','Wnt9a','Wnt9b','Wnt16','Chat','Thbs1','Fam46a','Kcna4','Fam129b','Tdpoz4','Ero1l','Cd244','A430078G23Rik','D1Ertd622e','Fam49a')

signatures <- map(.x = signatures, .f = ~ .x[!(.x %in% blacklist_elements)])


#
#
#
#
#

seurat <- AddModuleScore_UCell(seurat, features = signatures, name = NULL)

# Scale signatures (z-score) for comparison and visualization purposes
seurat@meta.data  <-  seurat@meta.data |>
    mutate(across(names(signatures), ~ scale(.x, center = T, scale = T), .names = '{col}_scaled'))
#
#
#
#
#
#
#| layout-ncol: 3

# Visualization
signature_violin_plot <- function (signature) {
    plot1 <- VlnPlot(seurat, features = paste0(signature, '_scaled'), group.by = 'Groups', pt.size = 0) + labs(title = signature) + NoLegend()
    print(plot1)
    ggsave(plot = plot1, filename = paste0('VlnPlot_', signature, '_by_group_', object_annotations, '.pdf'), path = figures_path, width = 5, height = 4)
}
for (signature in names(signatures)) {
    signature_violin_plot(signature) 
}


#
#
#
#
#
#| layout-ncol: 3

# Visualization
signature_violin_plot <- function (signature) {
    plot1 <- VlnPlot(seurat, features = paste0(signature, '_scaled'), group.by = 'seurat_clusters', pt.size = 0) + NoLegend() + labs(title = signature)
    print(plot1)
    ggsave(plot = plot1, filename = paste0('VlnPlot_', signature, '_by_seurat_clusters_', object_annotations, '.pdf'), path = figures_path, width = 8, height = 4)
}
for (signature in names(signatures)) {
    signature_violin_plot(signature) 
}

#
#
#
#
#
#
#| fig-width: 4
#| fig-height: 12
#| layout-ncol: 4 

signature_dot_plot <- function (signature) {
    plot1 <- DotPlot_scCustom(seurat, features = signatures[[signature]] |> rev() |> unique() , group.by = 'Groups', colors_use = sequential_palette_dotplot, flip_axes = T, scale = F, dot.scale = 6) + labs(title = signature)
    print(plot1)
    ggsave(plot = plot1, filename = paste0('DotPlot_', signature, '_by_groups_', object_annotations, '.pdf'), path = figures_path, width = 5, height = length(signatures[[signature]])/4 + 1)
}
for (signature in names(signatures)) {
    signature_dot_plot(signature) 
}

#
#
#
#
#
#
#| layout-ncol: 2

diverging_palette_2 <- hcl.colors(n = 20,'RdBu',rev = T)
signature_feature_plot <- function (signature) {
    min_value  <-  min(seurat@meta.data |> pull(paste0(signature, '_scaled')))
    max_value  <-  max(seurat@meta.data |> pull(paste0(signature, '_scaled')))
    cut_off_value <- min(abs(min_value), abs(max_value))
    plot1 <- FeaturePlot_scCustom(seurat, features = paste0(signature, '_scaled'), colors_use = diverging_palette_2, na_cutoff=NULL, min.cutoff = -1*cut_off_value, max.cutoff = cut_off_value) + labs(title = signature)
    print(plot1)
    ggsave(plot = plot1, filename = paste0('FeaturePlot Underscore', signature, '_', object_annotations, '.pdf'), path = figures_path, width = 5, height = 5)
}
for (signature in names(signatures)) {
    signature_feature_plot(signature) 
}

FeaturePlot_scCustom(seurat, features = c('Sell', 'Lef1', 'S1pr1', 'Cd44'), colors_use = sequential_palette)
DotPlot_scCustom(seurat, features = c('Sell', 'Lef1', 'S1pr1', 'Cd44') , group.by = 'seurat_clusters', colors_use = sequential_palette_dotplot, flip_axes = T, scale = T, dot.scale = 6) + labs(title = 'LN MAIT cell markers')
#
#
#
#
#
#
#
#
#
saveRDS(seurat, file = here('data/MAIT_Abau_Ftula_CT_Lungs.rds'))
#
#
#
#
#
#
#
#
#
seurat <- readRDS(here('data/MAIT_Abau_Ftula_CT_Lungs.rds'))
#
#
#
#| fig-width: 7
#| fig-height: 7
#| layout-ncol: 2
#| results: hide

gProfiler2 <- TRUE 

path <- paste0(results_path, 'DEG_pseudobulk_analysis_Lymphocytes_Abau_vs_Ftula_seurat_clusters')
dir.create(path)
Idents(seurat) <- 'seurat_clusters'

DEG_counts <- data.frame(matrix(ncol=3, nrow=0))
colnames(DEG_counts) <- c('DEG_count', 'DEG_UP_count', 'DEG_DOWN_count')
rnames <- c()

DEG_counts <- pseudobulk(seurat, comparison='Samples', group1='Ftula', group2='Abau', cluster= 'all_MAITs' , path=path, label_threshold = 10000, max_overlaps = 20,gene_lists_to_plot = NA, FC_threshold = 1,p_value_threshold = 0.05, run_gProfiler2 = gProfiler2)
rnames <- c(rnames, 'all_MAITs')

for (x in levels(seurat$seurat_clusters)) {
    # cat(paste0("Cell type: ", x, "\n"))
    seurat_small <- subset(seurat, subset = seurat_clusters == x )
    counts <- pseudobulk(seurat_small, comparison='Samples', group1='Ftula', group2='Abau', cluster= x , path=path, label_threshold = 10000, max_overlaps = 15,gene_lists_to_plot = NA, FC_threshold = 1,p_value_threshold = 0.05, run_gProfiler2 = gProfiler2)
    rm(seurat_small)
    DEG_counts <- rbind(DEG_counts, counts)
    rnames <- c(rnames, x)
}
rownames(DEG_counts) <- rnames
write.csv(DEG_counts, file = file.path(path, "DEG_counts.csv"))
DEG_counts
#
#
#
#
#
#
#
#
#
#
#
#
#
seurat <- readRDS(here('data/MAIT_Abau_Ftula_CT_Lungs.rds'))
#
#
#
#
#
#
#
#
#
#
library(Matrix)
local_path <- here('data', 'cNMF')
unlink(local_path, recursive = T)
dir.create(local_path)
local_path_input <- here('data', 'cNMF','input')
unlink(local_path_input, recursive = T)
dir.create(local_path_input)
local_path_input

filtered_dir  <- local_path_input
data_dir <- local_path

counts <- seurat@assays$RNA$counts
barcodes <- colnames(counts)
gene_names <- rownames(counts)
counts[1:5, 1:5]

# Output counts matrix
writeMM(counts, here(filtered_dir, 'matrix.mtx'))

# Output cell barcodes
write.table(as.data.frame(barcodes), here(filtered_dir, 'barcodes.tsv'),
           col.names = FALSE, row.names = FALSE, sep = "\t")


# Output feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names,"gene_name" = gene_names,type = "Gene Expression")
write.table(as.data.frame(features), sep = "\t", here(filtered_dir, 'genes.tsv'),
           col.names = FALSE, row.names = FALSE)
           
#
#
#
#
#
#
#
#
#
#
#

runname = "cNMF_run"
cmd = paste("cnmf prepare --output-dir", data_dir,
            "--name", runname,
            "-c", here(filtered_dir, 'matrix.mtx'),
            "--max-nmf-iter 2000", 
            "-k 5 6 7 8 9 10 --n-iter 20", sep=" ")
print(cmd)
system(cmd)
#
#
#
#
#
#
#
#
#

cmd = paste("cnmf factorize --output-dir", data_dir,
            "--name", runname,
            "--worker-index 0 --total-workers 4", sep=" ")
print(cmd)
system(cmd)

cmd = paste("cnmf factorize --output-dir", data_dir,
            "--name", runname,
            "--worker-index 1 --total-workers 4", sep=" ")
print(cmd)
system(cmd)

cmd = paste("cnmf factorize --output-dir", data_dir,
            "--name", runname,
            "--worker-index 2 --total-workers 4", sep=" ")
print(cmd)
system(cmd)

cmd = paste("cnmf factorize --output-dir", data_dir,
            "--name", runname,
            "--worker-index 3 --total-workers 4", sep=" ")
print(cmd)
system(cmd)

#
#
#
#
#
#
#
#
#
#
cmd = paste("cnmf combine --output-dir", data_dir,
            "--name", runname, sep=" ")
print(cmd)
system(cmd)
#
#
#
#
#
#
#
#
#
cmd = paste("cnmf k_selection_plot --output-dir", data_dir,
            "--name", runname, sep=" ")
print(cmd)
system(cmd)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cmd = paste("cnmf consensus --output-dir", data_dir,
            "--name", runname,
            '--components', 7,
            '--local-density-threshold', 0.1,
            '--show-clustering', sep=" ")
print(cmd)
system(cmd)
#
#
#
#
#
#
#
#
usage_file <- paste(data_dir[1:length(data_dir)], runname, paste(runname, "usages", "k_7.dt_0_1", 'consensus', 'txt', sep="."), sep="/")
spectra_score_file <- paste(data_dir[1:length(data_dir)], runname, paste(runname, "gene_spectra_score", "k_7.dt_0_1", 'txt', sep="."), sep="/")
spectra_tpm_file <- paste(data_dir[1:length(data_dir)], runname, paste(runname, "gene_spectra_tpm", "k_7.dt_0_1", 'txt', sep="."), sep="/")

usage <- read.table(usage_file, sep='\t', row.names=1, header=TRUE)
spectra_score <- read.table(spectra_score_file, sep='\t', row.names=1, header=TRUE)
spectra_tpm <- read.table(spectra_tpm_file, sep='\t', row.names=1, header=TRUE)
head(usage)
#
#
#
#
#
usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))
#
#
#
#
#
#


# library(SeuratDisk)
# seurat_local <- seurat
# # DefaultAssay(seurat_local) <- 'RNA'
# # seurat_local <- DietSeurat(seurat_local, assays = c("RNA", 'SCT'), dimreducs = c("pca", "umap"), layers = c('counts'))
# SaveH5Seurat(seurat_local, filename = here('data', "seurat.h5Seurat"), overwrite = T)
# SeuratDisk::Convert(source = here('data', "seurat.h5Seurat"), dest = "h5ad", overwrite = T)

#
#
#
#
#
#
#
#
#
#
#
#| fig.width: 8
#| fig.height: 8
table(seurat$cloneSize)
#UMAP clonotype frequency
slot(seurat, "meta.data")$cloneSize <- factor(slot(seurat, "meta.data")$cloneSize, 
                levels = c(
                           "Small (1e-04 < X <= 0.001)",
                            'Medium (0.001 < X <= 0.01)', 
                            'Large (0.01 < X <= 0.1)',
                            'Hyperexpanded (0.1 < X <= 1)',
                                                        NA)) 
DimPlot_scCustom(seurat, group.by = "cloneSize", pt.size = 0.5, order=T, colors_use = viridis(length(levels(seurat$cloneSize)))) 
ggsave(filename = 'UMAP_VDJ_clone_frequencies.pdf', path = figures_path, width = 8, height = 8)
DimPlot_scCustom(seurat, group.by = "cloneSize", pt.size = 0.5, order=T, colors_use = viridis(length(levels(seurat$cloneSize))), split.by = 'Groups')
ggsave(filename = 'UMAP_VDJ_clone_frequencies_by_group.pdf', path = figures_path, width = 13, height = 8)  
#
#
#
#| fig-width: 9
#| fig-height: 6
DimPlot_scCustom(seurat, group.by = "cloneSize", pt.size = 0.5, order=T, colors_use = viridis(length(levels(seurat$cloneSize))), split.by = 'Samples', num_columns = 4)& theme(text = element_text(size = 7))
ggsave(filename = 'UMAP_VDJ_clone_frequencies_by_sample.pdf', path = figures_path, width = 9, height = 6)
#
#
#
#
#| fig-width: 10
#| fig-height: 5
# #| layout-ncol: 2
# Clonotype size table for visualization
clonotype_size_table <- seurat@meta.data |>
    dplyr::select(c('Samples', 'cloneSize')) |>
    table() |>
    as_tibble() |>
    mutate(Groups = case_when(
        str_detect(Samples, 'Abau') ~ 'Abau',
        str_detect(Samples, 'Ftula') ~ 'Ftula',
        str_detect(Samples, 'CT') ~ 'CT',
    ))  |> 
    mutate(
        Groups = factor(Groups, levels = c('CT', 'Abau', 'Ftula')),
        Samples = factor(Samples, levels = levels(seurat$Samples))
        ) |>
    arrange((Samples))

clonotype_size_table <- clonotype_size_table |>
    group_by(Samples) |>
    mutate(frequency_within_sample = n/sum(n)*100) |>
    ungroup()
head(clonotype_size_table)

clonotype_size_table |> 
    ggplot(aes(x=Groups, y = n, fill=cloneSize))  +
    geom_col() +
    theme_classic() +
    labs(x = 'Group', y = 'Cell count', title = 'Clonotype size by group')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_viridis(discrete = T, na.value = 'grey90', direction = -1)
ggsave(filename = paste0('clonotype_size_by_group_', object_annotations, '.pdf'), path = figures_path, width = 6, height = 6)

clonotype_size_table |> 
    ggplot(aes(x=Samples, y = n, fill=cloneSize))  +
    geom_col() +
    theme_classic() +
    labs(x = 'Group', y = 'Cell count', title = 'Clonotype size by group')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_viridis(discrete = T, na.value = 'grey90', direction = -1)
ggsave(filename = paste0('clonotype_size_by_sample_', object_annotations, '.pdf'), path = figures_path, width = 10, height = 6)


write.csv(clonotype_size_table, file=here(tables_path, 'clonotype_size_table_results.csv'), row.names=T)
#
#
#
#
#
#
#
#| fig.width: 8
#| fig.height: 8

# Checking for full TCR sequences
seurat@meta.data <- seurat@meta.data |>
    mutate(TCR_exists = ifelse(is.na(CTaa), 'NO', 'YES'),
        TCR_full = ifelse(str_detect(CTaa, '_NA') | str_detect(CTaa, 'NA_'), 'NO', 'YES'))
#scRNAseq_test <- subset(seurat, subset = TCR_full == 'YES')
print('Number of total cells: ')
seurat@meta.data  |> nrow()
print('Number of cells with either TCR sequence: ')
seurat@meta.data |> filter(TCR_exists == 'YES') |> nrow()
print('Number of cells with full TCR sequences detected: ')
seurat@meta.data |> filter(TCR_full == 'YES') |> nrow()

# Extracting the TCR alpha chain
local_dataframe <- seurat@meta.data |> dplyr::select(starts_with('CT')) 
local_dataframe <- local_dataframe |>
                            separate_wider_delim(CTgene,delim = '_',names = c('CT_alpha', 'CT_beta'), cols_remove = FALSE)  |>
                            separate_wider_delim(CT_alpha,delim = '.',names = c('CT_V_alpha', 'CT_J_alpha', NA), cols_remove = FALSE, too_many = 'drop')  |>
                            mutate(CT_alpha_final = str_c(CT_V_alpha, '_', CT_J_alpha))
seurat$CT_alpha_final <- pull(local_dataframe, CT_alpha_final)


# Extracting the TCR beta chain
local_dataframe <- seurat@meta.data |> dplyr::select(starts_with('CT')) 
local_dataframe <- local_dataframe |>
                            separate_wider_delim(CTgene,delim = '_',names = c('CT_alpha', 'CT_beta'), cols_remove = FALSE)  |>
                            separate_wider_delim(CT_beta,delim = '.',names = c('CT_V_beta', 'CT_J_beta', NA), cols_remove = FALSE, too_many = 'drop')  |>
                            mutate(CT_beta_final = str_c(CT_V_beta, '_', CT_J_beta))
seurat$CT_beta_final <- pull(local_dataframe, CT_beta_final)


# for (column in seurat@meta.data |> dplyr::select(starts_with('CT'))) {
#     print(colnames(column))
#     print(length(column))
# }
#
#
#
#
#
#
#
#| layout-ncol: 2
#| fig.width: 8
#| fig.height: 8

Idents(seurat) <- 'CT_alpha_final'
seurat@meta.data <- mutate(seurat@meta.data, highlight=ifelse(CT_alpha_final == "TRAV11_TRAJ18", 'TRAV11_TRAJ18', 'Other'))

seurat$highlight <- factor(seurat$highlight, levels=c('Other','TRAV11_TRAJ18'))

DimPlot_scCustom(seurat, group.by = "highlight", order = T, pt.size = 1, colors_use = hcl.colors(n = 2, palette = 'ag_GrnYl')) + ggtitle('iNKTs (TRAV11_TRAJ18)')
ggsave(filename = 'UMAP_VDJ_iNKTs.pdf', width = 6, height = 5, path = figures_path)
DimPlot_scCustom(seurat, group.by = "highlight", order = T, pt.size = 1, colors_use = hcl.colors(n = 2, palette = 'ag_GrnYl'), split.by = 'Groups') & labs(color = 'iNKT TCR')
ggsave(filename = 'UMAP_VDJ_iNKTs_by_group.pdf', width = 8, height = 8, path = figures_path)
#
#
#
#
#
#| layout-ncol: 2
#| fig.width: 8
#| fig.height: 8

#MAITs
Idents(seurat) <- 'CT_alpha_final'
seurat@meta.data <- mutate(seurat@meta.data, highlight_MAIT=ifelse(CT_alpha_final == "TRAV1_TRAJ33", 'TRAV1_TRAJ33', 'Other'))

seurat$highlight_MAIT <- factor(seurat$highlight_MAIT, levels=c('Other','TRAV1_TRAJ33'))
DimPlot_scCustom(seurat, group.by = "highlight_MAIT", order = T, pt.size = 1, colors_use = hcl.colors(n = 2, palette = 'ag_GrnYl')) + ggtitle('MAIT TCR (TRAV1_TRAJ33)') + labs(color = 'MAIT TCR')
ggsave(filename = 'UMAP_VDJ_MAIT.pdf', width = 6, height = 5, path = figures_path)
DimPlot_scCustom(seurat, group.by = "highlight_MAIT", order = T, pt.size = 1, colors_use = hcl.colors(n = 2, palette = 'ag_GrnYl'), split.by = 'Groups') & labs(color = 'MAIT TCR')
ggsave(filename = 'UMAP_VDJ_MAIT_by_group.pdf', width = 8, height = 8, path = figures_path)  
#
#
#
DimPlot_scCustom(seurat, group.by = "highlight_MAIT", order = T, pt.size = 1, colors_use = hcl.colors(n = 2, palette = 'ag_GrnYl'), split.by = 'Samples', num_columns = 4) & labs(color = 'MAIT TCR') & theme(text = element_text(size = 7))
ggsave(filename = 'UMAP_VDJ_MAIT_by_sample.pdf', width = 8, height = 8, path = figures_path)
#
#
#
#
#
#
#
#| fig.width: 8
#| fig.height: 8

#QFL
Idents(seurat) <- 'CT_alpha_final'
seurat@meta.data <- mutate(seurat@meta.data, highlight_QFL=ifelse(CT_alpha_final == "TRAV9D-3_TRAJ21", 'TRAV9D-3_TRAJ21', 'Other'))

seurat$highlight_QFL <- factor(seurat$highlight_QFL, levels=c('Other','TRAV9D-3_TRAJ21'))
DimPlot_scCustom(seurat, group.by = "highlight_QFL", order = T, pt.size = 1, colors_use = hcl.colors(n = 2, palette = 'ag_GrnYl')) + ggtitle('QFL TCR (TRAV9D-3_TRAJ21)')
ggsave(filename = 'UMAP_VDJ_QFL.pdf', width = 8, height = 5, path = figures_path)
DimPlot_scCustom(seurat, group.by = "highlight_QFL", order = T, pt.size = 1, colors_use = hcl.colors(n = 2, palette = 'ag_GrnYl'), split.by = 'Groups') & labs(color = 'QFL TCR')

ggsave(filename = 'UMAP_VDJ_QFL_by_group.pdf', width = 8, height = 8, path = figures_path)
#
#
#
#
#
#
#
#| output: false
calculate_D50(seurat, cell_grouping_var = seurat_clusters, replicate_var = Samples, replicate_group_var = Groups, results_path = results_path, figures_path = figures_path) 
#
#
#
#| include: false
#| vscode: {languageId: r}
# #Extracting TCR data for clusters of interest
# Idents(scRNAseq) <- 'cell_types'
# combined2 <- scRepertoire:::.expression2List(scRNAseq, split.by ='ident')[

# diversity_measure_results <- clonalDiversity(combined2, cloneCall = 'aa', exportTable = T, n.boots = 100)
# #diversity_measure_results <- diversity_measure_results |>
#     #separate_wider_delim(Group, '_', names = c('cell_type', 'Sample'))
# write.csv(diversity_measure_results, file=paste0(path, 'diversity_measure_results.csv'), row.names=T)
#
#
#
#
#
#
#
#
#
#| layout-ncol: 2
#| fig-width: 10
#| fig-height: 10

# By sample
immunarch_table <- exportClones(seurat, format = 'immunarch', write.file = FALSE, group.by = 'Samples')
overlap <- repOverlap(immunarch_table$data, .method = 'morisita', .verbose=FALSE)
vis(overlap) + ggtitle('Morisita overlap index by samples')
ggsave(filename = paste0(figures_path, 'morisita_overlap_samples.pdf'), width = 10, height = 10)
combined2 <- scRepertoire:::.expression2List(seurat, split.by ='Samples')
morisita_table <- clonalOverlap(combined2, cloneCall = 'aa', chain = 'both',exportTable = T,method = 'morisita')
write.csv(morisita_table, file = paste0(results_path, 'morisita_table_samples.csv'))

# By Cluster
immunarch_table <- exportClones(seurat, format = 'immunarch', write.file = FALSE, group.by = 'seurat_clusters')
overlap <- repOverlap(immunarch_table$data, .method = 'morisita', .verbose=FALSE)
vis(overlap) + ggtitle('Morisita overlap index') + labs(x = NULL, y = NULL)
ggsave(filename = paste0(figures_path, 'morisita_overlap_seurat_clusters.pdf'), width = 10, height = 8)
combined2 <- scRepertoire:::.expression2List(seurat, split.by ='seurat_clusters')
morisita_table <- clonalOverlap(combined2, cloneCall = 'aa', chain = 'both',exportTable = T,method = 'morisita')
write.csv(morisita_table, file = paste0(results_path, 'morisita_table_seurat_clusters.csv'))

#
#
#
#
#
#| layout-ncol: 2
#| fig-width: 10
#| fig-height: 10
#| 
# by Sample
immunarch_table <- exportClones(seurat, format = 'immunarch', write.file = FALSE, group.by = 'Samples')
overlap <- repOverlap(immunarch_table$data, .method = 'overlap', .verbose=FALSE)
vis(overlap) + ggtitle('Overlap Coefficient')
ggsave(filename = paste0(figures_path, 'overlap_coefficient.pdf'))
combined2 <- scRepertoire:::.expression2List(seurat, split.by ='Samples')
overlap_coefficient <- clonalOverlap(combined2, cloneCall = 'aa', chain = 'both',exportTable = T,method = 'overlap')
write.csv(overlap_coefficient, file = paste0(results_path, 'overlap_coefficient.csv'))

# By Cluster
immunarch_table <- exportClones(seurat, format = 'immunarch', write.file = FALSE, group.by = 'seurat_clusters')
overlap <- repOverlap(immunarch_table$data, .method = 'overlap', .verbose=FALSE)
vis(overlap) + ggtitle('Overlap coefficient')
ggsave(filename = paste0(figures_path, 'overlap_coefficient_seurat_clusters.pdf'))
combined2 <- scRepertoire:::.expression2List(seurat, split.by ='seurat_clusters')
overlap_coefficient <- clonalOverlap(combined2, cloneCall = 'aa', chain = 'both',exportTable = T,method = 'overlap')
write.csv(overlap_coefficient, file = paste0(results_path, 'overlap_coefficient_overlap_coefficient_seurat_clusters.csv'))

#
#
#
#
#
#
#
#| results: hide
color_palette <-c(
        "CT_Mouse-1"  =  'skyblue1', 
        "CT_Mouse-2"  =  'skyblue2',
        "CT_Mouse-3"  =  'skyblue3',
        "CT_Mouse-4"  =  'skyblue4',
        "Ftula_Mouse-1"  =  'salmon1',
        "Ftula_Mouse-2"  =  'salmon2',
        "Ftula_Mouse-3"  =  'salmon3',
        "Ftula_Mouse-4"  =  'salmon4',
        "Ftula_Mouse-5"  =  'salmon',
        "Abau_Mouse-1"  =  'palegreen1',
        "Abau_Mouse-2"  =  'palegreen2',
        "Abau_Mouse-3"  =  'palegreen3',
        "Abau_Mouse-4"  =  'palegreen4')
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

grouping_variable <- 'Samples'
variables_to_color_by <- NULL

overlap_circos_and_tables(
    seurat,
    grouping_variable = 'Samples',
    results_path = results_path,
    variables_to_color_by = NULL,
    cell_types_column = 'seurat_clusters',
    write_table = TRUE,
    figures_path = figures_path,
    circle_margin = 1,
    major_ticks = 200,
    cex = 0.6,
    alpha_col = 0.25,
    samples = levels(seurat$Samples),
    groups = levels(seurat$Groups),
    color_palette = color_palette)
#
#
#
#| fig-height: 6
path <- paste0(figures_path, 'Circos_clonotypes_per_', grouping_variable, '.png')
knitr::include_graphics(path, dpi = 100)
#
#
#
#
#
#| results: hide 
# color_palette <-c(
#         "MAIT-B2m-Cre--uninfected"  =  'skyblue1', 
#         "MAIT-B2m-Cre+-uninfected"  =  'skyblue3', 
#         "MAIT-B2m-Cre--Abau"  =  'slateblue1',       
#         "MAIT-B2m-Cre+-Abau"  =  'slateblue3',       
#         "abDN-B2m-Cre--uninfected"  =  'wheat1', 
#         "abDN-B2m-Cre+-uninfected"  =  'wheat3', 
#         "abDN-B2m-Cre--Abau"  =  'tan1',       
#         "abDN-B2m-Cre+-Abau"  =  'tan3')        
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

grouping_variable <- 'seurat_clusters'
variables_to_color_by <- NULL

TCR_data <- overlap_circos_and_tables(
    seurat,
    grouping_variable = grouping_variable,
    results_path = results_path,
    variables_to_color_by = NULL,
    cell_types_column = 'seurat_clusters',
    write_table = TRUE,
    figures_path = figures_path,
    circle_margin = 1,
    major_ticks = 200,
    cex = 0.6,
    alpha_col = 0.1,
    samples = levels(seurat$seurat_clusters),
    groups = levels(seurat$seurat_clusters),
    color_palette = NULL)
#
#
#
#| fig-height: 6
path <- paste0(figures_path, 'Circos_clonotypes_per_', grouping_variable, '.png') 
knitr::include_graphics(path, dpi = 100)
#
#
#
#
#
#
#
#| fig-width: 10
#| fig-height: 20
combined2 <- scRepertoire:::.expression2List(seurat, split.by ='seurat_clusters')

vizGenes(combined2,
         x.axis = "TRBV",
         y.axis = NULL, # No specific y-axis variable, will group all samples
         plot = "barplot",
         summary.fun = "percent") + theme(text = element_text(size = 10))

ggsave(filename = paste0(figures_path, 'TRBV_gene_usage_per_seurat_clusters.pdf'), width = 12, height = 18)

combined4 <- scRepertoire:::.expression2List(seurat, split.by ='Samples')

vizGenes(combined4,
         x.axis = "TRBV",
         y.axis = NULL, # No specific y-axis variable, will group all samples
         plot = "barplot",
         summary.fun = "percent") + theme(text = element_text(size = 10
         ))

ggsave(filename = paste0(figures_path, 'TRBV_gene_usage_per_samples.pdf'), width = 12, height = 18)
#
#
#
#
#
#| fig-width: 12
#| fig-height: 16
vizGenes(combined4,
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",
         summary.fun = "percent") # Display percentages

ggsave(filename = paste0(figures_path, 'TRBV_TRBJ_gene_pairings.pdf'), width = 12, height = 16)
#
#
#
#
#
#| fig-width: 10
#| fig-height: 20


vizGenes(combined4,
         x.axis = "TRAV",
         y.axis = NULL, # No specific y-axis variable, will group all samples
         plot = "barplot",
         summary.fun = "percent") 

ggsave(filename = paste0(figures_path, 'TRAV_gene_usage_barplot.pdf'), width = 22, height = 18)
#
#
#
#
#
#| fig-width: 20
#| fig-height: 30

vizGenes(combined4,
         x.axis = "TRAV",
         y.axis = "TRAJ",
         plot = "heatmap",
         summary.fun = "percent")

ggsave(filename = paste0(figures_path, 'TRAV_TRAJ_pairings_heatmap.pdf'), width = 20, height = 30)
#
#
#
#
#
#
#
#
#
#| vscode: {languageId: r}
head(TCR_data, n = 30) |> knitr::kable()
#
#
#
#
#
#
