## This code is based on the tutorial - https://docs.scvi-tools.org/en/1.0.1/tutorials/notebooks/scvi_in_R.html
## Download scvi tools with using conda forge 
# download scanpy  
# Calling python from R - please read https://rstudio.github.io/reticulate/articles/calling_python.html


rm(list = ls())

library(reticulate)
library(anndata)
library(sceasy)
library(Seurat)

## it is very important to download scvi and all dependencies  (use - mamba and conda forge if needed) 
# first step was downloading mamba itself from https://github.com/conda-forge/miniforge#mambaforge
# install scvi tools and all other dependencies from mamba via condaforge
# install scanpy via mamba condaforge conda install -c conda-forge scanpy  from https://scanpy.readthedocs.io/en/stable/installation.html)
# usepython("path to where you have currently downloaded mamba needs to be set and then rerun library(reticulate) if it still cannot find the package)

reticulate::conda_list()

Sys.setenv(RETICULATE_PYTHON = "/path_to_/python")
path_to_python <- "path_to_/python"
use_python(path_to_python)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

merged_data <- readRDS("path_to_merged_seurat_object/2023_filt_merged_all_RA.RDS") 

merged_data <- NormalizeData(merged_data, normalization.method = "LogNormalize", scale.factor = 10000)
merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 3000)
top3000 <- head(VariableFeatures(merged_data), 3000)
merged_data <- merged_data[top3000]
saveRDS(merged_data, "path_to_folder/merged_data.rds")

adata <- convertFormat(merged_data, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata)
scvi$model$SCVI$setup_anndata(adata,  batch_key = "orig.ident")
model = scvi$model$SCVI(adata)
model$train()
model$train()
latent = model$get_latent_representation()

latent <- as.matrix(latent)
write.table(latent, 'path_to_folder_where_to_store_latents_from_scvi/latents_orig_idents_batch_key.txt')
rm(merged_data)
merged_data <- readRDS("path_to_folder/merged_seurat_object.rds") 
rownames(latent) = colnames(merged_data)
merged_data[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(merged_data))


# Find clusters, then run UMAP, and visualize
merged_data <- FindNeighbors(merged_data, dims = 1:10, reduction = "scvi")
merged_data <- FindClusters(merged_data, resolution =.3)
merged_data <- RunUMAP(merged_data, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(merged_data, reduction = "umap")

merged_data <- NormalizeData(merged_data, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(merged_RA_HRUH, "path_to_folder/seurat_orig_idents_batch_key.rds")

## Visualise the data 
scvi_origident <-  readRDS("path_to_folder/seurat_orig_idents_batch_key.rds")
DimPlot(scvi_origident, label = TRUE,  raster = TRUE ,  pt.size = 0.8,repel = TRUE, cols = colorRampPalette(wes_palette("Darjeeling1",5))(15)) + theme_void() + NoLegend()
ggsave("path_to_folder/scvi_origident.pdf", width = 4, height = 3) 

DimPlot(scvi_origident, label = FALSE, group.by = "Tissue", raster = TRUE ,  pt.size = 0.8,repel = TRUE,cols = wes_palette("Darjeeling1",2)) + theme_void() 
ggsave("path_to_folder/scvi_origident_group.by_tissue.pdf",width = 4, height = 3)

