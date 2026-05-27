library(tidyverse)
library(Seurat)
library(reticulate)

#Loading the file via Python
ad <- import("anndata", convert = FALSE)
scipy <- import("scipy.sparse", convert = FALSE)
adata <- ad$read_h5ad("zf_atlas_10dpf_v1_release.h5ad")
print(adata)

#extraction of metadata to R
metadata <- py_to_r(adata$obs)
dim(metadata)        #name rows and columns 
head(metadata[, 1:4])

#Matrix Extraction (Counts)
X_transposed <- adata$X$transpose()
counts_matrix <- py_to_r(X_transposed)
# request to Python to transform the index into a pure list (.to_list)
gene_names <- py_to_r(adata$var_names$to_list())
cell_names <- py_to_r(adata$obs_names$to_list())

#We assign the proper names to the matrix
rownames(counts_matrix) <- gene_names
colnames(counts_matrix) <- cell_names

#The final assembly in Seurat
seurat_10dpf <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata)

#Final Verification Line
print(seurat_10dpf)

counts_matrix[1:10, 1:5]
colnames(metadata)
head(rownames(metadata), 10)

#Look for the details names
colnames(seurat_10dpf@meta.data)
#How many sample + typecell for each one 
table(seurat_10dpf$fish)

#list all the tissues present and count the number of cells in each one
table(seurat_10dpf$zebrafish_anatomy_ontology_class)

#Crosses the fish (rows) with the cell types (columns)
table(seurat_10dpf$fish, seurat_10dpf$zebrafish_anatomy_ontology_class)

#Check the computer format of the matrix 
class(seurat_10dpf[["RNA"]]$counts)

#Presence of raw accounts (UMI accounts) and technical noise.
seurat_10dpf[["RNA"]]$counts["actb1", 1:5]
#Biological specificity and the zero rate (sparsity)
seurat_10dpf[["RNA"]]$counts["rho", 1:5]
