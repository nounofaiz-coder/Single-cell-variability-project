# ==============================================================================
# DATA IMPORT AND SEURAT OBJECT INITIALIZATION
# ==============================================================================

# Load required libraries for data manipulation, single-cell analysis, and Python interoperability
library(tidyverse)
library(Seurat)
library(reticulate)

# Import required Python modules using reticulate (without automatic conversion to preserve complex structures)
ad <- import("anndata", convert = FALSE)
scipy <- import("scipy.sparse", convert = FALSE)

# Load the raw transcriptomic data stored in AnnData format (.h5ad)
adata <- ad$read_h5ad("zf_atlas_10dpf_v1_release.h5ad")
print(adata)

# Extract metadata (cell annotations) from the AnnData object and convert it to an R data frame
metadata <- py_to_r(adata$obs)
dim(metadata)        # Verify the dimensions (rows and columns) of the metadata
head(metadata[, 1:4])

# Extract the expression count matrix. 
# Transposition is required because Python/AnnData stores matrices as (cells x genes), 
# whereas R/Seurat strictly requires a (genes x cells) structural format.
X_transposed <- adata$X$transpose()
counts_matrix <- py_to_r(X_transposed)

# Request Python to transform the indices into pure lists, then convert them to R vectors
gene_names <- py_to_r(adata$var_names$to_list())
cell_names <- py_to_r(adata$obs_names$to_list())

# Assign the extracted gene and cell names to the R count matrix
rownames(counts_matrix) <- gene_names
colnames(counts_matrix) <- cell_names

# Construct the Seurat object using the transposed count matrix and the extracted metadata
seurat_10dpf <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata)

# ==============================================================================
# DATA VERIFICATION AND EXPLORATION
# ==============================================================================

# Final verification of the newly created Seurat object dimensions and size
print(seurat_10dpf)

# Inspect the first few rows and columns of the matrix and metadata for quality control
counts_matrix[1:10, 1:5]
colnames(metadata)
head(rownames(metadata), 10)

# Explore the available metadata categories provided by the dataset authors
colnames(seurat_10dpf@meta.data)

# Quantify the number of cells originating from each individual biological replicate (fish)
table(seurat_10dpf$fish)

# List all annotated tissues present in the dataset and quantify the cellular population for each
table(seurat_10dpf$zebrafish_anatomy_ontology_class)

# Generate a cross-tabulation to evaluate the distribution of cell types across the different biological replicates
table(seurat_10dpf$fish, seurat_10dpf$zebrafish_anatomy_ontology_class)

# Verify the computational format of the count matrix (ensuring it is a sparse matrix structure, e.g., dgCMatrix)
class(seurat_10dpf[["RNA"]]$counts)

# Inspect raw Unique Molecular Identifier (UMI) counts and potential technical noise for a broadly expressed gene (actb1)
seurat_10dpf[["RNA"]]$counts["actb1", 1:5]

# Inspect counts for a highly specific gene (rho) to evaluate biological specificity and the high zero rate (sparsity)
seurat_10dpf[["RNA"]]$counts["rho", 1:5]
