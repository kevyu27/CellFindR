#initialize libraries
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
dev.off()

# please run functions in the CellfindR script.

#load RDS file from created Seurat object


file_loc <- 'location to place data'
proj_name <- 'name of project'
tenx <- readRDS('locate dataset')

# if you want to create from the scratch dataset run
load_tenx(file_loc)


# view UMAP
DimPlot(tenx, reduction = "umap",  label = TRUE)

######################################
# cluster layer 1:

# finds the resolution (this is seurat clustering)
res <- find_res(tenx, initial_res = 0.1, jump = 0.1)

# sets resolution of the subclustering as "res"
# you can modify res if you don't want to run the above res to set initial clustering resolution
# as the above find_res does take a lot of time.

# sets resolution
tenx <- FindClusters(tenx, resolution = res)
DimPlot(tenx, reduction = "umap",  label = TRUE)

###################
# running CellFindR
# create output folder for results
output_folder <- paste(file_loc, '/', Sys.time(), sep = '')
dir.create(output_folder)

#get subclusters based on initial clustering
tenx_labeled <- sub_clustering(tenx, output_folder)
saveRDS(tenx_labeled, file = paste(file_loc,'/', proj_name, ".rds", sep = ''))

# Dimplot
DimPlot(tenx_labeled, reduction = "umap",  label = TRUE)

################################################################################################3
# create matrices

z <- get_matrix(tenx_labeled)
write.csv(z, file = paste(output_folder, '/', 'matrix_cellfindr.csv', sep = ''))

a <- get_stats(tenx_labeled)
write.csv(a, file = paste(file_loc, '/', 'all_stats.csv', sep = ''))

tenx_labeled <-SetIdent(tenx_labeled, value = 'seurat_clusters')
z <- get_matrix(tenx_labeled)
write.csv(z, file = paste(file_loc, '/', 'matrix_big_groups.csv', sep = ''))



