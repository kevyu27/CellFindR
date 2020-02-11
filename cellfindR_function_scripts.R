#test code for iteration from raw file
#02072020 CellFindR2 Kevin Yu

#######################################


file_loc <- '/Users/kyu/Desktop/CellFindR_2.0/TM_Epi/'
proj_name <- 'TM_epi_test'

# if from Raw: load_tenx uses standard seurat, point to folder with
# matrix.mtx barcode genes (may need to change features to barcode)
load_tenx(file_loc)

# if RDS file is already created, the above load_tenx will create file.
tenx <- readRDS('/Users/kyu/Desktop/CellFindR_2.0/TM_Epi/TM_epi.rds')
file_loc <- '/Users/kyu/Desktop/CellFindR_2.0/TM_Epi/'
proj_name <- 'TM_epi_test'

setwd(file_loc)
DimPlot(tenx, reduction = "umap",  label = TRUE)

# cluster first layer
res <- find_res(tenx, initial_res = 0.1, jump = 0.1)
tenx <- FindClusters(tenx, resolution = res)
DimPlot(tenx, reduction = "umap",  label = TRUE)

# create output folder
output_folder <- paste(file_loc, '/', Sys.time(), sep = '')
dir.create(output_folder)

#get subclusters based on initial clustering
tenx_labeled <- sub_clustering(tenx, output_folder)

# save file
saveRDS(tenx_labeled, file = paste(output_folder,'/', proj_name, "_labeled.rds", sep = ''))

# load updated file
tenx_labeled <- readRDS(paste(output_folder,'/', proj_name, "_labeled.rds", sep = ''))

DimPlot(tenx_labeled, reduction = "umap",  label = TRUE)

tenx_labeled <-SetIdent(tenx_labeled, value = 'CellfindR')
levels(tenx_labeled) <-str_sort(levels(tenx_labeled), numeric = TRUE)

z <- get_matrix(tenx_labeled)
write.csv(z, file = paste(output_folder, '/', 'matrix_cellfindr.csv', sep = ''))

a <- get_stats(tenx)
write.csv(a, file = paste(output_folder, '/', 'all_stats.csv', sep = ''))

tenx_labeled <-SetIdent(tenx_labeled, value = 'seurat_clusters')
y <- get_matrix(tenx_labeled)
write.csv(y, file = paste(output_folder, '/', 'matrix_big_groups.csv', sep = ''))

#markers <-FindAllMarkers(tenx_labeled,only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25)
