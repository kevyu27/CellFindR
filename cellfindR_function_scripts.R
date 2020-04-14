
#test code for iteration
#02062020 CellFindR2 Kevin Yu

################################
## load data from processed RDS
tenx <- readRDS(file = "/Users/kyu/Desktop/Embryonic_New/Data/E14.5/E14.5.rds")
file_loc <- '/Users/kyu/Desktop/CellFindR_2.0/E14.5_test/'
proj_name <- 'E14.5_test'

# if you need to load from matrix files, use load_tenx
################################

setwd(file_loc)
DimPlot(tenx, reduction = "umap",  label = TRUE)
################################
# cluster first layer
#manual set res for first layer
res <- 0.6
#or find a resolution
res <- find_res(tenx, initial_res = 0.1, jump = 0.1)

tenx <- FindClusters(tenx, resolution = res) 
DimPlot(tenx, reduction = "umap",  label = TRUE)

################################
# create output folder
output_folder <- paste(file_loc, '/', Sys.time(), sep = '')
dir.create(output_folder)

#get subclusters based on initial clustering
tenx_labeled <- sub_clustering(tenx, output_folder, proj_name)

# save file
saveRDS(tenx_labeled, file = paste(output_folder,'/', proj_name, "_labeled.rds", sep = ''))

# load updated file
tenx_labeled <- readRDS(paste(output_folder,'/', proj_name, "_labeled.rds", sep = ''))

################################
# output matrice and dataQ files

get_analysis(tenx_labeled, output_folder, proj_name)
