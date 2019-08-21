###########################################################################
# running CellFindR script:
# Script for running CellFindR for a 10x dataset, from the matrix folders
# Version 1.0.0
# 5/19/2019
# written by Kevin Yu
###########################################################################
#
# INITIAL SETTINGS

# please entire your own here!
working_directory <- './data/E14.5/'
project_name <- 'E14.5_cochlea'

max_genes_per_cell_cut_off <- 7000
mito_filter <- FALSE


#INITIALIIZE DEPENDENT LIBRARIES
library(Seurat)
library(dplyr)
library(Matrix)
library(CellFindR)

##############################################################################
setwd(working_directory)
# loading 10x data and outputing as tenx/Seurat object
res <- 0.7
load_tenx('.', res = res, proj_name = project_name, cutoff= max_genes_per_cell_cut_off <- 7000, mito = mito_filter )

load(paste(project_name, '.Robj', sep = ''))
TSNEPlot(tenx)

# create directory and destination for CellFindR output:
# they will appear in time-marked output folder in the working directory
dir_add <- paste( 'output', Sys.time(),'/', sep ='')
dir.create(dir_add)
file_dest <- paste(dir_add, project_name,'_', sep = '')

# you can run find_res for the overall cluster but we found it sufficient to use the original
# resolution to continue: if you want that validation, run this next commented lines
#res <- find_res(tenx)
#print(res)
#tenx <- FindClusters(tenx,pc.use = 1:10,resolution = res ,print.output = 0,save.SNN = T)

#############################################################################
# runs CellFindR Clusters and the outputs:
run_clustering(tenx, file_dest, dir_add)

# sets the new CellFindR groups as a new meta.data -finalcluster, reading from the
# assign_test.csv that has CellFindR groups for each cell.
tenx@data.info <- tenx@data.info[order(row.names(tenx@data.info)),]

# if this doesnt work, please read.csv the assign_test file that is generated to locate cluster destination for cells
label_cells <- read.csv(paste(dir_add, "/assign_test.csv", sep = ''), header = TRUE)
label_cells <- label_cells[order(label_cells$"cellnames"),]

tenx@data.info$finalcluster <- label_cells$clusterid
tenx<-SetAllIdent(tenx, id = 'finalcluster')

# outputs TSNE
pdf('./TSNE_all_groups.pdf', width = 15, height = 15)
TSNEPlot(tenx, do.label = TRUE)
dev.off()

# generates differential gene matrix with respect to each CellFindR group
getmatrix(tenx, file_dest, 'CellFindR')
get_stats(tenx, file_dest, project_name)







