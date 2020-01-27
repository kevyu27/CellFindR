library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(gplot2)
library(RColorBrewer)
dev.off()
# load data

# loading tenx data
# input: file_loc = location of file to load
# res = resolution to run once
# proj_name = name of project
# cutoff = high end cut off of num genes
# output saves to file_loc a rdata file.
load_tenx <- function(file_loc, res = 0.1, proj_name = 'tenx_data', cutoff= 10000, mito = FALSE){
  tenx.data <- Read10X(data.dir = file_loc)
  tenx <- CreateSeuratObject(counts = tenx.data, project = proj_name, min.cells = 3, min.features = 200)
  
  #processing
  tenx[["percent.mt"]] <- PercentageFeatureSet(tenx, pattern = "^mt-")
  tenx[["percent.MT"]] <- PercentageFeatureSet(tenx, pattern = "^MT-")
  
  VlnPlot(tenx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  ###############################
  tenx <- NormalizeData(tenx, normalization.method = "LogNormalize", scale.factor = 10000)
  tenx <- FindVariableFeatures(tenx, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(tenx), 10)
  plot1 <- VariableFeaturePlot(tenx)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  CombinePlots(plots = list(plot1, plot2))
  
  all.genes <- rownames(tenx)
  tenx <- ScaleData(tenx, features = all.genes)
  
  #####################################
  tenx <- RunPCA(tenx, features = VariableFeatures(object = tenx))
  
  tenx <- FindNeighbors(tenx, dims = 1:10)
  tenx <- FindClusters(tenx, resolution = 0.1)
  tenx <- RunUMAP(tenx, dims = 1:10)
  DimPlot(tenx, reduction = "umap")
  
  saveRDS(tenx, file = paste(file_loc,'/', proj_name, ".rds", sep = ''))
}

##############
# asking if the grouping is a cluster
# input: tenx = tenx object
# thresh_genes = threshold of genes at thresh_val
# thresh_val = value of the threshold in log space
# pval = cut off of pval for the signficance
is_cluster <- function(tenx, thresh_genes = 10, thresh_val = log(2), pval = 1e-4){
  val = 0 # groups that does not satisfy threshold genes
  counter = 0 # groups that satisfy threshold genes 
  # loop through the identitiy
  matrix_output <- data.frame(row.names = row.names(tenx))
  
  for (j in sort(unique(tenx@active.ident))){
    markers <- FindMarkers(tenx, ident.1 = j, min.pct = 0.25)
    markers <- markers[markers$p_val_adj < pval,]
    #find if the 10th biggest is less than log2, sum 
    print(sort(markers$avg_logFC, decreasing = TRUE)[thresh_genes])
    # if less than 10 significant genes
    if (sum(tenx@active.ident == j) < 5){
      return(FALSE)
    }
    if (length((markers$avg_logFC)) < 10){
      val <- val + 1
    } else if (sort(markers$avg_logFC, decreasing = TRUE)[thresh_genes] < thresh_val){
      #print(val)
      val <- val + 1
    } else{
      counter = counter + 1
    }
  }
  if (val > 1){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

# finds resolution that satisfy
# input: tenx object
# initial resolution of starting clustering
# how much to increment up 
# threshold of genes
# value of the threshold 
find_res <- function(tenx, initial_res = 0.1, jump = 0.1, thresh_genes = 10, thresh_val = log(2)) {
  RES_POST <- initial_res # keeping
  RES_IT <- initial_res # iterative
  
  while(TRUE){
    print(paste('Trying',RES_IT, sep = ' '))
    tenx <- FindNeighbors(tenx, dims = 1:10)
    tenx <- FindClusters(tenx, resolution = RES_IT)

    # also check if theres only 1 cluster/ then can go up higher es
    # Find number of clusters
    length_group <- length(unique(tenx@active.ident))
    # if only one group then need to look deeper
    if (length_group == 1){
      print(paste('noclusterfound',RES_IT, sep = ' '))
      # still not groups at 0.7 res stop and just step as 1
      if (RES_IT == 0.7){
        print(paste('noclusterfound',RES_IT, sep = ' '))
        break
      }
    } else{
      testing <- is_cluster(tenx)
      if (testing == FALSE){ # if not real group
        print(paste('broke', RES_IT, sep = ' '))
        RES_IT <- RES_IT - jump
        RES_POST <- RES_IT
        print(RES_POST)
        break
      } else{ # valid groups
        RES_POST <- RES_IT
        print(paste('ok',RES_IT, sep = ' '))
      }
    }
    RES_IT <- RES_IT + jump
  }
  return(RES_POST)
}

# getsubclustering
# input: tenx object
# location of output folder
# project_name
sub_clustering <- function(tenx, output_folder = '.', proj_name = 'proj_name', thresh_genes = 10, thresh_val = log(2)){
  # set cell directory
  print('running through clusters to find subclusters')
  res_keep <- data.frame('cluster'= NA,'res'= NA, 'num_clusters' =NA)
  celltocluster <- data.frame(row.names = colnames(tenx))
  celltocluster$cellnames <- colnames(tenx)
  res_counter <- 1
  
  for (j in sort(unique(tenx@active.ident))){
    sub_tenx <- subset(tenx, idents = toString(j))
    print(paste('clustering ', j, sep = ''))

    sub_tenx <- FindVariableFeatures(sub_tenx, selection.method = "vst", nfeatures = 2000)
    sub_tenx <- FindNeighbors(sub_tenx, dims = 1:10)
    sub_tenx <- RunUMAP(sub_tenx, dims = 1:10, n.neighbors = 10)
    
    # get resolution
    set_res <- find_res(sub_tenx)
    if (set_res > 0){
      #### label everything
      sub_tenx <-FindClusters(sub_tenx,pc.use = 1:10, resolution = set_res)
    }
    
    # only one group
    num_groups <- length(unique(sub_tenx@active.ident))
    if (num_groups == 1){
      print(paste(j, ' has no subclusters', sep = ''))
      framer <- data.frame(row.names = colnames(sub_tenx))
      framer$cellnames <- row.names(framer)
      framer$clusterid <- j
      celltocluster <- merge(celltocluster, framer, all = TRUE)
    } else{
      print(paste(j, ' has ', num_groups, ' subgroups', sep = ''))
      
      # plot data
      gen_matrix_plot(sub_tenx, output_folder, j)
      
      # layer 2: X.X
      for (k in sort(unique(sub_tenx@active.ident))){
        print(paste('clustering ', j,'.',k, sep = ''))
        sub2_tenx <- subset(sub_tenx, idents = toString(k))
        
        sub2_tenx <- FindVariableFeatures(sub2_tenx, selection.method = "vst", nfeatures = 2000)
        sub2_tenx <- FindNeighbors(sub2_tenx, dims = 1:10)
        sub2_tenx <- RunUMAP(sub2_tenx, dims = 1:10, n.neighbors = 10)
        
        # get resolution
        set_res <- find_res(sub2_tenx)
        if (set_res > 0){
          #### label everything
          sub2_tenx <-FindClusters(sub2_tenx,pc.use = 1:10, resolution = set_res)
        }
        
        num_groups <- length(unique(sub2_tenx@active.ident))
        if (num_groups == 1){
          print(paste(j,'.', k, ' has no subclusters', sep = ''))
          framer <- data.frame(row.names = colnames(sub2_tenx))
          framer$cellnames <- row.names(framer)
          framer$clusterid <- paste(j,k, sep = '.')
          celltocluster <- merge(celltocluster, framer, all = TRUE)
        } else{
          print(paste(j,'.', k, ' has ', num_groups, ' subgroups', sep = ''))
          
          # plot data
          gen_matrix_plot(sub2_tenx, output_folder, paste(j, k, sep = '.'))
          
          # layer 3: X.X.X
          for (l in sort(unique(sub2_tenx@active.ident))){
            sub3_tenx <- subset(sub2_tenx, idents = toString(l))
            framer <- data.frame(row.names = colnames(sub3_tenx))
            framer$cellnames <- row.names(framer)
            framer$clusterid <- paste(j,k,l, sep = '.')
            celltocluster <- merge(celltocluster, framer, all = TRUE)
          }
        }
      }
    }
  }
  celltocluster <- na.omit(celltocluster)
  write.csv(celltocluster, paste(output_folder, "/cell_labels.csv", sep = ''), row.names = TRUE)
  celltocluster <- celltocluster[order(celltocluster$"cellnames"),]
  tenx@meta.data <- tenx@meta.data[order(rownames(tenx@meta.data)),]
  tenx@meta.data$CellfindR <- celltocluster$clusterid
  tenx <-SetIdent(tenx,cells = rownames(tenx@meta.data), value = tenx@meta.data$CellfindR)
  levels(tenx) <-str_sort(levels(tenx), numeric = TRUE)
  
  ggsave(paste(output_folder, '/', 'UMAP_CellfindR.pdf', sep = ''),    
         DimPlot(tenx, label = TRUE),width = 8, height = 8)
  return(tenx)
}

# generate matrix and plots
gen_matrix_plot <-function(tenx, output_folder = '.', proj_name = 'proj_name'){
  file_create <-paste(output_folder,'/', proj_name,sep = '')
  print(file_create)
  ### Output data destinations
  # create folder
  ggsave(paste(file_create, '_umap.pdf', sep = ''), DimPlot(tenx, label = TRUE))
  markers <-FindAllMarkers(tenx,only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25)
  markers_filtered <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) 
  genes <- unique(markers_filtered$gene)
  write.csv(markers,paste(file_create, '_matrix.csv', sep = ''), row.names = TRUE)
  
  #create subdirectories
  file_create2 <- paste(output_folder, '/',proj_name, sep ='')
  dir.create(file_create2)
  dir.create(paste(file_create2, 'Cluster', sep = '/'))
  dir.create(paste(file_create2, 'Violin', sep = '/'))
  for (i in genes) {
    # cluster maps
    ggsave(paste(paste(file_create, 'Cluster', '', sep = '/'),i,'.pdf',sep = ''),
           FeaturePlot(tenx, features = c(i),pt.size = 2), width = 8, height = 8)
    #violin plot
    ggsave(paste(paste(file_create, 'Violin', '', sep = '/'), i, '.pdf', sep = ''),    
           VlnPlot(tenx, c(i)),width = 6, height = 4)
  }
}

# generate matrix
get_matrix <- function(tenx){
  print("getting matrix")
  avg_expression <- AverageExpression(tenx)
  matrix_all <- data.frame(row.names = rownames(avg_expression$RNA))
  for (i in levels(tenx@active.ident)){
    print(i)
    markers <- FindMarkers(tenx, ident.1 = i, logfc.threshold = 0.1)
    avg_val <- avg_expression$RNA[i]
    avg_diff <- markers[rownames(avg_expression$RNA),]$avg_logFC
    avg_diff[is.na(avg_diff)] <-0
    p_val <- markers[rownames(avg_expression$RNA),]$p_val_adj
    p_val[is.na(p_val)] <-1
    
    matrix_all <- cbind(matrix_all, avg_val)
    matrix_all <- cbind(matrix_all, avg_diff)
    matrix_all <- cbind(matrix_all, p_val)
  }
  name_col <- c()
  for (k in levels(tenx@active.ident)){
    print(k)
    name_col <- c(name_col,(c(paste(k,'Mean', sep = '_'),paste(k,'Avg_diff', sep = '_') , paste(k,'Pval', sep = '_'))))
  }
  colnames(matrix_all) <- name_col
  matrix
  return(matrix_all)
}

# get stats
get_stats <- function(tenx, num_genes = 50){
  aoe <- c("Group", "cell_number", "avg_read", "avg_umi")
  for (i in 1:num_genes){
    aoe <- c(aoe, paste('top_',i, sep = ""))
  }
  df <- data.frame(aoe)
  #initialize matrix
  for (groups in levels(tenx@active.ident)){
    subgroup <-subset(tenx, idents = groups)
    # group name
    aod <- c(groups)
    # cell number
    aod <- c(aod, length(subgroup@meta.data$nCount_RNA))
    # avg_read
    aod <- c(aod, mean(subgroup@meta.data$nCount_RNA))
    # avg_umi
    aod <- c(aod, mean(subgroup@meta.data$nFeature_RNA))
    # top 10 diff genes
    markers <- FindMarkers(tenx, groups)
    top_markers <- row.names(markers)[1:num_genes]
    for (topm in top_markers){
      aod <- c(aod, topm)
    }
    df[groups] <-aod
  }
  return(df)
}

# getting plots
get_plots<- function(tenx, output_folder = '.'){
  dir_creater <- paste(output_folder, '/plots', sep = '')
  dir.create(dir_creater)
  for (groups in levels(tenx@active.ident)){
    markers <-FindMarkers(tenx, groups)
    maxer <- min(50, length(rownames(markers)))
    for (gene in rownames(markers)[1:maxer]){
      ggsave(paste(dir_creater,'/', gene, '_cluster.pdf', sep = ''),    
             FeaturePlot(tenx, features = gene), width = 6, height = 6)
      ggsave(paste(dir_creater,'/', gene, '_violin.pdf', sep = ''),    
             VlnPlot(tenx, features = gene, slot = "counts", log = TRUE), width = length(levels(tenx@active.ident)), height = 5)
    }
  }
}

# output metrics
metrics_output <- function(tenx, output_folder = '.', species = 'mouse'){
  tenx@active.assay
  # mito percent
  if (species == 'mouse'){
    ggsave(paste(output_folder, '/', 'percent_mito.pdf', sep = ''),    
           VlnPlot(tenx, features = 'percent.mt'), width =  length(levels(tenx@active.ident)), height = 5)
  }
  if (species == 'human'){
    ggsave(paste(output_folder, '/', 'percent_mito.pdf', sep = ''),    
           VlnPlot(tenx, features = 'percent.MT'), width =  length(levels(tenx@active.ident)), height = 5)
  }
  # umi
  print('test')
  ggsave(paste(output_folder, '/', 'uMI.pdf', sep = ''),    
           VlnPlot(tenx, features = 'nFeature_RNA'), width =length(levels(tenx@active.ident)), height = 5)
}


# get the top 100 genes from the values
get_top100 <- function(a, s){
  sorted <- a[order(a[s], decreasing = TRUE), ]
  return(row.names(sorted)[1:100])
}

intersection_top100 <- function(new_expression_table, old_expression_table){
  mat_cor <- 0
  names <- (colnames(new_expression_table))
  mat_cor <-data.frame(row.names = names)
  
  for (i in colnames(old_expression_table)){
    print(paste('Subset', i))
    #sort based on column one and extract list
    tester <- get_top100(old_expression_table, i)
    list_of_i <- c()
    for (j in colnames(new_expression_table)){
      # new tables for correlation
      o <-get_top100(new_expression_table, j)
      list_of_i <- c(list_of_i, length(intersect(tester,o)))
      print(length(intersect(tester,o)))
    }
    df <- data.frame(j=list_of_i)
    mat_cor <- cbind(mat_cor, df)
  }
  colnames(mat_cor) <- colnames(old_expression_table)
  return(mat_cor)
}

##############################
#mouse human comparison
# get intersection of mouse to human 
# gets the index of the reference in mouse or in human, 

ortholog_index <- function(list_genes, id = "Mouse"){
  table_ref <- read.csv('/Users/kyu/Desktop/Project_Cochlea/look_up_table/look_up_table.csv')
  
  index_list <- c()
  if (id == "Mouse"){
    for (gene in list_genes){
      index_list <- c(index_list, which(table_ref$Mouse == gene)[1])
    }
  }
  if (id =="Human"){
    for (gene in list_genes){
      index_list <- c(index_list, which(table_ref$Human == gene)[1])
    }
  }
  return(index_list)
}

# creates intersection based on the 
intersection_h_m <- function(human_matrix, mouse_matrix){
  # set up ending coorelation matrix
  mat_cor <- 0
  names <- (colnames(mouse_matrix))
  mat_cor <-data.frame(row.names = names)
  
  # iterate through each subgroup in the human matrix
  for (i in colnames(human_matrix)){
    print(paste('Subset', i))
    
    #set up all the human genes and gets the index with respect to the hash table
    human_genes<- get_top100(human_matrix, i)
    human_index <- ortholog_index(human_genes, id= 'Human')
    list_of_i <- c()
    
    for (j in colnames(mouse_matrix)){
      #convert to human
      mouse_genes <- get_top100(mouse_matrix, j)
      mouse_index <- ortholog_index(mouse_genes)
      
      # new tables for correlation
      # could do it by homology index (returns might be easier)
      
      #get intersect of the values:
      intersect_hm <- intersect(human_index, mouse_index)
      val_hm <- length(intersect_hm)
      # if has a NA remove 1
      if (is.element(NA, intersect_hm)){
        val_hm <- val_hm -1
      }
      list_of_i <- c(list_of_i, val_hm)
      print(paste('Intersect', val_hm))
    }
    df <- data.frame(j=list_of_i)
    mat_cor <- cbind(mat_cor, df)
  }
  colnames(mat_cor) <- colnames(human_matrix)
  return(mat_cor)
}


# get deafness gene plots need to reupdate 
get_deafness_plot <- function(tenx, id = 'Mouse', name = 'deafness'){
  if (id == "Mouse"){
    df_deafness <- read.csv("/Users/kyu/Desktop/Project_Cochlea/look_up_table/Deafness_gene_list.csv", header = FALSE)
  }
  else{
    df_deafness <- read.csv("/Users/kyu/Desktop/Project_Cochlea/look_up_table/Deafness_gene_list_human.csv", header = FALSE)
  }
  
  Genes <- 0
  # loop through all the different subclusters
  # initialize data.frame 
  names <- rownames(tenx@scale.data)
  Genes <-data.frame(row.names = names)
  
  for (j in levels(tenx@ident)){
    # subset each one and then make each column that way. 
    test <- SubsetData(tenx, ident.use = j)
    mat_simp <- test@scale.data # scaled data from the values
    mat_simp <-as.matrix(mat_simp)
    row_means <- rowMeans(mat_simp, na.rm=TRUE)
    df <- data.frame( j=row_means)
    colnames(df) <-j
    Genes <- cbind(Genes, df)
    print(j)
  }
  
  # remove genes column
  # set upper limit as 1
  Genes[Genes <0 ] <- 0
  Genes[Genes >1] <- 1
  # write to heatmap (find the file for the deafness genes)
  list_deafness <- c()
  #get only the unique ones
  genenames <- row.names(tenx@scale.data)
  for (gene in df_deafness[,1]){
    if (any(genenames == gene)){
      list_deafness <- c(list_deafness, gene)
    }
  }
  list_deafness <- unique(list_deafness)
  
  ###
  # get rid of it 
  
  heatmap_test <- Genes[list_deafness,]
  heatmap_test <- heatmap_test[1:length(heatmap_test)]
  pdf(paste('./', name,'.pdf', sep = ''), width = 12, height = 20)
  pheatmap(heatmap_test, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
  dev.off()
  return(heatmap_test)
}

