# CellFindR Version 1.0.0

library(Seurat)
library(RColorBrewer)
library(Matrix)
library(gplots)
library(pheatmap)


# loads 10x data from 10x files
load_tenx <- function(file_loc, res = 1, proj_name = 'tenx_data', cutoff= 7000, mito = FALSE){
  #load file
  tenx.data <- Read10X(file_loc)

  # initialize object
  tenx <- new("seurat", raw.data = tenx.data)

  # Setup
  tenx <-
    Setup(tenx,min.cells = 3,min.genes = 200,do.logNormalize = T,
          total.expr = 1e4,project = proj_name)

  #check upper bound
  pdf(paste(file_loc,"/",proj_name, 'hist_nGenes.pdf',sep = ''))
  hist(tenx@data.info$nGene)
  dev.off()
  #filter out
  tenx <-
    SubsetData(tenx, subset.name = "nGene", accept.high = cutoff)

  # mitochondria check, all the human/mouse ones
  mito.genes <- c(grep("^mt-", rownames(tenx@data), value = T),  grep("^MT-", rownames(tenx@data), value = T))
  percent.mito <- colSums(as.matrix(expm1(tenx@data[mito.genes, ])))/colSums(as.matrix(expm1(tenx@data)))
  # add metadata
  tenx <- AddMetaData(tenx, percent.mito, "percent.mito")
  percent_keep <- sum(tenx@data.info$percent.mito <0.05)/length(tenx@data.info$nGene)
  # plot it out
  pdf(paste(file_loc,"/",proj_name, 'hist_percent_mito.pdf',sep = ''))
  hist(tenx@data.info$percent.mito)
  dev.off()
  # no dots, since confusing
  p <- ggplot(tenx@data.info,
              aes(x=proj_name, y=percent.mito, fill = proj_name)) +
    geom_violin()+ geom_hline(yintercept = 0.05, linetype = "dashed", color = 'black', size =2) +
    annotate("text", x=1, y=-0.03, label= paste('Percent Cells: ', signif(percent_keep, digits = 4)))
  pdf(paste(file_loc,"/",proj_name, 'vln_percent_mito.pdf',sep = ''), width = 4, height = 12)
  print(p)
  dev.off()

  # sort out mito if mito = TRUE
  if(mito == TRUE){
    tenx <- SubsetData(tenx, subset.name = "percent.mito", accept.high = 0.05)
  }

  # regress out
  tenx <- RegressOut(tenx, latent.vars = c("nUMI"))

  # plot gene variance
  tenx <-
    MeanVarPlot(tenx ,fxn.x = expMean,fxn.y = logVarDivMean,x.low.cutoff = 0.0125,
                x.high.cutoff = 3,y.cutoff = 0.5,do.contour = F)
  # pca
  tenx <-
    PCA(tenx,pc.genes = tenx@var.genes,do.print = TRUE,pcs.print = 5,genes.print = 5)
  tenx <- ProjectPCA(tenx)

  # find clusters
  tenx <-
    FindClusters(tenx,pc.use = 1:10,resolution = res,print.output = 0,save.SNN = T)

  # run TSNE
  tenx <- RunTSNE(tenx, dims.use = 1:10, do.fast = T, perplexity = 20)
  #save(tenx, file = paste(paste(file_loc, proj_name, sep = '/'), '.Robj', sep = ''))

  TSNEPlot(tenx, do.label = TRUE)
}

#plots 2 different violin plots for the mitochondria distribution
plot_mito <- function(cochobj, file_loc){
  mito.genes <- c(grep("^mt-", rownames(cochobj@data), value = T),  grep("^MT-", rownames(cochobj@data), value = T))
  percent.mito <- colSums(as.matrix(expm1(cochobj@data[mito.genes, ])))/colSums(as.matrix(expm1(cochobj@data)))

  # add metadata
  cochobj <- AddMetaData(cochobj, percent.mito, "percent.mito")
  percent_keep <- sum(cochobj@data.info$percent.mito <0.05)/length(cochobj@data.info$nGene)

  cochobj@data.info$finalcluster <- as.character(cochobj@data.info$finalcluster)
  pdf(paste(file_loc,"/",cochobj@project.name, '_vln_percent_mito.pdf',sep = ''), width = 40, height = 5)
  VlnPlot(cochobj, 'percent.mito')
  dev.off()
  #no dots, since confusing
  p <- ggplot(cochobj@data.info,  aes(x=finalcluster, y=percent.mito, fill = cochobj@project.name)) +
    geom_violin()+ geom_hline(yintercept = 0.05, linetype = "dashed", color = 'black', size =2) + theme(legend.position="none")
  pdf(paste(file_loc,"/",cochobj@project.name, '_vln_percent_mito2.pdf',sep = ''), width = 40, height = 5)
  print(p)
  dev.off()
}

#plots the cluster maps for given gene
just_plot <- function(tenx, gene, size = 1 ){
  FeaturePlot(tenx,
              gene,
              cols.use = c("grey", "blue"),
              pt.size = size)
}

#generates cluster and vln plots of given genes in the 10x object
plot_gen <- function(tenx, list_genes, ex_loc = '.'){
  # load file
  library(Seurat)
  dir_add <- paste(ex_loc, '/output', Sys.time(), sep ='')
  dir.create(dir_add)
  for (i in list_genes)
  {

    print(i)

    # cluster maps
    pdf(paste(dir_add,'/cluster_', i, '.pdf',sep = ''))
    FeaturePlot(tenx, c(i),cols.use = c("grey","blue"), pt.size = 1)
    dev.off()

    #violin plot
    pdf(paste(dir_add,'/violin_', i, '.pdf',sep = ''))
    VlnPlot(tenx, c(i))
    dev.off()

  }

}

# finds resolution that satisfy
find_res <- function(tenx, exact = FALSE) {
  RES_POST <- 0.1 # keeping
  RES_IT <- 0.1 # iterative

  gap <- 0.2
  if (exact){
    g
  }
  for (i in 1:15) {
    tenx <-
      FindClusters(
        tenx,
        pc.use = 1:10,
        resolution = RES_IT,
        print.output = 0,
        save.SNN = T
      )
    # also check if theres only 1 cluster/ then can go up higher es
    # Find number of clusters
    a <- tenx@ident
    # get iteration
    maxer <- length(levels(tenx@ident))
    # check cluster here
    testing <- is_cluster(tenx)

    # check if doesnt work and has more than 1 cluster
    if (testing == FALSE & maxer > 1) {
      print(paste('broke', RES_IT, sep = ' '))
      RES_IT <- RES_IT - gap
      RES_POST <- RES_IT
      print(RES_POST)
      break
    }
    else{
      # valid cluster
      if (testing == TRUE & maxer > 1){
        RES_POST <- RES_IT

        print(paste('ok',RES_IT, sep = ' '))
      }
      # invalid cluster stop at 0.7
      if (maxer == 1 & RES_IT == 0.7){
        print(paste('noclusterfound',RES_IT, sep = ' '))
        break
      }
    }
    RES_IT <- RES_IT + gap
  }
  return(RES_POST)
}

# asks if the given cluster definition satisfies our conditions
is_cluster <- function(tenx) {
  thresh <- 10 # genes
  mat <- gen_matrix_small(tenx)
  b_list <- c()
  p_list <- c()
  #get the 2x and 10x difference
  for (j in 1:(length(mat[1, ]))) {
    #check odd is the avg_diff
    if (j %% 2 != 0){
      b_list <- c(b_list, sum(mat[, j] > log(2)))
    }
    # even number need to check pvalue
    else{
      p_list <- c(p_list, sum(mat[, j] < 10e-5))
    }
  }
  # test sum
  # both conditions must be kept
  return(sum(b_list > thresh) >= (length(b_list) - 1) &&  sum(p_list > thresh) >= (length(p_list) - 1) )
}

# generates the shorthand matrix for the clustering alg
gen_matrix_small <- function(tenx) {
  names <- rownames(tenx@scale.data)
  # initialize data.frame
  Genes <- data.frame(row.names = names)

  looper <- sort(unique(tenx@ident))

  # loop through all the different subclusters
  for (j in looper)
  {
    print(j)
    # subset each one and then make each column that way.
    test <- SubsetData(tenx, ident.use = toString(j))
    mat_simp <- test@scale.data
    #add cluster markers if they are listed in high differential
    markers <- FindMarkers(tenx, ident.1 = j, min.pct = 0.25)
    mat_markers <- markers[names, ]
    if (length(markers$p_val) >10){
      print('2')
      holder <- data.frame(row.names = names)
      holder$avgdiff <- 0
      holder$pval <- 0
      # gotta actualy loop and then append
      for (g in row.names(markers)) {
        holder[g, 'avgdiff'] <- markers[2][g, ]
        holder[g, 'pval'] <- markers[1][g, ]
      }
      avg_diff <- holder$avgdiff
      pval <- holder$pval
      pval_adj <-
        p.adjust(holder$pval,
                 method = "hochberg",
                 n = length(holder$pval))
    }else {
      avg_diff <- 0
      pval_adj <- 1
    }
    Genes <- cbind(Genes, avg_diff, row.names = NULL)
    Genes <- cbind(Genes, pval_adj, row.names = NULL)
  }
  # rename columns

  Genes[Genes == 0] <- 0.5 # makes it 0.5
  return(Genes)
}

# gets the matrix for the dataset
getmatrix <- function(tenx, file_dest,d){
  names <- rownames(tenx@scale.data)
  # initialize data.frame
  Genes <- data.frame(row.names = names)

  looper <- sort(unique(tenx@ident))
  for (j in looper){
    # subset each one and then make each column that way.
    test <- SubsetData(tenx, ident.use = toString(j))
    mat_simp <- test@scale.data
    rowmeans <- rowMeans(mat_simp)
    rowsd <- apply(mat_simp, 1, sd)
    Genes <- cbind(Genes, data.frame(Mean = rowmeans), row.names = NULL)
    Genes <- cbind(Genes, data.frame(SD = rowsd), row.names = NULL)

    #add cluster markers if they are listed in high differential
    markers <- FindMarkers(tenx, ident.1 = j, min.pct = 0.25)
    mat_markers <- markers[names, ]

    holder <- data.frame(row.names = names)
    holder$avgdiff <- 0
    holder$pval <- 0
    # gotta actualy loop and then append
    for (g in row.names(markers)) {
      holder[g, 'avgdiff'] <- markers[2][g, ]
      holder[g, 'pval'] <- markers[1][g, ]
    }
    avg_diff <- holder$avgdiff
    pval <- holder$pval
    Genes <- cbind(Genes, avg_diff)
    Genes <- cbind(Genes, pval)
  }

  # rename columns
  name_col <- c()
  for (k in looper){
    name_col <- c(name_col,(c(paste(k,'Mean', sep = '_'),paste(k,'SD', sep = '_') ,paste(k,'Avg_diff', sep = '_') , paste(k,'Pval', sep = '_'))))
  }
  colnames(Genes) <- name_col
  # write to file
  row.names(Genes) <- names
  Genes[Genes == 0] <- ""

  write.csv(Genes,
            paste(file_dest, d, '_matrix.csv', sep = ''),
            row.names = TRUE)
  print(paste(file_dest, d, '_matrix.csv', sep = ''))
}

# plots top 100 genes of subcluster
plottingj <-function(sub_obj, file_dest, d){
  file_create <-paste(paste(file_dest, d,sep = ''))
  markers <-
    FindAllMarkers(
      sub_obj,
      only.pos = TRUE,
      min.pct = 0.25,
      thresh.use = 0.1
    )
  markers <-
    markers[order(markers$avg_diff, decreasing = TRUE), ]
  markers <- na.omit(markers[1:100,])
  a <- (markers[, 6])

  ### Output data destinations
  # create folder
  dir.create(file_create)
  pdf(paste(file_create, '/TSNE.pdf', sep = ''))
  TSNEPlot(sub_obj)
  dev.off()

  dir.create(paste(file_create, 'Cluster', sep = '/'))
  dir.create(paste(file_create, 'Violin', sep = '/'))
  for (i in a) {
    # cluster maps
    pdf(paste(paste(file_create, 'Cluster', '', sep = '/'),
              i,
              '.pdf',
              sep = ''))
    FeaturePlot(sub_obj,
                c(i),
                cols.use = c("grey", "blue"),
                pt.size = 2)
    dev.off()

    #violin plot
    pdf(paste(paste(file_create, 'Violin', '', sep = '/'), i, '.pdf', sep = ''))
    VlnPlot(sub_obj, c(i))
    dev.off()
  }
}

# runs clustering algorithm
run_clustering <-function(cochobj, file_dest = '.', dir_add){
  celltocluster <- data.frame(row.names = cochobj@cell.names)
  celltocluster$cellnames <- cochobj@cell.names

  #save resolutions
  res_keep <- data.frame('cluster'= NA,'res'= NA, 'num_clusters' =NA)
  res_counter <- 1
   # Initial Clusters need to have already been initiated
  ### get number of clusters given as maxer
  a <- cochobj@ident
  # get iteration
  class(a) <- 'numeric'
  it_max <- max(a) - 1

  for (j in 0:it_max) {
    print(j)
    ### get all the cell names to put in the cluster definitions

    ###
    # Function to iterate through resolutions

    # get the subset to test
    sub <- SubsetData(cochobj, ident.use = toString(j))

    if (length(sub@data.info[,1]) >20){

      sub <-PCA(sub,pc.genes = sub@var.genes,do.print = TRUE,pcs.print = 5,genes.print = 5)

      sub <- ProjectPCA(sub)
      sub <- RunTSNE(sub,dims.use = 1:10,do.fast = T, perplexity = 6)
      ###############################################

      #find clusters

      set_res <- find_res(sub)
      print(set_res)
      if (set_res >0){
        print('Subclustering!')
        #### label everything
        sub <-FindClusters(sub,pc.use = 1:10, resolution = set_res, print.output = 0, save.SNN = T)
      }
        # Find number of clusters
        a <- sub@ident
        # get iteration
        class(a) <- 'numeric'
        maxer <- max(a)

      #save(sub, file = paste(paste(paste(file_dest, j,sep = ''),'.Robj', sep = '')))

      #terminate and set cells if no subgroups
      if (maxer == 1 || set_res < 0){
        framer <- data.frame(row.names = sub@cell.names)
        framer$cellnames <- row.names(framer)
        framer$clusterid <- j
        celltocluster <- merge(celltocluster, framer, all = TRUE)

        #####
        # document no clusters
        res_keep[res_counter, 1] <- j
        res_keep[res_counter, 2] <- 0
        res_keep[res_counter, 3] <- NA

        res_counter <- res_counter + 1
      }
      else{
        # for each group
        # Save resolutions
        print('Doing Subgroups!')
        res_keep[res_counter, 1] <- j
        res_keep[res_counter, 2] <- set_res
        res_keep[res_counter, 3] <- maxer

        res_counter <- res_counter + 1

        # plot graph and then matrix
        pdf(paste(paste(file_dest,j,sep = ''), '.pdf', sep = ''))
        TSNEPlot(sub)
        dev.off()

        # matrix
        getmatrix(sub,file_dest, j)

        plottingj(sub, file_dest,j)

        for (k in 0:(maxer - 1)) {
          print(paste(j, k, sep = '  '))
          # get the cell names and put it into celltocluster
          sub2 <- SubsetData(sub, ident.use = toString(k))

          if (length(sub2@data.info[,1]) >20){
            ############################################################################
            # reiterate what we did earlier: set up new subcluster
            sub2 <-PCA(sub2, pc.genes = sub2@var.genes, do.print = TRUE,pcs.print = 5,genes.print = 5)
            sub2 <- ProjectPCA(sub2)
            sub2 <- RunTSNE(sub2,dims.use = 1:10,do.fast = T,perplexity = 4)

            ################################################
            # find clusters

            set_res <- find_res(sub2)
            print(set_res)
            # if the subcluster has more than 20 cells do the next part and positive result
            if (set_res > 0){
              print('Subclustering!')

              sub2 <-FindClusters(sub2,pc.use = 1:10,resolution = set_res, print.output = 0,save.SNN = T)
            }
              # Find number of clusters
              a <- sub2@ident
              # get iteration
              class(a) <- 'numeric'
              maxer2 <- max(a)


            #save obj
            #save(sub2, file = paste(paste(file_dest, paste(j,k,sep='.'),sep = ''),'.Robj', sep = ''))

            #terminate and set cells if no subgroups
            if (maxer2 == 1 | set_res < 0) {
              framer <- data.frame(row.names = sub2@cell.names)
              framer$cellnames <- row.names(framer)
              framer$clusterid <- paste(j, k, sep = '.')
              celltocluster <- merge(celltocluster, framer, all = TRUE)

              # document no clusters
              res_keep[res_counter, 1] <- paste(j,k, sep = '.')
              res_keep[res_counter, 2] <- 0
              res_keep[res_counter, 3] <- NA

              res_counter <- res_counter + 1

            }
            else{

              #res counter 2 deep
              res_keep[res_counter, 1] <- paste(j,k, sep = '.')
              res_keep[res_counter, 2] <- set_res
              res_keep[res_counter, 3] <- maxer2

              res_counter <- res_counter + 1

              # get matrix and plot tsne
              pdf(paste(paste(file_dest,paste(j,k,sep = '.'),sep = ''), '.pdf', sep = ''))
              TSNEPlot(sub2)
              dev.off()

              getmatrix(sub2,file_dest, paste(j,k,sep = '.'))

              plottingj(sub2, file_dest,paste(j,k,sep = '.'))


              # iterate through and assign values

              for (l in 0:(maxer2 - 1)) {
                sub3 <- SubsetData(sub2, ident.use = toString(l))
                framer2 <- data.frame(row.names = sub3@cell.names)
                framer2$cellnames <- row.names(framer2)
                framer2$clusterid <- paste(j, k,l, sep = '.')
                celltocluster <- merge(celltocluster, framer2, all=TRUE)

                #save object
                #save(sub3, file = paste(paste(file_dest, paste(j,k,l,sep='.'),sep = ''),'.Robj', sep = ''))

              }
            }
          }
          # if less than 20 cells. just output the data
          else{
            # add to celltocluster
            framer <- data.frame(row.names = sub2@cell.names)
            framer$cellnames <- row.names(framer)
            framer$clusterid <- paste(j, k, sep = '.')
            celltocluster <- merge(celltocluster, framer, all = TRUE)
            #save(sub2, file = paste(paste(file_dest, paste(j,k,sep='.'),sep = ''),'.Robj', sep = ''))
            # document no clusters
            res_keep[res_counter, 1] <- paste(j,k, sep = '.')
            res_keep[res_counter, 2] <- 0
            res_keep[res_counter, 3] <- NA

            res_counter <- res_counter + 1
          }
        }
      }
    }
    # if less than 20 cells. just output the data
    else{
      #add to celltocluster
      framer <- data.frame(row.names = sub@cell.names)
      framer$cellnames <- row.names(framer)
      framer$clusterid <- j
      celltocluster <- merge(celltocluster, framer, all = TRUE)
      #save(sub, file = paste(paste(paste(file_dest, j,sep = ''),'.Robj', sep = '')))

      #####
      # document no clusters
      res_keep[res_counter, 1] <- j
      res_keep[res_counter, 2] <- 0
      res_keep[res_counter, 3] <- NA

      res_counter <- res_counter + 1
    }
  }

  write.csv(na.omit(celltocluster), paste(dir_add, "/assign_test.csv", sep = ''), row.names = TRUE)
  write.csv(res_keep, paste(dir_add, "/res_test.csv", sep = ''),row.names = TRUE)

}

# gets the stats for each subgroup provided
get_stats <- function(tenx, file_dest,d){
  num_find <- 50
  aoe <- c("Group", "cell_number", "avg_gene", "avg_read", "avg_umi")
  for (i in 1:num_find){
    aoe <- c(aoe, paste('top_',i, sep = ""))
  }
  df <- data.frame(aoe)
  #initialize matrix
    for (groups in levels(tenx@ident)){
      subgroup <-SubsetData(tenx, ident.use = groups)
      # group name
      aod <- c(groups)
      # cell number
      aod <- c(aod, length(subgroup@data.info$nGene))
      # avg_gene
      aod <- c(aod, mean(subgroup@data.info$nGene))
      # avg_read
      aod <- c(aod, sum(subgroup@raw.data)/length(subgroup@data.info$nGene))
      # avg_umi
      aod <- c(aod, mean(subgroup@data.info$nUMI))
      # top 10 diff genes
      markers <- FindMarkers(tenx, groups)
      top_markers <- row.names(markers)[1:num_find]
      for (topm in top_markers){
        aod <- c(aod, topm)
      }
      print(aod)
      df[groups] <-aod

      print(groups)
    }

  write.csv(df,
            paste(file_dest, d, '_stats.csv', sep = ''),
            row.names = TRUE)
}
