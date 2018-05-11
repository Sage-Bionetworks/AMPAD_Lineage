
RunSlicerTraj <- function(Dat, m, kmin = 5, kmax = 50, by = 5, reg = 2){
  
  # m is instrinsic dimension of data
  # kmin is minimum number of nearest neighbors
  
  library(SLICER)
  
  Dat <- t(Dat)
  
  k = select_k(Dat, kmin, kmax, by)
  print(k)
  traj_lle = lle(Dat, m, k, reg = reg)$Y
  traj_graph = conn_knn_graph(traj_lle,5)
  ends = find_extreme_cells(traj_graph, traj_lle)
  
  l <- list()
  l$traj_lle <- traj_lle
  l$traj_graph <- traj_graph
  l$traj_graph2 <- as.matrix(as_adj(traj_graph))
  l$ends <- ends 
  
  return(l)
  
}

RunSlicerTrajK <- function(Dat, m, k = 5, reg = 2){
  
  # m is instrinsic dimension of data
  # kmin is minimum number of nearest neighbors
  
  library(SLICER)
  
  Dat <- t(Dat)
  
  traj_lle = lle(Dat, m, k, reg = reg)$Y
  
  if (m == 1){
    l <- length(traj_lle)
    traj_lle = rep(traj_lle,2)
    dim(traj_lle) <- c(l,2)
  }
    
  
  
  traj_graph = conn_knn_graph(traj_lle,5)
  ends = find_extreme_cells(traj_graph, traj_lle)
  
  l <- list()
  l$traj_lle <- traj_lle
  l$traj_graph <- traj_graph
  l$traj_graph2 <- as.matrix(as_adj(traj_graph))
  l$ends <- ends 
  
  return(l)
  
}

SlicerOrderCells <- function(traj_graph, start){
  
  #Given a starting cell, orders the other cells
  library(SLICER)
  
  cells_ordered = cell_order(traj_graph, start)
  branches = assign_branches(traj_graph,start)
  
  l <- list()
  
  l$cells_ordered <- cells_ordered
  l$branches <- branches 
  
  return(l)
  
}


WriteToMat <- function(l, fileName){
  
  library(R.matlab)
  writeMat(fileName, l = l)
  
  return(NULL)
}


RunMonocleDefault <- function(Dat, Labels, max_components = 2, meth = 'DDRTree'){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  Genes <- c(1:dim(Dat)[1])
  Genes <- data.frame(Genes)
  Labels <- data.frame(Labels)
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = Genes)
  

  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=negbinomial.size())
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  
  #HSMM <- setOrderingFilter(HSMM, c(1:length(Genes)))
  
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  HSMM <- orderCells(HSMM)
  plot_cell_trajectory(HSMM, color_by="Labels")
  
  return(HSMM)
  
  }


RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree'){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  Genes <- c(1:dim(Dat)[1])
  Genes <- data.frame(Genes)
  Labels <- data.frame(Labels)
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = Genes)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  HSMM <- orderCells(HSMM)
  plot_cell_trajectory(HSMM, color_by="Labels")
  
  return(HSMM)
  
}

Monocle.Plot <- function(MonRun, Labels, Discrete = T, 
                         titleVal = 'Lineage Scatter Plot'){
  
  library(ggplot2)
  library(reshape2)
  
  S <- MonRun@reducedDimS
  K <- MonRun@reducedDimK
  
  if(Discrete == T){
    Labels = factor(Labels)
  }
  
  df1 <- data.frame('x' = S[1,], 'y' = S[2,],
                   'Lab' = Labels)
  
  Pstime <- MonRun@phenoData@data$Pseudotime
  temp <- sort(Pstime, decreasing = F, index.return = T)
  
  
  #df2 <- data.frame('x2' = K[1,temp$ix], 'y2' = K[2,temp$ix])
  
  
  
  if(Discrete==T){
    g <- ggplot(df1, aes(x = x,y = y, color = Lab)) + geom_point(size=3)  + ggtitle(titleVal)
#    geom_point(data = na.omit(dat), aes(x = x2,y = y2))  
  } else{
    g <- ggplot(df1, aes(x = x, y = y, color = Lab)) + 
      geom_point(size=3) + 
      scale_color_gradient(low = "black", high = "red") + ggtitle(titleVal)
#      geom_line(data = df2, aes(x = x,y = y))  
  }
  
  return(g)

}


Mon.Plot.Genes <- function(MonRun, Dat, GeneNames, GeneID){
  
  In <- which(GeneNames == GeneID)
  Lab = as.numeric(Dat[In,])
  g <- Monocle.Plot(MonRun, Labels = Lab, Discrete = F, titleVal = GeneID)
  return(g)
  
}
