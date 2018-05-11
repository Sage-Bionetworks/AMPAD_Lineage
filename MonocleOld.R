RunMonocleOld <- function(Dat, Labels, num_paths=2){ 
  
  library('monocle')
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  Genes <- c(1:dim(Dat)[1])
  Genes <- data.frame(Genes)
  Labels <- data.frame(Labels)
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = Genes)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd)
  
  HSMM <- reduceDimension(HSMM, use_irlba=FALSE)
  HSMM <- orderCells(HSMM, num_paths=num_paths, reverse=TRUE)
  plot_spanning_tree(HSMM, color_by="Labels")
  
  return(HSMM)
  
}


ReduceDimension2 <- function (cds, max_components = 2, use_irlba = TRUE) 
{
  library(fastICA)
  FM <- exprs(cds)
  if (is.null(fData(cds)$use_for_ordering) == FALSE) 
    FM <- FM[fData(cds)$use_for_ordering, ]
  FM <- FM[rowSds(FM) > 0, ]
  FM <- log10(FM + 1)
  init_ICA <- ica_helper(t(FM), max_components, use_irlba = use_irlba)
  x_pca <- t(t(FM) %*% init_ICA$K)
  print('i was here')
  W <- t(init_ICA$W)
  weights <- W
  A <- t(solve(weights) %*% t(init_ICA$K))
  colnames(A) <- colnames(weights)
  rownames(A) <- rownames(FM)
  S <- weights %*% x_pca
  rownames(S) <- colnames(weights)
  colnames(S) <- colnames(FM)
  reducedDimW(cds) <- W
  reducedDimA(cds) <- A
  reducedDimS(cds) <- S
  cds
}
