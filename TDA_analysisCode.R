library(igraph)
library(TDAmapper)
library("RColorBrewer")

X <- MonRun@reducedDimS
d <- dist(t(X))

filter <- X[1,] # height projection
num_intervals <- 10
percent_overlap <- 50
num_bins_when_clustering <- 10
m3 <- mapper1D(
distance_matrix = d,
filter_values = filter,
# num_intervals = 10, # use default
# percent_overlap = 50, # use default
# num_bins_when_clustering = 10 # use default
)
g3 <- graph.adjacency(m3$adjacency, mode="undirected")
V(g3)$label <- NA
plot(g3, layout = layout.auto(g3) )

Diag <- (temp2$Tissue.Diagnosis == 'TCX.AD') + 0

Col <- rep(0,length(m3$points_in_vertex))

for ( i in 1:length(m3$points_in_vertex)){
  
  for (j in 1:length(m3$points_in_vertex[[i]])){
    
    Col[i] <- Col[i] + Diag[m3$points_in_vertex[[i]][[j]]]
    
  }
  
  Col[i] <- Col[i]/length(m3$points_in_vertex[[i]])
  
}

fine = 500
pal = colorRampPalette(c('red','green'))
graphCol = pal(fine)[as.numeric(cut(Col,breaks = fine))]

plot(g3, layout = layout.auto(g3), vertex.color=graphCol)

Diag <- (temp2$Tissue.Diagnosis == 'TCX.AD') + 0

Ape <- (temp2$Tissue.APOE4 == 'TCX.1') + 0

Col2 <- rep(0,length(m3$points_in_vertex))

for ( i in 1:length(m3$points_in_vertex)){
  
  for (j in 1:length(m3$points_in_vertex[[i]])){
    
    Col2[i] <- Col2[i] + Ape[m3$points_in_vertex[[i]][[j]]]
    
  }
  
  Col2[i] <- Col2[i]/length(m3$points_in_vertex[[i]])
  
}

V(g3)$size <- 5 + Col2/max(Col2)*5
plot(g3, layout = layout.auto(g3), vertex.color=graphCol)
