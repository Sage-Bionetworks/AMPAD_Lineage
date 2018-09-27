PlotMonocleContinuous <- function(pts, lines, G_name, clr = NULL, cp = NULL){
  library("tidyverse")
  library("igraph")
  
  pts <- as.data.frame(t(pts))
  #pts$clr <- as.factor(clr)
  lines <- as.data.frame(t(lines))
  
  #pts <- readRDS("Points.rds")
  #lines <- readRDS("Lines.rds")
  lines_ig <- read_graph(G_name)
  lines_df <- as.data.frame(get.edgelist(lines_ig))
  
  ## Add an id column to the lines (nodes) data frame
  lines <- rowid_to_column(lines, "id")
  
  ## Join with edgelist, then join back with nodes to get start/end points of
  ## edges
  lines_edges <- lines %>%
    inner_join(lines_df, by = c("id" = "V1")) %>%
    inner_join(lines, by = c("V2.y" = "id")) %>%
    ## Rename columns for clarity
    rename(x = V1.x, y = V2.x, xend = V1.y, yend = V2)
  
  ## Plot
  if(is.null(cp)){
    ggplot() +
      geom_point(data = pts, aes(x = V1, y = V2, color = clr)) +
      geom_segment(data = lines_edges, aes(x = x, y = y, xend = xend, yend = yend),size = 1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else{
    ggplot() +
      geom_point(data = pts, aes(x = V1, y = V2, color = clr)) +
      geom_segment(data = lines_edges, aes(x = x, y = y, xend = xend, yend = yend),size = 1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))  + 
      scale_fill_gradient(low = cp[1], high = cp[2])
  }
  

}



PlotMonocleDiscrete <- function(pts, lines, G_name, clr = NULL, cp = NULL){
  library("tidyverse")
  library("igraph")
  
  pts <- as.data.frame(t(pts))
  #pts$clr <- as.factor(clr)
  lines <- as.data.frame(t(lines))
  
  #pts <- readRDS("Points.rds")
  #lines <- readRDS("Lines.rds")
  lines_ig <- read_graph(G_name)
  lines_df <- as.data.frame(get.edgelist(lines_ig))
  
  ## Add an id column to the lines (nodes) data frame
  lines <- rowid_to_column(lines, "id")
  
  ## Join with edgelist, then join back with nodes to get start/end points of
  ## edges
  lines_edges <- lines %>%
    inner_join(lines_df, by = c("id" = "V1")) %>%
    inner_join(lines, by = c("V2.y" = "id")) %>%
    ## Rename columns for clarity
    rename(x = V1.x, y = V2.x, xend = V1.y, yend = V2)
  
  ## Plot
  if(is.null(cp)){
    ggplot() +
      geom_point(data = pts, aes(x = V1, y = V2, color = clr)) +
      geom_segment(data = lines_edges, aes(x = x, y = y, xend = xend, yend = yend),size = 1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else{
    ggplot() +
      geom_point(data = pts, aes(x = V1, y = V2, color = clr)) +
      geom_segment(data = lines_edges, aes(x = x, y = y, xend = xend, yend = yend),size = 1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))  + 
      scale_color_brewer(palette=cp)
  }
  
  
}


