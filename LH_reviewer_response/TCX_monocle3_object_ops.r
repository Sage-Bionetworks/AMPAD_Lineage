library(monocle3)

counts <- readRDS(file="TCX_processed_countmatrix_Mono3.rds")
gene_short_names <- readRDS(file="TCX_processed_gene_metadata_Mono3.rds")
temp2 <- readRDS(file="TCX_processed_cell_metadata_Mono3.rds")
cds_tcx <- readRDS(file="TCX_cds_Mono3.rds")
class(counts)
class(gene_short_names)
class(temp2)

cds_tcx <- new_cell_data_set(counts, cell_metadata=temp2, gene_metadata=gene_short_names)

cds_tcx <- preprocess_cds(cds_tcx, method="PCA", norm_method="size_only", num_dim=30)
cds_tcx <- reduce_dimension(cds_tcx)
plot_pc_variance_explained(cds_tcx)

cds_tcx <- cluster_cells(cds_tcx)
cds_tcx@clusters$UMAP$partitions[cds_tcx@clusters$UMAP$partitions == "2"] <- "1"


cds_tcx <- learn_graph(cds_tcx, use_partition=FALSE)
cds_tcx >- order_cells(cds_tcx, root_pr_nodes=get_earliest_principal_node(cds_tcx))

plot_cells(cds_tcx, xcolor_cells_by="pseudotime")
plot_cells(cds_tcx, color_cells_by="partition")







MonRun <- tobit(counts, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_names)
