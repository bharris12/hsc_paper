library(tidyverse)
library(monocle3)
library(SingleCellExperiment)

source('../beam_fit/read_loom.R')
source('monocle_utils.r')

## Code commented out and intermediate results saved for speed
sce_all <-
  get_sce_alt('processed_droplet_data.loom')
tm_ds <- sce_all$study_id %in% c('tm10x', 'tmSS')
sce_all <- sce_all[,!tm_ds]

sce_all$scNym %>% table() %>% sort()
sce_all$study_id %>% table() %>% sort()


##Saving intermediates because this takes a lot of RAM and time
saveRDS(sce_all,
        '~/pseudotime/data/monocle/integrated_unprocessed_sce.RDS')
sce_all <-
  readRDS('~/pseudotime/data/monocle/integrated_unprocessed_sce.RDS')

cds <-
  new_cell_data_set(
    counts(sce_all),
    cell_metadata = colData(sce_all),
    gene_metadata = rowData(sce_all)
  )



cds <-
  preprocess_cds(cds, verbose = T, alignment_group = 'study_id')
cds <- align_cds(cds, alignment_group = 'study_id', verbose = T)
saveRDS(cds, '~/pseudotime/data/monocle/integrated_aligned_cds.RDS')

cds <-
  readRDS('~/pseudotime/data/monocle/integrated_aligned_cds.RDS')
lfile <-
  loomR::connect(
    '~/pseudotime/data/processed_droplet_data.loom',
    mode = "r",
    skip.validate = T
  )
obs_names <- lfile[['col_attrs/obs_names']][]
umap1 <- lfile[['col_attrs/X_umap']][1,]
umap2 <- lfile[['col_attrs/X_umap']][2,]

umap <- cbind(umap1, umap2)
rownames(umap) <- obs_names
colnames(umap) <- c('UMAP1', 'UMAP2')
umap <- umap[colnames(cds),]
if (lfile$is_valid) {
  lfile$close_all()
}


cds <-
  reduce_dimension(cds, verbose = T, cores = 12)

cds <- cluster_cells(cds, cluster_method = 'louvain')

#cds@clusters@listData[["UMAP"]][["clusters"]] <-
#  as.numeric(as.factor(colData(sce)[, cluster_col]))
rowData(cds)$gene_short_name <- rownames(cds)
plot_cells(cds,
           genes = c('Procr', 'Cd34', 'Gata1'),
           show_trajectory_graph = F)
plot_cells(cds, color_cells_by = 'scNym', label_cell_groups = F)

plot_cells(cds, color_cells_by = 'study_id', label_cell_groups = F)

gc()
cds <- learn_graph(cds, close_loop = F, verbose = T)
gc()
max_procr <- get_root(cds)
root_node <- get_earliest_principal_node(cds, max_procr)
cds <- order_cells(cds, root_pr_nodes = root_node)
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_principal_points = T)

deg1_nodes <- igraph::degree(principal_graph(cds)[["UMAP"]]) ==  1
deg1_nodes <- names(deg1_nodes)[deg1_nodes]
#deg1_nodes <- as.numeric(sapply(strsplit(deg1_nodes, "_"), "[[", 2))

principle_leaves <-
  t(cds@principal_graph_aux$UMAP$dp_mst) %>% as.data.frame()
principle_leaves <-
  principle_leaves %>% filter(rownames(principle_leaves) %in% deg1_nodes)

plot_cells(
  cds,
  color_cells_by = 'pseudotime',
  show_trajectory_graph = F,
  label_cell_groups = F
) + facet_wrap(factor(colData(cds)$study_id))


path_vect <- get_eryth_path(cds, root_node  , 'Y_1192')

plot_cells(cds[, path_vect],
           color_cells_by = 'pseudotime',
           show_trajectory_graph = F)

plot_cells(
  cds[, path_vect],
  color_cells_by = 'pseudotime',
  show_trajectory_graph = F,
  label_cell_groups = F
) + facet_wrap(factor(colData(cds[, path_vect])$study_id))


plot_cells(cds[, path_vect],
           genes = c('Procr', 'Gata1','Klf1'),
           show_trajectory_graph = F)
ps <- cbind(pseudotime(cds), path_vect)
colnames(ps) <- c('pseudotime', 'eryth')
ps_fn <-
  '/home/bharris/pseudotime/data/monocle/integrated_erythroid_pseudotime.csv'
write.table(ps,
            ps_fn,
            sep = ',',
            col.names = T,
            quote = F)
saveRDS(cds, '~/pseudotime/data/monocle/integrated_monocle_cds.RDS')

umap <- umap %>% as.data.frame()
sum(path_vect) / dim(cds)[2]

plot_cells(cds, color_cells_by = 'pseudotime', label_branch_points = F, label_leaves = F, label_roots = T)
