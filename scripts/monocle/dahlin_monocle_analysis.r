library(tidyverse)
library(monocle3)
library(SingleCellExperiment)
library(tradeSeq)
library(BiocParallel)
source('../beam_fit/read_loom.R')
source('monocle_utils.r')

if(!file.exists('processed_droplet_obs_data.csv')){
  system('curl -O ftp://gillisdata.cshl.edu/data/HSC_atlas/processed_droplet_obs_data.csv')
}
if(!file.exists('dahlin_raw.loom')){
  system('curl -O ftp://gillisdata.cshl.edu/data/HSC_atlas/dahlin_raw.loom')
}
metadata <-
  read.csv(
    'processed_droplet_obs_data.csv',
    stringsAsFactors = F
  )
sce <-
  get_sce_alt('dahlin_raw.loom')
dahlin_meta <- metadata %>% filter(study_id == 'dahlin' &
                                     X %in% colnames(sce)) %>% column_to_rownames('X')
sce$scNym <- dahlin_meta[colnames(sce), 'scNym'] %>% as.vector()
sce$scNym %>% table() %>% sort()
cds <- build_cds(sce, use_clusterCol = F)
plot_principle_points(cds)

plot_cells(cds, color_cells_by = 'scNym')
cds <- monocle_graph(cds, minimal_branch_len = 100)
plot_principle_points(cds)

plot_cells(cds, color_cells_by = 'scNym')
plot_cells(
  cds,
  genes = c('Gata1', 'Procr', 'Elane', 'Irf8', 'Dntt', 'Cd19', 'Itga2b'),
  show_trajectory_graph = F,
  label_cell_groups = F
)

max_procr <- get_root(cds)
root_node <- get_earliest_principal_node(cds, max_procr)

#Root when using scNym Clusters
root_node <- 'Y_46'
cds <- order_cells(cds, root_pr_nodes = root_node)
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_principal_points = T)

## Node Y_1326 was identified by looking at Gata1 and Itga2b expression
#path_vect <- get_eryth_path(cds, root_node  , 'Y_1326')
path_vect <- get_eryth_path(cds, root_node, 'Y_264')

plot_cells(cds[, path_vect], color_cells_by = 'pseudotime')
plot_cells(cds[, path_vect],
           genes = c('Procr', 'Gata1'),
           show_trajectory_graph = F)
ps_fn <-
  '/home/bharris/pseudotime/data/monocle/dahlin_erythroid_psuedotime.csv'

compare_ps(cds, ps_fn)

colnames(ps) <- c('pseudotime', 'eryth')
# write.table(ps,
#             ps_fn,
#             sep = ',', col.names = T,
#             quote=F)
sum(path_vect) / dim(sce)[2]

saveRDS(cds,
        '~/pseudotime/data/monocle/dahlin_cds_reduced_complexity.RDS')


plot_cells(
  cds,
  color_cells_by = 'pseudotime',
  label_branch_points = F,
  label_leaves = F,
  label_roots = T
)
