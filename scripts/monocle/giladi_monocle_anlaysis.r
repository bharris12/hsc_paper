library(tidyverse)
library(monocle3)
library(SingleCellExperiment)

source('../beam_fit/read_loom.R')
source('monocle_utils.r')
if(!file.exists('processed_droplet_obs_data.csv')){
  system('curl -O ftp://gillisdata.cshl.edu/data/HSC_atlas/processed_droplet_obs_data.csv')
}
if(!file.exists('giladi_raw.loom')){
  system('curl -O ftp://gillisdata.cshl.edu/data/HSC_atlas/giladi_raw.loom')
}
metadata <-
  read.csv(
    'processed_droplet_obs_data.csv',
    stringsAsFactors = F
  )
sce <-
  get_sce_alt('giladi_raw.loom')
giladi_meta <- metadata %>% filter(study_id == 'giladi' &
                                     X %in% colnames(sce)) %>% column_to_rownames('X')

sce$scNym <- giladi_meta[colnames(sce), 'scNym'] %>% as.vector()
sce$scNym %>% table() %>% sort()
cds <- build_cds(sce)

cds <- monocle_graph(cds, minimal_branch_len = 75)
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
cds <- order_cells(cds, root_pr_nodes = root_node)
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_principal_points = T)
## Node Y_1326 was identified by looking at Gata1 and Itga2b expression
path_vect <- get_eryth_path(cds, root_node  , 'Y_368')

plot_cells(cds[, path_vect], color_cells_by = 'pseudotime')
plot_cells(cds[, path_vect],
           genes = c('Procr', 'Gata1'),
           show_trajectory_graph = F)
ps_fn <-
  '/home/bharris/pseudotime/data/monocle/giladi_erythroid_psuedotime.csv'
compare_ps(cds ,ps_fn)
saveRDS(cds,
        '~/pseudotime/data/monocle/giladi_cds_reduced_complexity.RDS')

ps<-cbind(pseudotime(cds),path_vect)
colnames(ps)<-c('pseudotime','eryth')

# write.table(ps,
#             ps_fn,
#             sep = ',', col.names = T,
#             quote=F)
sum(path_vect) / dim(sce)[2]
plot_cells(cds, color_cells_by = 'pseudotime', label_branch_points = F, label_leaves = F, label_roots = T)
