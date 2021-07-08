library(SingleCellExperiment)
library(monocle3)
library(tidyverse)
library(igraph)

datasets <-c(
  'cheng',
  'dahlin',
  'giladi',
  'rf_LARRY1',
  'tusi_batch1',
  'tusi_batch2',
  'weinreb16',
  'weinreb2',
  'weinreb9'
)

for (ds in datasets){
    message(ds)
    fn <- paste0('~/pseudotime/data/monocle/',ds,'_cds_reduced_complexity.RDS')
    cds<-readRDS(fn)
    
    adj<-as_adjacency_matrix(cds@principal_graph$UMAP)
    
    principle_points_loc<-cds@principal_graph_aux$UMAP$dp_mst %>% t() %>% as.data.frame()
    mst_vertex_info <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
    colnames(mst_vertex_info)<-c('nearest_vertex')
    expand_locs<-principle_points_loc[mst_vertex_info[,1],]
    
    rnames<-rownames(mst_vertex_info)
    mst_vertex_info <- mst_vertex_info %>% as.data.frame() %>% mutate(UMAP1=expand_locs[,1], UMAP2=expand_locs[,2])
    rownames(mst_vertex_info)<-rnames
    
    adjaceny_fn <- paste0('~/pseudotime/data/monocle/mst_adjacency_', ds, '.csv')
    mst_vertex_info_fn <- paste0('~/pseudotime/data/monocle/mst_vertex_info_', ds, '.csv')
    write.csv(as.matrix(adj), adjaceny_fn,quote=F)
    write.csv(mst_vertex_info, mst_vertex_info_fn,quote=F)
}
message('Done')


