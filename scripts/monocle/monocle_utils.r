library(monocle3)
build_cds <-
  function(sce,
           cluster_col = 'scNym',
           num_dim = 30,
           k = 20,
           scale = T,
           use_hvg = T,
           use_clusterCol = F,
           remove_clusters = c(),
           remove_genes = c()) {
    sce <- sce[,!colData(sce)[, cluster_col] %in% remove_clusters]
    sce <- sce[!rownames(sce) %in% remove_genes, ]
    if (use_hvg)
    {
      hvg <- scran::getTopHVGs(scran::modelGeneVar(sce))
    } else{
      hvg <- NULL
    }
    print(length(hvg))
    print(dim(sce))
    cds <-
      new_cell_data_set(counts(sce),
                        cell_metadata = colData(sce),
                        gene_metadata = rowData(sce))
    cds <-
      preprocess_cds(
        cds,
        use_genes = hvg,
        num_dim = num_dim,
        scaling = scale,
        verbose = T
      ) %>%
      reduce_dimension(cores = 12, verbose = T) %>%
      cluster_cells(cluster_method = 'louvain',
                    k = k,
                    verbose = T)
    
    if (use_clusterCol) {
      cds@clusters@listData[["UMAP"]][["clusters"]] <-
        as.factor(colData(sce)[, cluster_col])
    }
    
    rowData(cds)$gene_short_name <- rownames(cds)
    return(cds)
  }

monocle_graph <-
  function(cds,
           close_loop = F,
           minimal_branch_len = 10,
           rann.k = 25) {
    learn_graph_control <-
      list(minimal_branch_len = minimal_branch_len, rann.k = rann.k)#, ncenter=100)
    learn_graph(
      cds,
      verbose = T,
      close_loop = close_loop,
      learn_graph_control = learn_graph_control
    )
  }

get_eryth_path <- function(cds, n1, n2) {
  path_to_ertyh <-
    igraph::all_shortest_paths(cds@principal_graph[['UMAP']], from = n1, to =
                                 n2)$res[[1]]  %>% as.vector()
  path_filter <-
    cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex %in% path_to_ertyh
  path_filter
}
collapse_branch <- function(cds, n1, n2) {
  replace_node <- readr::parse_number(n1)
  fix_path <- get_eryth_path(cds, n1, n2)
  cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex[fix_path] <-
    replace_node
  cds
}

get_earliest_principal_node <- function(cds, cluster) {
  cell_ids <- which(clusters(cds) == cluster)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds),])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

get_root <- function(cds, gene = 'Procr') {
  max_procr <-
    sapply(split(counts(cds[gene]) %>%  as.vector(),
                 clusters(cds)), mean) %>%
    which.max() %>%
    as.vector()
  max_procr
}

plot_gene <-
  function(cds ,
           gene,
           groupby_study = F,
           study_name = 'study') {
    if ('study_id' %in% colnames(colData(cds))) {
      study_name <- factor(colData(cds)$study_id)
    }
    tmp <- data.frame(
      gene = assay(cds, 'counts')[gene,] %>% as.vector(),
      study_id = study_name,
      umap1 = reducedDims(cds)$UMAP[, 1] %>% as.vector(),
      umap2 = reducedDims(cds)$UMAP[, 2] %>% as.vector()
    )
    p1 <-
      ggplot(tmp %>% mutate(gene = log2(gene + 1)), aes(x = umap1,
                                                        y = umap2,
                                                        color = gene))  +
      geom_point(size = .15, alpha = .7)  +
      theme_classic() +
      scale_color_distiller(palette = 'Blues', direction = 1) + theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y =
          element_blank(),
        axis.ticks = element_blank(),
      ) + labs(color = gene)
    
    if (groupby_study) {
      p1 <- p1 + facet_wrap(vars(study_id))
    }
    p1
  }

get_mst_graph <- function(cds) {
  non_2_degree <- igraph::degree(cds@principal_graph[['UMAP']]) != 2
  mst_graph <-
    t(cds@principal_graph_aux$UMAP$dp_mst) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'nodename') %>%
    filter(non_2_degree) %>%
    mutate(degree = factor(igraph::degree(cds@principal_graph[['UMAP']])[nodename])) %>%
    mutate(leaf = degree == 1)
  
  if (!('UMAP1' %in% colnames(mst_graph))) {
    mst_graph %>% rename(V1 = 'UMAP1') -> mst_graph
  }
  if (!('UMAP2' %in% colnames(mst_graph))) {
    mst_graph %>% rename(V2 = 'UMAP2') -> mst_graph
  }
  mst_graph
}
get_closest_split <- function(cds, mst_graph) {
  mst_dist <- igraph::distances(cds@principal_graph[['UMAP']])
  
}
plot_principle_points <- function(cds) {
  non_2_degree <- igraph::degree(cds@principal_graph[['UMAP']]) != 2
  
  mst_dist <- igraph::distances(cds@principal_graph[['UMAP']])
  mst_graph <- get_mst_graph(cds)
  
  ct_df <- data.frame(
    UMAP1 = reducedDims(cds)$UMAP[, 1] %>% as.vector(),
    UMAP2 = reducedDims(cds)$UMAP[, 2] %>% as.vector(),
    scNym = factor(colData(cds)[, 'scNym'] %>% as.vector())
  )
  
  leaf_vs_bisect <- mst_dist[mst_graph %>%
                               filter(degree == 1) %>%
                               select(nodename) %>%
                               unlist() %>%
                               as.vector(),
                             mst_graph %>%
                               filter(degree != 1) %>%
                               select(nodename) %>%
                               unlist() %>%
                               as.vector()]
  
  closest_split <-
    colnames(leaf_vs_bisect)[apply(leaf_vs_bisect, 1, which.min)]
  names(closest_split) <- rownames(leaf_vs_bisect)
  ggplot(ct_df, aes(x = UMAP1, y = UMAP2, color = scNym)) +
    geom_point(alpha = .01) +
    ggrepel::geom_label_repel(data =
                                mst_graph,
                              aes(label = nodename, color = leaf),
                              max.overlaps = 100) +
    geom_point(data = mst_graph %>% filter(nodename %in%
                                             union(
                                               names(closest_split),
                                               unlist(closest_split)
                                             )),
               aes(label = nodename, color =
                     factor(leaf))) +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y =
        element_blank(),
      axis.ticks = element_blank()
    )
}

prune_tree <- function(cds , keep_nodes) {
  non_2_degree <- igraph::degree(cds@principal_graph[['UMAP']]) != 2
}
compare_ps <- function(cds, ps_fn) {
  old_ps <- read.table(ps_fn, sep = ',')$pseudotime
  ggplot(
    data.frame(old_pseudotime = old_ps, new_pseudotime = pseudotime(cds)),
    aes(x = old_pseudotime, y = new_pseudotime)
  ) + geom_point(size = .2) + theme_classic() + coord_fixed()
}

check_node_degree<-function(cds,node){
  mst <- principal_graph(cds)$UMAP
  igraph::degree(mst)[node]
}

plot_potential_figures<-function(cds, 
                                 ds_name, 
                                 path='~/pseudotime/figures/monocle_pseudotime/',
                                 extension='.pdf'){
    sce <- get_sce_alt('~/pseudotime/data/monocle//erythroid_and_monocyte_lineage_adata.loom')
    
    cds@principal_graph_aux$UMAP$pseudotime <- pseudotime(cds) / max(pseudotime(cds[,is.finite(pseudotime(cds))]))
    print(head(pseudotime(cds)))
    file_roots <- paste0(path, ds_name)
    plot_cells(cds,
               color_cells_by = 'scNym', 
               label_cell_groups = F, 
               label_branch_points = F, 
               label_roots = T, 
               label_leaves = F) +  
               theme_void() + 
               ggsave(paste0(file_roots,'_snNym_clusters',extension),
                      useDingbats=FALSE)                                                            
    plot_cells(cds,
               color_cells_by = 'pseudotime', 
               label_cell_groups = F, 
               label_branch_points = F, 
               label_roots = T, 
               label_leaves = F) +  
               theme_void()+ 
               ggsave(paste0(file_roots,'_pseudotime_all',extension),
                      useDingbats=FALSE) 
    plot_cells(cds[,colnames(sce)[sce$study_id==ds_name]],
               color_cells_by = 'pseudotime', 
               label_cell_groups = F, 
               label_branch_points = F, 
               label_roots = T, 
               label_leaves = F) +  
               theme_void()+ 
               ggsave(paste0(file_roots,'_pseudotime_erythroid_and_monocyte',extension),
                      useDingbats=FALSE) 
    plot_cells(cds[,colnames(sce)[(sce$study_id==ds_name & sce$lineage =='Erythroid')]],
               color_cells_by = 'pseudotime', 
               label_cell_groups = F, 
               label_branch_points = F, 
               label_roots = T, 
               label_leaves = F) +  
               theme_void()+ 
               ggsave(paste0(file_roots,'_pseudotime_erythroid',extension),
                      useDingbats=FALSE) 
        plot_cells(cds[,colnames(sce)[(sce$study_id==ds_name & sce$lineage =='Monocyte')]],
               color_cells_by = 'pseudotime', 
               label_cell_groups = F, 
               label_branch_points = F, 
               label_roots = T, 
               label_leaves = F) +  
               theme_void()+ 
               ggsave(paste0(file_roots,'_pseudotime_monocyte',extension),
                      useDingbats=FALSE) 
}