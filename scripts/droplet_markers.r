library(SingleCellExperiment)
library(MetaMarkers)
library(tidyverse)
library(ComplexHeatmap)
source('~/r_utils/read_loom.R')

get_ds_markers <- function(sce) {
    ds_name <- unique(sce$study_id)[[1]]
    message(ds_name)
    sce$scNym <- factor(sce$scNym)
    markers <- compute_markers(counts(sce), sce$scNym)
    export_markers(markers,
                   paste0(
                       '~/pseudotime/data/markers/markers_',
                       ds_name,
                       '_scNym.csv'
                   ))
    markers
}
if (!file.exists('processed_droplet_data_no_OBSM.loom')){
    system('curl -O ftp://gillisdata.cshl.edu/data/HSC_atlas/processed_droplet_data_no_OBSM.loom')
}
sce <-
    get_sce_alt('processed_droplet_data_no_OBSM.loom')
sce_TM <- sce[, sce$study_id %in% c('tm10x', 'tmSS')]
studies <- sce$study_id %>% unique()
exclude_studies <-
    c('tmSS', 'tm10x', 'tiknova', 'rf_LARRY2', 'rf_LARRY3')
studies <- studies[!studies %in% exclude_studies]
sce <- sce[, sce$study_id %in% studies]
#Compute Markers By Dataset
ds_markers <- lapply(studies, function(i) {
    get_ds_markers(sce[, sce$study_id == i])
})

#Compute MetaMarkers
meta_marks <- make_meta_markers(ds_markers, detailed_stats = T)

cts <- meta_marks$cell_type %>% unique()

plot_pareto_summary(meta_marks)
lapply(cts, function(x) {
    plot_pareto_markers(meta_marks, x, min_recurrence = 0)
} + ggtitle(x))

lapply(seq_along(studies), function(i) {
    export_markers(
        ds_markers[[i]],
        paste0(
            '~/pseudotime/data/markers/droplet_scNym_',
            studies[[i]],
            '.csv'
        )
    )
})
export_meta_markers(meta_marks,
                    '~/pseudotime/data/markers/droplet_scNym_metamarkers.csv',
                    studies)
plot_pareto_summary(meta_marks, min_recurrence = 3) + theme_classic()
ct_order <-
    meta_marks %>% filter(rank <= 3) %>% select(cell_type) %>% unlist() %>% unname()

top_3_markers <-
    meta_marks %>%
    filter(rank <= 3) %>%
    select(gene) %>%
    unlist() %>%
    unname()
dup <- !duplicated(top_3_markers)
top_3_markers <- top_3_markers[dup]
ct_order <- ct_order[dup]

expr <- counts(sce[top_3_markers,])

label_matrix <-
    design_matrix(MetaNeighbor::makeClusterName(sce$study_id, sce$scNym))
label_matrix <-
    scale(label_matrix,
          center = F,
          scale = colSums(label_matrix))

centroids <- scale(t(expr %*% label_matrix))
centroids <- centroids %>%
    as.tibble(rownames = 'cluster') %>%
    pivot_longer(cols = -cluster,
                 names_to = 'gene',
                 values_to = 'average_expression') %>%
    mutate(
        study = MetaNeighbor::getStudyId(cluster),
        cell_type = MetaNeighbor::getCellType(cluster)
    ) %>%
    group_by(gene, cell_type) %>%
    summarise(average_expression = mean(average_expression)) %>%
    pivot_wider(id_cols = everything(),
                names_from = cell_type,
                values_from = average_expression) %>%
    column_to_rownames('gene')

centroids_top <- centroids[top_3_markers, ]
label_order <- labels(as.dendrogram(hclust(dist(t(
    centroids_top
)))))
gene_order <-
    factor(ct_order, levels = label_order, ordered = T) %>% order()

pdf(
    '~/pseudotime/figures/scNym_cluster_figures/scNym_marker_heatmap.pdf',
    width = 7,
    height = 7
)
centroids_top[gene_order, label_order]  %>%
    as.matrix() %>%
    write.csv(.,
              '~/pseudotime/data/markers/marker_droplet_top_genes_heatmap.csv',
              quote = F)
centroids_top[gene_order, label_order] %>%
    as.matrix() %>%
    Heatmap(
        cluster_rows = F,
        cluster_columns = F,
        width = 1,
        height = 1,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        column_title = 'scNym Labeled Clusters',
        row_title = 'Marker Genes',
        show_row_names = T,
        name = 'Average \nExpression',
        col = circlize::colorRamp2(c(-4, 0, 4), c('blue', 'white', 'red'))
    )
dev.off()


exprs_TM <- counts(sce_TM[top_3_markers,])
label_matrix_TM <-
    design_matrix(MetaNeighbor::makeClusterName(sce_TM$study_id, sce_TM$scNym))
label_matrix_TM <-
    scale(label_matrix_TM,
          center = F,
          scale = colSums(label_matrix_TM))

centroids_TM <- scale(t(exprs_TM %*% label_matrix_TM))
centroids_TM <- centroids_TM %>%
    as.tibble(rownames = 'cluster') %>%
    pivot_longer(cols = -cluster,
                 names_to = 'gene',
                 values_to = 'average_expression') %>%
    mutate(
        study = MetaNeighbor::getStudyId(cluster),
        cell_type = MetaNeighbor::getCellType(cluster)
    ) %>%
    group_by(gene, cell_type) %>%
    summarise(average_expression = mean(average_expression)) %>%
    pivot_wider(id_cols = everything(),
                names_from = cell_type,
                values_from = average_expression) %>%
    column_to_rownames('gene')

centroids_TM[rownames(centroids_top[gene_order, ]), label_order] %>%
    as.matrix() %>%
    Heatmap(
        cluster_rows = F,
        cluster_columns = F,
        width = 1,
        height = 1,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5),
        column_title = 'scNym Labeled Clusters',
        row_title = 'Marker Genes',
        name = 'Tabula Muris \n Expression',
        col = circlize::colorRamp2(c(-4, 0, 4), c('blue', 'white', 'red'))
    )