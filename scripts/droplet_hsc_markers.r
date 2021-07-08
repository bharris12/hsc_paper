library(SingleCellExperiment)
library(MetaMarkers)
library(tidyverse)

source('./beam_fit/read_loom.R')

get_ds_markers <- function(sce) {
    ds_name <- unique(sce$study_id)[[1]]
    message(ds_name)
    sce$FACS_labels <- factor(sce$FACS_labels)
    markers <- compute_markers(counts(sce), sce$FACS_labels)
    export_markers(
        markers,
        paste0(
            '~/pseudotime/data/markers/markers_',
            ds_name,
            '_HSC_FACS_labels.csv'
        )
    )
    markers
}

if(!file.exists('mouse_hsc_labeled.loom')){
    system('curl -O ftp://gillisdata.cshl.edu//data/HSC_atlas/mouse_hsc_labeled.loom')
}
sce <-
    get_sce_alt('mouse_hsc_labeled.loom')
sce <- sce[, sce$FACS_labels != 'Unknown']

studies <- unique(sce$study_id)
sce <- sce[!rownames(sce) %in% c('Slamf1', 'Cd48', 'Flt3'),]
ds_markers <- lapply(studies, function(i) {
    get_ds_markers(sce[, sce$study_id == i])
})
names(ds_markers) <- studies


#Compute MetaMarkers
meta_marks <- make_meta_markers(ds_markers, detailed_stats = T)

plot_pareto_summary(meta_marks)
meta_marks %>% filter(rank < 50) %>% ggplot(aes(x = auroc, col = cell_type, fill =
                                                    cell_type)) + stat_density()
export_meta_markers(
    meta_marks,
    '~/pseudotime/data/markers/droplet_FACS_labels_metamarkers.csv',
    studies
)
meta_markers_loo = lapply(names(ds_markers), function(leave_out) {
    make_meta_markers(ds_markers[names(ds_markers) != leave_out])
})
names(meta_markers_loo) <- names(ds_markers)


n_genes = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
mean_pr = list()
all_pr = tibble()
for (n in n_genes) {
    pr = list()
    for (test_dataset in studies) {
        true_labels = sce[, sce$study_id == test_dataset]$FACS_labels
        top_markers = filter(meta_markers_loo[[test_dataset]], rank <= n)
        
        ct_scores = score_cells(counts(sce[, sce$study_id == test_dataset]), top_markers)
        ct_enrichment = compute_marker_enrichment(ct_scores)
        
        pr[[test_dataset]] = summarize_precision_recall(ct_enrichment, true_labels, seq(1, 5, 0.1)) %>%
            filter(get_cell_type(marker_set) == true_label)
    }
    pr = bind_rows(pr, .id = "test_dataset")
    pr$n_genes<-n
    all_pr<-rbind(all_pr, pr)
    mean_pr[[as.character(n)]] = pr %>%
        group_by(true_label, score_threshold) %>%
        summarize(f1 = mean(f1, na.rm = TRUE)) %>%
        group_by(true_label) %>%
        filter(f1 == max(f1, na.rm = TRUE)) %>%
        summarize(score_threshold = median(score_threshold),
                  f1 = first(f1))
}
mean_pr = bind_rows(mean_pr, .id = "n_genes") %>%
    mutate(n_genes = as.integer(n_genes))
write_csv(mean_pr,'~/pseudotime/data/markers/drop_FACS_labels_F1_scores.csv')
write_csv(all_pr,'~/pseudotime/data/markers/drop_FACS_labels_all_F1_scores.csv')
mean_pr %>%
    ggplot(aes(x = n_genes, y = f1, col = true_label)) +
    geom_line() +
    theme_bw() +
    scale_x_log10()


aurocs<- list()
n<-10
for (test_dataset in studies) {
    true_labels = sce[, sce$study_id == test_dataset]$FACS_labels
    top_markers = filter(meta_markers_loo[[test_dataset]], rank <= n)
    
    ct_scores = score_cells(counts(sce[, sce$study_id == test_dataset]), 
                            top_markers)
    ct_enrichment = compute_marker_enrichment(ct_scores)
    aurocs[[test_dataset]] <-summarize_auroc(ct_enrichment, true_labels)
}