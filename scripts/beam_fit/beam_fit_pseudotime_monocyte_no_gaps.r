library(tidyverse)
library(monocle)
library(SingleCellExperiment)
source('read_loom.R')

build_cds <- function(sce) {
    cds <-
        monocle::newCellDataSet(
            counts(sce),
            phenoData = AnnotatedDataFrame(colData(sce) %>% as.data.frame()),
            featureData = AnnotatedDataFrame(rowData(sce) %>% as.data.frame())
        ) %>% estimateSizeFactors()
}

if(!file.exists('erythroid_and_monocyte_lineage_adata_no_gaps.loom')){
    system('curl -O  ftp://gillisdata.cshl.edu//data/HSC_atlas/erythroid_and_monocyte_lineage_adata_no_gaps.loom')
}
sce <-
    get_sce_alt('erythroid_and_monocyte_lineage_adata_no_gaps.loom')
studies <- unique(sce$study_id)
results <- data.frame()
for (study in studies) {
    message(study)
    sce_s <- sce[, sce$study_id == study]
    test <- build_cds(sce_s)
    
    diff_gene <-
        monocle::differentialGeneTest(
            test,
            fullModelFormulaStr = '~sm.ns(monocle_ps,df=3)',
            cores = 1,
            verbose = F
        )
    diff_gene$study <- study
    write.csv(
        diff_gene,
        paste(
            '~/pseudotime/data/monocle/gene_models/',
            study,
            '_BEAM_pseudotime_test_monocyte_no_gaps.csv'
        )
    )
    results <-
        rbind(results, diff_gene %>% rownames_to_column(var = 'Gene'))
}
write.csv(
    results,
    '~/pseudotime/data/monocle/gene_models/all_studies_BEAM_pseudotime_test_monocyte_no_gaps.csv'
)
