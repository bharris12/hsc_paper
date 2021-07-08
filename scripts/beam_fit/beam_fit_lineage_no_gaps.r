library(tidyverse)
library(monocle)
library(SingleCellExperiment)
source('read_loom.R')


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
    test <-
        monocle::newCellDataSet(
            counts(sce_s),
            phenoData = AnnotatedDataFrame(colData(sce_s) %>% as.data.frame()),
            featureData = AnnotatedDataFrame(rowData(sce_s) %>% as.data.frame())
        )
    message('Estimate Size Factors')
    test <- estimateSizeFactors(test)
    
    
    diff_gene <-
        monocle::differentialGeneTest(
            test,
            fullModelFormulaStr = '~sm.ns(monocle_ps,df=3) * new_lineage',
            reducedModelFormulaStr = '~sm.ns(monocle_ps,df=3)',
            cores = 1,
            verbose = F
        )
    diff_gene$study <- study
    write.csv(
        diff_gene,
        paste0(
            '~/pseudotime/data/monocle/gene_models/',
            study,
            '_BEAM_lineage_test_no_gaps.csv'
        )
    )
    results <-
        rbind(results, diff_gene %>% rownames_to_column(var = 'Gene'))
}
write.csv(
    results,
    '~/pseudotime/data/monocle/gene_models/all_studies_BEAM_lineage_test_no_gaps.csv'
)
