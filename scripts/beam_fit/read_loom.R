#'Gets a SingleCellExperiment from a loom file.
#' Uses native Seurat and LoomR functions
#' Will natively get you counts and logcounts, but is slower that above version
#' @param filename The path to th file being read
#' @return SingleCellExperiment object 
#' @export
get_sce_alt <-function(filename){
    seurat_vers<- SeuratDisk::LoadLoom(filename,features='var_names',cells='obs_names')
    sce<-Seurat::as.SingleCellExperiment(seurat_vers)
    return(sce)
