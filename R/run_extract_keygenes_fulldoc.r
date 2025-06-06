#' @title Extract key genes of a group of spots/cells
#' @description Run Structformer trained model to extract key genes of a group spots/cells. 
#' 
#' @param seu_obj The seurat object which will be used for predicting.
#' @param save_checkpoint_path The save path of Structformer trained model.
#' @param envir_path The python env path.
#' @param target_name The colname of which you want to predict.
#' @param phenotype_num The number represents the phenotype in the target_name column.
#' @param metadata Whether use metadata of Seurat object
#' @param top_n The number of top gene you want to get.
#' @param with_gpu Defulat FALSE, TRUE/FALSE, use GPU or not.
#'
#' @return Key genes of a group spots/cells
#' @export 

run_Structformer_extract_keygenes <- function(seu_obj,save_checkpoint_path,envir_path,target_name,phenotype_num,metadata = TRUE,top_n=60,with_gpu = FALSE){
  
  reticulate::use_condaenv(envir_path, required = TRUE)
  reticulate::source_python(system.file("python", "extract_kg.py", package = "Structformer"))

  if(metadata){ 
    doc_input <- seu_obj[seu_obj[,target_name]%in%phenotype_num,]
  }else{
    doc_input <- seu_obj@meta.data[seu_obj@meta.data[,target_name]%in%phenotype_num,]
  }
  
  
  keygenes_df <- get_keywords(pretrainedmodel_path = save_checkpoint_path,doc_input = doc_input,top_n = as.integer(top_n),use_gpu=with_gpu)
  
  # if((system.file(package = "SCpubr")!="")&(system.file(package = "clusterProfiler")!="")&(system.file(package = "org.Hs.eg.db")!="")){
  #   library(SCpubr)
  #   library(clusterProfiler)
  #   library(org.Hs.eg.db)
  #   suppressMessages({
  #     options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  #     options(enrichR.live = TRUE)
  #     options(modEnrichR.use = TRUE)
  #     options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  #     options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))
      
  #     # Set the search to Human genes.
  #     enrichR::setEnrichrSite(site = "Enrichr")
      
  #     websiteLive <- TRUE
  #     dbs <- enrichR::listEnrichrDbs()
  #     # Get all the possible databases to query.
  #     dbs <- sort(dbs$libraryName)
  #   })
  #   dbs_use <- c("GO_Biological_Process_2023")
  #   enriched_terms <- enrichR::enrichr(keygenes_df$Gene, dbs_use)
  #   SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms)
  #   CirclePlot <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms)
  #   HeatmapPlot_BP <- SCpubr::do_FunctionalAnnotationPlot(genes = keygenes_df$Gene,
  #                                                         org.db = org.Hs.eg.db,
  #                                                         #database = "KEGG",
  #                                                         p.adjust.cutoff = adjp_pathway)
  #   HeatmapPlot_KEGG <- SCpubr::do_FunctionalAnnotationPlot(genes = keygenes_df$Gene,
  #                                                           org.db = org.Hs.eg.db,
  #                                                           database = "KEGG",
  #                                                           p.adjust.cutoff = adjp_pathway)
  #   keygenes_df <- list(keygenes_df,CirclePlot,HeatmapPlot_BP,HeatmapPlot_KEGG)
  #   names(keygenes_df) <- c("keygenes","CirclePlot","HeatmapPlot_BP","HeatmapPlot_KEGG")
  # }
  
  return(keygenes_df)
  
}