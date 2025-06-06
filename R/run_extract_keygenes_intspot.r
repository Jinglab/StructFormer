#' @title Extract key genes of each spot and intergrate genes by RRA
#' @description Run Structformer trained model to extract key genes of a each spots/cells. 
#' 
#' @param seu_obj The seurat object which will be used for predicting.
#' @param save_checkpoint_path The save path of Structformer trained model.
#' @param envir_path The python env path.
#' @param target_name The colname of which you want to predict
#' @param phenotype_num The number represents the phenotype in the target_name column
#' @param top_n The number of top gene you want to get
#' @param with_gpu Defulat FALSE, TRUE/FALSE, use GPU or not.
#'
#' @return The intergrated key genes of user input gene sentences of each spot/cell.
#' @export 

run_Structformer_extract_keygenes_rra <- function(seu_obj,save_checkpoint_path,envir_path,target_name,phenotype_num,top_n=30000,with_gpu = FALSE){
  
  reticulate::use_condaenv(envir_path, required = TRUE)
  reticulate::source_python(system.file("python", "extract_kg_foreach_point.py", package = "Structformer"))
  
  doc_input <- seu_obj@meta.data[seu_obj@meta.data[,target_name]%in%phenotype_num,]

  model_prepared <- model_init(save_checkpoint_path,use_gpu=with_gpu)

  # Initialise a progress bar
  # pb <- txtProgressBar(min = 1, max = n, style = 3)

  if((system.file(package = "RobustRankAggreg")!="")){
    library(RobustRankAggreg)
    for(n in 1:nrow(doc_input)){
        if(n == 1){
            words_rank_list <- list(rep(c(),nrow(doc_input)))
        }
        words_rank_tmp <- get_keywords(model_prepared,doc_input = doc_input$sentence[n],top_n = as.integer(top_n))
        words_rank_list[[n]] <- words_rank_tmp$Gene
        cat("\rFinished", n, "of", nrow(doc_input))
        # # Print progress
	      # setTxtProgressBar(pb, i)
    }
    #close(pb)
    glist <- words_rank_list
    freq=as.data.frame(table(unlist(glist)))
    # Aggregate the inputs
    ag=aggregateRanks(glist = glist)
    ag$Freq=freq[match(ag$Name,freq$Var1),2]
    key_gene_df <- as.data.frame(ag)
    colnames(key_gene_df) <- c("Gene","P value","Freq")
  }else{
    BiocManager::install("RobustRankAggreg",ask = F,update = F) 
    library(RobustRankAggreg)
    for(n in 1:nrow(doc_input)){
        if(n == 1){
            words_rank_list <- list(rep(c(),nrow(doc_input)))
        }
        words_rank_tmp <- get_keywords(model_prepared,doc_input = doc_input$sentence[n],top_n = as.integer(top_n))
        words_rank_list[[n]] <- words_rank_tmp$Gene
        cat("\rFinished", n, "of", nrow(doc_input))
        # # Print progress
        # setTxtProgressBar(pb, i)
    }
    #close(pb)
    glist <- words_rank_list
    freq=as.data.frame(table(unlist(glist)))
    # Aggregate the inputs
    ag=aggregateRanks(glist = glist)
    ag$Freq=freq[match(ag$Name,freq$Var1),2]
    key_gene_df <- as.data.frame(ag)
    colnames(key_gene_df) <- c("Gene","P value","Freq")
  }

  return(key_gene_df)
  
}
