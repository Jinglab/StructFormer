#' @title Calculate the embedding for spots/cells
#' @description Run the Structformer trained model or a pre-trained BERT model to extract the embeddings of the spots/cells. 
#' 
#' @param seu_obj The Seurat object contains the data for spots/cells.
#' @param save_checkpoint_path The save path for the Structformer trained model.
#' @param pretrained_model The default is Structformer_BERT. Options include Structformer_BERT or Geneformer. Use this function when the save_checkpoint_path parameter is NULL. 
#' @param pretrained_model_path The save path for the pre-trained BERT model. Use this function when the save_checkpoint_path parameter is NULL.
#' @param envir_path The python env path.
#' @param with_gpu Default: FALSE. Specify TRUE or FALSE to indicate whether to use a GPU.
#'
#' @return The seurat object
#' @export 

run_Structformer_extract_cs_embedding <- function(seu_obj,save_checkpoint_path,pretrained_model,envir_path,with_gpu = FALSE){
    reticulate::use_condaenv(envir_path, required = TRUE)
    reticulate::source_python(system.file("python", "sentence_embedding.py", package = "Structformer"))
    #seu_obj@meta.data$save <- ifelse(rownames(seu_obj@meta.data)%in%rownames(seu_obj@meta.data[seu_obj@meta.data[,target_name]%in%phenotype_num,]),1,0)
    
    metadata <- seu_obj@meta.data

    sentences_embedding_vec <- get_sentences_embedding(sentences = metadata, save_path_checkpoint = save_checkpoint_path,pretrained_model = pretrained_model,with_gpu=with_gpu)

    sentences_embedding_vec <- t(sentences_embedding_vec)

    sentences_embedding_vec <- Seurat::CreateSeuratObject(
                    sentences_embedding_vec,
                    project = "Embedding",
                    assay="Embedding",
                    min.cells = 0,
                    min.features =0
                    )
    sentences_embedding_vec <- ScaleData(sentences_embedding_vec, features = rownames(sentences_embedding_vec))
    sentences_embedding_vec <- RunPCA(sentences_embedding_vec, features = rownames(sentences_embedding_vec))
    sentences_embedding_vec <- FindNeighbors(sentences_embedding_vec, dims = 1:10)
    sentences_embedding_vec <- FindClusters(sentences_embedding_vec, resolution = 0.7)
    sentences_embedding_vec <- RunUMAP(sentences_embedding_vec, dims = 1:10)
    sentences_embedding_vec@meta.data$barcode <- rownames(sentences_embedding_vec@meta.data)

    seu_obj@graphs <- sentences_embedding_vec@graphs
    seu_obj@neighbors <- sentences_embedding_vec@neighbors
    seu_obj@reductions <- sentences_embedding_vec@reductions
    seu_obj@meta.data$embedding_clusters <- sentences_embedding_vec@meta.data$seurat_clusters

    # results <- list(gene_modules,sentences_embedding_vec)

    # names(results) <- c("gene_modules","gene_embedding")
    
    return(seu_obj)
}