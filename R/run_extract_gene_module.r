#' @title Extract gene modules
#' @description Run Structformer trained model to extract key genes of a group spots/cells. 
#' 
#' @param seu_obj The seurat object which will be used for predicting.
#' @param save_checkpoint_path The save path of Structformer trained model.
#' @param envir_path The python env path.
#' @param target_name The colname of which you want to predict
#' @param phenotype_num The number represents the phenotype in the target_name column
#' @param hvg_num The number of top gene you want to get
#' @param genes_representor The file path of genes were used in the pre-training the gene word encoder.
#' @param res The resoultion to distinguish clusters.
#' @param with_gpu Defulat FALSE, TRUE/FALSE, use GPU or not.
#'
#' @return The relative distance of predicted single cells or spots with TLS prototype and non-TLS prototype, the prediction label of whether a single cell or spot belong to TLS region.
#' @export 

run_Structformer_extract_gene_module <- function(seu_obj,save_checkpoint_path,envir_path,target_name,phenotype_num,hvg_num=3000,res,genes_representor,with_gpu = FALSE){
    reticulate::use_condaenv(envir_path, required = TRUE)
    reticulate::source_python(system.file("python", "gene_embedding.py", package = "Structformer"))
    seu_obj@meta.data$save <- ifelse(rownames(seu_obj@meta.data)%in%rownames(seu_obj@meta.data[seu_obj@meta.data[,target_name]%in%phenotype_num,]),1,0)
    doc_input <- subset(seu_obj,subset = save == 1)
    doc_input <- FindVariableFeatures(doc_input, selection.method = "vst", nfeatures = hvg_num,assay = names(doc_input)[1])

    hvg_doc_input <- VariableFeatures(doc_input,assay = names(doc_input)[1])

    hvg_doc_input <- as.data.frame(hvg_doc_input)

    colnames(hvg_doc_input)[1] <- "Genes"
    hvg_doc_input$id <- rownames(hvg_doc_input)

    gene_embedding_vec <- get_gene_embedding(genes = hvg_doc_input, save_path_checkpoint = save_checkpoint_path,with_gpu=with_gpu)

    gene_embedding_vec <- t(gene_embedding_vec)

    gene_embedding_vec <- Seurat::CreateSeuratObject(
                        gene_embedding_vec,
                        project = "Spatial",
                        assay="Spatial",
                        min.cells = 0,
                        min.features =0
                        )
    gene_embedding_vec <- ScaleData(gene_embedding_vec, features = rownames(gene_embedding_vec))
    gene_embedding_vec <- RunPCA(gene_embedding_vec, features = rownames(gene_embedding_vec))
    gene_embedding_vec <- FindNeighbors(gene_embedding_vec, dims = 1:10)
    gene_embedding_vec <- FindClusters(gene_embedding_vec, resolution = res)
    gene_embedding_vec <- RunUMAP(gene_embedding_vec, dims = 1:10)
    #gene_embedding_vec <- RunTSNE(gene_embedding_vec, dims = 1:10)
    gene_embedding_vec@meta.data$Genes <- rownames(gene_embedding_vec@meta.data)
    gene_modules <- gene_embedding_vec@meta.data[,c("Genes","seurat_clusters")]
    colnames(gene_modules) <- c("Genes","Clusters")
    results <- list(gene_modules,gene_embedding_vec)
    names(results) <- c("gene_modules","gene_embedding")
    return(results)
}