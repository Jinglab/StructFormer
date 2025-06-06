#' @title Predict single cell or spots whether belong to TLS region
#' @description Run Structformer trained model to calculate relative distance of single cells or spot with reference TLSs or Non-TLSs spots and infer whether ir belong to TLS region. 
#' 
#' @param seu_obj The seurat object which will be used for predicting.
#' @param pretrained_model_path The pre-trained gene word encoder model saved path.
#' @param save_checkpoint_path The save path of Structformer trained model.
#' @param envir_path The python env path.
#' @param pretrained_model Default is TLSPredictor_BERT. Structformer_BERT or geneformer.
#' @param sen_len Default is 260. The sentence length, the generated sentences length will be minus or equal this parameter. If the gene expression level is zero, the gene will not be invovled.
#' @param data_type Default is sc_st, the other is bulk.
#' @param metadata Use metadata of Seurat object, default is FALSE.
#' @param class_num The classification number, default is 2, it can support multiple classification.
#' @param with_gpu, TRUE/FALSE, use GPU or not.
#'
#' @return The relative distance of predicted single cells or spots with TLS prototype and non-TLS prototype, the prediction label of whether a single cell or spot belong to TLS region.
#' @export 

run_Structformer_pred <- function(seu_obj,pretrained_model_path,save_checkpoint_path,envir_path,
                                  pretrained_model = "Structformer_BERT", sen_len=260,data_type = "sc_st", metadata=FALSE, class_num = 2,with_gpu = TRUE){
    reticulate::use_condaenv(envir_path, required = TRUE)
    reticulate::source_python(system.file("python", "pred_singlecells_spots.py", package = "Structformer"))

    if(data_type=="sc_st"){
        if(metadata){
            dat_pred <- pred_tls(dat_pred = seu_obj,
                                pretrained_model = pretrained_model,
                                sen_len = as.integer(sen_len),
                                pretrained_model_path = pretrained_model_path,
                                save_checkpoint_path = save_checkpoint_path,
                                with_gpu = with_gpu)
            rownames(dat_pred) <- dat_pred$cell_barcode
            if(class_num == 2){
                dat_pred <- dat_pred[,c("Non-TLS Distance","TLS Distance","Non-TLS Cosine similarity","TLS Cosine similarity", "relative_distance","pred_label")]
                colnames(dat_pred) <- c("Non-Target Distance","Target Distance","Non-Target Cosine similarity","Target Cosine similarity", "relative_distance","pred_label")
            }else{
                dat_pred <- dat_pred[,(ncol(seu_obj)+1):ncol(dat_pred)]
            }
            seu_obj <- cbind(seu_obj,dat_pred)
        }else{
            dat_pred <- pred_tls(dat_pred = seu_obj@meta.data,
                                pretrained_model = pretrained_model,
                                sen_len = as.integer(sen_len),
                                pretrained_model_path = pretrained_model_path,
                                save_checkpoint_path = save_checkpoint_path,
                                with_gpu = with_gpu)
            rownames(dat_pred) <- dat_pred$cell_barcode
            if(class_num == 2){
                dat_pred <- dat_pred[,c("Non-TLS Distance","TLS Distance","Non-TLS Cosine similarity","TLS Cosine similarity", "relative_distance","pred_label")]
                colnames(dat_pred) <- c("Non-Target Distance","Target Distance","Non-Target Cosine similarity","Target Cosine similarity", "relative_distance","pred_label")
            }else{
                dat_pred <- dat_pred[,(ncol(seu_obj@meta.data)+1):ncol(dat_pred)]
            }
            seu_obj@meta.data <- cbind(seu_obj@meta.data,dat_pred)
       }
    }else if(data_type=="bulk"){
        seu_obj <- as.data.frame(data.frame(cell_barcode=unlist(colnames(seu_obj)),sentence = unlist(seu_obj["sentence",])))
        dat_pred <- pred_tls(dat_pred = seu_obj,
                                pretrained_model = pretrained_model,
                                sen_len = as.integer(sen_len),
                                pretrained_model_path = pretrained_model_path,
                                save_checkpoint_path = save_checkpoint_path,
                                with_gpu = with_gpu)
        rownames(dat_pred) <- dat_pred$cell_barcode
        if(class_num == 2){
            dat_pred <- dat_pred[,c("Non-TLS Distance","TLS Distance","Non-TLS Cosine similarity","TLS Cosine similarity", "relative_distance","pred_label")]
            colnames(dat_pred) <- c("Non-Target Distance","Target Distance","Non-Target Cosine similarity","Target Cosine similarity", "relative_distance","pred_label")
        }else{
            dat_pred <- dat_pred[,(ncol(dat_pred)-(class_num+1)):ncol(dat_pred)]
        }
        seu_obj <- dat_pred
    }
    return(seu_obj)
}