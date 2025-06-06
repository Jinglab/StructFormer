#' @title Train Structformer on spatial transcriptomics data
#' @description Train the Structformer that is based on metrics learing few-shot learning to deal with the unblance problem. 
#' Beacuse The number of TLSs spots for training is generally much smaller that non-TLSs spots.    
#' 
#' @param seu_obj, The seurat object which will be used for training.
#' @param pretrained_model, Default is Structformer_BERT. Structformer_BERT or geneformer.
#' @param sen_len, Default is 260. The sentence length, the generated sentences length will be minus or equal this parameter. If the gene expression level is zero, the gene will not be invovled.
#' @param pretrained_model_path, The pre-trained model saved path.
#' @param save_checkpoint_path, The save path of model.
#' @param batch_size, Batch size of training process.
#' @param train_K, Support set numbers.
#' @param train_Q, Query set numbers.
#' @param train_episodes, Training episodes.
#' @param val_episodes, Validation is performed after the set number of episodes. The default is 30.
#' @param val_steps, After val_steps, the accuracy of the current model is compared to the previous model. If the accuracy score improves, the current model will be saved..
#' @param metadata, Use metadata of Seurat object.
#' @param reproduce, Whether reproduce training process.
#' @param set_seed, Set a seed.
#' @param target_name, The colname of which you want to predict.
#' @param envir_path, The python env path.
#' @param with_gpu, TRUE/FALSE, use GPU or not.
#'
#' @return The trained Structformer model and support vector.
#' @export 

run_Structformer_train <- function(seu_obj,pretrained_model = "Structformer_BERT",
                               sen_len=260,pretrained_model_path,save_checkpoint_path,batch_size,
                               train_K,train_Q,train_episodes,val_episodes = 30,val_steps = 10,metadata = FALSE,reproduce = FALSE,set_seed = 42,target_name,envir_path,with_gpu = TRUE){
    reticulate::use_condaenv(envir_path, required = TRUE)
    reticulate::source_python(system.file("python", "train_prototype_vector.py", package = "Structformer"))
    
    if(metadata){
        input_dat <- seu_obj
    }else{
        input_dat <- seu_obj@meta.data
    }
    colnames(input_dat)[colnames(input_dat) %in% target_name] <- "region"
    input_dat$region <- as.integer(input_dat$region)
    dat_train <- train_prototype_hypersphere(dat_train = input_dat,
                                pretrained_model = pretrained_model,
                                sen_len = as.integer(sen_len),
                                pretrained_model_path = pretrained_model_path,
                                save_checkpoint_path = save_checkpoint_path,
                                batch_size = as.integer(batch_size),
                                train_K = as.integer(train_K),
                                train_Q = as.integer(train_Q), 
                                train_episodes = as.integer(train_episodes),
                                val_episodes = as.integer(val_episodes),
                                val_steps = as.integer(val_steps),
                                reproduce = reproduce,
                                set_seed = as.integer(set_seed),
                                with_gpu = with_gpu
                                )
    rownames(dat_train) <- dat_train$cell_barcode
    if(length(unique(input_dat$region))==2){
        dat_train <- dat_train[,c("Non-TLS Distance","TLS Distance","Non-TLS Cosine similarity","TLS Cosine similarity", "relative_distance","pred_label")]
        colnames(dat_train) <- c("Non-Target Distance","Target Distance","Non-Target Cosine similarity","Target Cosine similarity", "relative_distance","pred_label")
    }else{
        dat_train <- dat_train[,(ncol(input_dat)+1):ncol(dat_train)]
    }
    if(metadata){
        seu_obj <- cbind(seu_obj,dat_train)
    }else{
        seu_obj@meta.data <- cbind(seu_obj@meta.data,dat_train)
    }

    return(seu_obj)
}
