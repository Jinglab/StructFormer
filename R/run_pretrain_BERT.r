#' @title Calculate the embedding for spots/cells
#' @description Run the Structformer trained model or a pre-trained BERT model to extract the embeddings of the spots/cells. 
#' 
#' @param dat The Seurat object or dataframe which contains sentenences. 
#' @param gene_list The complete set of genes used to construct the sentences.
#' @param block_size The input sentences have a maximum word count for a single sentence.
#' @param vocab_size The total number of genes available.
#' @param hidden_size The length of the output vector in the last layer of the pre-trained model.
#' @param num_hidden_layers The hidden layer numbers.
#' @param num_attention_heads The anttention heads in model.
#' @param max_position_embeddings The defined maximum length of a sentence.
#' @param mlm_probability The probability of masked words in a sentence.
#' @param num_train_epochs Training epochs.
#' @param per_device_train_batch_size The batch size. 
#' @param save_steps After save_steps to automately save model. 
#' @param save_total_limit The limit on steps that can be used for saving the model. If this number is exceeded, the model will not be automatically saved when running the save_steps.
#' @param save_checkpoint_path The save path for the Structformer trained model. 
#' @param envir_path The python env path.
#' @param with_gpu Default: FALSE. Specify TRUE or FALSE to indicate whether to use a GPU.
#'
#' @return Runing the pre-training process
#' @export 

run_Structformer_pretrainBERT <- function(dat,gene_list,block_size = 260,vocab_size=50000,hidden_size=768,num_hidden_layers=6,
                                       num_attention_heads=12,max_position_embeddings=512,mlm_probability=0.15,
                                       num_train_epochs=10000,per_device_train_batch_size=16,save_steps=10000,save_total_limit=60000,
                                       save_checkpoint_path,envir_path){
    reticulate::use_condaenv(envir_path, required = TRUE)
    reticulate::source_python(system.file("python", "pretraining_BERT.py", package = "Structformer"))
    
    start_time <- Sys.time()

    if(class(dat)[1]=="Seurat"){
        metadata <- dat@meta.data
        if(dir.exists(tempdir())){
            write.table(metadata$sentence,file=paste0(tempdir(),"/temp_sentences.txt"),quote = F,col.names = F,row.names = F,sep = "\t")
            pretrain_BERT(file_path = paste0(tempdir(),"/temp_sentences.txt"),
                          vocab_file = system.file("python", "vocab.txt", package = "Structformer"),
                          gene_list = gene_list,
                          block_size = as.integer(block_size),vocab_size=as.integer(vocab_size+6),hidden_size=as.integer(hidden_size),
                          num_hidden_layers=as.integer(num_hidden_layers),num_attention_heads=as.integer(num_attention_heads),
                          max_position_embeddings=as.integer(max_position_embeddings),mlm_probability=mlm_probability,
                          output_dir=save_checkpoint_path,
                          num_train_epochs=as.integer(num_train_epochs),per_device_train_batch_size=as.integer(per_device_train_batch_size),
                          save_steps=as.integer(save_steps),save_total_limit=as.integer(save_total_limit))
        }else{
            dir.create(paste0(save_checkpoint_path,"/sentences/"))
            Sys.setenv(TEMP = paste0(save_checkpoint_path,"/sentences/"), TMPDIR = paste0(save_checkpoint_path,"/sentences/"))
            write.table(metadata$sentence,file=paste0(save_checkpoint_path,"/sentences","/temp_sentences.txt"),quote = F,col.names = F,row.names = F,sep = "\t")
            pretrain_BERT(file_path = paste0(tempdir(),"/temp_sentences.txt"),
                          vocab_file = system.file("python", "vocab.txt", package = "Structformer"),
                          gene_list = gene_list,
                          block_size = as.integer(block_size),vocab_size=as.integer(vocab_size),hidden_size=as.integer(hidden_size),
                          num_hidden_layers=as.integer(num_hidden_layers),num_attention_heads=as.integer(num_attention_heads),
                          max_position_embeddings=as.integer(max_position_embeddings),mlm_probability=mlm_probability,
                          output_dir=save_checkpoint_path,
                          num_train_epochs=as.integer(num_train_epochs),per_device_train_batch_size=as.integer(per_device_train_batch_size),
                          save_steps=as.integer(save_steps),save_total_limit=as.integer(save_total_limit))
        }
    }else{
        if(dir.exists(tempdir())){
            write.table(dat$sentence,file=paste0(tempdir(),"/temp_sentences.txt"),quote = F,col.names = F,row.names = F,sep = "\t")
            pretrain_BERT(file_path = paste0(tempdir(),"/temp_sentences.txt"),
                          vocab_file = system.file("python", "vocab.txt", package = "Structformer"),
                          gene_list = gene_list,
                          block_size = as.integer(block_size),vocab_size=as.integer(vocab_size),hidden_size=as.integer(hidden_size),
                          num_hidden_layers=as.integer(num_hidden_layers),num_attention_heads=as.integer(num_attention_heads),
                          max_position_embeddings=as.integer(max_position_embeddings),mlm_probability=mlm_probability,
                          output_dir=save_checkpoint_path,
                          num_train_epochs=as.integer(num_train_epochs),per_device_train_batch_size=as.integer(per_device_train_batch_size),
                          save_steps=as.integer(save_steps),save_total_limit=as.integer(save_total_limit))
        }else{
            dir.create(paste0(save_checkpoint_path,"/sentences/"))
            Sys.setenv(TEMP = paste0(save_checkpoint_path,"/sentences/"), TMPDIR = paste0(save_checkpoint_path,"/sentences/"))
            write.table(dat$sentence,file=paste0(save_checkpoint_path,"/sentences","/temp_sentences.txt"),quote = F,col.names = F,row.names = F,sep = "\t")
            pretrain_BERT(file_path = paste0(tempdir(),"/temp_sentences.txt"),
                          vocab_file = system.file("python", "vocab.txt", package = "Structformer"),
                          gene_list = gene_list,
                          block_size = as.integer(block_size),vocab_size=as.integer(vocab_size),hidden_size=as.integer(hidden_size),
                          num_hidden_layers=as.integer(num_hidden_layers),num_attention_heads=as.integer(num_attention_heads),
                          max_position_embeddings=as.integer(max_position_embeddings),mlm_probability=mlm_probability,
                          output_dir=save_checkpoint_path,
                          num_train_epochs=as.integer(num_train_epochs),per_device_train_batch_size=as.integer(per_device_train_batch_size),
                          save_steps=as.integer(save_steps),save_total_limit=as.integer(save_total_limit))
        }
    }

    end_time <- Sys.time()

    cost_time <- end_time - start_time
    
    print(paste0("Cost time:"))
    print(cost_time)
}


