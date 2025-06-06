#' @title Extract gene embedding
#' @description Run Structformer trained model to extract key genes of a group spots/cells. 
#' 
#' @param gene_list The gene list input.
#' @param save_checkpoint_path The save path of Structformer trained model.
#' @param pretrained_model_path, The pre-trained model saved path.
#' @param res The resoultion to distinguish clusters.
#' @param envir_path The python env path.
#' @param with_gpu Defulat FALSE, TRUE/FALSE, use GPU or not.
#'
#' @return The relative distance of predicted single cells or spots with TLS prototype and non-TLS prototype, the prediction label of whether a single cell or spot belong to TLS region.
#' @export 

run_Structformer_extract_gene_embedding <- function(gene_list,save_checkpoint_path,pretrained_model_path = NULL, res = 0.7, envir_path, with_gpu = FALSE){
  reticulate::use_condaenv(envir_path, required = TRUE)
  reticulate::source_python(system.file("python", "gene_embedding.py", package = "Structformer"))
  # seu_obj@meta.data$save <- ifelse(rownames(seu_obj@meta.data)%in%rownames(seu_obj@meta.data[seu_obj@meta.data[,target_name]%in%phenotype_num,]),1,0)
  # doc_input <- subset(seu_obj,subset = save == 1)
  # doc_input <- FindVariableFeatures(doc_input, selection.method = "vst", nfeatures = hvg_num,assay = names(doc_input)[1])
  # 
  # hvg_doc_input <- VariableFeatures(doc_input,assay = names(doc_input)[1])
  # 
  # hvg_doc_input <- as.data.frame(hvg_doc_input)
  # 
  # colnames(hvg_doc_input)[1] <- "Genes"
  #hvg_doc_input$id <- rownames(hvg_doc_input)
  
  gene_list <- as.data.frame(gene_list)
  colnames(gene_list) <- "Genes"
  rownames(gene_list) <- gene_list$Genes
  
  gene_embedding_vec <- get_gene_embedding(genes = gene_list, save_path_checkpoint = save_checkpoint_path,with_gpu=with_gpu)
  
  if(is.null(pretrained_model_path)){
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
  }else{
    if(system.file(package = "umap")!=""){
      library(umap)
    }else{
      install.packages("umap")
    }
    
    # pre-current gene embedding
    gene_embedding_raw_vec <- get_gene_embedding(genes = gene_list, save_path_checkpoint = pretrained_model_path,with_gpu=with_gpu)
    pre_gene_embedding <- as.data.frame(t(gene_embedding_raw_vec))
    colnames(pre_gene_embedding) <- paste0(colnames(pre_gene_embedding),"_pre")
    
    current_gene_embedding <- as.data.frame(t(gene_embedding_vec))
    colnames(current_gene_embedding) <- paste0(colnames(current_gene_embedding),"_current")
    
    gene_embedding_com <- cbind(pre_gene_embedding,current_gene_embedding)
    
    cosine_sim_net_pre_current <- cosine(as.matrix(gene_embedding_com)) 
    for(i in 1:nrow(cosine_sim_net_pre_current)){
      ncol_turn <- ncol(cosine_sim_net_pre_current)-i+1
      cosine_sim_net_pre_current[i,((ncol(cosine_sim_net_pre_current)-ncol_turn+1)):ncol_turn] <- 0
    }
    cosine_sim_net_pre_current <- melt(cosine_sim_net_pre_current)
    cosine_sim_net_pre_current <- cosine_sim_net_pre_current[!(cosine_sim_net_pre_current$Var1==cosine_sim_net_pre_current$Var2),]
    cosine_sim_net_pre_current <- cosine_sim_net_pre_current[str_split(cosine_sim_net_pre_current$Var1,"_",simplify = TRUE)[,1]==str_split(cosine_sim_net_pre_current$Var2,"_",simplify = TRUE)[,1],]
    cosine_sim_net_pre_current <- cosine_sim_net_pre_current[cosine_sim_net_pre_current$value>0,]
    
    cosine_sim_net_pre_current$Gene <- str_split(cosine_sim_net_pre_current$Var2,"_",simplify = TRUE)[,1]
    colnames(cosine_sim_net_pre_current)[3] <- "cosine_similarity"
    cosine_sim_net_pre_current <- arrange(cosine_sim_net_pre_current,cosine_similarity)
    cosine_sim_net_pre_current$Gene=factor(cosine_sim_net_pre_current$Gene,levels=cosine_sim_net_pre_current$Gene,order=T)
    
    # umap
    umap_gene_embedding <- umap::umap(gene_embedding_vec)
    umap_gene_embedding <- as.data.frame(umap_gene_embedding$layout)
    
    umap_gene_pre_embedding <- umap(gene_embedding_raw_vec)
    umap_gene_pre_embedding <- as.data.frame(umap_gene_pre_embedding$layout)
    
    umap_gene_embedding$gene <- rownames(umap_gene_embedding)
    umap_gene_embedding$state <- "current"
    umap_gene_pre_embedding$gene <- rownames(umap_gene_pre_embedding)
    umap_gene_pre_embedding$state <- "pre"
    
    umap_gene_embedding_com <- rbind(umap_gene_embedding,umap_gene_pre_embedding)
    
    colnames(umap_gene_embedding_com)[1:2] <- c("umap_1","umap_2")
    
    umap_gene_embedding_com <- as.data.frame(umap_gene_embedding_com)
    
    colnames(umap_gene_embedding) <- c("current_umap_1","current_umap_2","gene")
    colnames(umap_gene_pre_embedding) <- c("pre_umap_1","pre_umap_2","gene")
    umap_gene_pre_embedding <- umap_gene_pre_embedding[rownames(umap_gene_embedding),]
    umap_gene_embedding_line <- cbind(umap_gene_pre_embedding[,c(1:2)],umap_gene_embedding[,c(1:2)])
    umap_gene_embedding_line$gene <- rownames(umap_gene_embedding_line)
    
    umap_gene_embedding_line_selected_top <- umap_gene_embedding_line[umap_gene_embedding_line$gene%in%(c(cosine_sim_net_pre_current[order(cosine_sim_net_pre_current$cosine_similarity,decreasing = F),]$Gene[1:3])),]
    
    # plot umap gene embedding
    PointPlot <-ggplot(data = umap_gene_embedding_com,aes(x = umap_1,y = umap_2)) +
      geom_point(aes(color=state))+
      geom_segment(
        data = umap_gene_embedding_line_selected_top,
        aes(x = pre_umap_1,y = pre_umap_2, xend = current_umap_1,yend = current_umap_2),
        arrow = arrow(length = unit(0.01,units = "npc")),
        size = 0.7
      )+
      scale_color_manual(values = c('#8DA0CB','#66C2A5'))+
      geom_label( 
        data = umap_gene_embedding_line_selected_top,
        aes(x = pre_umap_1,y = pre_umap_2,label=gene)
      )+
      theme_test()

    # plot top cosine similarity gene
    CosineBarplot <- ggplot(cosine_sim_net_pre_current[1:10,],aes(x = Gene, y = cosine_similarity)) +
      geom_segment(aes(x = Gene, xend = Gene, y = (min(cosine_sim_net_pre_current$cosine_similarity)-0.05), yend = cosine_similarity),
                   linetype = "solid",size = 1,color = "gray40") +
      geom_point(
        color = "#ff630d",
        alpha = 0.8,
        size=3)+
      scale_y_continuous(limits = c(min(cosine_sim_net_pre_current$cosine_similarity)-0.05, 1.0)) +
      labs(y = "Cosine similarity score") +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(size = 7, angle = 45, hjust=1, face = 'bold'),
        axis.text.y = element_text(size = 7, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = 'bold'),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
      )
    
    results <- list(cosine_sim_net_pre_current[,c(4,3)],PointPlot,CosineBarplot)
    names(results) <- c("GeneState","PointPlot","CosineBarplot")
    
  }
  return(results)
}