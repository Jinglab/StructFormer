#' @title Gene correlation analysis and network construction
#' @description Run Structformer trained model to extract key genes of a group spots/cells. 
#' 
#' @param gene_list The seurat object which will be used for predicting.
#' @param save_checkpoint_path The save path of Structformer trained model.
#' @param res The resoultion to distinguish clusters.
#' @param link_thres The threshold of cosine similarity score of two genes' linking.
#' @param col_list The color list for gene modules.
#' @param out_dist The distance of label between circle.
#' @param envir_path The python env path.
#' @param with_gpu Defulat FALSE, TRUE/FALSE, use GPU or not.
#'
#' @return The relative distance of predicted single cells or spots with TLS prototype and non-TLS prototype, the prediction label of whether a single cell or spot belong to TLS region.
#' @export 

run_Structformer_gene_network <- function(gene_list,save_checkpoint_path, res = 0.7, link_thres = 0.6, col_list = NULL,out_dist=2.5,envir_path = NULL, with_gpu = FALSE){
  
  discrete_palette_tmp <- function (name, n = 1) {
    palettes <- list(stallion = c("#D51F26", "#272E6A", "#208A42", 
                                  "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", 
                                  "#E6C2DC", "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", 
                                  "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", 
                                  "#D8A767", "#3D3D3D"), calm = c("#7DD06F", "#844081", 
                                                                  "#688EC1", "#C17E73", "#484125", "#6CD3A7", "#597873", 
                                                                  "#7B6FD0", "#CF4A31", "#D0CD47", "#722A2D", "#CBC594", 
                                                                  "#D19EC4", "#5A7E36", "#D4477D", "#403552", "#76D73C", 
                                                                  "#96CED5", "#CE54D1", "#C48736"), kelly = c("#FFB300", 
                                                                                                              "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", 
                                                                                                              "#817066", "#007D34", "#F6768E", "#00538A", "#FF7A5C", 
                                                                                                              "#53377A", "#FF8E00", "#B32851", "#F4C800", "#7F180D", 
                                                                                                              "#93AA00", "#593315", "#F13A13", "#232C16"), bear = c("#faa818", 
                                                                                                                                                                    "#41a30d", "#fbdf72", "#367d7d", "#d33502", "#6ebcbc", 
                                                                                                                                                                    "#37526d", "#916848", "#f5b390", "#342739", "#bed678", 
                                                                                                                                                                    "#a6d9ee", "#0d74b6", "#60824f", "#725ca5", "#e0598b"), 
                     ironMan = c("#371377", "#7700FF", "#9E0142", "#FF0080", 
                                 "#DC494C", "#F88D51", "#FAD510", "#FFFF5F", "#88CFA4", 
                                 "#238B45", "#02401B", "#0AD7D3", "#046C9A", "#A2A475", 
                                 "grey35"), circus = c("#D52126", "#88CCEE", "#FEE52C", 
                                                       "#117733", "#CC61B0", "#99C945", "#2F8AC4", "#332288", 
                                                       "#E68316", "#661101", "#F97B72", "#DDCC77", "#11A579", 
                                                       "#89288F", "#E73F74"), paired = c("#A6CDE2", "#1E78B4", 
                                                                                         "#74C476", "#34A047", "#F59899", "#E11E26", "#FCBF6E", 
                                                                                         "#F47E1F", "#CAB2D6", "#6A3E98", "#FAF39B", "#B15928"), 
                     grove = c("#1a1334", "#01545a", "#017351", "#03c383", 
                               "#aad962", "#fbbf45", "#ef6a32", "#ed0345", "#a12a5e", 
                               "#710162", "#3B9AB2"), tableau = c("#4E79A7", "#E15759", 
                                                                  "#F28E2B", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", 
                                                                  "#FF9DA7", "#9C755F", "#BAB0AC"), summerNight = c("#2a7185", 
                                                                                                                    "#a64027", "#fbdf72", "#60824f", "#9cdff0", "#022336", 
                                                                                                                    "#725ca5"), captain = c("grey", "#A1CDE1", "#12477C", 
                                                                                                                                            "#EC9274", "#67001E"))
    if (!(name %in% names(palettes))) {
      stop(sprintf("palette name must be one of: %s", toString(names(palettes))))
    }
    palette <- palettes[[name]]
    if (n > length(palette)) {
      palette <- (grDevices::colorRampPalette(palette, space = "Lab"))(n)
    }
    return(palette)
  }
  
  tmp_results <- run_Structformer_extract_gene_embedding(
    gene_list,
    save_checkpoint_path,
    res = res,
    envir_path = envir_path,
    with_gpu = FALSE
  )
  
  gene_embedding_dat <- as.matrix(tmp_results$gene_embedding@assays$Spatial@counts)

  if(require(lsa)){
    "load lsa"
    library(lsa)
  }else{
    install.packages("lsa")
    library(lsa)
  }
  if(require(reshape2)){
    "load reshape2"
    library(reshape2)
  }else{
    install.packages("reshape2")
    library(reshape2)
  }
  if(require(igraph)){
    "load igraph"
    library(igraph)
  }else{
    install.packages("igraph")
    library(igraph)
  }
  
  cosine_sim_net <- cosine(gene_embedding_dat) 
  
  for(i in 1:nrow(cosine_sim_net)){
    ncol_turn <- ncol(cosine_sim_net)-i+1
    cosine_sim_net[i,((ncol(cosine_sim_net)-ncol_turn+1)):ncol_turn] <- 0
  }
  
  tmp_links <- melt(cosine_sim_net)
  colnames(tmp_links) <- c("source", "target", "importance")
  cos_tmp_links <- tmp_links[tmp_links$importance!=0,]
  tmp_links$new_name <- paste0(tmp_links$source,"_",tmp_links$target)
  tmp_links <- tmp_links[tmp_links$importance>=link_thres,]
  tmp_nodes <- tmp_results$gene_modules
  colnames(tmp_nodes) <- c("name","carac")
  tmp_nodes$carac <- paste0("Module ",tmp_nodes$carac)
  tmp_nodes <- tmp_nodes[(tmp_nodes$name)%in%c(unique(c(as.vector(tmp_links$source),as.vector(tmp_links$target)))),]
  
  # Turn it into igraph object
  network <- graph_from_data_frame(d=tmp_links, vertices=tmp_nodes, directed=F) 
  # color
  if(is.null(col_list)){
    coul  <- discrete_palette_tmp("stallion") 
  }else{
    coul = col_list
  }
  #layout <- layout_with_kk(network)
  # Create a vector of color
  my_color <- coul[as.numeric(as.factor(V(network)$carac))]

  plot(network, 
       #vertex.color = as.factor(V(g)$carac, 
       vertex.color=my_color,
       vertex.label.color = 'black', 
       vertex.label.family = "Arial",
       #layout = layout, 
       #layout = NULL, 
       vertex.label.dist = out_dist, 
       edge.width=E(network)$importance*2, 
       #edge.width = ifelse(E(network)$is_formal, 5, 1),
       edge.curved=.1,main = "The gene-gene interaction")
  
  legend("bottomright", 
         legend=levels(as.factor(V(network)$carac)), 
         col = coul , 
         bty = "n", pch=20 , pt.cex = 2, cex = 1, text.col=coul , horiz = FALSE, inset = c(0, 0))

  
  return(cos_tmp_links)
  
}