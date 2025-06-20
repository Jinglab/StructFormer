# StructFormer  </a>

# Cell-resolved analysis of spatially organized cellular architecture using transformer models

StructFormer is a computational biology method designed to trace back and identify cells within specific spatial structures in single-cell transcriptomics data by learning from spatial transcriptomics data. In this study, we focused primarily on tertiary lymphoid structures (TLSs) as the key spatial structure of interest. By leveraging StructFormer, we aim to enable precise identification and characterization of cells residing in TLS regions, providing deeper insights into their spatial organization and functional roles. Unlike previous methods, StructFormer leverages transformer model and meta-learning to acquire general knowledge from the limited available specific spatial structures information. This knowledge is then transferred to identify single cells or spots within specific spatial structures from scRNA-seq.

## Requirements and Installation
This toolkit is implemented in both R and Python programming languages. Therefore, users need firstly install both R and Python before using the tool.

### Installation of StructFormer

To use StructFormer, first use the [yaml file](https://github.com/Jinglab/StructFormer/blob/main/StructFormer_env_cuda12.4.yml) to create StructFormer Conda environment (The CUDA version must be between 12.0 and 12.6 (inclusive). If a version outside this range is used, users may need to perform manual installation according to the software versions specified in the documentation.).

    conda env create -f StructFormer_env.yml

After successfully creating the environment, the Python path of this environment can be found like this, which will be used in subsequent workflow steps.

    /home/xfan/miniconda3/envs/StructFormer_env/bin/python
    
Install StructFormer by devtools in R
[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/)

    devtools::install_github("Jinglab/StructFormer")
    
Alternatively, you can download the [StructFormer_1.0.tar.gz](https://github.com/Jinglab/StructFormer/blob/main/Structformer_1.0.tar.gz) package file from this GitHub repository and install it locally.

    install.packages("home/xfan/MLTLS_package/StructFormer_1.0.tar.gz")

## Quick Start

### Run StructFormer under a TLS knowledge transfer scenario
To use StructFormer, no complex preprocessing is needed; we only require the counts from a Seurat object as input.
Download the 10x Visium breast cancer pre-trained BERT and demo data from Google Cloud. The saved location of this pre-trained BERT will be utilized in the subsequent workflow steps.
- [breast cancer pre-trained gene word encoder](https://drive.google.com/drive/folders/1qLsl22T3IU2EEyXYM3z52_8MLNsFDyjO?usp=drive_link)
- [demo data](https://drive.google.com/drive/folders/1DZJ-f_RjpnRUszXNKm_KRGXpbHcwsEBK?usp=drive_link)


1.Load package and demo data

    library(Seurat)
    library(Structformer)
    library(reticulate)
    library(tidyverse)
    st_dat_train <- readRDS("home/xfan/MLTLS_package/demo_data/bc_st_demo_data.rds")
    st_dat_pred <- readRDS("home/xfan/MLTLS_package/demo_data/melanoma_st_demo_data.rds")

2.Set parameters

    sen_len = 260
    save_inseu = TRUE
    genes_representor = "home/xfan/MLTLS_package/demo_data/pretrained_models_rank260/genelist.txt"
    envir_path = "/home/xfan/miniconda3/envs/StructFormer_env/bin/python"
    pretrained_model = "StructFormer_BERT"
    pretrained_model_path = "home/xfan/MLTLS_package/demo_data/pretrained_models_rank260/"
    save_checkpoint_path = "home/xfan/MLTLS_package/demo_data/"
    batch_size = 1 # depend on your GPU memory to change
    train_K = 2 # depend on your GPU memory to change
    train_Q = 2 # depend on your GPU memory to change
    train_episodes = 600

3.Generate sentences
    
    # Training data
    st_dat_train <- generate_sentences(
      seu_obj = st_dat_train,
      sen_len = sen_len,
      region_info = st_dat_train@meta.data$region,
      save_inseu = save_inseu,
      genes_representor = genes_representor,
      envir_path = envir_path
    )
    
    # Prediction data
    st_dat_pred <- generate_sentences(
      seu_obj = st_dat_pred,
      sen_len = sen_len,
      region_info = st_dat_pred@meta.data$region,
      save_inseu = save_inseu,
      genes_representor = genes_representor,
      envir_path = envir_path
    )

4.Training StructFormer
    
    # Training
    st_dat_train <- run_structformer_train(
        seu_obj = st_dat_train,
        pretrained_model = pretrained_model,
        sen_len = sen_len,
        pretrained_model_path = pretrained_model_path,
        save_checkpoint_path = save_checkpoint_path,
        batch_size = batch_size,
        train_K = train_K,
        train_Q = train_Q,
        train_episodes = train_episodes,
        envir_path = envir_path
    )

5.Use trained StructFormer to predict

    # Prediction
    st_dat_pred <- run_structformer_pred(
                        seu_obj = st_dat_pred,
                        pretrained_model_path = pretrained_model_path,
                        save_checkpoint_path = save_checkpoint_path,
                        envir_path = envir_path,
                        pretrained_model = pretrained_model,
                        sen_len=sen_len)
    # Normalization -- 0-1 scale
    st_dat_pred$relative_distance <- 1- (st_dat_pred$relative_distance - min(st_dat_pred$relative_distance))/(max(st_dat_pred$relative_distance)-min(st_dat_pred$relative_distance))
    SpatialFeaturePlot(st_dat_pred,features = c("region","relative_distance"))

## Built With
  - [Python](https://www.python.org/)
  - [PyTorch](https://pytorch.org/)
  - [R](https://www.contributor-covenant.org/](https://www.r-project.org/about.html))

## Visitor Map

[![ClustrMaps](https://www.clustrmaps.com/map_v2.png?d=GofyOV-Ko5ZvOMmDh4yGs2s5SYgmT7ig2NjwcxfvYD4&cl=ffffff)](https://clustrmaps.com/site/1c6i8)
