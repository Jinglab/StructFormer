import warnings
warnings.filterwarnings("ignore")
import os
import dill
import pandas as pd
import numpy as np
import torch
import tokenizers
from pathlib import Path
from tqdm.notebook import tqdm
from transformers import Trainer, TrainingArguments
from transformers import BertTokenizer, LineByLineTextDataset
from transformers import BertConfig, BertForMaskedLM, DataCollatorForLanguageModeling
from transformers import AutoModelForMaskedLM,AutoTokenizer

def pretrain_BERT(file_path,vocab_file,gene_list,block_size = 260,vocab_size=50000,hidden_size=768,num_hidden_layers=6,num_attention_heads=12,max_position_embeddings=512,mlm_probability=0.15,output_dir="./",num_train_epochs=1000,per_device_train_batch_size=16,save_steps=10000,save_total_limit=60000):
    
    tokenizer = BertTokenizer(vocab_file= vocab_file)
    tokenizer.add_tokens(gene_list)
    tokenizer.save_pretrained(output_dir)
    
    dataset= LineByLineTextDataset(
        tokenizer = tokenizer,
        file_path = file_path,
        block_size = block_size  # maximum sequence length
    )
    print('No. of lines: ', len(dataset)) # No of lines in your datset

    config = BertConfig(
        vocab_size=vocab_size,
        hidden_size=hidden_size, 
        num_hidden_layers=num_hidden_layers, 
        num_attention_heads=num_attention_heads,
        max_position_embeddings=max_position_embeddings
    )
    
    model = BertForMaskedLM(config)

    print('No of parameters: ', model.num_parameters())

    data_collator = DataCollatorForLanguageModeling(
        tokenizer=tokenizer, mlm=True, mlm_probability=mlm_probability
    )

    training_args = TrainingArguments(
        output_dir=output_dir,
        overwrite_output_dir=True,
        num_train_epochs=num_train_epochs,
        per_device_train_batch_size=per_device_train_batch_size,
        save_steps=save_steps,
        save_total_limit=save_total_limit
        # ,
        # use_cpu = with_gpu
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        data_collator=data_collator,
        train_dataset=dataset
    )

    trainer.train()

    trainer.save_model(output_dir)
    
    
