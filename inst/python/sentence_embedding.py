### import packages
import warnings
warnings.filterwarnings("ignore")
import os
import pandas as pd
import numpy as np
import torch
from pathlib import Path
from transformers import BertTokenizer, BertModel
from tqdm import tqdm
#import matplotlib.pyplot as plt
#from sklearn.metrics.pairwise import cosine_similarity
import operator
import dill

def get_sentences_embedding(sentences, save_path_checkpoint,pretrained_model,hidden_out = True,with_gpu=False):

    tokenizer = BertTokenizer.from_pretrained(save_path_checkpoint, use_fast=True)

    # Load pre-trained model (weights)
    model_bert = BertModel.from_pretrained(save_path_checkpoint,
                                    output_hidden_states = hidden_out, # Whether the model returns all hidden-states.
                                    )

    all_dat = sentences
    
    if pretrained_model=="Structformer_BERT":
        if with_gpu==False:
            model_bert.to("cpu")
            # Put the model in "evaluation" mode, meaning feed-forward operation.
            model_bert.eval()
            sentence_vec = torch.tensor(np.zeros(shape = (len(all_dat),768))).to("cpu")
            for row_num in  tqdm(range(0,len(all_dat))):
                #print(row_num)
                # word_vec = np.zeros((1,300))
                # if word in w2v_model.wv:
                # word_detected = [i for i in word if i in w2v_model.wv]
                text = all_dat.sentence[row_num]
                marked_text = text
                # Split the sentence into tokens.
                tokenized_text = tokenizer.tokenize(marked_text)
                # Map the token strings to their vocabulary indeces.
                indexed_tokens = tokenizer.convert_tokens_to_ids(tokenized_text)
                # Mark each of the 22 tokens as belonging to sentence "1".
                segments_ids = [1] * len(tokenized_text)
                # Convert inputs to PyTorch tensors
                tokens_tensor = torch.tensor([indexed_tokens]).to("cpu")
                segments_tensors = torch.tensor([segments_ids]).to("cpu")
                # Run the text through BERT, and collect all of the hidden states produced
                # from all 12 layers. 
                with torch.no_grad():
                    outputs = model_bert(tokens_tensor, segments_tensors)
                    hidden_states = outputs[2]
                # `token_vecs` is a tensor with shape [126 x 768]
                token_vecs = hidden_states[-2][0]
                # Calculate the average of all 126 token vectors.
                sentence_embedding = torch.mean(token_vecs, dim=0)
                sentence_embedding = sentence_embedding.view(1,768)
                sentence_vec[row_num] = sentence_embedding
        else:
            model_bert.to("cuda")
            # Put the model in "evaluation" mode, meaning feed-forward operation.
            model_bert.eval()
            sentence_vec = torch.tensor(np.zeros(shape = (len(all_dat),768))).to("cuda")
            for row_num in tqdm(range(0,len(all_dat))):
                #print(row_num)
                # word_vec = np.zeros((1,300))
                # if word in w2v_model.wv:
                # word_detected = [i for i in word if i in w2v_model.wv]
                text = all_dat.sentence[row_num]
                marked_text = text
                # Split the sentence into tokens.
                tokenized_text = tokenizer.tokenize(marked_text)
                # Map the token strings to their vocabulary indeces.
                indexed_tokens = tokenizer.convert_tokens_to_ids(tokenized_text)
                # Mark each of the 22 tokens as belonging to sentence "1".
                segments_ids = [1] * len(tokenized_text)
                # Convert inputs to PyTorch tensors
                tokens_tensor = torch.tensor([indexed_tokens]).to("cuda")
                segments_tensors = torch.tensor([segments_ids]).to("cuda")
                # Run the text through BERT, and collect all of the hidden states produced
                # from all 12 layers. 
                with torch.no_grad():
                    outputs = model_bert(tokens_tensor, segments_tensors)
                    hidden_states = outputs[2]
                # `token_vecs` is a tensor with shape [126 x 768]
                token_vecs = hidden_states[-2][0]
                # Calculate the average of all 126 token vectors.
                sentence_embedding = torch.mean(token_vecs, dim=0)
                sentence_embedding = sentence_embedding.view(1,768)
                sentence_vec[row_num] = sentence_embedding
    elif pretrained_model=="Geneformer":
        if with_gpu==False:
            model_bert.to("cpu")
            # Put the model in "evaluation" mode, meaning feed-forward operation.
            model_bert.eval()
            sentence_vec = torch.tensor(np.zeros(shape = (len(all_dat),512))).to("cpu")
            for row_num in  tqdm(range(0,len(all_dat))):
                #print(row_num)
                # word_vec = np.zeros((1,300))
                # if word in w2v_model.wv:
                # word_detected = [i for i in word if i in w2v_model.wv]
                text = all_dat.sentence[row_num]
                marked_text = text
                # Split the sentence into tokens.
                tokenized_text = tokenizer.tokenize(marked_text)
                # Map the token strings to their vocabulary indeces.
                indexed_tokens = tokenizer.convert_tokens_to_ids(tokenized_text)
                # Mark each of the 22 tokens as belonging to sentence "1".
                segments_ids = [1] * len(tokenized_text)
                # Convert inputs to PyTorch tensors
                tokens_tensor = torch.tensor([indexed_tokens]).to("cpu")
                segments_tensors = torch.tensor([segments_ids]).to("cpu")
                # Run the text through BERT, and collect all of the hidden states produced
                # from all 12 layers. 
                with torch.no_grad():
                    outputs = model_bert(tokens_tensor, segments_tensors)
                    hidden_states = outputs[2]
                # `token_vecs` is a tensor with shape [126 x 768]
                token_vecs = hidden_states[-2][0]
                # Calculate the average of all 126 token vectors.
                sentence_embedding = torch.mean(token_vecs, dim=0)
                sentence_embedding = sentence_embedding.view(1,512)
                sentence_vec[row_num] = sentence_embedding
        else:
            model_bert.to("cuda")
            # Put the model in "evaluation" mode, meaning feed-forward operation.
            model_bert.eval()
            sentence_vec = torch.tensor(np.zeros(shape = (len(all_dat),512))).to("cuda")
            for row_num in tqdm(range(0,len(all_dat))):
                #print(row_num)
                # word_vec = np.zeros((1,300))
                # if word in w2v_model.wv:
                # word_detected = [i for i in word if i in w2v_model.wv]
                text = all_dat.sentence[row_num]
                marked_text = text
                # Split the sentence into tokens.
                tokenized_text = tokenizer.tokenize(marked_text)
                # Map the token strings to their vocabulary indeces.
                indexed_tokens = tokenizer.convert_tokens_to_ids(tokenized_text)
                # Mark each of the 22 tokens as belonging to sentence "1".
                segments_ids = [1] * len(tokenized_text)
                # Convert inputs to PyTorch tensors
                tokens_tensor = torch.tensor([indexed_tokens]).to("cuda")
                segments_tensors = torch.tensor([segments_ids]).to("cuda")
                # Run the text through BERT, and collect all of the hidden states produced
                # from all 12 layers. 
                with torch.no_grad():
                    outputs = model_bert(tokens_tensor, segments_tensors)
                    hidden_states = outputs[2]
                # `token_vecs` is a tensor with shape [126 x 768]
                token_vecs = hidden_states[-2][0]
                # Calculate the average of all 126 token vectors.
                sentence_embedding = torch.mean(token_vecs, dim=0)
                sentence_embedding = sentence_embedding.view(1,512)
                sentence_vec[row_num] = sentence_embedding
    sentence_vec = pd.DataFrame(np.array(sentence_vec.to("cpu")))
    sentence_vec.index = all_dat.index   
    return(sentence_vec)
