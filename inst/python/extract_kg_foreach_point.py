### import packages
import warnings
warnings.filterwarnings("ignore")
import time
import pandas as pd
from keybert import KeyBERT
from sentence_transformers import SentenceTransformer

def model_init(pretrainedmodel_path,use_gpu=False):
    print("Running")
    start_time = time.time()
    if use_gpu==False:
        sentence_model = SentenceTransformer(pretrainedmodel_path,device = "cpu")
    else:
        sentence_model = SentenceTransformer(pretrainedmodel_path,device = "cuda")
    kw_model = KeyBERT(model=sentence_model)
    return kw_model

def get_keywords(kw_model,doc_input,top_n):
    doc = doc_input 
    keygenes = kw_model.extract_keywords(doc,top_n=top_n)
    keygenes = pd.DataFrame(keygenes,columns=["Gene","Score"])
    keygenes["Gene"] = keygenes["Gene"].str.upper()
    return keygenes
