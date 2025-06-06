### import packages
import warnings
warnings.filterwarnings("ignore")
import time
import pandas as pd
from keybert import KeyBERT
from sentence_transformers import SentenceTransformer

def get_keywords(pretrainedmodel_path,doc_input,top_n,use_gpu=False):
    print("Running")
    start_time = time.time()
    if use_gpu==False:
        sentence_model = SentenceTransformer(pretrainedmodel_path,device = "cpu")
    else:
        sentence_model = SentenceTransformer(pretrainedmodel_path,device = "cuda")
    kw_model = KeyBERT(model=sentence_model)
    doc = doc_input['sentence'].str.cat(sep=' ') 
    keygenes = kw_model.extract_keywords(doc,top_n=top_n)
    keygenes = pd.DataFrame(keygenes,columns=["Gene","Score"])
    keygenes["Gene"] = keygenes["Gene"].str.upper()
    end_time = time.time()
    print("Cost time: {:.2f}s".format(end_time - start_time))
    if len(keygenes)>=10:
        print("The top 10 genes are {}".format(" ".join(list(keygenes["Gene"][range(0,10)]))))
    return keygenes




