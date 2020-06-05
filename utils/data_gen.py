from __future__ import absolute_import, division
import torch
from torch.functional import F
from torch_geometric.data import Dataset
from torch_geometric.utils import from_networkx, degree
import pandas as pd
import numpy as np
import pickle
import networkx as nx

from sklearn.preprocessing import MaxAbsScaler, MinMaxScaler, Normalizer, RobustScaler, StandardScaler
from sklearn.model_selection import train_test_split

def to_categorical(y, num_classes):
    """ 1-hot encodes a tensor """
    return torch.eye(num_classes, dtype=torch.long)[y]

def wcsv2graph(fname, prot_dict, normalize=False):
    """
    Weighted Graph Creator
    """
    sample = pd.read_csv('../../snac_data/' + fname)
    
    G = nx.from_pandas_edgelist(sample, source='node1', target='node2', 
                            edge_attr=['sign','weight'], create_using=nx.DiGraph())
    n1a1 = sample[['node1','upact1','downact1']]
    n1a1.columns = ['node','upact','downact']

    n2a2 = sample[['node2','upact2','downact2']]
    n2a2.columns = ['node','upact','downact']
    na = pd.concat([n1a1,n2a2])
    na = na.drop_duplicates('node')
    na = na.set_index('node')
    na['acts'] = na[['upact','downact']].apply(lambda x: np.hstack(x), axis=1)
    na = na.drop(['upact', 'downact'], axis=1)['acts'].to_dict()

    nx.set_node_attributes(G, prot_dict, 'x')
    nx.set_node_attributes(G, na, 'acts')

    data = from_networkx(G)
    data.weight = data.weight.float()
    data.sign = data.sign.float()
    data.x = torch.cat((data.x, data.acts.double()), axis=1)
    data.edge_attr = data.weight
    
    data.weight = data.acts = None
    #Normalize data
    if normalize:
        data.x = F.normalize(data.x)
    
    return data

def ucsv2graph(fname, global_dict, normalize=False):
    """
    Unweighted Graph Creator
    """
    sample = pd.read_csv('../snac_data/' + fname)
    
    G = nx.from_pandas_edgelist(sample, source='node1', target='node2', 
                            edge_attr=['sign'], create_using=nx.DiGraph())

    n1a1 = sample[['node1','activity1']]
    n1a1.columns = ['node','act']

    n2a2 = sample[['node2','activity2']]
    n2a2.columns = ['node','act']
    na = pd.concat([n1a1,n2a2])
    na = na.drop_duplicates('node')
    na = na.set_index('node')
    na['acts'] = na[['act']].apply(lambda x: np.hstack(x), axis=1)
    na = na.drop(['act'], axis=1)['acts'].to_dict()
    
    nx.set_node_attributes(G, global_dict,'global_idx')
    nx.set_node_attributes(G, na, 'acts')
    
    data = from_networkx(G)
    
    deg = degree(data.edge_index.reshape(-1).sort()[0]).long()
    data.acts[data.acts < 0] = 0
    data.acts = to_categorical(data.acts, 2).reshape(-1,2).long()
    data.degree = deg.reshape(-1,1).long()
    
    data.sign[data.sign < 0] = 0
    data.sign = to_categorical(data.sign, 2).reshape(-1,2).float()
    
    #data.ad = torch.cat([data.acts, data.degree], axis=-1)
    
    return data

def csv2graph(fname, prot_dict):
    
    df = pd.read_csv('../../snac_data/' + fname)
    G = nx.from_pandas_edgelist(df, source='node1', target='node2', edge_attr='sign', create_using=nx.DiGraph)
    nx.set_node_attributes(G, prot_dict, 'x')
    data = from_networkx(G)
    
    # Custom Standard Scaler
    m = data.x.mean(0, keepdim=True)
    s = data.x.std(0, unbiased=False, keepdim=True)
    data.x -= m
    data.x /= s
    
    return data

def load_prot_embs(size, norm=False):
    prot_path = 'data/prot_embeddings/linear_corex_embeddings/lc_seqveq_embedding_dicts/lc_seq_embeddings_dict_{}.pkl'.format(size)
    prot_dict = {}
    global_dict = {}

    unique_prots = 'data/prot_embeddings/new_features/proteins.csv'
    unique_df = pd.read_csv(unique_prots)

    for idx, prot in enumerate(unique_df.proteins.to_numpy()):
        global_dict[prot] = idx

    with open(prot_path, 'rb') as f:
        prot_dict = pickle.load(f)

    prot_dict['NKX3~1'] = prot_dict.pop('NKX3-1') # due to key mismatch
    prot_embs = np.array([prot_dict[key] for key in global_dict.keys()])
    
    if norm and size != 1024:
        norm = Normalizer()
        prot_embs = norm.fit_transform(prot_embs)
    
    return prot_embs, global_dict

class SNDatasetAuto(Dataset):
    def __init__(self, fnames, global_dict):
        super(SNDatasetAuto, self).__init__()
        self.fnames = fnames
        self.gd = global_dict
        
    def len(self):
        return len(self.fnames)
        
    def get(self, idx):
        return ucsv2graph(self.fnames[idx], self.gd)
    
class W_SNDataset(Dataset):
    def __init__(self, fnames, prot_dict, norm=False, transform=None):
        super(W_SNDataset, self).__init__()
        self.fnames = fnames
        self.pd = prot_dict
        self.norm = norm
        self.transform = transform
        
    def len(self):
        return len(self.fnames)
        
    def get(self, idx):
        
        data = wcsv2graph(self.fnames[idx], self.pd, self.norm)

        return data