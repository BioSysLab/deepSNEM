from __future__ import absolute_import, division
import torch
from torch.functional import F
from torch_geometric.data import Dataset
from torch_geometric.utils import from_networkx
import pandas as pd
import numpy as np
import pickle
import networkx as nx

from sklearn.preprocessing import MaxAbsScaler, MinMaxScaler, Normalizer, RobustScaler, StandardScaler
from sklearn.model_selection import train_test_split

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

def load_prot_dict():
    prot_path = '../data/lc_embeddings_raw.pkl'
    prot_dict = {}
    with open(prot_path, 'rb') as f:
        prot_dict = pickle.load(f)
        
    prot_dict['Perturbation'] = np.zeros((384,))
    prot_dict['IKBKAP'] = np.zeros((384,))
    prot_dict['WHSC1'] = np.zeros((384,))
    prot_dict['ADRBK1'] = np.zeros((384,))
    prot_dict['FIGF'] = np.zeros((384,))
    prot_dict['FAM175A'] = np.zeros((384,))
    prot_dict['PARK2'] = np.zeros((384,))
    prot_dict['NKX3~1'] = np.zeros((384,))
    prot_dict['ERBB2IP'] = np.zeros((384,))

    return prot_dict

def load_go_emb_prot_dict(path_to_dict):
    prot_dict = {}
    with open(path_to_dict, 'rb') as f:
        prot_dict = pickle.load(f)
    
    return prot_dict

def train_val_paths(train_size, val_ratio, samples_csv, path_to_files=None):
    path_to_files = '../../snac_data/file_info.csv'

    assert path_to_files is not None
    
    finfo = pd.read_csv(path_to_files)

    samples = '../data/samples_all.csv'
    samples = pd.read_csv(samples)
    samples = samples.path_list.to_numpy()
    samples = set(samples)
    non_samples = set(finfo.files_combined).difference(samples)

    fnames = np.random.choice(list(samples), train_size)
    X, val = train_test_split(fnames, test_size=val_ratio)

    return X, val, finfo, np.array(list(samples)) 

class SNDatasetAuto(Dataset):
    def __init__(self, fnames, prot_dict, transforms=None):
        super(SNDatasetAuto, self).__init__()
        self.fnames = fnames
        self.pd = prot_dict
        self.transforms = transforms
        
    def len(self):
        return len(self.fnames)
        
    def get(self, idx):
        
        data = csv2graph(self.fnames[idx], self.pd)
        data.train_mask = data.val_mask = data.test_mask = data.y = None
        
        return data
    
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