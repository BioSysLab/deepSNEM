from __future__ import absolute_import, division
import torch
from torch.functional import F
from torch_geometric.data import Dataset
from torch_geometric.utils import from_networkx, degree, to_networkx, negative_sampling, add_self_loops
import pandas as pd
import numpy as np
import pickle
import networkx as nx
from tqdm.auto import tqdm

from sklearn.preprocessing import MaxAbsScaler, MinMaxScaler, Normalizer, RobustScaler, StandardScaler
from sklearn.model_selection import train_test_split

def to_categorical(y, num_classes):
    """ 1-hot encodes a tensor """
    return torch.eye(num_classes, dtype=torch.long)[y]

def wcsv2graph(fname, global_dict, y):
    """
    Weighted Graph Creator
    """
    sample = pd.read_csv('../snac_data/' + fname)
    
    G = nx.from_pandas_edgelist(sample, source='node1', target='node2', 
                            edge_attr=['sign', 'weight'], create_using=nx.DiGraph())

    n1a1d = sample[['node1','downact1']]
    n1a1u = sample[['node1','upact1']]
    n1a1d.columns = ['node','downact']
    n1a1u.columns = ['node', 'upact']

    n2a2d = sample[['node2','downact2']]
    n2a2u = sample[['node2','upact2']]
    n2a2d.columns = ['node','downact']
    n2a2u.columns = ['node', 'upact']
    
    nad = pd.concat([n1a1d,n2a2d])
    nad = nad.drop_duplicates('node')
    nad = nad.set_index('node')
    nad['downacts'] = nad[['downact']].apply(lambda x: np.hstack(x), axis=1)
    nad = nad.drop(['downact'], axis=1)['downacts'].to_dict()
    
    nau = pd.concat([n1a1u,n2a2u])
    nau = nau.drop_duplicates('node')
    nau = nau.set_index('node')
    nau['upacts'] = nau[['upact']].apply(lambda x: np.hstack(x), axis=1)
    nau = nau.drop(['upact'], axis=1)['upacts'].to_dict()
    
    nx.set_node_attributes(G, global_dict,'global_idx')
    nx.set_node_attributes(G, nad, 'downacts')
    nx.set_node_attributes(G, nau, 'upacts')
    
    data = from_networkx(G)
    
    G = to_networkx(data)
    data.weight = data.weight.float()
    data.downacts, data.upacts = data.downacts.double(), data.upacts.float()
    data.acts = torch.cat([data.downacts.float(), data.upacts.float()], dim=-1).double()
    data.downacts = data.upacts = None
    
    data.sign[data.sign < 0] = 0
    data.sign = data.sign.long()
    data.sign = to_categorical(data.sign, 2).reshape(-1,2).float()
    
    data.y = torch.tensor(y)
    data.label = torch.tensor(np.argmax(y)).view(-1).long()
    
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
    
    G = to_networkx(data)
    G = G.to_undirected()
    
    neighbors = {}
    for g in G.nodes():
        neighbors[g] = [n for n in G.neighbors(g)]
    
    nx.set_node_attributes(G, neighbors, 'pos_childs')
    
    
    data_2 = from_networkx(G)
    data_2.edge_index, _ = add_self_loops(data_2.edge_index)
    data_2.edge_index = negative_sampling(data_2.edge_index, num_neg_samples=700)
    
    G = to_networkx(data_2)
    G = G.to_undirected()
    
    neighbors = {}
    for g in G.nodes():
        neighbors[g] = [n for n in G.neighbors(g)]
    
    nx.set_node_attributes(G, neighbors, 'neg_childs')
    
    data_3 = from_networkx(G)

    data.acts[data.acts < 0] = 0
    data.acts = to_categorical(data.acts, 2).reshape(-1,2).long()
    
    data.sign[data.sign < 0] = 0
    data.sign = to_categorical(data.sign, 2).reshape(-1,2).float()
    
    data.pos_childs = data_2.pos_childs
    data.neg_childs = data_3.neg_childs
    
    return data

def ucsv2graph_infomax(fname, global_dict):
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
    
    G = to_networkx(data)
    data.acts[data.acts < 0] = 0
    data.acts = to_categorical(data.acts, 2).reshape(-1,2).long()
    
    data.sign[data.sign < 0] = 0
    data.sign = to_categorical(data.sign, 2).reshape(-1,2).float()
    
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

class SNDatasetInfomax(Dataset):
    def __init__(self, fnames, global_dict):
        super(SNDatasetInfomax, self).__init__()
        self.fnames = fnames
        self.gd = global_dict
        
    def len(self):
        return len(self.fnames)
        
    def get(self, idx):
        return ucsv2graph_infomax(self.fnames[idx], self.gd)

class SNLDataset(Dataset):
    def __init__(self, fnames, y, global_dict):
        super(SNLDataset, self).__init__()
        self.fnames = fnames
        self.gd = global_dict
        self.y = y
        
    def len(self):
        return len(self.fnames)
        
    def get(self, idx):
        return wcsv2graph(self.fnames[idx], self.gd, self.y[idx])

def deg_distr(train_data):
    deg = torch.zeros(58, dtype=torch.long)
    for data in tqdm(train_data):
        d = degree(data.edge_index[0], num_nodes=data.num_nodes, dtype=torch.long)
        deg += torch.bincount(d, minlength=58)
    return deg