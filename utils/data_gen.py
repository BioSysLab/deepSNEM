from __future__ import absolute_import, division
import torch
import torch.nn as nn
from torch.functional import F
from torch_geometric.data import Dataset, Data
from torch_geometric.utils import from_networkx, degree, to_networkx, negative_sampling, add_self_loops, to_undirected

from torch_sparse import coalesce

import pandas as pd
import numpy as np
import pickle
import networkx as nx
from tqdm.auto import tqdm

from sklearn.preprocessing import MaxAbsScaler, MinMaxScaler, Normalizer, RobustScaler, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from sklearn.utils.class_weight import compute_class_weight

from karateclub.utils.treefeatures import WeisfeilerLehmanHashing
from karateclub.node_embedding import Role2Vec

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
    aps = nx.floyd_warshall_numpy(G)
    aps[aps == np.inf] = 0.
    aps = torch.FloatTensor(aps)

    G1 = G.reverse()
    aps_r = nx.floyd_warshall_numpy(G1)
    aps_r[aps_r == np.inf] = 0.
    aps_r = torch.FloatTensor(aps_r)

    G = G.to_undirected()

    neighbors = {}
    for g in G.nodes():
        neighbors[g] = [n for n in G.neighbors(g)]
    
    nx.set_node_attributes(G, neighbors, 'pos_childs')

    data_2 = from_networkx(G)
    data_2.edge_index, _ = add_self_loops(data_2.edge_index)
    data_2.edge_index = negative_sampling(data_2.edge_index, num_neg_samples=800)
    
    G = to_networkx(data_2)
    G = G.to_undirected()
    
    neighbors = {}
    for g in G.nodes():
        neighbors[g] = [n for n in G.neighbors(g)]
    
    nx.set_node_attributes(G, neighbors, 'neg_childs')

    data_3 = from_networkx(G)

    data.weight = data.weight.float()
    data.downacts, data.upacts = data.downacts.double(), data.upacts.float()
    data.acts = torch.cat([data.downacts.float(), data.upacts.float()], dim=-1).double()
    data.downacts = data.upacts = None
    
    data.sign[data.sign < 0] = 0
    data.sign = data.sign.long()
    data.sign = to_categorical(data.sign, 2).reshape(-1,2).float()

    data.pos_childs = data_2.pos_childs
    data.neg_childs = data_3.neg_childs
    
    data.y = torch.tensor(y)
    data.label = torch.tensor(np.argmax(y)).view(-1).long()

    data.seq_mat = torch.add(aps, aps_r)
    
    return data

def wcsv2graph_infomax(fname, global_dict):
    """
    Weighted Graph Creator
    """
    sample = pd.read_csv('../snac_data/' + fname)
    
    G = nx.from_pandas_edgelist(sample, source='node1', target='node2', 
                            edge_attr=['sign', 'weight'], create_using=nx.DiGraph())

    data1 = from_networkx(G)
    G1 = to_networkx(data1)
    rv = Role2Vec(dimensions=512, workers=-1, seed=69)
    rv.fit(G1.to_undirected())

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
    rv = Role2Vec(dimensions=512, workers=-1, seed=69)
    rv.fit(G.to_undirected())
    #aps = nx.floyd_warshall_numpy(G)
    #aps[aps == np.inf] = 0.
    #aps = torch.FloatTensor(aps)

    #G1 = G.reverse()
    #aps_r = nx.floyd_warshall_numpy(G1)
    #aps_r[aps_r == np.inf] = 0.
    #aps_r = torch.FloatTensor(aps_r)

    data.weight = data.weight.float()
    data.downacts, data.upacts = data.downacts.double(), data.upacts.float()
    data.acts = torch.cat([data.downacts.float(), data.upacts.float()], dim=-1).double()
    data.downacts = data.upacts = None
    
    data.sign[data.sign < 0] = 0
    data.sign = data.sign.long()
    data.sign = to_categorical(data.sign, 2).reshape(-1,2).float()

    data.r2v = torch.tensor(rv.get_embedding())

    #data.seq_mat = torch.add(aps, aps_r)
    
    return data

def ucsv2graph(fname, global_dict):
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
    aps = nx.floyd_warshall_numpy(G)
    aps[aps == np.inf] = 0.
    aps = torch.FloatTensor(aps)

    G1 = G.reverse()
    aps_r = nx.floyd_warshall_numpy(G1)
    aps_r[aps_r == np.inf] = 0.
    aps_r = torch.FloatTensor(aps_r)
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

    data.seq_mat = torch.add(aps, aps_r)
    
    return data

def ucsv2graph_infomax(fname, global_dict, sig_one_hot=None, y=None):
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

    data.weight = torch.ones(data.num_edges)

    if sig_one_hot is not None:
        data.sig = torch.tensor(sig_one_hot)
        data.node_sig = data.sig.view(1,-1).repeat(data.num_nodes,1)

    if y is not None:
        data.y = y
        data.node_y = data.y.view(-1,1).repeat(data.num_nodes,1)

    #data.seq_mat = torch.add(aps, aps_r)

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

def load_prot_embs_go(size, norm=False):
    prot_path = 'data/prot_embeddings/linear_corex_embeddings/lc_go_term_embeddings/lc_go_term_embeddings_dict_{}.pkl'.format(size)
    prot_dict = {}
    global_dict = {}

    unique_prots = 'data/prot_embeddings/new_features/proteins.csv'
    unique_df = pd.read_csv(unique_prots)

    for idx, prot in enumerate(unique_df.proteins.to_numpy()):
        global_dict[prot] = idx

    with open(prot_path, 'rb') as f:
        prot_dict = pickle.load(f)

    #prot_dict['NKX3~1'] = prot_dict.pop('NKX3-1') # due to key mismatch
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

class WSNDataset(Dataset):
    def __init__(self, fnames, global_dict):
        super(WSNDataset, self).__init__()
        self.fnames = fnames
        self.gd = global_dict
        
    def len(self):
        return len(self.fnames)
        
    def get(self, idx):
        return wcsv2graph_infomax(self.fnames[idx], self.gd)

class SNDatasetInfomax(Dataset):
    def __init__(self, fnames, global_dict, sig_id_one_hot=None):
        super(SNDatasetInfomax, self).__init__()
        self.fnames = fnames
        self.gd = global_dict
        self.sig_id = sig_id_one_hot

    def len(self):
        return len(self.fnames)
        
    def get(self, idx):
        if self.sig_id is not None:
            data = ucsv2graph_infomax(self.fnames[idx], self.gd, self.sig_id[idx])
        else:
            data = ucsv2graph_infomax(self.fnames[idx], self.gd)
        return data

class SNDatasetInfomaxSemi(Dataset):
    def __init__(self, fnames, global_dict, sig_one_hot, y):
        super(SNDatasetInfomaxSemi, self).__init__()
        self.fnames = fnames
        self.sig_id = sig_one_hot
        self.gd = global_dict
        self.y = y

    def len(self):
        return len(self.fnames)
        
    def get(self, idx):
        data = ucsv2graph_infomax(self.fnames[idx], self.gd, self.sig_id[idx], self.y[idx])
        return data

class PositionalEmbedding(nn.Module):
    def __init__(self, d):
        super().__init__()
        self.d = d
        inv_freq = 1 / (10000 ** (torch.arange(0.0, d, 2.0) / d))
        self.register_buffer("inv_freq", inv_freq)
        
    def forward(self, positions: torch.LongTensor, # (seq, )
               ):
        # outer product
        sinusoid_inp = torch.einsum("i,j->ij", positions.float(), self.inv_freq)
        pos_emb = torch.cat([sinusoid_inp.sin(), sinusoid_inp.cos()], dim=-1)
        return pos_emb[:,None,:]

    def __init__(self, seq_mat):
        super(SigNetData, self).__init__()
        self.seq_mat = seq_mat

    def none(self):
        pass