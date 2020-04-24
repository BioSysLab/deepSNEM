from __future__ import absolute_import, division
import torch
from torch_geometric.data import Dataset
import pandas as pd
import numpy as np
import pickle

from sklearn.preprocessing import MaxAbsScaler, MinMaxScaler, Normalizer, RobustScaler, StandardScaler
from sklearn.model_selection import train_test_split

def csv2edge(path, prot_dict, shapes_dict, x_scaler=None, b_scaler=MinMaxScaler()):
    """
        Reads the csv file and after converting it to a DataFrame, creates the edge matrix in the appropriate
        COO format requested by torch geometric. Afterwards, creates the x tensor of the node features taken
        from the prot_dict dictionary created.
        Inputs 
            path: path to graph csv file
            prot_dict: protein features dictionary
            shapes_dict: dictionary specifying the shapes of each graph with keys
                        'max_prots', 'num_prot_features', 'num_edge_features'
        Outputs 
            edge_index: edge matrix tensor in COO format : (2, num_edges)
            x: node features tensor with shape (num_nodes, num_features)
    """
    npf = 'num_prot_features' in shapes_dict
    nef = 'num_edge_features' in shapes_dict
    assert npf and nef
    
    graph = pd.read_csv('../../snac_data/' + path)
    node1 = graph.node1.to_numpy()
    node2 = graph.node2.to_numpy()
    sign = graph.sign.to_numpy()
    
    unique = np.unique(np.concatenate((node1, node2)))
    
    node_dict = {}
    node_act = {}
    
    n1a1 = graph[['node1','activity1']]
    n1a1.columns = ['node', 'act']
    n2a2 = graph[['node2','activity2']]
    n2a2.columns = ['node', 'act']
    na = pd.concat([n1a1,n2a2])
    na = na.drop_duplicates('node').sort_values(by='node')
    
    for idx, node in enumerate(unique):
            node_dict[node] = idx
            
    node1 = [node_dict[n] for n in node1]
    node2 = [node_dict[n] for n in node2]
    n1 = np.reshape(node1, (len(node1), 1))
    n2 = np.reshape(node2, (len(node2), 1))
    sign = np.reshape(sign, (len(sign), 1))
    
    # For the edge index matrix in the appropriate coo format
    pairs_1 = np.hstack((n1, n2))
    pairs_ws = np.hstack((pairs_1, sign))
    
    pairs_1 = pairs_1[np.argsort(pairs_1[:,0])]
    pairs_ws = pairs_ws[np.argsort(pairs_ws[:,0])]
    sign = pairs_ws[:, 2]
    sign = np.reshape(sign, (len(sign), 1))
    
    pairs_flipped = np.flip(pairs_1, axis=1)
    pairs_flipped_ws = np.hstack((pairs_flipped, sign))
    
    pairs = np.concatenate((pairs_ws, pairs_flipped_ws))
    pairs = pairs[np.argsort(pairs[:, 0])]
    edge_index = pairs[:, [0,1]]
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    
    # Creating the node attributes tensor using the prot_dict and shapes_dict files
    num_prot_features, num_edge_features = [p for p in shapes_dict.values()]
    
    max_prots = len(unique)
    x = np.zeros((max_prots, num_prot_features))
    
    for idx, node in enumerate(unique):
        x[idx] = np.hstack((prot_dict[node], na['act'].to_numpy()[idx]))

    x_scale = x_scaler
    if x_scale is not None:
        x = x_scale.fit_transform(x)
    
    x = torch.tensor(x, dtype=torch.float) 
    
    # Creating the edge_features matrix for the edges present in the graph
    edge_feats = pairs[:, 2]
    edge_feats = edge_feats.reshape(-1,1)

    scaler = b_scaler
    edge_feats = scaler.fit_transform(edge_feats)
    edge_feats = edge_feats.reshape(-1)
    edge_feats = torch.tensor(edge_feats, dtype=torch.float)
    
    return edge_index, x, edge_feats

def load_prot_dict(pkl_file):
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

def train_val_paths(train_size, val_ratio, samples_csv, path_to_files=None):
    path_to_files = '../../snac_data/file_info.csv'

    assert path_to_files is not None
    
    finfo = pd.read_csv(path_to_files)

    samples = '../data/samples_all.csv'
    samples = pd.read_csv(samples)
    samples = samples.path_list.to_numpy()
    samples = set(samples)
    non_samples = set(finfo.files_combined).difference(samples)

    fnames = np.random.choice(list(non_samples), train_size)
    X, val = train_test_split(fnames, test_size=val_ratio)

    return X, val, finfo, np.array(list(samples)) 

class SNDataset(Dataset):
    def __init__(self, fnames, prot_dict, shape_dict):
        super(SNDataset, self).__init__()
        self.fnames = fnames
        self.pd = prot_dict
        self.sd = shape_dict 
        
    def len(self):
        return len(self.fnames)
        
    def get(self, idx):
        edge_index, x, edge_attr = csv2edge(self.fnames[idx], self.pd, self.sd)
        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
        return data