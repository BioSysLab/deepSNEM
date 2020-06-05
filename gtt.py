from __future__ import absolute_import, division

import argparse

import numpy as np
import pandas as pd
import math
import sys
sys.path.append('..')

import torch
import torch.nn as nn
from torch.functional import F

from torch_geometric.data import Data, DataLoader, Dataset
from torch_geometric.utils import from_networkx, to_networkx

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import Normalizer, StandardScaler
from sklearn.metrics import accuracy_score
from tqdm.auto import tqdm

from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder, PostEncoding
from models.graph_transformer.autoencoder_base import DeepSNEM, LinearDecoder, FermiDiracDecoder
from utils.data_gen import ucsv2graph, SNDatasetAuto, load_prot_embs

torch.backends.cudnn.benchmark = True

# Arguments Parser
parser = argparse.ArgumentParser(description='Transformer training hyperparameters.')

parser.add_argument('--emb_dim', type=int, help='Embedding dimension')
parser.add_argument('--n_heads', type=int, help='Number of attention heads', default=4)
parser.add_argument('--n_layers', type=int, help='Number of transformer layers', default=1)

parser.add_argument('--epochs', type=int, help='Number of epochs to train', default=1)
parser.add_argument('--batch_size', metavar='BS', type=int, help='Batch size', default=32)

args = parser.parse_args()

# Load Global Dict and prot embeddings
SIZE = args.emb_dim
EMB_DIM = args.emb_dim

prot_embs, global_dict = load_prot_embs(SIZE, norm=False)

# Load Unweighted Dataset
unweighted_fnames = 'data/graph_info_df/samples_all.csv'
u_fnames = pd.read_csv(unweighted_fnames)
u_path_list = u_fnames.path_list.to_numpy()

# Create train val splits
X, val = train_test_split(u_path_list, test_size=0.3)
train_data = SNDatasetAuto(X, global_dict)
val_data = SNDatasetAuto(val, global_dict)

train_loader = DataLoader(train_data, batch_size=args.batch_size, num_workers=12)
val_loader = DataLoader(val_data, batch_size=args.batch_size, num_workers=6)


# Network
dev = torch.device('cuda')

encoder = GraphTransformerEncoder(n_layers=args.n_layers, n_heads=args.n_heads, n_hid=EMB_DIM, 
                            pretrained_weights=prot_embs).to(dev)
decoder = LinearDecoder(emb_dim=EMB_DIM, original_dim=SIZE).to(dev)
autoenc = DeepSNEM(encoder, decoder).to(dev)

# Training Hyperparameters
lr = 5.0
optimizer = torch.optim.SGD(autoenc.parameters(), lr=lr)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1, gamma=0.9)

EPOCHS = args.epochs

def train(epoch):
        autoenc.train()
        for tb in train_data_iterator:
            tb.to(dev)
            optimizer.zero_grad()
            
            loss = F.mse_loss(autoenc(tb), autoenc.encoder.emb_layer(tb.global_idx))
            train_data_iterator.set_postfix(Epoch=epoch,mse ='%.4f' % float(loss.item()))
            
            loss.backward()
            torch.nn.utils.clip_grad_norm_(autoenc.parameters(), 0.01)
            optimizer.step()

            del tb
        
        val_data_iterator = tqdm(val_loader,leave=True,unit='batch',postfix={'Epoch': epoch+1,'val_loss': '%.4f' % 0.0,})
        autoenc.eval()
        for vb in val_data_iterator:
            vb = vb.to(dev)
            with torch.no_grad():
                v_rec = autoenc(vb)
                v_original = autoenc.encoder.emb_layer(vb.global_idx)
                
                v_loss = F.mse_loss(v_rec, v_original).item()
                val_data_iterator.set_postfix(Epoch=epoch,val_loss='%.4f' % float(v_loss))
                
        scheduler.step() 

for epoch in range(EPOCHS):
    print('-' * 100)
    print('-' * 100)
    train_data_iterator = tqdm(train_loader,leave=True,unit='batch',postfix={'Epoch': epoch+1,'mse': '%.4f' % 0.0,})
    train(epoch) 

<<<<<<< master
torch.save(autoenc.state_dict(), f'../gt_{args.emb_dim}_tl_{args.n_layers}_leaky_relu.pt')
=======
torch.save(autoenc.state_dict(), f'../gt_{args.emb_dim}_tl_{args.n_layers}_relu.pt')
>>>>>>> various changes and fixes

# Save the embeddings
