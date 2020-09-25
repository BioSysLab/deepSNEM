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
from torch_geometric.nn import Set2Set

from tqdm.auto import tqdm

from models.deep_graph_infomax.infomax import DeepGraphInfomax
from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder
from utils.data_gen import SNDatasetInfomax, load_prot_embs

from sklearn.model_selection import train_test_split

torch.backends.cudnn.benchmark = True

# Arguments Parser
parser = argparse.ArgumentParser(description='Transformer training hyperparameters.')

parser.add_argument('--emb_dim', type=int, help='Embedding dimension')
parser.add_argument('--n_heads', type=int, help='Number of attention heads', default=4)
parser.add_argument('--n_layers', type=int, help='Number of transformer layers', default=1)

parser.add_argument('--epochs', type=int, help='Number of epochs to train', default=1)
parser.add_argument('--batch_size', metavar='BS', type=int, help='Batch size', default=1)

args = parser.parse_args()

# Load Global Dict and prot embeddings
SIZE = args.emb_dim
EMB_DIM = args.emb_dim
EPOCHS = args.epochs

prot_embs, global_dict = load_prot_embs(SIZE, norm=False)

# Load Unweighted Dataset
unweighted_fnames = 'data/graph_info_df/samples_all.csv'
u_fnames = pd.read_csv(unweighted_fnames)
u_path_list = u_fnames.path_list.to_numpy()

# Create train val splits
X, val = train_test_split(u_path_list, test_size=0.3)
train_data = SNDatasetInfomax(u_path_list, global_dict)
val_data = SNDatasetInfomax(val, global_dict)

train_loader = DataLoader(train_data, batch_size=args.batch_size, num_workers=12)
val_loader = DataLoader(val_data, batch_size=args.batch_size, num_workers=6)

# Network
dev = torch.device('cuda')

encoder = GraphTransformerEncoder(n_layers=args.n_layers, n_heads=args.n_heads, n_hid=EMB_DIM, 
                            pretrained_weights=prot_embs).to(dev)

model = DeepGraphInfomax(hidden_channels=args.emb_dim, encoder=encoder,
                                     summary=lambda z, *args, **kwargs: z.mean(dim=0))).to(dev)

# Training Hyperparameters
lr = 1e-3
optimizer = torch.optim.SGD(model.parameters(), lr=lr)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1, gamma=0.98)

def train(epoch):
    model.train()
    for tb in train_data_iterator:
        tb = tb.to(dev)
        optimizer.zero_grad()

        pos_z, neg_z_x, neg_z_edges, summary = model(tb)

        loss = model.loss(pos_z, neg_z_x, neg_z_edges, summary)
        train_data_iterator.set_postfix(Epoch=epoch+1, MI ='%.4f' % float(loss.item()))

        loss.backward(retain_graph=True)
        #torch.nn.utils.clip_grad_norm_(model.parameters(), 0.05)

        optimizer.step()

        del tb

def eval(epoch):
    model.eval()
    for vb in val_data_iterator:
        vb = vb.to(dev)
        with torch.no_grad():
            pos_z, neg_z, summary = model(vb)

            loss = model.loss(pos_z, neg_z, summary).item()

            val_data_iterator.set_postfix(Epoch=epoch+1, MI_val ='%.4f' % float(loss))

        del vb

    scheduler.step()

for epoch in range(EPOCHS):
    print('-' * 100)
    print('-' * 100)
    train_data_iterator = tqdm(train_loader, leave=True, unit='batch', postfix={'Epoch': epoch+1,'MI': '%.4f' % 0.0,})
    train(epoch)
    #val_data_iterator = tqdm(val_loader, leave=True, unit='batch', postfix={'Epoch': epoch+1,'MI_val': '%.4f' % 0.0,})
    #eval(epoch)

torch.save(model.state_dict(), f'embeddings/deep_graph_infomax/dgi_{args.emb_dim}_tl_{args.n_layers}.pt')