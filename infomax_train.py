from __future__ import absolute_import, division

import argparse

import numpy as np
import pandas as pd
import math
import re
import sys
sys.path.append('..')

import torch
import torch.nn as nn
from torch.functional import F

from torch_geometric.data import Data, DataLoader, Dataset
from torch_geometric.nn import Set2Set

from tqdm.auto import tqdm

from models.deep_graph_infomax.infomax import DeepGraphInfomax
from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder, MultipleOptimizer
from utils.data_gen import SNDatasetInfomax, load_prot_embs, calc_pos_mat, cid, load_prot_embs_go

from sklearn.model_selection import train_test_split

torch.backends.cudnn.benchmark = True

# Arguments Parser
parser = argparse.ArgumentParser(description='Transformer training hyperparameters.')

parser.add_argument('--emb_type', type=str, help='Embedding type (GO, seq, mixed)')
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
TYPE = args.emb_type

prot_embs, global_dict = load_prot_embs(SIZE, norm=False) if TYPE=='seqveq' else load_prot_embs_go(SIZE, norm=False)
pos_mat = calc_pos_mat(512)

# Load Unweighted Dataset
unweighted_fnames = 'data/graph_info_df/samples_all.csv'
u_fnames = pd.read_csv(unweighted_fnames)
u_path_list = u_fnames.path_list.to_numpy()
us_cellid = []
for us in u_path_list:
    x = re.split('_', us)
    us_cellid.append(x[2])
u_fnames['cell_id'] = us_cellid
cellid = np.array(us_cellid)

# Network
dev = torch.device('cuda')

encoder = GraphTransformerEncoder(n_layers=args.n_layers, n_heads=args.n_heads, n_hid=EMB_DIM, 
                            pretrained_weights=prot_embs).to(dev)

model = DeepGraphInfomax(hidden_channels=args.emb_dim, encoder=encoder,
                                     summary=lambda z, *args, **kwargs: z.mean(dim=0)).to(dev)

# Training Hyperparameters
lr_sparse = 1e-5
lr = 1e-5
optimizer_sparse = torch.optim.SparseAdam(model.encoder.emb_layer.parameters(), lr=lr_sparse)
optimizer_dense = torch.optim.AdamW(list(model.encoder.parameters())[1:], lr=lr)
optimizer = MultipleOptimizer(optimizer_sparse, optimizer_dense)
#optimizer = torch.optim.SGD(model.parameters(), lr=lr)

def train(epoch):
    model.train()
    for tb, neg in train_data_iterator:
        tb = tb.to(dev)
        neg = neg.to(dev)
        optimizer.zero_grad()

        pos_z, _, _, summary = model(tb)
        neg_sig, _, _, _ = model(neg)

        pos_loss = -torch.log(
            model.discriminate(pos_z, summary, sigmoid=True) + 1e-5).mean()
        neg_loss_sigs = -torch.log(
            1 - model.discriminate(neg_sig, summary, sigmoid=True) + 1e-5).mean()
        loss = pos_loss + neg_loss_sigs

        train_data_iterator.set_postfix(Epoch=epoch+1, MI ='%.4f' % float(loss.item()))

        loss.backward(retain_graph=True)
        #torch.nn.utils.clip_grad_norm_(model.encoder.emb_layer.parameters(), 0.05)

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

for epoch in range(EPOCHS):
    print('-' * 100)
    print('Finding Negative Signatures')
    negs = cid(u_path_list, cellid)
    print('Done!')
    print('-' * 100)
    print('Starting training...')
    pn = SNDatasetInfomax(u_path_list, negs, global_dict, pos_mat)
    train_loader = DataLoader(pn, batch_size=1, num_workers=12)
    train_data_iterator = tqdm(train_loader, leave=True, unit='batch', postfix={'Epoch': epoch+1,'MI': '%.4f' % 0.0,})
    train(epoch)
    #val_data_iterator = tqdm(val_loader, leave=True, unit='batch', postfix={'Epoch': epoch+1,'MI_val': '%.4f' % 0.0,})
    #eval(epoch)

torch.save(model.state_dict(), f'embeddings/deep_graph_infomax/dgi_{args.emb_dim}_tl_{args.n_layers}.pt')