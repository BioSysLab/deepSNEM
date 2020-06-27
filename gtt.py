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
from torch_geometric.utils import from_networkx, to_networkx, batched_negative_sampling, to_dense_adj, add_self_loops
from torch_geometric.nn import Set2Set

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import Normalizer, StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, average_precision_score
from tqdm.auto import tqdm

from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder, PostEncoding
from models.graph_transformer.autoencoder_base import DeepSNEM, LinearDecoder, FermiDiracDecoder, pdist, get_edge_dists, reshape_probs
from utils.data_gen import ucsv2graph, SNDatasetAuto, load_prot_embs

torch.backends.cudnn.benchmark = True

# Arguments Parser
parser = argparse.ArgumentParser(description='Transformer training hyperparameters.')

parser.add_argument('--emb_dim', type=int, help='Embedding dimension')
parser.add_argument('--n_heads', type=int, help='Number of attention heads', default=4)
parser.add_argument('--n_layers', type=int, help='Number of transformer layers', default=1)

parser.add_argument('--epochs', type=int, help='Number of epochs to train', default=1)
parser.add_argument('--batch_size', metavar='BS', type=int, help='Batch size', default=32)
parser.add_argument('--task', type=str, help='Either lp or nr (link prediction or node reconstruction)', default='nr')
parser.add_argument('--ps', type=int, help='Set2Set processing steps', default=2)

args = parser.parse_args()

# Load Global Dict and prot embeddings
SIZE = args.emb_dim
EMB_DIM = args.emb_dim

prot_embs, global_dict = load_prot_embs(SIZE, norm=False)

# Load Unweighted Dataset
unweighted_fnames = 'data/graph_info_df/samples_all.csv'
u_fnames = pd.read_csv(unweighted_fnames)
u_path_list = u_fnames.path_list.to_numpy()

X_path = 'data/graph_info_df/all_pairs3_train_graphs.csv'
X = pd.read_csv(X_path)
X = X.x.to_numpy()

val_path = 'data/graph_info_df/val_set_1.csv'
val = pd.read_csv(val_path)
val = val.graphs.to_numpy()


# Create train val splits

X, val = train_test_split(u_path_list, test_size=0.3)
train_data = SNDatasetAuto(X, global_dict)
val_data = SNDatasetAuto(val, global_dict)

train_loader = DataLoader(train_data, batch_size=args.batch_size, num_workers=12)
val_loader = DataLoader(val_data, batch_size=args.batch_size, num_workers=6)


# Network
dev = torch.device('cuda')

summarizer = Set2Set(args.emb_dim, args.ps).to(dev)
encoder = GraphTransformerEncoder(n_layers=args.n_layers, n_heads=args.n_heads, n_hid=EMB_DIM, 
                            pretrained_weights=prot_embs, summarizer=summarizer).to(dev)
decoder = LinearDecoder(emb_dim=EMB_DIM, original_dim=SIZE).to(dev) if args.task=='nr' else FermiDiracDecoder(1.0).to(dev)

autoenc = DeepSNEM(encoder, decoder).to(dev)

# Training Hyperparameters
lr = 5.0
optimizer = torch.optim.SGD(autoenc.parameters(), lr=lr)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1, gamma=0.9)

EPOCHS = args.epochs

def train_nr(epoch):
    """Training function for node reconstruction."""
    autoenc.train()
    for tb in train_data_iterator:
        tb.to(dev)
        optimizer.zero_grad()

        node_embs = encoder(tb)
        summary = encoder.summarize(node_embs, tb.batch)
        dec = node_embs * summary
            
        loss = F.mse_loss(dec, autoenc.encoder.emb_layer(tb.global_idx))
        train_data_iterator.set_postfix(Epoch=epoch,mse ='%.4f' % float(loss.item()))
            
        loss.backward(retain_graph=True)
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

def train_lp(epochs):
    """Training function for link prediction."""
    autoenc.train()
    for tb in train_data_iterator:
        tb.to(dev)

        optimizer.zero_grad()
            
        embs = autoenc.encode(tb)
        summary = autoenc.encoder.summarize(embs, tb.batch)

        embs = embs * summary

        pos_dists = get_edge_dists(embs, tb.pos_childs)
        neg_dists = get_edge_dists(embs, tb.neg_childs)

        pos_dists = reshape_probs(pos_dists)
        neg_dists = reshape_probs(neg_dists)

        pos_probs = autoenc.decode(pos_dists)
        neg_probs = autoenc.decode(neg_dists)

        pred = torch.cat([pos_probs, neg_probs], dim=-1).detach().cpu()
        pos_y = torch.ones_like(pos_probs).detach().cpu()
        neg_y = torch.zeros_like(neg_probs).detach().cpu()
        y = torch.cat([pos_y, neg_y], dim=-1)

        roc = roc_auc_score(y, pred)
        avg_p = average_precision_score(y, pred)

        pos_loss = F.binary_cross_entropy(pos_probs, torch.ones_like(pos_probs))
        neg_loss = F.binary_cross_entropy(neg_probs, torch.zeros_like(neg_probs))
        loss = pos_loss + neg_loss

        train_data_iterator.set_postfix(Epoch=epoch+1,
                                        loss ='%.4f' % float(loss.item()),
                                        pos_loss ='%.4f' % float(pos_loss.item()),
                                        neg_loss ='%.4f' % float(neg_loss.item()), 
                                        roc_auc = '%.4f' % float(roc.item()),
                                        ap = '%.4f' % float(avg_p.item()))
            
        loss.backward(retain_graph=True)
        torch.nn.utils.clip_grad_norm_(autoenc.parameters(), 0.05)
        optimizer.step()

        del tb
        
    val_data_iterator = tqdm(val_loader,leave=True,unit='batch',
                        postfix={'Epoch': epoch+1,'val_loss': '%.4f' % 0.0, 'roc_auc': '%.4f' % 0.0, 'ap': '%.4f' % 0.0})
    autoenc.eval()
    for tb in val_data_iterator:
        tb = tb.to(dev)
        with torch.no_grad():
            embs = autoenc.encode(tb)
            summary = autoenc.encoder.summarize(embs, tb.batch)

            embs = embs * summary

            pos_dists = get_edge_dists(embs, tb.pos_childs)
            neg_dists = get_edge_dists(embs, tb.neg_childs)

            pos_dists = reshape_probs(pos_dists)
            neg_dists = reshape_probs(neg_dists)

            pos_probs = autoenc.decode(pos_dists)
            neg_probs = autoenc.decode(neg_dists)

            pred = torch.cat([pos_probs, neg_probs], dim=-1).detach().cpu()
            pos_y = torch.ones_like(pos_probs).detach().cpu()
            neg_y = torch.zeros_like(neg_probs).detach().cpu()
            y = torch.cat([pos_y, neg_y], dim=-1)
            roc = roc_auc_score(y, pred)
            avg_p = average_precision_score(y, pred)

            pos_loss = F.binary_cross_entropy(pos_probs, torch.ones_like(pos_probs)).mean()
            neg_loss = F.binary_cross_entropy(neg_probs, torch.zeros_like(neg_probs)).mean()
            loss = pos_loss + neg_loss
            val_data_iterator.set_postfix(Epoch=epoch+1,val_loss='%.4f' % float(loss.item()),
                                        roc_auc = '%.4f' % float(roc.item()),
                                        ap = '%.4f' % float(avg_p.item()))
                
    scheduler.step() 

train_func = train_nr if args.task=='nr' else train_lp

for epoch in range(EPOCHS):
    print('-' * 100)
    print('-' * 100)
    if args.task=='nr':
        train_data_iterator = tqdm(train_loader,leave=True,unit='batch',postfix={'Epoch': epoch+1,'mse': '%.4f' % 0.0,})
    else:
        train_data_iterator = tqdm(train_loader,leave=True,unit='batch',
                        postfix={'Epoch': epoch+1,'loss': '%.4f' % 0.0,
                        'loss': '%.4f' % 0.0,
                        'pos_loss': '%.4f' % 0.0,
                        'neg_loss': '%.4f' % 0.0,
                        'roc_auc': '%.4f' % 0.0, 
                        'ap': '%.4f' % 0.0})
    train_func(epoch) 

torch.save(autoenc.state_dict(), f'embeddings/autoencoder_graph/gt_{args.emb_dim}_tl_{args.n_layers}_{args.task}.pt')

# Save the embeddings
