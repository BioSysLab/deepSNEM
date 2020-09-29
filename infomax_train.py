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

from tqdm.auto import tqdm

from models.deep_graph_infomax.infomax import DeepGraphInfomax
from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder, MultipleOptimizer
from utils.data_gen import WSNDatasetInfomaxUn, load_prot_embs, load_prot_embs_go, SNDatasetInfomax, WSNDatasetInfomaxSemi
from models.infograph_semi_supervised.infomax import get_positive_expectation, get_negative_expectation

torch.backends.cudnn.benchmark = True

# Arguments Parser
parser = argparse.ArgumentParser(description='Transformer training hyperparameters.')

parser.add_argument('--emb_type', type=str, help='Embedding type (GO, seq, mixed)')
parser.add_argument('--emb_dim', type=int, help='Embedding dimension')
parser.add_argument('--n_heads', type=int, help='Number of attention heads', default=4)
parser.add_argument('--n_layers', type=int, help='Number of transformer layers', default=1)

parser.add_argument('--epochs', type=int, help='Number of epochs to train', default=1)
parser.add_argument('--batch_size', metavar='BS', type=int, help='Batch size', default=1)
parser.add_argument('--graph_data_type', type=str, help='W, U or None', default='U')
parser.add_argument('--task', type=str, help='semi or un', default='un')

args = parser.parse_args()

# Load Global Dict and prot embeddings
SIZE = args.emb_dim
EMB_DIM = args.emb_dim
EPOCHS = args.epochs
TYPE = args.emb_type
GT = args.graph_data_type
TASK = args.task

prot_embs, global_dict = load_prot_embs(SIZE, norm=False) if TYPE=='seqveq' else load_prot_embs_go(SIZE, norm=False)

# Load Weighted, Unweighted and SemiSupervision Dataset
unweighted_fnames = 'data/graph_info_df/samples_all.csv'
u_fnames = pd.read_csv(unweighted_fnames)
u_path_list = u_fnames.path_list.to_numpy()

w_path_list = 'data/graph_info_df/file_info_weighted.csv'
w_fnames = pd.read_csv(w_path_list)
w_path_list = w_fnames.files_weighted.to_numpy()

semi_path_list = 'data/graph_info_df/semi_infomax_pos_neg_triads.csv'
semi_fnames = pd.read_csv(semi_path_list)
sig = semi_fnames.sig.to_numpy()
neg = semi_fnames.neg.to_numpy()

pn = WSNDatasetInfomaxUn(w_path_list, global_dict) if GT=='W' else SNDatasetInfomax(u_path_list, global_dict)
if TASK=='semi':
    pn = WSNDatasetInfomaxSemi(sig, neg, global_dict)
train_loader = DataLoader(pn, batch_size=1, num_workers=12, shuffle=True)

# Network
dev = torch.device('cuda')

encoder = GraphTransformerEncoder(n_layers=args.n_layers, n_heads=args.n_heads, n_hid=EMB_DIM, 
                            pretrained_weights=prot_embs).to(dev)

model = DeepGraphInfomax(hidden_channels=args.emb_dim, encoder=encoder,
                                     summary=lambda z, *args, **kwargs: z.mean(dim=0)).to(dev)


# Training Hyperparameters
lr_sparse = 1e-4
lr = 1e-4
optimizer_sparse = torch.optim.SparseAdam(model.encoder.emb_layer.parameters(), lr=lr_sparse)
optimizer_dense = torch.optim.AdamW(list(model.encoder.parameters())[1:], lr=lr)

if TASK=='semi':
    encoder2 = GraphTransformerEncoder(n_layers=args.n_layers, n_heads=args.n_heads, n_hid=EMB_DIM, 
                            pretrained_weights=prot_embs).to(dev)

    model2 = DeepGraphInfomax(hidden_channels=args.emb_dim, encoder=encoder,
                                     summary=lambda z, *args, **kwargs: z.mean(dim=0)).to(dev)
    model.load_state_dict(torch.load('embeddings/deep_graph_infomax/GO_dgi_512_tl_1_un.pt'))
    model2.load_state_dict(torch.load('embeddings/deep_graph_infomax/GO_dgi_512_tl_1_un.pt'))
    print('-' * 100)
    print('Keys Loaded Succesfully!')
    print('-' * 100)
    optimizer_sparse2 = torch.optim.SparseAdam(model2.encoder.emb_layer.parameters(), lr=lr_sparse)
    optimizer_dense2 = torch.optim.AdamW(list(model2.encoder.parameters())[1:], lr=lr)

    optimizer = MultipleOptimizer(optimizer_sparse, optimizer_dense, optimizer_sparse2, optimizer_dense2)

scheduler_dense = torch.optim.lr_scheduler.StepLR(optimizer_dense, step_size=1, gamma=0.05)
scheduler_sparse = torch.optim.lr_scheduler.StepLR(optimizer_sparse, step_size=1, gamma=0.05)
if TASK != 'semi':
    optimizer = MultipleOptimizer(optimizer_sparse, optimizer_dense)
#optimizer = torch.optim.SGD(model.parameters(), lr=lr)

def train(epoch):
    model.train()
    for tb in train_data_iterator:
        tb = tb.to(dev)
        optimizer.zero_grad()

        pos_z, neg_z, summary = model(tb)

        #pos_loss = -torch.log(
            #model.discriminate(pos_z, summary, sigmoid=True) + 1e-5).mean()
        #neg_loss = -torch.log(
            #1 - model.discriminate(neg_z, summary, sigmoid=True) + 1e-5).mean()
        
        pos_loss = get_positive_expectation(model.discriminate(pos_z, summary, fd=False), 'JSD')
        neg_loss = get_negative_expectation(model.discriminate(neg_z, summary, fd=False), 'JSD')

        loss = neg_loss - pos_loss

        train_data_iterator.set_postfix(Epoch=epoch+1, MI ='%.4f' % float(loss.item()))

        loss.backward(retain_graph=True)
        #torch.nn.utils.clip_grad_norm_(model.encoder.emb_layer.parameters(), 0.05)

        optimizer.step()

        del tb
    
    scheduler_dense.step()
    scheduler_sparse.step()

def train_semi(epoch):
    model.train()
    for tb, nb in train_data_iterator:
        tb = tb.to(dev)
        nb = nb.to(dev)
        optimizer.zero_grad()

        sig_z_un, neg_z_un, summary_un = model2(tb)
        sig_z, _, summary = model(tb)
        neg_z, _, _ = model(nb)
        l = 1e-3

        # Supervised Losses
        sig_loss = get_positive_expectation(model.discriminate(sig_z, summary, fd=False), 'JSD')
        neg_sig_loss = get_negative_expectation(model.discriminate(neg_z, summary, fd=False), 'JSD')
        # Unsupervised Losses
        pos_un_loss = get_positive_expectation(model2.discriminate(sig_z_un, summary_un, fd=False), 'JSD')
        neg_un_loss = get_negative_expectation(model2.discriminate(neg_z_un, summary_un, fd=False), 'JSD')
        # Global Embedding Losses
        mul = torch.matmul(summary, summary_un)
        pos_glob_loss = get_positive_expectation(mul, 'JSD', average=False)
        neg_glob_loss = get_negative_expectation(mul, 'JSD', average=False)

        loss = neg_sig_loss - sig_loss + neg_un_loss - pos_un_loss + l*(neg_glob_loss - pos_glob_loss)

        train_data_iterator.set_postfix(Epoch=epoch+1, MI ='%.4f' % float(loss.item()))

        loss.backward(retain_graph=True)
        #torch.nn.utils.clip_grad_norm_(model.encoder.emb_layer.parameters(), 0.05)

        optimizer.step()

        del tb
    
    scheduler_dense.step()
    scheduler_sparse.step()

for epoch in range(EPOCHS):
    print('-' * 100)
    print(f'Training with task: {TASK}')
    #negs = cid(u_path_list, cellid)
    #print('Done!')
    print('-' * 100)
    print('Starting training...')
    train_data_iterator = tqdm(train_loader, leave=True, unit='batch', postfix={'Epoch': epoch+1,'MI': '%.4f' % 0.0,})
    if TASK=='semi':
        train_semi(epoch)
    else:
        train(epoch)
    #torch.save(model.state_dict(), f'embeddings/deep_graph_infomax/dgi_{args.emb_dim}_tl_{args.n_layers}_epoch_{epoch+1}.pt')
    #val_data_iterator = tqdm(val_loader, leave=True, unit='batch', postfix={'Epoch': epoch+1,'MI_val': '%.4f' % 0.0,})
    #eval(epoch)

torch.save(model.state_dict(), f'embeddings/deep_graph_infomax/{GT}_dgi_{args.emb_dim}_tl_{args.n_layers}_{TASK}.pt')

def eval(epoch):
    model.eval()
    for vb in val_data_iterator:
        vb = vb.to(dev)
        with torch.no_grad():
            pos_z, neg_z, summary = model(vb)

            loss = model.loss(pos_z, neg_z, summary).item()

            val_data_iterator.set_postfix(Epoch=epoch+1, MI_val ='%.4f' % float(loss))

        del vb