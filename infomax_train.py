from __future__ import absolute_import, division
from utils.data_gen import load_prot_embs, load_prot_embs_go, SNDatasetInfomax
from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder, MultipleOptimizer, MultipleScheduler
from models.deep_graph_infomax.infomax import SNInfomax
from sklearn.preprocessing import OneHotEncoder
from tqdm.auto import tqdm
from torch_geometric.nn import Set2Set
from torch_geometric.data import DataLoader
import torch

import argparse

import numpy as np
import pandas as pd
import sys
import random
import os
sys.path.append('..')

torch.backends.cudnn.benchmark = True


def seed_everything(seed=42):
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.benchmark = False


# Arguments Parser
parser = argparse.ArgumentParser(
    description='Transformer training hyperparameters.')

parser.add_argument('--emb_type',
                    type=str,
                    help='Embedding type (GO, seq, mixed)')
parser.add_argument('--emb_dim', type=int, help='Embedding dimension')
parser.add_argument('--n_heads',
                    type=int,
                    help='Number of attention heads',
                    default=4)
parser.add_argument('--n_layers',
                    type=int,
                    help='Number of transformer layers',
                    default=1)

parser.add_argument('--epochs',
                    type=int,
                    help='Number of epochs to train',
                    default=1)
parser.add_argument('--batch_size',
                    metavar='BS',
                    type=int,
                    help='Batch size',
                    default=1)
parser.add_argument('--task', type=str, help='un or semi', default='un')

args = parser.parse_args()

# Load Global Dict and prot embeddings
SIZE = args.emb_dim
BS = args.batch_size
EMB_DIM = args.emb_dim
EPOCHS = args.epochs
TYPE = args.emb_type
TASK = args.task
SEMI = False if TASK == 'un' else True

prot_embs, global_dict = load_prot_embs(
    SIZE, norm=False) if TYPE == 'seqveq' else load_prot_embs_go(SIZE,
                                                                 norm=False)
if TYPE == 'random':
    prot_embs = None

# Load Weighted, Unweighted and SemiSupervision Dataset
oh = OneHotEncoder()

u_fnames = pd.read_csv('data/graph_info_df/full_dataset.csv')
u_path_list = u_fnames.files_combined.to_numpy()
labels = u_fnames.sigs_g.to_numpy().reshape(-1, 1)
labels = oh.fit_transform(labels).toarray()

train_data = SNDatasetInfomax(u_path_list, global_dict, labels)

train_loader = DataLoader(train_data,
                          batch_size=BS,
                          num_workers=12,
                          shuffle=True)

# Network
dev = torch.device('cuda')

summarizer = Set2Set(args.emb_dim, 3)
encoder = GraphTransformerEncoder(n_layers=args.n_layers,
                                  n_heads=args.n_heads,
                                  n_hid=EMB_DIM,
                                  pretrained_weights=prot_embs).to(dev)

model = SNInfomax(hidden_channels=args.emb_dim,
                  encoder=encoder,
                  summary=summarizer,
                  semi=SEMI).to(dev)

# Training Hyperparams
lr_sparse = 1e-4
lr = 1e-4
optimizer_sparse = torch.optim.SparseAdam(model.encoder.emb_layer.parameters(),
                                          lr=lr_sparse)
optimizer_dense = torch.optim.AdamW(list(model.parameters())[1:], lr=lr)

scheduler_dense = torch.optim.lr_scheduler.StepLR(optimizer_dense,
                                                  step_size=4,
                                                  gamma=0.05)
scheduler_sparse = torch.optim.lr_scheduler.StepLR(optimizer_sparse,
                                                   step_size=4,
                                                   gamma=0.05)

optimizer = MultipleOptimizer(optimizer_sparse, optimizer_dense)
scheduler = MultipleScheduler(scheduler_sparse, scheduler_dense)


def train(epoch):
    model.train()
    for tb in train_data_iterator:
        tb = tb.to(dev)
        optimizer.zero_grad()

        res, summary = model(tb)
        local, prior = model.loss(res, summary)
        loss = local + prior

        train_data_iterator.set_postfix(
            Epoch=epoch + 1,
            MI='%.4f' % float(loss.item()),
            local_loss='%.4f' % float(local.item()),
            prior_loss='%.4f' % float(prior.item()))

        loss.backward(retain_graph=True)
        # torch.nn.utils.clip_grad_norm_(model.encoder.emb_layer.parameters(), 0.05)
        optimizer.step()

        del tb

    scheduler.step()


print(f'Training with task: {TASK}')

for epoch in range(EPOCHS):
    print('-' * 100)
    print(f'SparseAdam Learning Rate: {lr_sparse}')
    print(f'AdamW Learning Rate: {lr}')
    print('-' * 100)
    print(f'Now on Epoch {epoch}')
    print('-' * 100)

    train_data_iterator = tqdm(train_loader,
                               leave=True,
                               unit='batch',
                               postfix={
                                   'Epoch': epoch + 1,
                                   'MI': '%.4f' % 0.0,
                               })
    train(epoch)

torch.save(
    model.state_dict(),
    f'embeddings/deep_graph_infomax/unsupervised/DGI_JSD_{args.emb_dim}_{TYPE}_uniform_{TASK}.pt'
)
