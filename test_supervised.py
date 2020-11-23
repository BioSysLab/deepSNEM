#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import argparse

import torch
from torch.nn.modules.loss import _WeightedLoss
from torch.functional import F
from torch_geometric.data import DataLoader
from torch_geometric.nn import Set2Set

import pytorch_lightning as pl
from pytorch_lightning.callbacks import EarlyStopping

from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder
from utils.data_gen import load_prot_embs, SNDatasetSup

from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import accuracy_score, average_precision_score

dev = torch.device('cuda')

# Parse input arguments
parser = argparse.ArgumentParser(
    description='Transformer training hyperparameters.')

parser.add_argument('--emb_type',
                    type=str,
                    help='Embedding type (GO, seq)',
                    default='seq')
parser.add_argument('--emb_dim',
                    type=int,
                    help='Embedding dimension',
                    default=512)
parser.add_argument('--n_heads',
                    type=int,
                    help='Number of attention heads',
                    default=4)
parser.add_argument('--n_layers',
                    type=int,
                    help='Number of transformer layers',
                    default=2)
parser.add_argument('--epochs',
                    type=int,
                    help='Number of epochs to train',
                    default=1)
parser.add_argument('--batch_size',
                    metavar='BS',
                    type=int,
                    help='Batch size',
                    default=32)

# Load protein data
prot_embs, global_dict = load_prot_embs(512, norm=False)

# Load graph files
u_fnames = '../snac_data/graph_classification_all.csv'
samples_all = '../snac_data/samples_all.csv'

uf = pd.read_csv(u_fnames)[['sig_id', 'files_combined', 'moa_v1']]
sa = pd.read_csv(samples_all).drop('Unnamed: 0', axis=1)
sa.columns = ['files_combined']
sa_l = pd.merge(uf, sa, on='files_combined')

# Get paths and moas
moas = sa_l.moa_v1.values.reshape(-1, 1)

# One hot encode
oh = OneHotEncoder()
moas = oh.fit_transform(moas).toarray()

# Validation and Test splits
val_set = 'data/graph_info_df/val_set_1.csv'
vs = pd.read_csv(val_set).drop('Unnamed: 0', axis=1)
vs = vs[['graphs', 'moa_v1']]
vs.columns = ['files_combined', 'moa_v1']
vs = pd.merge(uf, vs,
              on='files_combined')[['sig_id', 'moa_v1_x', 'files_combined']]

test_set = 'data/graph_info_df/test_set.csv'
ts = pd.read_csv(test_set)[['graphs', 'moa_v1']]
ts.columns = ['files_combined', 'moa_v1']
ts = pd.merge(ts, uf,
              on='files_combined')[['sig_id', 'moa_v1_x', 'files_combined']]

uf_v = uf[uf['sig_id'].isin(vs.sig_id.values)]
uf = pd.concat([uf, uf_v]).drop_duplicates(keep=False)

uf_t = uf[uf['sig_id'].isin(ts.sig_id.values)]
uf = pd.concat([uf, uf_t]).drop_duplicates(keep=False)

# Datasets and DataLoaders
X_train, y_train = uf.files_combined.values, oh.fit_transform(
    uf.moa_v1.values.reshape(-1, 1)).toarray()
X_val, y_val = vs.files_combined.values, oh.transform(
    vs.moa_v1_x.values.reshape(-1, 1)).toarray()
X_test, y_test = ts.files_combined.values, oh.transform(
    ts.moa_v1_x.values.reshape(-1, 1)).toarray()

train_data = SNDatasetSup(X_train, global_dict, y_train)
val_data = SNDatasetSup(X_val, global_dict, y_val)
test_data = SNDatasetSup(X_test, global_dict, y_test)

# Model Creation

if __name__ == "__main__":
    args = parser.parse_args()

    # Create DataLoader
    train_loader = DataLoader(train_data,
                              batch_size=args.batch_size,
                              num_workers=12,
                              shuffle=True,
                              pin_memory=True)

    val_loader = DataLoader(val_data,
                            batch_size=args.batch_size,
                            num_workers=12,
                            shuffle=False)

    test_loader = DataLoader(test_data,
                             batch_size=args.batch_size,
                             num_workers=12)

    class SmoothBCEwLogits(_WeightedLoss):
        def __init__(self, weight=None, reduction='mean', smoothing=0.0):
            super().__init__(weight=weight, reduction=reduction)
            self.smoothing = smoothing
            self.weight = weight
            self.reduction = reduction

        @staticmethod
        def _smooth(targets: torch.Tensor, n_labels: int, smoothing=0.0):
            assert 0 <= smoothing < 1
            with torch.no_grad():
                targets = targets * (1.0 - smoothing) + 0.5 * smoothing
            return targets

        def forward(self, inputs, targets):
            targets = SmoothBCEwLogits._smooth(targets, inputs.size(-1),
                                               self.smoothing)
            loss = F.binary_cross_entropy_with_logits(inputs, targets,
                                                      self.weight)

            if self.reduction == 'sum':
                loss = loss.sum()
            elif self.reduction == 'mean':
                loss = loss.mean()

            return loss

    class SupervisedTransformer(pl.LightningModule):
        """A Lighting Module that holds the graph transfomer along with
        training options."""

        def __init__(self):
            super().__init__()

            self.encoder = GraphTransformerEncoder(
                n_layers=args.n_layers,
                n_heads=args.n_heads,
                n_hid=args.emb_dim,
                summarizer=Set2Set(args.emb_dim, processing_steps=3),
                pretrained_weights=prot_embs)
            self.criterion = SmoothBCEwLogits(smoothing=0.001)
            self.decoder = torch.nn.Sequential(
                torch.nn.Linear(2 * args.emb_dim, 4 * args.emb_dim),
                torch.nn.ReLU(), torch.nn.Linear(4 * args.emb_dim, 255))

            self.new_lr = 1e-1
            self.accuracy = pl.metrics.Accuracy()

        def forward(self, data):
            x = self.encoder(data)
            x = self.encoder.summarize(x, data.batch)
            x = self.decoder(x)
            return x

        def configure_optimizers(self):
            opt = torch.optim.SGD(self.parameters(), lr=self.new_lr)
            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                opt, 'min', verbose=True, patience=2)

            opt_dict = {
                'optimizer': opt,
                'lr_scheduler': scheduler,
                'monitor': 'val_log_loss'
            }
            return opt_dict

        def training_step(self, train_batch, batch_idx):

            x = self.forward(train_batch)
            y_train = train_batch.y.view(args.batch_size, -1)

            t_loss = self.criterion(x, y_train)
            t_acc = self.accuracy(x, y_train)

            self.log('train_loss', t_loss, on_step=True, prog_bar=True)
            self.log('train_acc',
                     t_acc,
                     on_step=True,
                     on_epoch=False,
                     prog_bar=True)

            return t_loss

        def validation_step(self, val_batch, batch_idx):

            x = self.forward(val_batch)
            y_val = val_batch.y.view(args.batch_size, -1)

            v_loss = self.criterion(x, y_val)
            v_acc = self.accuracy(x, y_val)
            self.log('val_loss', v_loss, on_step=True, prog_bar=True)
            self.log('val_acc', v_acc, on_step=True, prog_bar=True)

            return v_loss

        def backward(self, loss, optimizer, optimizer_idx):

            loss.backward(retain_graph=True)

    model = SupervisedTransformer().to(dev)
    es = EarlyStopping(monitor='val_loss', patience=4)
    trainer = pl.Trainer(gpus=1, auto_lr_find='new_lr', callbacks=[es])

    # trainer.tune(model, train_loader, val_loader)

    # Fit Model
    trainer.fit(model, train_loader, val_loader)
