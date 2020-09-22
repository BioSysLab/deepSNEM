import numpy as np
import pandas as pd
import sys
sys.path.append('..')

import torch
from torch.functional import F
import torch.nn as nn

from torch_geometric.data import Data, DataLoader, Dataset
from torch_geometric.utils import from_networkx, to_networkx, degree
from torch_geometric.nn import GATConv, GCNConv, global_add_pool, PNAConv, BatchNorm, CGConv, global_max_pool

class PostEncoding(nn.Module):
    def __init__(self, emb_dim):
        super(PostEncoding,self).__init__()
        self.emb_dim = emb_dim
        
        self.act_emb = nn.Linear(2, emb_dim)
        
    def forward(self, data):
        return self.act_emb(data.acts.float())

class Net(torch.nn.Module):
    def __init__(self, n_classes, deg, pretrained_weights=None, train_embs=True):
        super(Net, self).__init__()
        self.pretrained_weights = pretrained_weights
        if pretrained_weights is not None:
            self.n_prots, self.in_channels = pretrained_weights.shape
        else:
            self.n_prots = 919
            self.in_channels = 512

        self.node_emb = nn.Embedding(self.n_prots, self.in_channels, sparse=True)
        self.node_emb.weight.requires_grad = train_embs
        
        self.edge_emb = nn.Embedding(2, 50)
        self.pe = PostEncoding(self.in_channels)

        aggregators = ['mean', 'min', 'max', 'std']
        scalers = ['identity', 'amplification', 'attenuation']

        self.convs = nn.ModuleList()
        self.batch_norms = nn.ModuleList()
        for _ in range(2):
            conv = PNAConv(in_channels=self.in_channels, out_channels=self.in_channels,
                           aggregators=aggregators, scalers=scalers, deg=deg,
                           edge_dim=50, towers=4, pre_layers=1, post_layers=1,
                           divide_input=False)
            self.convs.append(conv)
            self.batch_norms.append(BatchNorm(self.in_channels))
            
        self.fc1 = nn.Linear(self.in_channels, 2 * self.in_channels)
        self.fc_out = nn.Linear(2 * self.in_channels, n_classes)
        self.act = nn.PReLU()    
        
        self.init_weights()
        
    def init_weights(self):
        if self.pretrained_weights is not None:
            initrange = 0.1
            self.node_emb.weight.data.copy_(torch.from_numpy(self.pretrained_weights))
            self.edge_emb.weight.data.uniform_(-initrange, initrange)
        else:
            initrange = 0.1
            emb_layer.weight.data.uniform_(-initrange, initrange)
            
        for m in self.modules():
            if isinstance(m, nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight.data)
                if m.bias is not None:
                      m.bias.data.fill_(0.0)

    def forward(self, data):
        x = self.node_emb(data.global_idx)
        edge_attr = self.edge_emb(data.sign.long())
        act = self.pe(data)
        x = torch.add(x, act)

        for conv, batch_norm in zip(self.convs, self.batch_norms):
            x = F.relu(batch_norm(conv(x, data.edge_index, edge_attr)))

        x = global_add_pool(x, data.batch)
        x = self.act(self.fc1(x))
        return F.log_softmax(self.fc_out(x), dim=-1)

class BigNet(nn.Module):
    def __init__(self, out_channels, pretrained_weights=None):
        super(BigNet, self).__init__()
        self.pretrained_weights = pretrained_weights
        
        assert pretrained_weights is not None
        
        n_prots, in_channels = pretrained_weights.shape
            
        self.emb_layer = nn.Embedding(n_prots, in_channels, sparse=True)
        self.emb_layer.weight.requires_grad=True
        self.pe = PostEncoding(in_channels)
        
        self.conv1 = CGConv(in_channels, dim=2)
        self.conv2 = GATConv(in_channels, in_channels)
        self.conv3 = GATConv(in_channels, in_channels)
        self.conv4 = GATConv(in_channels, in_channels)
        self.conv5 = GATConv(in_channels, in_channels)
        
        self.bn1 = nn.BatchNorm1d(in_channels)
        self.bn2 = nn.BatchNorm1d(in_channels)
        self.bn3 = nn.BatchNorm1d(in_channels)
        self.bn4 = nn.BatchNorm1d(in_channels)
        self.bn5 = nn.BatchNorm1d(in_channels)
        
        self.act1 = nn.PReLU()
        self.act2 = nn.PReLU()
        self.act3 = nn.PReLU()
        self.act4 = nn.PReLU()
        self.act5 = nn.PReLU()
            
        self.fc1 = nn.Linear(in_channels, 2*in_channels)
        self.fc2 = nn.Linear(2*in_channels , out_channels)
        
    def init_weights(self, emb_layer):
        emb_layer.weight.data.copy_(torch.from_numpy(self.pretrained_weights))
        for m in self.modules():
            if isinstance(m, nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight.data)
                if m.bias is not None:
                      m.bias.data.fill_(0.0)
        
    def forward(self, data):
        x = self.emb_layer(data.global_idx)
        act = self.pe(data)
        x = torch.add(x, act)
        
        x_1 = self.conv1(x, data.edge_index, data.sign)
        x_1 = self.bn1(x_1)
        x_1 = self.act1(torch.add(x_1, x))
        x_1 = F.dropout(x_1, training=self.training, p=0.1)
        
        x_2 = self.conv2(x_1, data.edge_index)
        x_2 = self.bn2(x_2)
        x_2 = self.act2(torch.add(x_2, x_1))
        x_2 = F.dropout(x_2, training=self.training, p=0.1)
        
        x_3 = self.conv3(x_2, data.edge_index)
        x_3 = self.bn3(x_3)
        x_3 = self.act3(torch.add(x_3, x_2))
        x_3 = F.dropout(x_3, training=self.training, p=0.1)
        
        x_4 = self.conv4(x_3, data.edge_index)
        x_4 = self.bn4(x_4)
        x_4 = self.act4(torch.add(x_4, x_3))
        x_4 = F.dropout(x_4, training=self.training, p=0.1)
        
        x_5 = self.conv5(x_4, data.edge_index)
        x_5 = self.bn5(x_5)
        x_5 = self.act5(torch.add(x_5, x_4))
        x_5 = F.dropout(x_5, training=self.training, p=0.1)
        
        x = global_add_pool(x, data.batch).squeeze()
        x = F.relu(self.fc1(x))
        x = F.log_softmax(self.fc2(x), dim=-1)
        return x    
    