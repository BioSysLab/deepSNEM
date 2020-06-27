import numpy as np
import math

import torch
import torch.nn as nn
from torch.functional import F

import torch_geometric.transforms as T
from torch_geometric.utils import to_dense_adj
from torch_geometric.nn import Set2Set

dev = torch.device('cuda:0')
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

class MultiHeadGraphAttention(nn.Module):
    def __init__(self, in_channels:int, n_heads: int) -> None:
        super(MultiHeadGraphAttention, self).__init__()
        self.n_heads = n_heads
        self.depth = in_channels // n_heads
        self.in_channels = in_channels

        self.wq = nn.Linear(in_channels, in_channels, bias=False)
        self.wk = nn.Linear(in_channels, in_channels, bias=False)
        self.wv = nn.Linear(in_channels, in_channels, bias=False)

        self.lin_out = nn.Linear(in_channels, in_channels, bias=False)
        
        self.conv = nn.Conv2d(2, 1, kernel_size=1, stride=1, bias=True)
        nn.init.uniform_(self.conv.weight.data)
        nn.init.zeros_(self.conv.bias)

        self.param = nn.Parameter(torch.tensor(10.)).clamp_min(1.).to(dev)
        
    def forward(self, x: torch.Tensor, data: torch.Tensor) -> torch.Tensor:
        q, k, v = self.create_qkv(x)

        adj_w_sign = to_dense_adj(data.edge_index, None, edge_attr=data.sign)
        
        #adj_after_ppr = to_dense_adj(data_2.edge_index, None, edge_attr=data_2.edge_attr)
        
        #Scaled Dot Product Attention of Heads
        self.attn_logits = ((q**2).sum(-1).view(4,-1,1) + (k**2).sum(-1).view(4,1,-1))  - 2.* torch.bmm(q, k.transpose(-2,-1))
        self.attn_logits = self.attn_logits/math.sqrt(self.in_channels)

        conv_s = adj_w_sign.view(1,2,adj_w_sign.size(-2),-1)
        #conv_ppr = adj_after_ppr.view(1,1,adj_after_ppr.size(-2), -1)
        #conv_f = torch.cat([conv_s, conv_ppr], dim=1)

        self.attn_logits = torch.add(self.attn_logits, self.param * self.conv(conv_s))
        
        self.attn_logits = F.softmax(self.attn_logits, dim=-1)

        attn_out = torch.matmul(self.attn_logits, v)

        attn_out = attn_out.reshape(x.size(0), self.in_channels)
        attn_out = self.lin_out(attn_out)
        attn_out = F.dropout(attn_out, training=self.training, p=0.1)

        return attn_out
        
    def split_heads(self, mat : torch.Tensor) -> torch.Tensor:
        return mat.view(self.n_heads, -1, self.depth)
    
    def create_qkv(self, x):

        q = self.split_heads(self.wq(x))
        k = self.split_heads(self.wk(x))
        v = self.split_heads(self.wv(x))
        
        return q, k, v

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

class FeedForward(nn.Module):
    def __init__(self, in_channels, emb_dim):
        super(FeedForward, self).__init__()

        self.ff1 = nn.Linear(in_channels, 4096)
        self.ff2 = nn.Linear(4096, emb_dim)
        
        self.init_weights(FeedForward)

    def init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.xavier_normal_(m.weight)
            nn.init.zeros_(m.bias)

    def forward(self, x):
        # Feedforward then add and norm
        x = F.relu(self.ff1(x))
        x = F.dropout(self.ff2(x), training=self.training, p=0.05)
        return x # Node embeddings before add and norm

#-----------------------------------------------------------------------------------------------    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

class GraphTransformerEncoderLayer(nn.Module):
    def __init__(self, in_channels, n_heads, n_hid):
        super(GraphTransformerEncoderLayer, self).__init__()
        self.n_heads = n_heads
        self.in_channels = in_channels
        self.n_hid = n_hid
        
        if in_channels % n_heads != 0:
            raise AssertionError('Number of heads and input_channels need to be perfectly divided.')

        # Multi Head Attention Layer
        self.mhga = MultiHeadGraphAttention(in_channels, n_heads)

        self.lnorm1 = nn.LayerNorm(in_channels)
        self.lnorm2 = nn.LayerNorm(n_hid)

        self.pwff = FeedForward(in_channels, n_hid)

        self.emb_act = nn.PReLU(self.in_channels)

    def forward(self, x, data):
        
        # Perform Multi Head Attention
        attn_out = self.mhga(x, data)

        # Add and Norm
        out1 = self.lnorm1(torch.add(x,attn_out))

        # Feedforward then add and norm
        ff = self.pwff(out1)
        ff = self.lnorm2(torch.add(ff, out1))
        
        ff = self.emb_act(ff)

        return ff

#-----------------------------------------------------------------------------------------------    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

class PositionalEmbedding(nn.Module):
    def __init__(self, demb):
        super(PositionalEmbedding, self).__init__()

        self.demb = demb

        self.inv_freq = (1 / ((10000) ** (torch.arange(0.0, demb, 2.0) / demb))).cuda()

    def forward(self, pos_seq):
        sinusoid_inp = torch.ger(pos_seq, self.inv_freq)
        pos_emb = torch.cat([sinusoid_inp.sin(), sinusoid_inp.cos()], dim=-1)
        return pos_emb[:,None,:]

class PostEncoding(nn.Module):
    def __init__(self, emb_dim):
        super(PostEncoding,self).__init__()
        self.emb_dim = emb_dim
        
        self.pos_enc = PositionalEmbedding(emb_dim)
        #self.degree_enc = torch.nn.Embedding(21, emb_dim, padding_idx=0)
        #self.degree_enc.weight.data = position_encoding_init(data.degree, emb_dim)
        
        self.act_emb = nn.Linear(2, emb_dim)
        
    def forward(self, data):
        return self.act_emb(data.acts.float())

#-----------------------------------------------------------------------------------------------    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

class GraphTransformerEncoder(nn.Module):
    def __init__(self, n_layers, n_heads, n_hid, summarizer=None, pretrained_weights=None, train_embs = True):
        super(GraphTransformerEncoder, self).__init__()
        self.pretrained_weights = pretrained_weights
        self.summarizer = summarizer
        if pretrained_weights is not None:
            self.n_prots, self.in_channels = pretrained_weights.shape
        else:
            self.n_prots = 919
            self.in_channels = emb_dim

        # Create Embedding Layer and initialize using pretrained weights
        self.emb_layer = nn.Embedding(self.n_prots, self.in_channels, sparse=True)
        self.emb_layer.weight.requires_grad = train_embs

        if train_embs:
            self.init_weights(self.emb_layer)
        
        self.pe = PostEncoding(self.in_channels)
    
        self.transformer_1 = GraphTransformerEncoderLayer(self.in_channels, n_heads, n_hid)
        self.transformers = nn.ModuleList([GraphTransformerEncoderLayer(self.in_channels, n_heads, n_hid) for _ in range(n_layers)])

    def init_weights(self, emb_layer):
        if self.pretrained_weights is not None:
            emb_layer.weight.data.copy_(torch.from_numpy(self.pretrained_weights))
        else:
            initrange = 0.1
            emb_layer.weight.data.uniform_(-initrange, initrange)

    def forward(self, data):
        global_idx = data.global_idx

        x = self.emb_layer(global_idx)
        act = self.pe(data)
        x = torch.add(x, act)

        #x = self.transformer_1(x, data)
        for t in self.transformers:
            x = t(x, data)

        #print(self.emb_layer(global_idx)[:5])
        #print(x[:5])

        return x

    def corrupt_forward_edges(self, data):
        global_idx = data.global_idx

        x = self.emb_layer(global_idx)
        act = self.pe(data)
        x = torch.add(x, act)

        data.edge_index[1] = data.edge_index[1][torch.randperm(data.edge_index[1].size(0))]

        for t in self.transformers:
            x = t(x, data)

        return x

    def corrupt_forward_x(self, data):
        global_idx = data.global_idx

        x = self.emb_layer(global_idx)
        act = self.pe(data)

        x = torch.add(x, act)
        x = x[torch.randperm(x.size(0))]

        #data.edge_index[1] = data.edge_index[1][torch.randperm(data.edge_index[1].size(0))]

        for t in self.transformers:
            x = t(x, data)

        return x

    def summarize(self, x, batch):
        summary = self.summarizer(x, batch)
        return summary
#-----------------------------------------------------------------------------------------------    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------