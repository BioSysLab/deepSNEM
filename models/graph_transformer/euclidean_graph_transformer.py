import numpy as np
import math

import torch
import torch.nn as nn
from torch.functional import F

import torch_geometric.transforms as T
from torch_geometric.utils import to_dense_adj

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
        
        #self.gdc = T.GDC(self_loop_weight=1, diffusion_kwargs={'alpha':0.15, 'method': 'ppr'})
        self.conv = nn.Conv2d(2, 1, kernel_size=1, stride=1, bias=True)
        
    def forward(self, x: torch.Tensor, data: torch.Tensor) -> torch.Tensor:
        q, k, v = self.create_qkv(x)

        #data = self.gdc(data)
        adj_w_sign = to_dense_adj(data.edge_index, None, edge_attr=data.edge_attr)
        #adj_w_sign[adj_w_sign < 0.0009] = -100000000
        
        #Scaled Dot Product Attention of Heads
        self.attn_logits = torch.matmul(q, k.transpose(-2,-1)) / math.sqrt(k.shape[-1])
        #self.attn_logits = torch.add(self.attn_logits, self.conv(adj_w_sign.view(1,2,adj_w_sign.size(-2),-1)))
        
        self.attn_logits = F.softmax(self.attn_logits, dim=-1)
        #self.attn_logits = self.attn_logits * adj_w_sign

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

    def forward(self, x, data):
        
        # Perform Multi Head Attention
        attn_out = self.mhga(x, data)

        # Add and Norm
        out1 = self.lnorm1(torch.add(x,attn_out))

        # Feedforward then add and norm
        ff = self.pwff(out1)
        ff = self.lnorm2(torch.add(ff, out1))
        
        ff = F.leaky_relu(ff, negative_slope=0.2)

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
        return self.act_emb(data.acts.float()), self.pos_enc(data.degree.view(-1).float())

#-----------------------------------------------------------------------------------------------    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

class GraphTransformerEncoder(nn.Module):
    def __init__(self, n_layers, n_heads, n_hid, pretrained_weights=None):
        super(GraphTransformerEncoder, self).__init__()
        self.pretrained_weights = pretrained_weights
        if pretrained_weights is not None:
            self.n_prots, self.in_channels = pretrained_weights.shape
        else:
            self.n_prots = 919
            self.in_channels = emb_dim

        # Create Embedding Layer and initialize using pretrained weights
        self.emb_layer = nn.Embedding(self.n_prots, self.in_channels, sparse=True)
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
        act, pos = self.pe(data)
        x = torch.add(x, act)
        x = torch.add(x, pos.squeeze())

        #x = self.transformer_1(x, data)
        for t in self.transformers:
            x = t(x, data)

        print(self.emb_layer(global_idx)[:5])
        print(x[:5])

        return x

#-----------------------------------------------------------------------------------------------    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------