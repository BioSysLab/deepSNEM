import numpy as np

import torch
import torch.nn as nn
from torch.functional import F

class PathTransformer(nn.Module):
    def __init__(self, n_layers, n_heads, pretrained_prots):
        super(PathTransformer, self).__init__()

        self.pretrained_prots = pretrained_prots
        n_prots, in_channels = pretrained_prots.shape

        self.emb_layer = nn.Embedding(n_prots, in_channels, sparse=True)

        self.transformer = nn.TransformerEncoderLayer(in_channels, n_heads,
                                                      activation='gelu')

    def forward(self, x):

        x = self.transformer(src=x, src_mask=None)

        return x
