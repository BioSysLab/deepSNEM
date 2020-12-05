import torch
import torch.nn as nn
from torch.functional import F


class PathTransformer(nn.Module):
    def __init__(self, n_layers, n_heads, pretrained_prots):
        super(PathTransformer, self).__init__()

        self.pretrained_prots = pretrained_prots
        if pretrained_prots is not None:
            n_prots, in_channels = pretrained_prots.shape
        else:
            n_prots = 919
            in_channels = 512

        self.emb_layer = nn.Embedding(n_prots, in_channels, sparse=True)

        self.transformer = nn.TransformerEncoderLayer(in_channels, n_heads,
                                                      activation='gelu')

    def forward(self, x):

        x = self.transformer(src=x, src_mask=None)

        return x
