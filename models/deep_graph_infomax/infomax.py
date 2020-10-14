import torch
import numpy as np
from torch.nn import Parameter
import torch.nn.functional as F
import math

from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder

from .inits import reset, uniform

EPS = 1e-15

def log_sum_exp(x, axis=None):
    """Log sum exp function
    Args:
        x: Input.
        axis: Axis over which to perform sum.
    Returns:
        torch.Tensor: log sum exp
    """
    x_max = torch.max(x, axis)[0]
    y = torch.log((torch.exp(x - x_max)).sum(axis)) + x_max
    return y

class LocalDiscriminator(torch.nn.Module):
    def __init__(self,input_dim, dim):
        super().__init__()
        self.block = torch.nn.Sequential(
            torch.nn.Linear(input_dim, input_dim),
            torch.nn.ReLU(),
            torch.nn.Linear(input_dim, input_dim),
            torch.nn.ReLU(),
            torch.nn.Linear(input_dim, dim),
            torch.nn.ReLU()
        )
        self.linear_shortcut = torch.nn.Linear(input_dim, dim)

    def forward(self, x):
        return self.block(x) + self.linear_shortcut(x)

class PriorDiscriminator(torch.nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.l0 = torch.nn.Linear(input_dim, input_dim)
        self.l1 = torch.nn.Linear(input_dim, input_dim)
        self.l2 = torch.nn.Linear(input_dim, 1)

    def forward(self, x):
        h = F.relu(self.l0(x))
        h = F.relu(self.l1(h))
        return torch.sigmoid(self.l2(h))

class SNInfomax(torch.nn.Module):

    def __init__(self, hidden_channels, encoder, summary, semi=False, beta=1.0):
        super(SNInfomax, self).__init__()
        self.hidden_channels = hidden_channels
        self.encoder = encoder
        self.summary = summary
        self.semi = semi

        self.local_d = LocalDiscriminator(hidden_channels, hidden_channels)
        self.global_d = LocalDiscriminator(2*hidden_channels, hidden_channels)

        self.prior_d = PriorDiscriminator(2*hidden_channels)
        self.beta = beta

        self.reset_parameters()
        self.init_emb()

    def reset_parameters(self):
        reset(self.encoder)
        reset(self.summary)

    def init_emb(self):
      initrange = -1.5 / self.hidden_channels
      for m in self.modules():
          if isinstance(m, torch.nn.Linear):
              torch.nn.init.xavier_uniform_(m.weight.data)
              if m.bias is not None:
                  m.bias.data.fill_(0.0)

    def forward(self, data):
        #neg_z = self.encoder.corrupt_forward(data)
        z = self.encoder(data)
        summary = self.summary(z, data.batch) if hasattr(data, 'batch') else self.summary(pos_z)
        mask = torch.matmul(data.node_sig, data.sig.view(data.num_graphs, -1).t()).bool()
        self.pos_mask = mask.float()
        self.neg_mask = (~mask).float()

        z_un = self.local_d(z)
        summary_un = self.global_d(summary)

        res_un = torch.matmul(z_un, summary_un.t()) 

        return res_un, summary

    def loss(self, res_un, summary):
        r"""Computes the mutual information maximization objective."""
        log_2 = math.log(2.)
        p_samples = res_un * self.pos_mask
        q_samples = res_un * self.neg_mask

        Ep = log_2 - F.softplus(- p_samples)
        Eq = F.softplus(-q_samples) + q_samples - log_2
        
        Ep = (Ep * self.pos_mask).sum() / self.pos_mask.sum()
        Eq = (Eq * self.neg_mask).sum() / self.neg_mask.sum()
        LOCAL = Eq - Ep

        prior = torch.rand_like(summary)

        term_a = torch.log(self.prior_d(prior)).mean()
        term_b = torch.log(1.0 - self.prior_d(summary)).mean()
        PRIOR = -(term_a + term_b) * self.beta

        return LOCAL, PRIOR

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.hidden_channels)
