import torch
import torch.nn as nn
from torch.functional import F

class DeepSNEM(nn.Module):
    """
    Abstract signaling networks autoencoder class.
    Input: 
        encoder: (nn.Module, specifying the encoder architecture, returns the network node embeddings)
        decoder: (nn.Module, receives the node embeddings and performs either reconstruction or link prediction)
    """
    def __init__(self, encoder, decoder):
        super(DeepSNEM, self).__init__()

        assert isinstance(encoder, nn.Module) and isinstance(decoder, nn.Module)

        self.encoder = encoder
        self.decoder = decoder

    def encode(self, *args, **kwargs):
        return self.encoder(*args, **kwargs)

    def decode(self, *args, **kwargs):
        return self.decoder(*args, **kwargs)

    def forward(self, *args, **kwargs):
        return self.decode(self.encode(*args, **kwargs))
    
class LinearDecoder(nn.Module):
    def __init__(self, emb_dim, original_dim):
        super(LinearDecoder, self).__init__()
        #self.fc_r = nn.Linear(emb_dim, emb_dim, bias=True)
        self.fc_out = nn.Linear(emb_dim, original_dim)
        self.init_weights()
        
    def init_weights(self):
        initrange = 0.1 
        #self.fc_r.bias.data.zero_()
        #self.fc_r.weight.data.uniform_(-initrange, initrange)
        
        self.fc_out.bias.data.zero_()
        self.fc_out.weight.data.uniform_(-initrange, initrange)
        
    def forward(self, z):
        return self.fc_out(z)

class FermiDiracDecoder(nn.Module):
    """Fermi Dirac to compute edge probabilities based on distances."""

    def __init__(self, r, t):
        super(FermiDiracDecoder, self).__init__()
        self.r = r
        self.t = t

    def forward(self, dist):
        probs = 1. / (torch.exp((dist - self.r) / self.t) + 1.0)
        return probs