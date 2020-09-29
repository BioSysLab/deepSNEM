import torch
import torch.nn as nn
from torch.functional import F

dev = torch.device('cuda:0')

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
        #self.fc_r = nn.Linear(emb_dim, emb_dim, bias=False)
        self.fc_out = nn.Linear(emb_dim, original_dim, bias=False)
        self.init_weights()
        
    def init_weights(self):
        initrange = 0.1 
        #self.fc_r.bias.data.zero_()
        #self.fc_r.weight.data.uniform_(-initrange, initrange)
        
        #self.fc_out.bias.data.zero_()
        self.fc_out.weight.data.uniform_(-initrange, initrange)
        
    def forward(self, z):
        return self.fc_out(z)

class FermiDiracDecoder(nn.Module):
    """Fermi Dirac to compute edge probabilities based on distances."""

    def __init__(self, t):
        super(FermiDiracDecoder, self).__init__()
        #self.r = nn.Parameter(torch.tensor(6.)).cuda()
        self.r = 5.0
        self.t = t

    def forward(self, dist):
        probs = 1. / (torch.exp((dist - self.r) / self.t) + 1.0)
        return probs

class InnerProductDecoder(torch.nn.Module):
    def forward(self, z, edge_index, sigmoid=True):
        r"""Decodes the latent variables :obj:`z` into edge probabilities for
        the given node-pairs :obj:`edge_index`.
        Args:
            z (Tensor): The latent space :math:`\mathbf{Z}`.
            sigmoid (bool, optional): If set to :obj:`False`, does not apply
                the logistic sigmoid function to the output.
                (default: :obj:`True`)
        """
        value = (z[edge_index[0]] * z[edge_index[1]]).sum(dim=1)
        return torch.sigmoid(value) if sigmoid else value

    def forward_all(self, z, sigmoid=True):
        r"""Decodes the latent variables :obj:`z` into a probabilistic dense
        adjacency matrix.
        Args:
            z (Tensor): The latent space :math:`\mathbf{Z}`.
            sigmoid (bool, optional): If set to :obj:`False`, does not apply
                the logistic sigmoid function to the output.
                (default: :obj:`True`)
        """
        adj = torch.matmul(z, z.t())
        return torch.sigmoid(adj) if sigmoid else adj


def pdist(sample_1, sample_2, norm=2, eps=1e-5):
    r"""Compute the matrix of all squared pairwise distances.
    Arguments
    ---------
    sample_1 : torch.Tensor or Variable
        The first sample, should be of shape ``(n_1, d)``.
    sample_2 : torch.Tensor or Variable
        The second sample, should be of shape ``(n_2, d)``.
    norm : float
        The l_p norm to be used.
    Returns
    -------
    torch.Tensor or Variable
        Matrix of shape (n_1, n_2). The [i, j]-th entry is equal to
        ``|| sample_1[i, :] - sample_2[j, :] ||_p``."""
    n_1, n_2 = sample_1.size(0), sample_2.size(0)
    norm = float(norm)
    if norm == 2.:
        norms_1 = torch.sum(sample_1**2, dim=1, keepdim=True)
        norms_2 = torch.sum(sample_2**2, dim=1, keepdim=True)
        norms = (norms_1.expand(n_1, n_2) +
                 norms_2.transpose(0, 1).expand(n_1, n_2))
        distances_squared = norms - 2 * sample_1.mm(sample_2.t())
        return torch.sqrt(eps + torch.abs(distances_squared))
    else:
        dim = sample_1.size(1)
        expanded_1 = sample_1.unsqueeze(1).expand(n_1, n_2, dim)
        expanded_2 = sample_2.unsqueeze(0).expand(n_1, n_2, dim)
        differences = torch.abs(expanded_1 - expanded_2) ** norm
        inner = torch.sum(differences, dim=2, keepdim=False)
        return (eps + inner) ** (1. / norm)

def get_link_labels(pos_edge_index, neg_edge_index):
    link_labels = torch.zeros(pos_edge_index.size(1) +
                              neg_edge_index.size(1)).float().to(device)
    link_labels[:pos_edge_index.size(1)] = 1.
    return link_labels

def get_edge_dists(x, edges):
    probs = []
    for node_id, child in enumerate(edges[0]):
        if child is not None:
            x_child = x[child]
            probs.append(pdist(x[node_id].unsqueeze(0), x_child)[0])

    return probs

def reshape_probs(probs):
    probs_all = probs[0]
    for prob in probs[1:]:
        probs_all = torch.cat((prob, probs_all), dim=-1)

    return probs_all