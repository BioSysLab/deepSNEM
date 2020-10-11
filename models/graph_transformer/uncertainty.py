import torch
import torch.nn as nn
from torch.functional import F

dev = torch.device('cuda:0')

def kl_divergence(alpha, num_classes, device=None):
    if not device:
        device = get_device()
    beta = torch.ones([1, num_classes], dtype=torch.float32, device=device)
    S_alpha = torch.sum(alpha, dim=1, keepdim=True)
    S_beta = torch.sum(beta, dim=1, keepdim=True)
    lnB = torch.lgamma(S_alpha) - \
        torch.sum(torch.lgamma(alpha), dim=1, keepdim=True)
    lnB_uni = torch.sum(torch.lgamma(beta), dim=1,
                        keepdim=True) - torch.lgamma(S_beta)

    dg0 = torch.digamma(S_alpha)
    dg1 = torch.digamma(alpha)

    kl = torch.sum((alpha - beta) * (dg1 - dg0), dim=1,
                   keepdim=True) + lnB + lnB_uni
    return kl

class HeteroFF(nn.Module):
    def ___init__(self, hidden_channels, N_CLASSES):
        super(HeteroFF, self).__init__()

        self.fc1 = nn.Linear(2*hidden_channels, hidden_channels)
        self.fc_out = nn.Linear(hidden_channels, N_CLASSES)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.log_softmax(self.fc_out(x))
        return x