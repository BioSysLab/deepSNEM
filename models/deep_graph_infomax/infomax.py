import torch
from torch.nn import Parameter
from sklearn.linear_model import LogisticRegression

from .inits import reset, uniform

EPS = 1e-15


class DeepGraphInfomax(torch.nn.Module):

    def __init__(self, hidden_channels, encoder, summary):
        super(DeepGraphInfomax, self).__init__()
        self.hidden_channels = hidden_channels
        self.encoder = encoder
        self.summary = summary

        self.weight = Parameter(torch.Tensor(hidden_channels, hidden_channels))

        self.reset_parameters()

    def reset_parameters(self):
        reset(self.encoder)
        reset(self.summary)
        uniform(self.hidden_channels, self.weight)

    def forward(self, data):
        """Returns the latent space for the input arguments, their
        corruptions and their summary representation."""
        pos_z = self.encoder(data)

        neg_z_x = self.encoder.corrupt_forward_x(data)
        neg_z_edges = self.encoder.corrupt_forward_edges(data)

        summary = self.summary(pos_z)

        return pos_z, neg_z_x, neg_z_edges, summary

    def discriminate(self, z, summary, sigmoid=True):
        r"""Given the patch-summary pair :obj:`z` and :obj:`summary`, computes
        the probability scores assigned to this patch-summary pair.
        Args:
            z (Tensor): The latent space.
            sigmoid (bool, optional): If set to :obj:`False`, does not apply
                the logistic sigmoid function to the output.
                (default: :obj:`True`)
        """
        value = torch.matmul(z, torch.matmul(self.weight, summary))
        return torch.sigmoid(value) if sigmoid else value

    def loss(self, pos_z, neg_z_x, neg_z_edges, summary):
        r"""Computes the mutual information maximization objective."""
        pos_loss = -torch.log(
            self.discriminate(pos_z, summary, sigmoid=True) + 1e-5).mean()
        pos_loss_dups = -torch.log(
            self.discriminate(pos_z, summary, sigmoid=True) + 1e-5).mean()
        neg_loss_x = -torch.log(
            1 - self.discriminate(neg_z_x, summary, sigmoid=True) + 1e-5).mean()
        neg_loss_edges = -torch.log(
            1 - self.discriminate(neg_z_edges, summary, sigmoid=True) + 1e-5).mean()

        return pos_loss + neg_loss_x + neg_loss_edges

    def test(self, train_z, train_y, test_z, test_y, solver='lbfgs',
             multi_class='auto', *args, **kwargs):
        r"""Evaluates latent space quality via a logistic regression downstream
        task."""
        clf = LogisticRegression(solver=solver, multi_class=multi_class, *args,
                                 **kwargs).fit(train_z.detach().cpu().numpy(),
                                               train_y.detach().cpu().numpy())
        return clf.score(test_z.detach().cpu().numpy(),
                         test_y.detach().cpu().numpy())

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.hidden_channels)