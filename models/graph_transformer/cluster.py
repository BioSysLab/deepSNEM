import torch
import torch.nn as nn

class ClusterAssignment(nn.Module):
    def __init__(self, emb_dim : int, n_clusters : int, n_heads : int, alpha: int = 1.0,
                 cluster_centers : Optional[torch.Tensor] = None) -> None:
        super(ClusterAssignment, self).__init__()
        
        self.n_clusters = n_clusters
        self.n_heads = n_heads
        self.alpha = alpha
        
        if cluster_centers is None:
            initial_cluster_centers = torch.zeros(n_heads, n_clusters, emb_dim, dtype=torch.float) 
            nn.init.uniform_(initial_cluster_centers)
        else:
            initial_cluster_centers = cluster_centers
            
        self.cluster_centers = nn.Parameter(initial_cluster_centers)

    def forward(self, queries : torch.Tensor) -> torch.Tensor:
        
        # Compute soft assignment
        norm_squared = torch.sum((queries - self.cluster_centers[:,:,None]).pow(2), 2)
        numerator = 1.0 / (1.0 + (norm_squared / self.alpha))
        power = float(self.alpha + 1) / 2
        numerator = numerator.pow(power)
        
        return numerator / torch.sum(numerator, dim=1, keepdim=True)