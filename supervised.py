import numpy as np
import pandas as pd
import sys
from tqdm.auto import tqdm
import networkx as nx
sys.path.append('..')

import torch
from torch.functional import F
import torch.nn as nn

from torch_geometric.data import Data, DataLoader, Dataset
from torch_geometric.utils import from_networkx, to_networkx, degree
from torch_geometric.nn import GATConv, GCNConv, global_add_pool, PNAConv, BatchNorm, CGConv, global_max_pool
from torch_geometric.utils.metric import accuracy, precision, f1_score
import torch_geometric.transforms as T

from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder, PostEncoding
from models.graph_pna.model import Net, BigNet

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import roc_auc_score, accuracy_score, average_precision_score, precision_score

from utils.data_gen import load_prot_embs, to_categorical, wcsv2graph, SNLDataset, deg_distr
torch.backends.cudnn.benchmark = True

dev = torch.device('cuda:0')

N_CLASSES = 255
EMB_DIM = 512
EPOCHS = 30

class FF(nn.Module):
    def __init__(self, encoder, in_channels, n_classes):
        super(FF, self).__init__()
        
        self.encoder = encoder
        self.fc_out = nn.Linear(in_channels, n_classes)
        
    def forward(self, data):
        x = self.encoder(data)
        
        x = global_max_pool(x, data.batch)
        x = F.log_softmax(self.fc_out(x), dim=1)
        return x

class LabelSmoothingLoss(nn.Module):
    def __init__(self, classes, smoothing=0.0, dim=-1):
        super(LabelSmoothingLoss, self).__init__()
        self.confidence = 1.0 - smoothing
        self.smoothing = smoothing
        self.cls = classes
        self.dim = dim

    def forward(self, pred, target):
        #pred = pred.log_softmax(dim=self.dim)
        with torch.no_grad():
            true_dist = torch.zeros_like(pred)
            true_dist.fill_(self.smoothing / (self.cls - 1))
            true_dist.scatter_(1, target.data.unsqueeze(1), self.confidence)
        return torch.mean(torch.sum(-true_dist * pred, dim=self.dim))

prot_embs, global_dict = load_prot_embs(512, norm=False)

labelled_ugraphs = pd.read_csv('../snac_data/graph_classification_all.csv')
weighted_df = pd.read_csv('../snac_data/file_info_weighted.csv')

val_set_1 = pd.read_csv('../snac_data/splits/val_set_1.csv')
val_set_2 = pd.read_csv('../snac_data/splits/val_set_2.csv')
val_set_3 = pd.read_csv('../snac_data/splits/val_set_3.csv')
val_set_4 = pd.read_csv('../snac_data/splits/val_set_4.csv')
test_set = pd.read_csv('../snac_data/splits/test_set.csv')

val_sets = [val_set_1, val_set_2, val_set_3, val_set_4]

ls = LabelSmoothingLoss(N_CLASSES, smoothing=0.1)

def train(epoch):
    model.train()
    total_train_acc = []
    for tb in train_data_iterator:
        tb = tb.to(dev)
        optimizer.zero_grad()
        
        pred = model(tb).view(tb.num_graphs, -1)
        y = tb.y.reshape(tb.num_graphs, -1)
        y_l = tb.label
        
        #loss = F.nll_loss(pred, y_l)
        loss = ls(pred, y_l)
        
        acc_t = accuracy(torch.argmax(pred, dim=1), y_l)
        roc_t = roc_auc_score(y.long().detach().cpu().numpy(), 
                                pred.detach().cpu().numpy(), average='samples')
        
        train_data_iterator.set_postfix(Epoch=epoch+1,
                                        loss ='%.4f' % float(loss.item()),
                                        acc = '%.4f' % float(acc_t),
                                        roc= '%.4f' % float(roc_t))
        
        loss.backward(retain_graph=True)
        torch.nn.utils.clip_grad_norm_(model.parameters(), 0.01)
        optimizer.step()
        total_train_acc.append(acc_t)

        del tb
    
    print(f'Total Training Accuracy: {np.mean(total_train_acc)}')
    
    val_data_iterator = tqdm(val_loader,leave=True,unit='batch',
                        postfix={'Epoch': epoch+1,'val_loss': '%.4f' % 0.0, 
                                 'acc': '%.4f' % 0.0,
                                 'roc': '%.4f' % 0.0})
    total_val_acc = []
    total_val_loss = 0.0
    with torch.no_grad():
        for tb in val_data_iterator:
            tb = tb.to(dev)

            pred = model(tb)
            y = tb.y.reshape(tb.num_graphs, -1)
            y_l = tb.label

            #loss = F.nll_loss(pred, y_l)
            val_loss = ls(pred, y_l)
            
            acc_v = accuracy(torch.argmax(pred, dim=1), y_l)
            roc_v = roc_auc_score(y.long().detach().cpu().numpy(), 
                                        pred.detach().cpu().numpy(), average='samples')

            val_data_iterator.set_postfix(Epoch=epoch+1,
                                            val_loss ='%.4f' % float(val_loss.item()),
                                            acc = '%.4f' % float(acc_v),
                                            roc= '%.4f' % float(roc_v))
            total_val_acc.append(acc_v)
            total_val_loss += val_loss
        
    print(f'Total Validation Accuracy: {np.mean(total_val_acc)}')

    scheduler.step(total_val_loss / len(val_df))
    return total_val_acc

for idx, val_set in enumerate(val_sets):
    print('-' * 100)
    print(f'VALIDATION DATASET {idx + 1}')
    print('-' * 100)

    BEST_VAL_ACC = 0.0

    #encoder = GraphTransformerEncoder(n_layers=2, n_heads=4, n_hid=512, 
                            #pretrained_weights=prot_embs, summarizer=None).to(dev)
    #model = FF(encoder, EMB_DIM, N_CLASSES).to(dev)

    usm = pd.DataFrame(labelled_ugraphs.groupby('sig_id').moa_v1.unique()).reset_index()
    usm_corr = np.array([np.array(i) for i in usm.moa_v1.to_numpy()]).reshape(-1)
    usm['moa_v1'] = usm_corr

    X_df = pd.merge(weighted_df, usm, on='sig_id')
    val_df =  pd.merge(X_df, val_set, on='sig_id')
    test_df = pd.merge(X_df, test_set, on='sig_id')

    for sig in tqdm(val_set.sig_id):
        X_df = X_df[X_df['sig_id'] != sig]
        
    for sig in tqdm(test_set.sig_id):
        X_df = X_df[X_df['sig_id'] != sig]

    X_train, y_train = X_df.files_weighted.to_numpy(), X_df.moa_v1.to_numpy()
    X_val, y_val = val_df.files_weighted.to_numpy(), val_df.moa_v1_x.to_numpy()
    X_test, y_test = test_df.files_weighted.to_numpy(), test_df.moa_v1_x.to_numpy()

    le = OneHotEncoder()
    y = np.concatenate([y_train, y_val, y_test])
    le = le.fit(y.reshape(-1,1))
    y_train = le.transform(y_train.reshape(len(y_train),-1)).toarray()
    y_val = le.transform(y_val.reshape(len(y_val),-1)).toarray()
    y_test = le.transform(y_test.reshape(len(y_test), -1)).toarray()

    train_data = SNLDataset(X_train, y_train, global_dict)
    val_data = SNLDataset(X_val, y_val, global_dict)
    test_data = SNLDataset(X_test, y_test, global_dict)

    #deg = deg_distr(train_data)
    #print(deg)

    model = BigNet(N_CLASSES, pretrained_weights=prot_embs).to(dev)
    lr = 5.0
    optimizer = torch.optim.SGD(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, factor=0.5, patience=2, verbose=True)

    train_loader = DataLoader(train_data, batch_size=64, num_workers=12, shuffle=True)
    val_loader = DataLoader(val_data, batch_size=16, num_workers=12, shuffle=True)
    test_loader = DataLoader(test_data, batch_size=1, num_workers=12)

    for epoch in range(EPOCHS):
        print('-' * 100)
        print('-' * 100)
        train_data_iterator = tqdm(train_loader,leave=True,unit='batch',
                            postfix={'Epoch': epoch+1,'loss': '%.4f' % 0.0,
                            'loss': '%.4f' % 0.0,
                            'acc': '%.4f' % 0.0,
                            'roc': '%.4f' % 0.0})
        total_val_acc = train(epoch)
        is_best = np.mean(total_val_acc) > BEST_VAL_ACC
        BEST_VAL_ACC = np.mean(total_val_acc) if is_best else BEST_VAL_ACC
        if is_best and epoch > 10:
            torch.save(model.state_dict(),
                       f'../snac_data/checks_conv/val_set_{idx+1}/moa_vs_{idx+1}_val_acc_{np.mean(total_val_acc)}.pt')
