import torch
import numpy as np
import pandas as pd

def test_eval(models, tloader):
    preds = []
    for model in models:
        model.eval()
        with torch.no_grad():
            for tb in tqdm(tloader):
                tb = tb.to(dev)
                y_l = tb.label
                pred = model(tb).view(1,-1).detach().cpu().numpy()

                preds.append(pred)
            
    return preds

def per_drug_acc(df):
    unique_drugs = df['rdkit.x'].unique()
    s = 0
    for drug in unique_drugs:
        filt = df[df['rdkit.x']==drug]
        score = accuracy_score(filt['true'], filt['predicted'])
        nunique_moa = filt['predicted'].nunique()
        if score >= (1/nunique_moa):
            s = s + 1
    return(s/len(unique_drugs))

def per_sig_acc(df):
    unique_sigs = df['sig_id'].unique()
    s = 0
    for sig in unique_sigs:
        filt = df[df['sig_id']==sig]
        score = accuracy_score(filt['true'], filt['predicted'])
        nunique_moa = filt['predicted'].nunique()
        if score >= (1/nunique_moa):
            s = s + 1
    return(s/len(unique_sigs))