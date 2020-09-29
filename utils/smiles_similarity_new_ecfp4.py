import urllib
import os
import numpy as np
import pandas as pd
import sys
from sys import getsizeof
import scipy
from scipy import stats

import csv
#import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.rdmolfiles import SDMolSupplier

from rdkit.Chem.Fingerprints import FingerprintMols
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 11:32:13 2019

@author: sofia
"""

def sim_two_serial():
#Load Data-----------------------------------------------------------------------
    path1 = input("Path for list 1: ")
    path2= input("Path for list 2: ")
    Ind_ch=input("Choose indexes to use (1) or run all (default=press any number): ")
    if (Ind_ch==1):
     lines2_begin=input("Indexes of rows to start reading for big: ")
     lines2_begin=int(lines2_begin)
     lines2_end=input("Indexes of rows to end reading for big: ")
     lines2_end=int(lines2_end)
     smis1=pd.read_csv(path1)
     smis1=smis1["smiles"]
     smis2=pd.read_csv(path2)
     smis2=smis2["smiles"]
     smis2=smis2[lines2_begin:lines2_end]
    else:
     smis1=pd.read_csv(path1)
     smis1=smis1["smiles"]
     smis2=pd.read_csv(path2)
     smis2=smis2["smiles"]
    l1=len(smis1)
    l2=len(smis2)
    l=l1*l2
    lp=round(l/20)

#Get molecules from smiles-----------------------------------------------------------------------
    bad1=[]
    molecules1=[]
    for i,smi in enumerate(smis1):
        m=Chem.MolFromSmiles(smi)
        if m is None:
            print('smile with number:',i,'in list 1 could not be converted to molecule')
            bad1.append(i)
            continue
        molecules1.append(m)
    
    bad2=[]
    molecules2=[]
    for i,smi in enumerate(smis2):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            print('smile with number:',i,'in list 2 could not be converted to molecule')
            bad2.append(i)
            continue
        molecules2.append(m)
			
    #can1=[Chem.MolToSmiles(x) for x in molecules1]
    #can2=[Chem.MolToSmiles(x) for x in molecules2]
    #for j in bad1:
        #can1.insert(j,"bad1")
    #for j in bad2:
        #can2.insert(j,"bad2")
    smis1=[]
    smis2=[]  


#Final output matrix-----------------------------------------------------------------------
    similarity=np.zeros(shape=(l1,l2),dtype=np.float32)

    from rdkit.Chem import MACCSkeys
    from rdkit.Chem.AtomPairs import Pairs
    from rdkit.Chem.AtomPairs import Torsions
    from rdkit.Chem import AllChem

    print('Begining fingerprint calculation...wait')
    fps_ecfp4_1 = [AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=1024) for x in molecules1]
    print('Begining fingerprint calculation...50%')
    fps_ecfp4_2 = [AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=1024) for x in molecules2]
    print('Begining fingerprint calculation...done\n')

    for j in bad1:
        fps_ecfp4_1.insert(j,1)


    for j in bad2:
        fps_ecfp4_2.insert(j,1)
    
    print('Begining of fingerprints similarity calculation\n')
    molecules1=[]
    molecules2=[]


    k=0
    for i in range(l1):
        for j in range(l2):
            if not((i in bad1) or (j in bad2)):
                similarities_ecfp4=DataStructs.FingerprintSimilarity(fps_ecfp4_1[i],fps_ecfp4_2[j]) 
                similarity[i][j]=similarities_ecfp4
            k=k+1
            if k%lp==0:
                print('running:',(k/l)*100,'%')
        #for other similarity metrics use for example DataStructs.FingerprintSimilarity(fps[0],fps[1], metric=DataStructs.DiceSimilarity)

    similarity[bad1,:]=10
    similarity[:,bad2]=10

    print('End of fingerprints similarity calculation')
    bad1=[]
    bad2=[]

    df_similarity = pd.DataFrame(similarity)
    similarity=[]
    return df_similarity


def sim_one_serial():
#Load Data-----------------------------------------------------------------------
    path=input("Path for list : ")
    smis=pd.read_csv(path)
    smis=smis["smiles"]
    l=len(smis)
    lp=round(l*l/20)
#Get molecules from smiles-----------------------------------------------------------------------
    bad=[]
    molecules=[]
    for i,smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            print('smile with number:',i,'in list could not be converted to molecule')
            bad.append(i)
            continue
        molecules.append(m)
    #can=[Chem.MolToSmiles(x) for x in molecules]
    #for j in bad:
        #can.insert(j,"bad")
    #smis=[]
#Final output matrix-----------------------------------------------------------------------
    similarity=np.zeros(shape=(l,l),dtype=np.float32)
	
    from rdkit.Chem import MACCSkeys
    from rdkit.Chem.AtomPairs import Pairs
    from rdkit.Chem.AtomPairs import Torsions
    from rdkit.Chem import AllChem
	
    print('Begining fingerprint calculation...wait')
    fps_ecfp4 = [AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=1024) for x in molecules]
    print('Begining fingerprint calculation...done\n')

    for j in bad:
        fps_ecfp4.insert(j,1)
			
    #molecules=[]
		
    print('Begining of fingerprints similarity calculation\n')
    k=0
    for i in range(l):
        for j in range(l):
            if i>=j:
                if not((i in bad) or (j in bad)):
                    similarities_ecfp4=DataStructs.FingerprintSimilarity(fps_ecfp4[i],fps_ecfp4[j])
                    similarity[i][j]=similarities_ecfp4
                    similarity[j][i]=similarity[i][j]
                k=k+1
                if k%lp==0:
                    print('running:',(k/(l*l/2))*100,'%')
        #for other similarity metrics use for example DataStructs.FingerprintSimilarity(fps[0],fps[1], metric=DataStructs.DiceSimilarity)

    similarity[bad,:]=10
    similarity[:,bad]=10

    print('End of fingerprints similarity calculation')
    bad=[]
    
    df_similarity = pd.DataFrame(similarity,columns = smis, index=smis)
    similarity=[]
    return df_similarity
  
print('Welcome to NTUA BIOSYS LAB compound similarity script')
user="wrong"
print('Enter user customization: Advanced or Default\n')
while ((user!='A') and (user!='D')):
	user=input("Enter A or D: ")
	if ((user!='A') and (user!='D')):
		print('Fatal error input,enter again')
nol=3
print('Choose to compare:')
print('[1].All compounds in one list')
print('[2].Two lists of compounds')
while ((nol!=1) and (nol!=2)):
    nol=int(input("Enter 1 or 2: "),10)
    if ((nol!=1) and (nol!=2)):
        print('Fatal error input,enter again')
out=input("Enter path/filename for outpout file [name.csv]: ")
if (user=='A'):
    print('Choose running setup:')
    run=int(input("Enter 1 for Serial or 2 for Parallel :"),10)

    if (nol==1):
        df_similarity=sim_one_serial()
#Write in csv--------------------------------------
        print('Saving results in csv...')
        df_similarity.to_csv(out)
        print('END')
    else:
        if (run==1):
            df_similarity=sim_two_serial()
#Write in csv--------------------------------------
            print('Saving results in csv...')
            df_similarity.to_csv(out)
            print('END')
        elif (run==2):
            import multiprocessing as mp
            from rdkit.Chem import MACCSkeys
            from rdkit.Chem.AtomPairs import Pairs
            from rdkit.Chem.AtomPairs import Torsions
            from rdkit.Chem import AllChem
            ncpu=mp.cpu_count()
            print('Number of CPUs core: ',ncpu)
            ncpu=int(input('Choose number of CPUs-workers: '),10)
            
# Step 1: Init multiprocessing.Pool()----------------------
            #pool = mp.Pool(ncpu)
            while (ncpu>mp.cpu_count()):
                print('Too many workers!! The max is: ',mp.cpu_count())
                ncpu=int(input('Choose number of CPUs-workers: '),10)
            #if (ncpu!=8):
                #sys.exit('sorry for know it works only with 8 cores!!')
# Step 2: `pool.apply` the `howmany_within_range()`
            #results = [pool.apply(sim_two_parallel, args=(sim_ind)) for sim_ind in range(1,9)]

# Step 3: Don't forget to close
            #pool.close()
            
            print('No parallel section yet!Sorry coming soon!')
            print('END')
else:
    run=1
    if (nol==1):
        df_similarity=sim_one_serial()
#Write in csv--------------------------------------
        print('Saving results in csv...')
        df_similarity.to_csv(out)
        print('END')
    else:
        if (run==1):
            df_similarity=sim_two_serial()
#Write in csv--------------------------------------
            print('Saving results in csv...')
            df_similarity.to_csv(out)
            print('END')
        elif (run==2):
            import multiprocessing as mp
            ncpu=mp.cpu_count()
            print('Number of CPUs core: ',ncpu)
            ncpu=int(input('Choose number of CPUs-workers: '),10)
            print('No parallel section yet!Sorry coming soon!')
            print('END')
