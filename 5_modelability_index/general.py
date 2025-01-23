import pandas as pd 
import numpy as np
import copy
from tqdm.auto import tqdm 
from rdkit import Chem 
tqdm.pandas()


def check_smiles(moldf):
    from rdkit.Chem import rdMolDescriptors,rdDepictor,Draw,MolFromSmiles,MolFromSmarts,Descriptors
    smiles_list = []    
    for i,row in moldf.iterrows():
        try:
            smiles_list.append( MolFromSmiles(row["SMILES"]) )    
        except:
            smiles_list.append(None) 
    return smiles_list


def fp_calculation(df , descriptor = "maccs"):
    from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect 
    from rdkit.Chem import MACCSkeys
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdmolops import LayeredFingerprint 
    moldf = copy.deepcopy(df)
    #verify valid smiles
    smiles_list =  check_smiles(df)
    test_smiles = None in smiles_list
    if test_smiles:
        print("SMILES with issues, fix it before you go")
    else:
        if (descriptor=='maccs'):
            maccs_fp = lambda mol: MACCSkeys.GenMACCSKeys(mol)  
            moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
            moldf['Descriptors'] = moldf.Mol.progress_apply(maccs_fp)
        elif (descriptor=='morgan'):  
            nBits = 2 * 1024
            radius= 2
            calcfp = lambda mol: GetMorganFingerprintAsBitVect(mol,nBits=nBits,radius=radius, useFeatures=False, useChirality=False)
            moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
            moldf['Descriptors'] = moldf.Mol.progress_apply(calcfp)      
        elif (descriptor=='morgan4'):    
            nBits = 2 * 1024
            radius= 4
            calcfp = lambda mol: GetMorganFingerprintAsBitVect(mol,nBits=nBits,radius=radius, useFeatures=False, useChirality=False)
            moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
            moldf['Descriptors'] = moldf.Mol.progress_apply(calcfp) 
        elif (descriptor=='morgan_chiral2'):
            nBits = 2 * 1024
            radius=2
            calcfp = lambda mol: GetMorganFingerprintAsBitVect(mol,nBits=nBits,radius=radius, useFeatures=False, useChirality=True)
            moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
            moldf['Descriptors'] = moldf.Mol.progress_apply(calcfp)  
        elif (descriptor=='morgan_chiral4'):
            nBits = 2 * 1024
            radius=4
            calcfp = lambda mol: GetMorganFingerprintAsBitVect(mol,nBits=nBits,radius=radius, useFeatures=False, useChirality=True)
            moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
            moldf['Descriptors'] = moldf.Mol.progress_apply(calcfp)          
        elif (descriptor=='rdkit'):
            calcfp = lambda mol: LayeredFingerprint(mol)  
            moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
            moldf['Descriptors'] = moldf.Mol.progress_apply(calcfp)                         
        else:
            print("The descriptor is not implemented")
    return moldf

def modi(df , outcome , descriptor = "maccs"):
    moldf = copy.deepcopy(df)
    import pandas as pd 
    import numpy as np
    from scipy.spatial import distance
    from funcs.sirms import calculate_sirms
    if (descriptor=='sirms'):
        data = calculate_sirms(moldf)
        data.fillna(0)
        data = data.drop(columns=['Compounds'])
        a = data.to_numpy()
    else:
        moldf = fp_calculation(moldf, descriptor)
        a = np.array(list(moldf['Descriptors']))
    dm = pd.DataFrame(distance.cdist(a, a, 'euclidean'))
    dm["nn_dist"] = dm.replace(0, np.nan).min()
    dm2 = moldf.merge(dm, left_index=True, right_index=True)    
    outcomeList = []
    for x in range(dm2.shape[0]):
        sel_index= dm2[x].index[dm2[x] == dm2["nn_dist"][x]].values
        outcomeList.append(dm2[outcome][dm2.index == sel_index[0]].to_list()[0])
    dm2["nn_outcome"] = outcomeList
    tot_agree = dm2[dm2[outcome] == dm2["nn_outcome"]].shape[0]
    tot = dm.shape[0]
    modi = tot_agree/tot
    return modi 

  # adapted from https://github.com/rcbraga/modi
  
  # A. Golbraikh, E. Muratov, D. Fourches and A. Tropsha, J. Chem. Inf. Model., 2014, 54, 1–4. (DOI: 10.1021/ci400572x).
