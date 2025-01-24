import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import rdMolDescriptors,rdDepictor,Draw,MolFromSmiles,MolFromSmarts,Descriptors
import time
import uuid
import os

def calculate_sirms(moldf):
  import sirms, os
  #generate UUIDS (avoid contamination)
  ids = str(uuid.uuid1())
  idhash = str(ids)
  idpath = idhash+'.sdf'
  idpathout = idhash+'.txt'
  idpathoutcsv = '../data/sirms_descs.json'
  sirmsexe = 'sirms -i '+''+idpath+' -o '+idpathout
  moldf['mol'] = moldf.SMILES.apply(MolFromSmiles)
  PandasTools.WriteSDF(moldf , idpath, molColName='mol', properties = list(moldf.columns))
  #_start = time.time()
  os.system(sirmsexe)
  #print(f"Execution time: {  (time.time() - _start )/60 }")
  data = pd.read_csv(idpathout, delimiter = "\t")
  os.remove(idpathout)
  os.remove(idpath)  
  return data
  
  # adapted from https://github.com/rcbraga/modi
  
  # A. Golbraikh, E. Muratov, D. Fourches and A. Tropsha, J. Chem. Inf. Model., 2014, 54, 1–4. (DOI: 10.1021/ci400572x).
