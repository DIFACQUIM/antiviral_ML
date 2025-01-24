import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import cdist

# Function to get quartile based on pairsim
def get_quartile(dtm, percentiles):
    if dtm <= percentiles["25%"]:
        return "Q1"
    elif dtm <= percentiles["50%"]:
        return "Q2"
    elif dtm <= percentiles["75%"]:
        return "Q3"
    elif dtm <= percentiles["max"]:
        return "Q4"
    else:
        return "Out"

# Specify paths to the datasets
data_dir = "/home/fsaldivar/Documents/Distance-to-Model-Virus/input/Quartiles"
fps_dir = "/home/fsaldivar/Documents/Distance-to-Model-Virus/input/Morgan"
external_compounds_path = "/home/fsaldivar/Documents/Distance-to-Model-Virus/input/molecules-to-screen.csv"

# List all CSV files in data and fps directories
data_files = [f for f in os.listdir(data_dir) if f.endswith('.csv')]
fps_files = [f for f in os.listdir(fps_dir) if f.endswith('.csv')]

# Initialize an empty DataFrame to store results
results = pd.DataFrame()

# Create a mapping for data files based on their initial parts
file_mapping = {}
for data_file in data_files:
    key = data_file.split('_')[0]  # Take the initial part before the first underscore
    if key not in file_mapping:
        file_mapping[key] = []
    file_mapping[key].append(data_file)

# Loop through each data file and corresponding fps file
for key in file_mapping:
    corresponding_fps_files = [f for f in fps_files if f.startswith(key)]
    for data_file in file_mapping[key]:
        for fps_file in corresponding_fps_files:
            try:
                # Read data and Morgan fingerprints
                data = pd.read_csv(os.path.join(data_dir, data_file))
                morgan_fp = pd.read_csv(os.path.join(fps_dir, fps_file), header=0, encoding='utf-8')

                # Calculate percentiles
                percentiles = data["pairsim"].describe(percentiles=[0.25, 0.5, 0.75])

                # Build external compounds DataFrame
                ext_compounds = pd.read_csv(external_compounds_path, sep=',')
                ext_compounds["mol"] = ext_compounds["smiles"].apply(lambda x: Chem.MolFromSmiles(x))

                # Compute fingerprints for external compounds
                ext_morgan_fp = pd.DataFrame([fp for fp in ext_compounds["mol"].apply(lambda x: [int(y) for y in AllChem.GetMorganFingerprintAsBitVect(x, radius=2, nBits=2048).ToBitString()])])

                # Calculate mean Jaccard distances
                ext_dtm = pd.DataFrame(cdist(ext_morgan_fp, morgan_fp, metric="jaccard")).mean(axis=1)
                ext_compounds["Pairsim"] = ext_dtm
                ext_compounds["Quartile"] = ext_compounds["Pairsim"].apply(lambda x: get_quartile(x, percentiles))
                ext_compounds["target"] = data["unique_target"].unique()[0]

                # Append results to the results DataFrame
                results = pd.concat([results, ext_compounds[['smiles', 'ID', 'Pairsim', 'Quartile', 'target']]], ignore_index=True)

            except Exception as e:
                print(f"Error en el archivo: {data_file} o {fps_file}. Detalles: {e}")

# Print paired files for reference
for key in file_mapping:
    corresponding_fps_files = [f for f in fps_files if f.startswith(key)]
    for data_file in file_mapping[key]:
        for fps_file in corresponding_fps_files:
            print(f"Processed data file: {data_file} with fps file: {fps_file}")

# Create a matrix from results
matrix = results.pivot_table(index='ID', columns='target', values='Quartile', aggfunc='first')
matrix.columns = [col.replace('target_', '') if isinstance(col, str) else col for col in matrix.columns]

# Reset the index to make 'ID' a column
matrix.reset_index(inplace=True)

# Save the results to CSV files
results.to_csv("/home/fsaldivar/Documents/Distance-to-Model-Virus/output/external_results_summary.csv", index=False)
matrix.to_csv("/home/fsaldivar/Documents/Distance-to-Model-Virus/output/matrix.csv", index=False)
