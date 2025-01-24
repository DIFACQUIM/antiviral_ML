import pandas as pd
import numpy as np
import os
from rdkit import Chem
from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler
import datamol as dm

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
data_dir = "/home/fernanda/Documents/avellaneda/4_virus_ml/4_distance_to_model/Distance-to-Model-Virus_properties_final/Distance-to-Model-Virus_properties/input/Quartiles"
properties_dir = "/home/fernanda/Documents/avellaneda/4_virus_ml/4_distance_to_model/Distance-to-Model-Virus_properties_final/Distance-to-Model-Virus_properties/input/Properties"
external_compounds_path = "/home/fernanda/Documents/avellaneda/4_virus_ml/4_distance_to_model/Distance-to-Model-Virus_properties_final/Distance-to-Model-Virus_properties/input/df_to_screen_complete.csv"



# List all CSV files in data and properties directories
data_files = [f for f in os.listdir(data_dir) if f.endswith('.csv')]
properties_files = [f for f in os.listdir(properties_dir) if f.endswith('.csv')]

# Initialize an empty DataFrame to store results
results = pd.DataFrame()

# Create a mapping for data files based on their initial parts
file_mapping = {}
for data_file in data_files:
    key = data_file.split('_')[0]  # Tomar la parte inicial antes del primer guion bajo
    if key not in file_mapping:
        file_mapping[key] = []
    file_mapping[key].append(data_file)

# Loop through each data file and corresponding fps file
for key in file_mapping:
    corresponding_properties_files = [f for f in properties_files if f.startswith(key)]
    for data_file in file_mapping[key]:
        for properties_file in corresponding_properties_files:
            try:
                data = pd.read_csv(os.path.join(data_dir, data_file))
                properties = pd.read_csv(os.path.join(properties_dir, properties_file))
                
                # Calculate percentiles
                percentiles = data["Distance"].describe(percentiles=[0.25, 0.5, 0.75])
                
                # Read external compounds DataFrame
                ext_compounds = pd.read_csv(external_compounds_path, sep=',')
                ext_compounds["mol"] = ext_compounds["smiles"].apply(lambda x: Chem.MolFromSmiles(x))
                
                # Compute molecular descriptors 
                dm_descriptors_df = dm.descriptors.batch_compute_many_descriptors(ext_compounds["smiles"].apply(dm.to_mol))
                
                # Filter only the first 22 columns of descriptors
                dm_descriptors_df_temp = dm_descriptors_df.iloc[:, 0:22].values
                
                # Normalize descriptors
                x_scaled = StandardScaler().fit_transform(dm_descriptors_df_temp)
                
                # Create a DataFrame with the normalized descriptors.
                dm_descriptors_df_norm = pd.DataFrame(x_scaled, index=dm_descriptors_df.index, columns=dm_descriptors_df.columns[0:22])
                
                # Make sure that the selected columns match the ones you have in `properties`.
                headers_list = properties.columns.tolist()  
                df_filtered = dm_descriptors_df_norm[headers_list]
                
                # Calculate Euclidean distances
                ext_dtm = pd.DataFrame(cdist(df_filtered, properties, metric="euclidean")).mean(axis=1)
                
                # Assign the calculated distances and quartiles
                ext_compounds["Distance"] = ext_dtm
                ext_compounds["Quartile"] = ext_compounds["Distance"].apply(lambda x: get_quartile(x, percentiles))
                ext_compounds["target"] = data["unique_target"].unique()[0]
                
                # Agregar los resultados al DataFrame final
                results = pd.concat([results, ext_compounds[['smiles', 'ID', 'Distance', 'Quartile', 'target']]], ignore_index=True)
                
            except Exception as e:
                print(f"Error en el archivo: {data_file} o {properties_file}. Detalles: {e}")


# Print paired files for reference
for key in file_mapping:
    corresponding_properties_files = [f for f in properties_files if f.startswith(key)]
    for data_file in file_mapping[key]:
        for properties_file in corresponding_properties_files:
            print(f"Processed data file: {data_file} with fps file: {properties_file}")

# Create a matrix of results
matrix = results.pivot_table(index='ID', columns='target', values='Distance', aggfunc='first')

# Clear column names
matrix.columns = [col.replace('target_', '') if isinstance(col, str) else col for col in matrix.columns]


matrix.reset_index(inplace=True)

# Save results in csv files
results.to_csv("/home/fernanda/Documents/avellaneda/4_virus_ml/4_distance_to_model/Distance-to-Model-Virus_properties_final/Distance-to-Model-Virus_properties/output/external_results_summary.csv", index=False)
matrix.to_csv("/home/fernanda/Documents/avellaneda/4_virus_ml/4_distance_to_model/Distance-to-Model-Virus_properties_final/Distance-to-Model-Virus_properties/output/matrix.csv", index=False)

print("Process completed successfully.")
