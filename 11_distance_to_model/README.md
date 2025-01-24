
# Distance to model

This script calculates similarity scores between external compounds and a training dataset using Morgan fingerprints and molecular properties. It categorizes the results into quartiles based on two distance metrics: Jaccard distance for structural similarity (based on Morgan fingerprints) and Euclidean distance for molecular property similarity. The output includes the quartile assignment (Q1–Q4 or "Out") for each compound, providing an indication of its confidence level in terms of similarity to the training dataset.

## Overview

- **Similarity Measurement**: The Jaccard distance is used to quantify the structural similarity between the input ligand and compounds in the training dataset. For evaluating molecular properties similarity, the Euclidean distance is applied to compare features derived from the molecular descriptors.
  
- **Confidence Categorization**: The results are divided into four quartiles (Q1–Q4) based on the **average Jaccard distance** or **Euclidian distance**. Compounds are categorized into:
  - **Q1–Q4**: According to this classification, compounds in the first quartile (Q1) have smaller distances from the training 183set, indicating higher similarity. In contrast, values placed in the fourth quartile (Q4) reflect that the predicted 184compound is farther away from the training set, suggesting lower similarity. 
  - **OUT**: Low confidence predictions due to significant dissimilarity with the training data.

- **Fingerprint**: The compounds are represented using **Morgan fingerprints** with a radius of 2 (2048 bits) for structural comparison.

- **Molecular properties** The molecular properties used to build the ML models are incorporated into the similarity assessment to further refine the comparison and confidence scoring.

## Usage

You can use this code locally (download this repository and install dependencies).

## Environment Setup

We provide an `environment.yml` file containing all the required dependencies for this code.

1. **Create the Conda environment**:

    ```bash
    conda env create -f environment.yml
    ```

    > **Note**: If you prefer, you can create the environment manually. Check "Alternatively create conda environment manually" for instructions.

2. **Activate the Conda environment**:

    ```bash
    conda activate Distance
    ```

Once the environment is activated, you are ready to work in it and start your calculations.


1. **Input folder**:
   - Training Data: CSV files containing training compounds and their similarity scores (pairsim column). The files are located in the Quartiles folder.
   - Morgan Fingerprints: CSV files containing the corresponding Morgan fingerprints for the training compounds, located in the Morgan folder.
   - Molecular properties: CSV files containing the corresponding Morgan fingerprints for the training compounds, located in the Properties folder.
  - External Compounds: A CSV file (molecules-to-screen.csv) containing the external compounds (ligands) to be screened. This file must contain a column called smiles.
5. **Otput **: The output of the script includes:

  **Results DataFrame**: A DataFrame containing the following columns for each external compound:
   - `smiles`: The SMILES representation of the compound.
   - `ID`: Unique identifier for the compound.
   - `Pairsim`: The mean Jaccard distance score.
   - `Quartile`: The quartile (Q1-Q4) or "Out" category based on the similarity to the training dataset.
   - `target`: The target type for the training data.
  **Quartile Matrix**: A pivot table that displays the quartile assignments of external compounds across different targets.

## Example

```bash
# Example to run the script
python distance_to_model2.py
```

After running the script, the resulting `matrix` DataFrame will be available, showing the quartile assignments of compounds across different targets.

## Error Handling

- The script includes error handling to catch issues during file reading or processing.
- Errors are logged with the specific file causing the issue.



---
## References

- Sánchez-Cruz, N.; Medina-Franco, J. L. [Epigenetic TargetProfiler: A Web Server to Predict Epigenetic Targets of Small Molecules](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00045). *J. Chem. Inf. Model.* 2021, 61 (4), 1550−1554.

