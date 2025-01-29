# IAV Polymerase (PA) Model Training

## Project Overview
This project focuses on developing supervised binary classification models using PyCaret (v3.3.2) to predict the activity of compounds against the Influenza A Virus (IAV) Polymerase (PA) target. The models are trained on ChEMBL data with binary activity labels (active/inactive) using molecular representations such as Morgan Chiral fingerprints and physicochemical descriptors.

## Dataset and Features
- **Source**: ChEMBL database
- **Molecular Representations**:
  - Morgan Chiral fingerprints (radius 2, 2048-bits)
  - Physicochemical descriptors
- **Preprocessing Steps**:
  - Z-score normalization
  - Stratified k-fold cross-validation
  - Imbalance correction via ADASYN

## Methodology
1. **Data Preprocessing**: Cleaning, feature engineering, and handling class imbalance.
2. **Model Training**: Various classification models are trained using PyCaret.
3. **Validation**:
   - Internal validation via cross-validation.
   - External validation using an unseen test set.
4. **Evaluation Metrics**:
   - Accuracy, AUC, Recall, Precision, F1 score, Kappa, MCC, and Balanced Accuracy.
   - MCC (Matthews Correlation Coefficient) is emphasized for robust classification performance.

## Requirements
To run the notebook, install the following dependencies:
```bash
pip install pycaret deepchem scikit-learn pandas numpy matplotlib seaborn
```

## Usage
1. Open the Jupyter Notebook:
   ```bash
   jupyter notebook 7_IAV_Polymerase_(PA)_model_training.ipynb
   ```
2. Run the notebook cells sequentially.
3. Review model performance metrics and visualization outputs.

## Results
The models are evaluated based on MCC and other classification metrics, ensuring reliable predictions for antiviral activity.

## License
This project is for research purposes only. Please cite appropriately if used in publications.
