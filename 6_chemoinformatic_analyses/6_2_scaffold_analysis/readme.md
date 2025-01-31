# Scaffold Analysis 

## Overview
This repository contains a Google Colab notebook for scaffold analysis in cheminformatics, focusing on antiviral compounds. Scaffold-based analysis helps identify core molecular frameworks, which are essential for understanding structure-activity relationships (SAR) and optimizing drug discovery pipelines.

## Features
- **Dataset Loading:** Reads antiviral compound datasets from various sources.
- **Scaffold Extraction:** Uses RDKit to generate Murcko scaffolds.
- **Visualization:** Displays molecular structures and their core scaffolds.
- **Frequency Analysis:** Identifies prevalent scaffolds within datasets.

## Requirements
The notebook uses Python and requires the following libraries:
- `rdkit`
- `pandas`
- `numpy`
- `google.colab`

To install RDKit, run:
```python
!pip install rdkit
```

## Usage
1. Open the notebook in Google Colab.
2. Load datasets.
3. Run the cells sequentially to perform scaffold extraction and analysis.
4. Interpret the results through visualizations and frequency tables.

## Data
The analysis uses antiviral compound datasets stored in CSV format. These datasets contain molecular structures and their associated properties.

## Contributions
Feel free to contribute by:
- Adding new datasets.
- Improving scaffold analysis methods.
- Enhancing visualizations.

## License
This project is for research and educational purposes. Please cite appropriate sources if using the notebook in publications.


