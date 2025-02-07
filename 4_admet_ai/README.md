# ADMET Properties Calculation (ADMET-AI)

This notebook calculates Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET) predictions. This was accomplished using ADMET-AI, with the open-source Python library for local predictions.

## Methodology

ADMET properties were computed using Python 3 and the following libraries:
- **pandas:** For data manipulation and organization.
- **ADMET-AI:** For the calculation of ADMET descriptors.

Predictions were calculated with a Jupyter notebook as described by the authors of ADMET-AI, with a conda environment:

```bash
conda create -y -n admet_ai python=3.10
conda activate admet_ai
```

Install ADMET-AI via pip.

```bash
pip install admet-ai
conda install jupyter 
```

## Authors
- **Diana L. Prado-Romero**: Adapted this section.
- **Pedro A. Laurel-García**: First example for the adaptation.

The methodology and implementation were inspired by the repository [ADMET-AI by swansonk14](https://github.com/swansonk14/admet_ai).

## Reference

The development of ADMET-AI models is detailed in:
- K. Swanson, P. Walther, J. Leitz, S. Mukherjee, J. C. Wu, R. V. Shivnaraine, J. Zou, Bioinformatics, 2024, 40, btae416.
<br />DOI: [10.1093/bioinformatics/btae416](https://doi.org/10.1093/bioinformatics/btae416)  
