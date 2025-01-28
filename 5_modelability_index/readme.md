# Modelability Index (MODI) Analysis for ChEMBL Data Sets

This notebook calculates the Modelability Index (MODI), which measures the proportion of compounds in a data set whose nearest neighbor belongs to the same class within a defined feature space. The analysis is conducted for ChEMBL data sets using two types of molecular fingerprints:

1. **MACCS keys** (166-bits)
2. **Morgan Chiral fingerprints** (radius 2, 2048-bits)

## Methodology

The MODI values were computed using Python 3 and the following libraries:
- **RDKit**: For molecular representation and fingerprint calculation.
- **NumPy**: For efficient numerical computations.
- **pandas**: For data manipulation and organization.
- **SciPy**: For spatial analysis and nearest-neighbor calculations.

MODI was calculated for each target in the ChEMBL data sets using two approaches:
1. **Including "Mixed" Compounds**: Compounds classified as "Mixed" were considered in the overall classification.
2. **Excluding "Mixed" Compounds**: Compounds classified as "Mixed" were excluded from the classification.

## Authors

- **Gabriela Valle-Núñez**: Adapted this section.
- **Fernanda I. Saldívar-González**: Supervised the adaptation.

The methodology and implementation were inspired by the repository [MODI by rcbraga](https://github.com/rcbraga/modi).

## Reference

The original MODI methodology is detailed in:
- A. Golbraikh, E. Muratov, D. Fourches, and A. Tropsha, *Journal of Chemical Information and Modeling*, 2014, 54, 1–4. DOI: [10.1021/ci400572x](https://doi.org/10.1021/ci400572x)

---

For further details, refer to the notebook and its step-by-step implementation.

