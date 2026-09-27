# p53/p21 dynamics and DNA-damage response heterogeneity

Analysis code supporting the study **“p53 and p21 dynamics encode single-cell DNA damage levels, fine-tuning proliferation and shaping population heterogeneity”** (*Communications Biology*, 2023).

**Publication:** https://doi.org/10.1038/s42003-023-05585-5

## Scientific question

Cells exposed to the same genotoxic stress can show markedly different signaling and proliferation responses. This project investigates how heterogeneous levels of endogenous and exogenous DNA damage are encoded in the dynamics of the tumor suppressor **p53** and its downstream target **p21**, and how these signaling dynamics shape proliferation behavior at the single-cell and population levels.

The analyses combine long-term live-cell measurements, time-series analysis, signal-feature extraction, wavelet analysis, and statistical comparisons. In particular, the work examines how p53/p21 signal properties vary with DNA-damage level and how changes in p53 pulse timing relate to transitions between low- and high-proliferative states.

## Repository contents

The repository is organized primarily by manuscript figure.

- `Fig1/` — long-term proliferation and signaling analyses.
- `Fig1S/` — supplementary analyses related to treatment sensitivity and proliferation.
- `Fig2/` — radiation-dose dependence, p53/p21 signal features, correlations, and proliferation measures.
- `Fig2S/` — supplementary analyses of arrested cells, signaling amplitudes, and division behavior.
- `Fig3/` — p53 pulse-period dynamics, proliferation-state comparisons, and wavelet-based analyses.
- `Fig3S/` and `Fig4S/` — additional supplementary analyses.
- `Box/` — supporting analyses and visualizations of intermitotic-time and proliferation behavior.

The scripts are research-analysis code corresponding to individual analyses and figures rather than a packaged software library.

## Methods represented in the code

The repository includes analyses involving:

- single-cell time-series processing;
- p53 and p21 signaling features;
- proliferation and intermitotic-time measurements;
- dose-response comparisons;
- correlation and distribution analyses;
- wavelet-based estimation of dynamic p53 pulse periods using `pyBOAT`;
- statistical visualization and figure generation.

## Requirements

The scripts are written in Python. Core packages used across the repository include:

```text
numpy
pandas
matplotlib
seaborn
pyboat
```

Individual scripts may require additional standard scientific-Python dependencies.

## Data and reproducibility

The raw experimental datasets are **not included in this repository**. Several scripts were written for the original analysis environment and contain local input/output paths. To rerun an analysis, the corresponding source data must be available and the paths in the relevant script should be adjusted.

Because the repository mirrors the analysis underlying specific manuscript figures, the most direct way to navigate it is to identify the figure or analysis of interest and use the corresponding directory.

## Citation

If you use this code or build on the analyses, please cite:

> Gutu, N. et al. **p53 and p21 dynamics encode single-cell DNA damage levels, fine-tuning proliferation and shaping population heterogeneity.** *Communications Biology* 6, 1196 (2023). https://doi.org/10.1038/s42003-023-05585-5

## Contact

**Nica Gutu**  
Computational Biology / Data Science  
https://nica-gutu.github.io/website/
