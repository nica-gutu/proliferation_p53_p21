# p53/p21 dynamics and DNA-damage response heterogeneity

Analysis code supporting:

**Gutu, N. et al. “p53 and p21 dynamics encode single-cell DNA damage levels, fine-tuning proliferation and shaping population heterogeneity.” _Communications Biology_ 6, 1196 (2023).**  
https://doi.org/10.1038/s42003-023-05585-5

## Overview

This project investigates how heterogeneous levels of DNA damage are encoded by the **p53–p21 signaling network** and how these dynamics contribute to cell-to-cell differences in proliferation.

The study combines two previously published single-cell datasets representing endogenous and exogenous DNA-damage conditions. The analyses connect long-term signaling dynamics to proliferation at the single-cell and population levels, with particular emphasis on whether DNA-damage strength is represented by changes in p53/p21 amplitude, pulse timing, and cell-division behavior.

## Study design

Two complementary experimental settings were analyzed:

1. **Endogenous DNA-damage heterogeneity**  
   Individual cells were tracked for 52 h to reconstruct division histories and were subsequently fixed and stained for the DNA-damage marker γH2AX. This links each cell's proliferation trajectory directly to its measured DNA-damage level.

2. **Radiation-induced DNA damage**  
   Human retinal pigment epithelial (RPE) cells exposed to **0, 2, 4, or 10 Gy** gamma radiation were followed for more than 5 days. Single-cell p53 and p21 dynamics were recorded together with division events, allowing signaling properties to be related to radiation dose and long-term proliferative outcome.

## Main findings represented by the analyses

- Endogenous DNA-damage levels are associated with heterogeneous proliferation behavior across individual cells.
- Differences in overall proliferation are not explained simply by progressive changes in intermitotic time; instead, the data are consistent with cells entering non-proliferative states after different numbers of divisions.
- Long-term **p53 amplitude and p21 levels change gradually with DNA-damage strength** and are strongly associated with proliferation.
- Radiation dose is encoded not only by the number of p53 pulses but also by changes in long-term p53/p21 amplitude.
- Time-resolved analysis reveals **dose-dependent prolongation of the p53 pulse period**.
- A subset of cells undergoes a temporal switch in p53 oscillatory behavior that is associated with transition from a low- to a higher-proliferative state.

These results support a quantitative view of the p53–p21 network in which DNA-damage strength is encoded continuously rather than as a simple binary arrest/proliferation decision.

## Repository structure

The repository is organized primarily by manuscript figure:

- `Fig1/` — proliferation trajectories and endogenous DNA-damage analyses.
- `Fig1S/` — supplementary proliferation and treatment-response analyses.
- `Fig2/` — radiation-dose dependence, p53/p21 signal features, correlations, and proliferation metrics.
- `Fig2S/` — supplementary signaling-amplitude, arrest, and division analyses.
- `Fig3/` — time-dependent p53 pulse-period analysis, proliferation-state comparisons, and wavelet analyses.
- `Fig3S/` and `Fig4S/` — additional supplementary analyses.
- `Box/` — supporting intermitotic-time and proliferation visualizations.

The scripts correspond to specific analyses and manuscript figures rather than forming a packaged software library.

## Computational methods

Methods represented in the repository include:

- single-cell time-series processing;
- p53 and p21 amplitude and trend quantification;
- proliferation metrics including number of divisions, cell age, and intermitotic time;
- correlation and distribution analyses;
- continuous-wavelet analysis of p53 oscillations with `pyBOAT`;
- peak-based pulse characterization;
- time-dependent classification of p53 period trajectories;
- stationarity testing;
- statistical comparison and scientific visualization.

## Requirements

The scripts are written in Python. Packages used across the analyses include:

```text
numpy
pandas
matplotlib
seaborn
pyboat
statsmodels
```

Individual scripts may require additional standard scientific-Python dependencies.

## Data and reproducibility

The manuscript analyzes previously published single-cell datasets; according to the paper, all supporting data are available through the original data publications.

Several scripts retain paths from the original research environment. To reproduce a specific analysis, obtain the corresponding source dataset and update the input/output paths in the relevant script.

The published code associated with the paper is also available through the Granada Lab repository:

https://github.com/Granada-Lab/proliferation-p53-p21

## Citation

> Gutu, N., Binish, N., Keilholz, U. et al. **p53 and p21 dynamics encode single-cell DNA damage levels, fine-tuning proliferation and shaping population heterogeneity.** _Communications Biology_ 6, 1196 (2023). https://doi.org/10.1038/s42003-023-05585-5

## Contact

**Nica Gutu**  
Computational Biology / Data Science  
https://nica-gutu.github.io/website/
