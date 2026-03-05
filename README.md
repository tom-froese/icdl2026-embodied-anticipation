# Embodied Temporal Anticipation in Perceptual Crossing

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b+-blue.svg)](https://www.mathworks.com/products/matlab.html)

Companion code and preprocessed data for the IEEE ICDL 2026 paper:

> **Embodied Temporal Anticipation Without Clocks: Poisson Process Dynamics in Perceptual Crossing**
>
> Tom Froese
>
> Embodied Cognitive Science Unit, Okinawa Institute of Science and Technology (OIST), Japan

## Overview

This repository contains the preprocessing scripts, analysis code, and preprocessed data needed to reproduce all figures and statistical results reported in the paper. The analysis demonstrates that physiological signals (electrodermal activity and respiration) and behavioral responses (click times) during a perceptual crossing experiment follow Poisson process dynamics with λ = *e*, providing an embodied mechanism for optimal stopping without explicit time representation.

## Repository Structure

```
.
├── code/
│   ├── preprocessing/          # Raw data → preprocessed CSV + JSON metadata
│   │   ├── preprocessEDA.m     # EDA preprocessing pipeline
│   │   ├── preprocessARU.m     # Respiration (ARU) preprocessing pipeline
│   │   └── preprocessClicks.m  # Click response-time extraction
│   └── analysis/               # Preprocessed data → figures + statistics
│       ├── plotEDAFigures.m    # Figures 5, 6 + SI onset-lag figure
│       ├── plotARUFigures.m    # Figure 7 (respiration)
│       └── plotClickTimesFigure.m  # Figure 4 (click response times)
├── data/
│   ├── EDA/                    # Preprocessed electrodermal activity
│   │   ├── EDA_Task_Preprocessed.csv
│   │   ├── EDA_Task_Preprocessed.json
│   │   ├── EDA_Rest_Preprocessed.csv
│   │   └── EDA_Rest_Preprocessed.json
│   ├── ARU/                    # Preprocessed respiration
│   │   ├── ARU_Task_Preprocessed.csv
│   │   ├── ARU_Task_Preprocessed.json
│   │   ├── ARU_Rest_Preprocessed.csv
│   │   └── ARU_Rest_Preprocessed.json
│   └── ClickTimes/             # Click response times
│       ├── ClickResponseTimes.csv
│       └── ClickResponseTimes.json
├── results/                    # Generated figures (600 dpi PNG)
│   ├── Fig4_ClickTimes.png
│   ├── Fig5_EDA_Rest.png
│   ├── Fig6_EDA_Task.png
│   ├── Fig7_ARU.png
│   └── FigSI_EDA_OnsetLag.png
├── LICENSE
└── README.md
```

## Requirements

- **MATLAB** R2023b or later (Signal Processing Toolbox, Statistics and Machine Learning Toolbox)

No additional toolboxes or third-party packages are required.

## Reproducing the Figures

All analysis scripts are designed to be run from `code/analysis/`:

```matlab
cd code/analysis

% Figure 4: Click response-time distribution with P(1) fit
plotClickTimesFigure

% Figures 5 & 6: EDA rest and task conditions with P(0) fits
% Also generates supplementary onset-lag figure
plotEDAFigures

% Figure 7: Respiration (ARU) with P(1) fit
plotARUFigures
```

Output figures are saved to `results/` at 600 dpi (EDA) or 300 dpi (click times, ARU).

## Data Description

The preprocessed data come from a perceptual crossing experiment with 32 dyads (64 participants). Each JSON sidecar file documents the preprocessing pipeline, column definitions, and summary statistics in BIDS-inspired format.

| Dataset | Condition | Samples | Duration | Sampling Rate |
|---------|-----------|---------|----------|---------------|
| EDA     | Task (18 trials × 60 s) | ~43 M rows | 60 s/trial | 25 Hz |
| EDA     | Rest (4 periods × 180 s) | ~28 M rows | 180 s/period | 25 Hz |
| ARU     | Task | ~43 M rows | 60 s/trial | 25 Hz |
| ARU     | Rest | ~28 M rows | 180 s/period | 25 Hz |
| Clicks  | Task | 1,152 rows | — | Event-based |

### Preprocessing Pipelines

The `code/preprocessing/` scripts document the full pipeline from raw experimental recordings to the preprocessed CSV files. Raw data are not included in this repository due to size and participant privacy, but the preprocessing scripts are provided for transparency and methodological reference.

**EDA**: 4th-order Butterworth low-pass filter (5 Hz cutoff, zero-phase), downsampling from 1000 Hz to 25 Hz, artifact flagging via range and gradient checks (Boucsein et al., 2012).

**ARU (Respiration)**: Analogous pipeline adapted for respiratory signal characteristics.

**Click Times**: First button-press detection per trial per participant from raw trial CSV recordings.

## Raw Data

The raw experimental data were collected using the open-source Perceptual Crossing Experiment (PCE) platform described in:

- Lenay, C. et al. (2021). Perceptual crossing device.
- Lerique, S. et al. (2024). The ECS-UPCE Dataset.

Access to the raw dataset can be arranged by contacting the corresponding author.

## Citation

If you use this code or data, please cite:

```bibtex
@inproceedings{froese2026embodied,
  title     = {Embodied Temporal Anticipation Without Clocks: Poisson Process Dynamics in Perceptual Crossing},
  author    = {Froese, Tom},
  booktitle = {Proceedings of the 2026 IEEE International Conference on Development and Learning (ICDL)},
  year      = {2026},
  publisher = {IEEE}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
