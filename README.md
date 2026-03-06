# Embodied Temporal Anticipation in Perceptual Crossing

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b+-blue.svg)](https://www.mathworks.com/products/matlab.html)

Companion code and preprocessed data for:

> **Conference paper (IEEE ICDL 2026, under review):**
> Embodied temporal anticipation: Stochastic processes approximate optimal waiting strategies in a social decision-making task.
> Tom Froese. Embodied Cognitive Science Unit (ECSU), OIST, Japan.

> **Journal paper (IEEE TCDS, in preparation):**
> [Title TBD — extends the conference paper with haptic feedback, PAS ratings, and EEG coherence analyses.]
> Tom Froese, Dini Putri. ECSU, OIST, Japan.

## Overview

This repository contains preprocessing scripts, analysis code, and preprocessed data for reproducing all figures and statistical results from both papers. The core dataset comes from a perceptual crossing experiment (PCE) with 32 dyads (64 participants).

**Conference paper (Figs. 4–7, ICDL):** Physiological signals (electrodermal activity, respiration) and behavioral responses (click times) follow Poisson process dynamics with λ = *e*, providing an embodied mechanism for optimal stopping.

**Journal extension (Figs. 5–7, TCDS):** Adds haptic feedback time series, Perceptual Awareness Scale (PAS) ratings with logistic regression crossover analysis, and EEG beta-band (13–30 Hz) magnitude-squared coherence comparing within-brain functional connectivity between sampling (PRE, 17–21 s) and selection (POST, 34–38 s) epochs.

## Repository Structure

```
.
├── code/
│   ├── preprocessing/                  # Raw data → preprocessed CSV + JSON
│   │   ├── preprocessClicks.m          # Click response-time extraction
│   │   ├── preprocessEDA.m             # EDA preprocessing pipeline
│   │   ├── preprocessARU.m             # Respiration (ARU) preprocessing
│   │   ├── preprocessHaptics.m         # Haptic feedback extraction        [TCDS]
│   │   ├── preprocessPAS.m             # PAS rating extraction             [TCDS]
│   │   └── preprocessEEG.m            # EEG beta-band coherence (MSC)     [TCDS]
│   └── analysis/                       # Preprocessed data → figures + stats
│       ├── plotClickTimesFigure.m      # Figure 4: Click times + P(1) fit
│       ├── plotEDAFigures.m            # Figures 5–6: EDA rest/task + P(0) fits
│       ├── plotARUFigures.m            # Figure 7: Respiration + P(1) fit
│       ├── plotHapticFeedbackFigure.m  # Figure 5: Haptic grand average    [TCDS]
│       ├── plotPASFigures.m            # Figure 6: PAS percentages + crossover [TCDS]
│       └── plotEEGFigures.m            # Figure 7: EEG coherence (3 panels) [TCDS]
├── data/
│   ├── ClickTimes/                     # Click response times
│   │   ├── ClickResponseTimes.csv
│   │   └── ClickResponseTimes.json
│   ├── EDA/                            # Preprocessed electrodermal activity
│   │   ├── EDA_Task_Preprocessed.csv
│   │   ├── EDA_Task_Preprocessed.json
│   │   ├── EDA_Rest_Preprocessed.csv
│   │   └── EDA_Rest_Preprocessed.json
│   ├── ARU/                            # Preprocessed respiration
│   │   ├── ARU_Task_Preprocessed.csv
│   │   ├── ARU_Task_Preprocessed.json
│   │   ├── ARU_Rest_Preprocessed.csv
│   │   └── ARU_Rest_Preprocessed.json
│   ├── Haptics/                        # Preprocessed haptic feedback       [TCDS]
│   │   ├── HapticFeedback.csv
│   │   └── HapticFeedback.json
│   ├── PAS/                            # Perceptual Awareness Scale ratings [TCDS]
│   │   ├── PASRatings.csv
│   │   └── PASRatings.json
│   └── EEG/                            # EEG beta-band coherence           [TCDS]
│       ├── BetaCoherence.csv
│       ├── BetaCoherence_ROI.csv
│       └── BetaCoherence.json
├── results/                            # Generated figures (300–600 dpi PNG)
│   ├── Fig4_ClickTimes.png             # ICDL Fig. 4
│   ├── Fig5_EDA_Rest.png               # ICDL Fig. 5
│   ├── Fig6_EDA_Task.png               # ICDL Fig. 6
│   ├── Fig7_ARU.png                    # ICDL Fig. 7
│   ├── FigSI_EDA_OnsetLag.png          # ICDL supplementary
│   ├── Fig_HapticFeedback.png          # TCDS Fig. 5                       [TCDS]
│   ├── Fig_PAS.png                     # TCDS Fig. 6 (two panels)          [TCDS]
│   └── Fig_EEG.png                     # TCDS Fig. 7 (three panels)        [TCDS]
├── CITATION.cff
├── LICENSE
└── README.md
```

Items marked `[TCDS]` were added for the journal extension.

## Requirements

- **MATLAB** R2023b or later (Signal Processing Toolbox, Statistics and Machine Learning Toolbox)

No additional toolboxes or third-party packages are required. EEG preprocessing uses `hann()` (not the deprecated `hanning()`).

## Reproducing the Figures

All analysis scripts run from `code/analysis/` and use relative paths (`../../data/`, `../../results/`):

```matlab
cd code/analysis

%% --- ICDL 2026 conference figures ---

% Figure 4: Click response-time distribution with P(1) fit
plotClickTimesFigure

% Figures 5 & 6: EDA rest and task conditions with P(0) fits
plotEDAFigures

% Figure 7: Respiration (ARU) with P(1) fit
plotARUFigures

%% --- IEEE TCDS journal figures ---

% Figure 5: Haptic feedback grand-average time series
plotHapticFeedbackFigure

% Figure 6: PAS relative percentages (A) + logistic crossover (B)
plotPASFigures

% Figure 7: EEG beta coherence — topography (A), ROI bars (B), histogram (C)
plotEEGFigures
```

Output figures are saved to `results/` at 300 dpi.

## Data Description

Preprocessed data from 32 dyads (64 participants). Dyad 31 is absent from the EEG data (no exclusion logic — simply not present in the source recordings). Each JSON sidecar documents preprocessing, column definitions, and summary statistics in BIDS-inspired format.

| Dataset | Source | Rows | Key Variables |
| --- | --- | --- | --- |
| ClickTimes | Behavioral | 1,152 | DyadID, ParticipantID, TrialNum, ClickTime_s, Clicked |
| EDA (Task) | Physiological | ~43 M | Dyad, Participant, Trial, Time, EDA_Clean |
| EDA (Rest) | Physiological | ~28 M | Dyad, Participant, Period, Time, EDA_Clean |
| ARU (Task) | Physiological | ~43 M | Dyad, Participant, Trial, Time, ARU_Clean |
| ARU (Rest) | Physiological | ~28 M | Dyad, Participant, Period, Time, ARU_Clean |
| Haptics | Behavioral | — | DyadID, ParticipantID, TrialNum, Time, HapticFeedback |
| PAS | Behavioral | — | DyadID, ParticipantID, TrialNum, PASRating |
| EEG (pairs) | Neural | 351 × 62 | ParticipantID, ChannelPair, PRE_Coherence, POST_Coherence |
| EEG (ROI) | Neural | 6 × 62 | ParticipantID, ROI, PRE_Coherence, POST_Coherence |

### Preprocessing Pipelines

The `code/preprocessing/` scripts document the full pipeline from raw experimental recordings to preprocessed CSV files. Raw data are not included due to size and participant privacy.

**EDA:** 4th-order Butterworth low-pass filter (5 Hz cutoff, zero-phase), downsampled from 1000 Hz to 25 Hz, artifact flagging via range and gradient checks (Boucsein et al., 2012).

**ARU (Respiration):** Analogous pipeline adapted for respiratory signal characteristics.

**Click Times:** First button-press detection per trial per participant.

**Haptic Feedback:** Haptic vibration motor state extracted from raw trial CSVs and aligned to a common time axis for grand averaging.

**PAS Ratings:** Perceptual Awareness Scale responses (1–4) extracted per trial per participant.

**EEG Coherence:** Magnitude-squared coherence (MSC) computed in the beta band (13–30 Hz) between all 351 unique pairs of 27 EEG channels, for PRE (17–21 s, sampling epoch) and POST (34–38 s, selection epoch) windows. Hann-windowed, 4 s segments. Individual participants (N = 62) as unit of analysis. ROI aggregation over six anatomical regions.

### Key Behavioral Landmarks

The PRE/POST epoch definitions are motivated by converging behavioral landmarks:

- **P(1) click-time peak:** ~24.4 s
- **PAS 4/3 crossover:** ~27.7 s (logistic regression)
- **Haptic stabilization:** ~4 s post-crossover

## Raw Data

Raw data were collected using the open-source PCE device described in [Estelle et al. (2024)](https://doi.org/10.1371/journal.pone.0305283). The full raw dataset is available at [Lerique et al. (2024)](https://osf.io/preprints/osf/6hjfy_v1).

The EEG coherence analyses build on methods from:

> Putri, D. & Froese, T. (2025). Rapid recognition of haptic feedback as mediated social touch: EEG evidence from the perceptual crossing paradigm. *[Journal TBD]*.

## Citation

If you use this code or data, please cite:

```bibtex
@inproceedings{froese2026embodied,
  title     = {Embodied temporal anticipation: Stochastic processes approximate
               optimal waiting strategies in a social decision-making task},
  author    = {Froese, Tom},
  booktitle = {Proc. IEEE Int. Conf. Development and Learning (ICDL)},
  year      = {2026},
  publisher = {IEEE},
  note      = {Under review}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
