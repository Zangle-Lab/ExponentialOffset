[README.md](https://github.com/user-attachments/files/34933033/README.md)
# Offset exponential fitting of primary B cell mass

Code and data for:

> Remick S, Mazumder B, Roy K, Zangle TA. Exponential offset model explains mass accumulation during B cell activation. Preprint (2026). DOI: [10.64898/2026.05.22.727310](https://doi.org/10.64898/2026.05.22.727310)

The script fits single-cell mass-over-time tracks from quantitative phase imaging (QPI) of naïve primary mouse B cells. Each track is fit to two models:

- **Pure exponential:** m(t) = A·e^(kt)
- **Offset exponential:** m(t) = A·e^(kt) + C

The script then extracts per-cell parameters: growth rate, initial mass, mass at division, and the offset fraction C/(A+C).

## Upstream tracking

Cell segmentation, mass measurement and tracking were performed with the QPI tracking code from [MultiparametricQPI](https://github.com/eddiethebiochemicalengineer/MultiparametricQPI). That code accompanies the manuscript *Multiparametric quantitative phase imaging for real-time, single cell, drug screening in breast cancer*. Its output is the `tracks` matrix described in [Data](#data). This repository starts from that output and does not include the tracking code.

## Contents

| File | Description |
|---|---|
| `ExponentialOffsetFitting.m` | Analysis script: filters tracks, fits both models, applies quality filters, and builds the per-cell results struct |
| `Founder_tracks.mat` | Single-cell mass tracks (see [Data](#data)). [Hosted on Zenodo: DOI] |

## Requirements

- MATLAB R2019b or later (uses `matches`)
- Signal Processing Toolbox (`medfilt1`)
- Optimization Toolbox (`lsqcurvefit`)
- Statistics and Machine Learning Toolbox (`mad`, `corr`)
- [Violinplot-Matlab](https://github.com/bastibe/Violinplot-Matlab) (`violinplot`), used only for the summary plot of initial mass

## Usage

1. Download `Founder_tracks.mat` and put it in the same folder as the script.
2. In MATLAB, load the data and run the script:

```matlab
load('Founder_tracks.mat')   % loads tracks, tracksColHeaders, drugLabels
ExponentialOffsetFitting
```

The per-cell results are stored in the struct array `valid_cells`.

## Data

`tracks` is a numeric matrix with one row per cell per timepoint, produced by the [MultiparametricQPI](https://github.com/eddiethebiochemicalengineer/MultiparametricQPI) tracking code. The column order is given in `tracksColHeaders`:

| Column | Name | Units / notes |
|---|---|---|
| 1 | X | Centroid x-position (pixels) |
| 2 | Y | Centroid y-position (pixels) |
| 3 | Mass (pg) | Dry mass from QPI |
| 4 | Time (h) | Time since start of imaging |
| 5 | id | Unique cell ID (`cell + 1000 × location`) |
| 6 | Area | Segmented cell area |
| 7 | shape factor | Segmentation shape factor |
| 8 | Location ID | Imaging location |
| 9 | Frame ID | Frame number |
| 10 | segmentLabel | Segmentation label |
| 11 | condition_drugID | Index into `drugLabels` (1 = MZ CpG, 2 = FO CpG, 3 = FO LPS, 4 = MZ LPS) |
| 12 | condition_concentration (um) | Stimulus concentration (µM) |

MZ: marginal zone B cells. FO: follicular B cells. CpG and LPS are the stimuli.

## Analysis steps

1. **Track selection.** Keep founder cells (cells present in the first frame) that have at least 10 timepoints (`minp = 10`).
2. **Smoothing.** Smooth each mass trace with a 5-point median filter (`medfilt1(..., 5, 'truncate')`). Sort points by time and shift time so each track starts at t = 0.
3. **Pure exponential fit.** Fit A and k by linear least squares on log(mass).
4. **Offset exponential fit.** Fit A, k and C by nonlinear least squares on mass (`lsqcurvefit`), with bounds A ≥ 0 and C ≥ 0. Starting values are A = max(m) − min(m), k = 0.05 h⁻¹ and C = min(m).
5. **Goodness of fit.** Compute the residual sum of squares (RSS) and AIC for both models, and R² for the offset exponential fit.
6. **Quality filter.** Keep cells whose offset-exponential fit has R² ≥ 0.75 and whose growth rate k is within 5 median absolute deviations of the median across all cells.
7. **Per-cell output.** Build `valid_cells` with these fields:
   - `id`, `condition`, `conditionNum`
   - `init_mass`: smoothed mass at the first timepoint (m₀)
   - `mass_at_div` and `time_of_div`: mass and time at the last tracked timepoint, taken as division (m_f)
   - `a`, `b`, `c`: offset exponential parameters A, k, C
   - `y_intercept`: fitted initial mass, A + C
   - `r2`, `RSS`

