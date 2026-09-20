# MATLAB analysis pipeline

All result figures of the manuscript are generated here as vector PDFs
(Times New Roman, colour-blind-safe palette) into `../figures/`, which
the LaTeX manuscript reads directly.

## Requirements

Base MATLAB only (developed on R2026a). No toolboxes are needed; the
regression, boosted-tree and Gaussian-process models are implemented
natively in the stage scripts.

## Workflow

| Stage | Script | Output |
|---|---|---|
| Style helpers | `fig_defaults.m`, `export_vector.m` | (helpers) |
| Source constants | `source_data.m` | every published number, in one place |
| 0. Setup montage | `s0_setup_montage.m` | `fig_setup.pdf` (Fig. 1, from the published photographs) |
| 1. Paris baseline | `s1_paris_baseline.m` | `fig_paris_published.pdf` (Fig. 2), `data/derived/paris_refit.csv` (Table 4) |
| 2. Arrhenius cross-check | `s2_arrhenius_crosscheck.m` | `fig_arrhenius_crosscheck.pdf` (Fig. 4) |
| 3. Inhibitor efficiency | `s3_inhibitor_efficiency.m` | `fig_inhibitor_efficiency.pdf` (Fig. 5), `data/derived/eta_surface.csv` |
| 4. ML layer | `s4_ml_hybrid.m` | `fig_parity.pdf`, `fig_diagnostics.pdf`, `fig_loto.pdf` (Figs. 6-8), `data/derived/model_metrics.csv` (Table 5) |
| 5. Life integration | `s5_life_integration.m` | `fig_aN_validation.pdf`, `fig_predictions_45_65.pdf` (Figs. 9-10), life tables |
| 6. Probabilistic layer | `s6_probabilistic.m` | `fig_design_maps.pdf` (Fig. 11), `data/derived/design_map.csv` |

Run the whole pipeline with `run_all`. Stage order matters only in that
s1 must run before s4-s6 (they read `data/derived/paris_refit.csv` and
`data/derived/geometry_calibration.csv`).

The consolidated numerical dataset of the two source campaigns lives in
`../data/digitized/` (see the README there for file naming and units).

`digitize/` contains internal tooling used during dataset assembly; it
is not part of the analysis pipeline and is not needed to reproduce any
number in the manuscript.
