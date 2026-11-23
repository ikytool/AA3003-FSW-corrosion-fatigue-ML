# Hybrid physics-GP modelling of fatigue crack growth from small datasets: pre-corroded AA3003 friction stir welds

Code and dataset for the manuscript

> Hybrid physics-Gaussian process modelling of fatigue crack growth
> from small datasets: pre-corroded AA3003 friction stir welds
> (Theoretical and Applied Fracture Mechanics, revised version, 2026;
> previous title: Physics-informed machine learning prediction of
> fatigue crack growth in pre-corroded AA3003 friction stir welded
> joints: coupled temperature and inhibitor effects)

This repository holds the version that accompanies the second revision
of the manuscript. It replaces the files of the earlier versions.

## Summary

A condition-wise Paris law forms the physics baseline, the temperature
dependence of the crack-growth rate is interpolated segment by segment
on an Arrhenius scale, and a Gaussian process learns the residual.
The models are assessed by grouped cross-validation in which every
model element (Paris parameters, temperature law, Gaussian process,
interval calibration) is estimated inside each fold. Under this
protocol the hybrid model improves on the physics baseline but is not
more accurate than a multilinear regression or a Gaussian process
alone; models that rely on the temperature law fail when the 25 or
85 degC condition is withheld. Fatigue lives at the untested 45 and
65 degC pre-corrosion temperatures are reported as exploratory
interpolations with parameter, rate-scale and model-form uncertainty.

## Data provenance

No new experiments were performed. All experimental data originate
from the authors' two previously published papers, which remain the
primary reference for the experiments:

- I. Chekalil et al., *Effect of corrosion environments on the
  mechanical properties of friction stir welded aluminum alloy
  AA3003*, J. Mater. Res. Technol. 33 (2024) 2353-2364.
  [10.1016/j.jmrt.2024.09.167](https://doi.org/10.1016/j.jmrt.2024.09.167)
- I. Chekalil et al., *Corrosion behavior of AA3003 friction stir
  welded joints*, Colloids Surf. A 680 (2024) 132673.
  [10.1016/j.colsurfa.2023.132673](https://doi.org/10.1016/j.colsurfa.2023.132673)

The numerical records were obtained by digitizing the published
figures. The dataset contains one curve per specimen (one specimen per
condition); the six pre-corroded specimens used for modelling give 141
crack-growth points. The raw automated digitization is kept unchanged
in `data/digitized/`, and every correction found by a marker-by-marker
audit is listed in `data/digitized/corrections.csv` and applied when
the data are loaded. See `data/digitized/README.md` for details,
including the units of the published rate axis.

## Requirements

Base MATLAB only (developed and tested on R2026a). No toolboxes are
needed; the regression, gradient-boosted-tree and Gaussian-process
models are implemented in the code.

## How to run

```matlab
cd matlab
run_all
```

This regenerates every figure into `figures/` (vector PDF, plus PNG
copies in `figures/png/`) and every derived table into
`data/derived/`. Total runtime is about 5 to 10 minutes; the
cross-validation and bootstrap of stage 4 take most of it.

## Repository structure

| Path | Content |
|---|---|
| `matlab/s0_setup_montage.m` | Fig. 1, setup photographs from the source papers |
| `matlab/s1_paris_baseline.m` | Fig. 2, Table 4: Paris refits, block-bootstrap intervals, effective sample size, digitization check |
| `matlab/s7_framework.m` | Fig. 3, framework diagram |
| `matlab/s2_arrhenius_crosscheck.m` | Fig. 4: global and segment apparent activation energies versus the gravimetric value |
| `matlab/s3_inhibitor_efficiency.m` | Fig. 5: inhibitor efficiency with confidence bands |
| `matlab/s4_ml_hybrid.m` | Figs. 6-8, Tables 6-7, Appendix B: grouped cross-validation, nested interval calibration, paired tests, Gaussian-process hyperparameters |
| `matlab/s5_life_integration.m` | Figs. 9-10, Tables 8-9: rate scale, reconstruction check, life-level blind test, endpoint sensitivity |
| `matlab/s6_probabilistic.m` | Figs. 11-12, Table 10: exploratory predictions at 45 and 65 degC |
| `matlab/run_all.m` | runs all stages in the required order |
| `matlab/README.md` | stage outputs and shared functions in detail |
| `matlab/digitize/` | digitization tool (not needed to reproduce the results) |
| `data/digitized/` | raw digitized records, one CSV per curve, and the audit log `corrections.csv` |
| `data/derived/` | derived tables as shipped, including the Gaussian-process hyperparameters of every fold (`gp_hyperparameters.csv`), the bootstrap draws of the Paris parameters (`paris_boot.mat`) and all predictions (`predictions.csv`, `predictions_45_65.csv`) |
| `data/source_pages/` | renders of the source-paper pages used for digitization and Fig. 1 |
| `figures/` | figures as shipped |

## Changes from the original submission

- Digitization audited marker by marker; corrections logged in
  `data/digitized/corrections.csv`.
- Replotted curves of the same specimens are no longer pooled: one
  curve per specimen (141 modelling points instead of 211).
- Every model element is estimated inside each cross-validation fold;
  prediction intervals are calibrated by nested cross-validation.
- Block bootstrap for the Paris parameters to account for serial
  correlation along each curve.
- Segment-wise Arrhenius temperature interpolation with alternative
  interpolants as a model-form check.
- Life-level blind test, endpoint sensitivity and rate-scale
  uncertainty added; the Monte Carlo design maps of the original
  submission were removed.

## Changes in the second revision

- `matlab/s5_life_integration.m`: the held-out life-level validation
  also reports, for every model, the 95% life band obtained by
  integrating the bounds of its calibrated empirical uncertainty
  interval (new columns `N_*_ei_lo`, `N_*_ei_hi` in `life_blind.csv`).
- `matlab/s6_probabilistic.m`: life-consistent anchors are the primary
  results of Figs. 11 and 12, with the range of the four temperature
  interpolants (new columns `N_lc_form_*`, `gain_lc_form_*`); the common
  rate scale is a sensitivity case. In Fig. 11 both bands are drawn
  before the lines, so the 65 degC band no longer hides the 45 degC
  prediction.
- All previously reported values in the derived tables are unchanged.

## Citation

If you use this code or dataset, please cite the manuscript above and
the two source papers. Citation metadata is in `CITATION.cff`.

## License

The code is released under the MIT License (see `LICENSE`). The
experimental data are described in the two source papers listed above
and are re-used here with attribution.

## Contact

Corresponding author: Sikandar Khan, sikandarkhan@kfupm.edu.sa
