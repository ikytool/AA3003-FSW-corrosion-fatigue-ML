# Uncertainty-aware modelling of fatigue crack growth in pre-corroded AA3003 FSW joints

Code and consolidated dataset for the manuscript

> Uncertainty-aware modelling of fatigue crack growth in pre-corroded
> AA3003 friction stir welded joints (2026)

A transparent secondary analysis of two published experimental
campaigns. The pre-corrosion temperature is treated strictly as a
conditioning variable of the exposure (the fatigue tests ran in air),
the unique specimen curve is the statistical unit (six pre-corroded
conditions, duplicated replots excluded), and every model element is
refit inside each cross-validation fold. All three temperature
hold-outs are reported; the physical rate unit is recovered from the
crack-length records; life integration is confined to the ASTM
E647-valid crack range; 45/65 degC scenarios carry parameter bands
plus a four-form model-form envelope.

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

## Requirements

Base MATLAB only (developed and tested on R2026a). No toolboxes are
needed; the regression, gradient-boosted-tree and Gaussian-process
models are implemented natively in the stage scripts.

## How to run

```matlab
cd matlab
run_all
```

This regenerates every result figure into `figures/` and every derived
table into `data/derived/`.

## What each script produces

| Script | Manuscript items |
|---|---|
| `s1_paris_baseline.m` | Fig. 2; Table 4 (`paris_refit.csv`, incl. duplicate-consistency check) |
| `s2_arrhenius_crosscheck.m` | Fig. 4; apparent sensitivity Q with bootstrap bands |
| `s3_inhibitor_efficiency.m` | Fig. 5; efficiency surface with joint bootstrap bands (`eta_surface.csv`) |
| `s4_ml_hybrid.m` | Figs. 6-8; Tables 6-7 (`model_metrics.csv`, `fold_contrasts.csv`); LOCO + all three LOTO folds, nested calibration |
| `s5_life_integration.m` | Figs. 9-10; unit recovery (`unit_recovery.csv`), E647 validity (`e647_validity.csv`), closure checks (`life_validation.csv`) |
| `s6_probabilistic.m` | Fig. 11; scenario map with model-form envelope (`design_map.csv`, `scenario_45_65.csv`) |

Helpers: `fig_defaults.m`, `export_vector.m`, `source_data.m` (every
constant taken from the two source papers, including the corrected
corrosion-rate unit).

## Data

`data/digitized/` holds the consolidated numerical dataset, one CSV
per measured curve (see its README for naming and units). The three
FSW.C curves appear in two source figures; the second appearance is
kept only for the extraction-consistency check and enters no
statistic. `data/derived/` holds the pipeline outputs as shipped.

## Citation

If you use this code or dataset, please cite the manuscript and the
two source papers. Citation metadata is in `CITATION.cff`.

## License

The code is released under the MIT License (see `LICENSE`). The
experimental data are described in the two source papers listed above
and are re-used here with attribution.

## Contact

Corresponding author: Sikandar Khan, sikandarkhan@kfupm.edu.sa
