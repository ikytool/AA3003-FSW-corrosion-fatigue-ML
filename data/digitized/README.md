# Consolidated dataset

Numerical dataset of the eleven curves of the two source papers
(Chekalil et al., J. Mater. Res. Technol. 33 (2024) 2353-2364 and
Colloids Surf. A 680 (2024) 132673), one CSV file per curve. The two
source papers remain the primary reference for the experiments.

## File naming

`fig<NN>_<condition>[_<temperature>c].csv`

`<NN>` is the figure number in the fatigue source paper:

- `fig12_*` / `fig14_*` - crack-growth curves (da/dN vs dK)
- `fig10_*` / `fig13_*` - crack-length curves (a vs N)

Conditions: `bm` (base metal), `fsw` (as-welded), `fswc`
(pre-corroded), `fswi` (pre-corroded with 10% ethylene glycol).

Example: `fig14_fswi_85c.csv` is the pre-corroded + inhibitor
crack-growth curve at 85 degC.

## Columns

- Crack growth: `dK` (MPa.m^1/2), `dadN` (plotted rate units; the
  factor to mm/cycle is calibrated in `matlab/s5_life_integration.m`)
- Crack length: `N` (cycles), `a` (mm)

Decimal point, not comma; one header row.

## Quality check

`matlab/s1_paris_baseline.m` refits every crack-growth curve and
compares the Paris parameters with the values printed in the source
paper (Table 3 of the manuscript); the refits agree with the published
exponents within about 1-2%. The per-curve residual scatter is carried
into the ML stage as input noise.
