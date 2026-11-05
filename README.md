# MC simulations-based companion for this repo (https://github.com/Art83/als_long)
 
Simulation based method selection for modeling plasma proteins change and functional decline in ALS (Answer ALS consortium). The workflow has three main stages. Simulations are calibrated using the parameters pulled from Answer ALS cohort and compare several methods on power, bias, false-positives etc. The chosen (best) method is then checked on the real cohort against the alternatives. Only after that are biological results produced from the chosen method (this part lives in `Art83/als_long`).

All scripts run from the project root on Windows or Linux-based OS with bash. Parallel work goes through `par_apply()` in `functions.R` (forked workers on Linux, a socket cluster on Windows; Linux is preferable).

Draw timing and protein keys match the paper pipeline in `Art83/als_long`. Draws are timed from the plasma collection date through `R/00_visit_time.R` (which both repos share), and proteins are keyed by SomaScan SeqId. The BLUP analysis in `03` reproduces the paper's 72 progression markers exactly.

## Run order

```
Rscript R/01_estimate_parameters.R 
Rscript R/02_run_simulation.R
Rscript R/02b_onestage.R 
Rscript R/03_validation.R 
Rscript R/04_blup_precision.R 
```

## Methods compared

| Method (code name) | Design |
|---|---|
| TwoStage_BLUP | Patient slopes for ALSFRS-R and protein from random-slope mixed models, then correlated. The chosen method |
| TwoStage_OLS | Per-patient OLS slopes, then correlated |
| Hybrid_BLUPsel_OLSest | BLUP p-value for selection, OLS coefficient for the estimate |
| RegCalibration | OLS estimate divided by the ALSFRS-R slope reliability |
| RepeatedMeasures_LMM | ALSFRS-R against protein at the same visit, participant random intercept |
| CrossSectional | First draw per patient against ALSFRS-R at that visit |
| CrudeRate_DFS | Baseline protein against Delta FS |
| Dichotomise | Baseline protein, fast against slow decliners |
| Baseline_Pred | ALSFRS-R mixed model with a baseline protein by time interaction |
| OneStage_Bivariate | ALSFRS-R and protein fitted jointly with correlated random slopes. Estimation only |
