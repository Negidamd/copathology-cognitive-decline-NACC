# Co-occurring neuropathologies and the rate of cognitive decline

Analysis code for: Negida A. Co-occurring Neuropathologies and the Rate of Cognitive Decline in Autopsy-Confirmed
Neurodegenerative Disease. *npj Dementia* (submitted).

Author: Ahmed Negida, Department of Neurology, Virginia Commonwealth University.

## Contents
| Folder | What it contains |
|---|---|
| scripts/ | Python pipeline: 00 cohort derivation from the raw NACC file; 10 all models and summaries; 20 figures; 30 runs the R and Stata re-analyses and compares all three engines |
| verify_R/ | Separate re-analysis in independently written R code (lme4): cohort rebuilt from the raw file, every reported model refitted; results/ holds the aggregate output |
| verify_StataMP/ | Separate re-analysis in independently written StataMP 18 code (mixed), built the same way; results/ holds the aggregate output and run log |
| tables/ | Aggregate Python results (results.json holds every estimate reported in the manuscript and supplement) |
| verification/ | Three-engine comparison: 17/17 cohort checks and 68/68 estimates agree (joint Wald χ² within 0.4%); largest R − Stata difference in coefficients and CIs 1.3e-07 |
| figures/ | Figures 1-4 (vector PDF) |

## Data
Individual-level NACC data are not distributed. Request the investigator file (UDS with linked FTLD and LBD
modules; December 2025 freeze, distributed March 2026) from https://naccdata.org. The analysis used
investigator_ftldlbd_nacc72.csv, SHA-256 7922a9e850d4449c39370b7c3972e2f34aef22cf5f6ac23f1935071388297ff5.

## Run order
    python scripts/00_cohort.py --raw /path/to/investigator_ftldlbd_nacc72.csv                   # writes data/
    python scripts/10_analysis.py                                                                # tables/results.json
    python scripts/20_figures.py                                                                 # figures/
    python scripts/30_verify_engines.py --raw /path/to/investigator_ftldlbd_nacc72.csv           # R + Stata, then verification/

Scripts 00-20 write participant-level files (data/, tables/individual_rates.csv, verify_*/results/*_primary_visits.csv);
keep these local under the NACC data use agreement.

## Software
Python 3.13.9 (pandas 2.3.3, numpy 2.3.5, scipy 1.16.3, statsmodels 0.14.5, matplotlib 3.10.6, seaborn 0.13.2);
R version 4.5.2 (2025-10-31), lme4 1.1.38, data.table 1.17.0; Stata 18 MP. Random seed 42.
Script 30 looks for Stata at the macOS default location; set STATA_BIN to point elsewhere, and NACC_RAW_FILE (or --raw)
to the investigator file.

## Model specification
Linear mixed-effects models, random intercept and slope for time (unstructured covariance), REML (ML for AIC
comparisons on a common sample). Baseline age, sex, education and APOE e4 enter as main effects and interactions with time.
Model 4 enters ADNC ABC score, Lewy body stage and LATE-NC stage as categories (reference stage 0) with FTLD-TDP and
FTLD-tauopathy as present/absent; joint Wald tests of each pathology's stage x time terms; a likelihood-ratio test compares
categorical with linear stage coding. Sensitivity analysis SA7 restricts to participants whose last visit was within
24 months of death (NACCINT <= 24).

## License
MIT (see LICENSE). NACC data are governed by the NACC data use agreement and are not covered by this license.
