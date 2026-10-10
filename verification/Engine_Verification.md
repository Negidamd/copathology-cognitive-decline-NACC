# Independent verification of the revised analysis in R and StataMP


## R

- Cohort checks: 17/17 identical (including the participant-visit set).
- Estimates: 68/68 agree (|Δ| < 0.001 for coefficients and CIs; < 0.1 for ΔAIC and the likelihood-ratio statistic; < 1% relative for joint Wald χ²; identical n and visits).
- Largest coefficient difference: 2.61e-06; largest CI-bound difference: 6.04e-04.
- Joint Wald χ² (Python vs R): ADNC 274.6 vs 275.7; LB 26.0 vs 26.1; LATE 1.8 vs 1.8. The small differences come from the fixed-effect covariance estimator (statsmodels vs lme4/Stata); all p-values lead to the same conclusion.
- Values as printed (2 decimals): 67/68 identical to the manuscript. Rounding-boundary cases: m2|1 (0.6853 vs 0.6853)

## StataMP

- Cohort checks: 17/17 identical (including the participant-visit set).
- Estimates: 68/68 agree (|Δ| < 0.001 for coefficients and CIs; < 0.1 for ΔAIC and the likelihood-ratio statistic; < 1% relative for joint Wald χ²; identical n and visits).
- Largest coefficient difference: 2.53e-06; largest CI-bound difference: 6.04e-04.
- Joint Wald χ² (Python vs StataMP): ADNC 274.6 vs 275.7; LB 26.0 vs 26.1; LATE 1.8 vs 1.8. The small differences come from the fixed-effect covariance estimator (statsmodels vs lme4/Stata); all p-values lead to the same conclusion.
- Values as printed (2 decimals): 67/68 identical to the manuscript. Rounding-boundary cases: m2|1 (0.6853 vs 0.6853)
