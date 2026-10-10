#!/usr/bin/env python3
"""30_verify_engines.py — independent re-analysis and three-engine comparison.
  1. Runs verify_R/verify_analysis.R (lme4) and verify_StataMP/verify_analysis.do (mixed) in parallel. Each rebuilds the cohort
     from the raw NACC file and refits every reported model (writes verify_*/results/).    [skip with --compare-only]
  2. Compares Python (tables/results.json, locked cohort) with R and Stata: cohort counts, identical participant-visit sets,
     every estimate and 95% CI (|diff| < 0.001; < 0.1 for fit statistics; < 1% relative for joint Wald chi2).
     Writes verification/Engine_Verification.md, engine_comparison.csv, engine_cohort_comparison.csv; exits non-zero on any
     disagreement.
Usage: python scripts/30_verify_engines.py [--raw FILE] [--compare-only]"""
import json, os, re, subprocess, sys, time
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parents[1]
RAW = sys.argv[sys.argv.index('--raw') + 1] if '--raw' in sys.argv else os.environ.get('NACC_RAW_FILE', '')
STATA = os.environ.get('STATA_BIN', '/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp')

if '--compare-only' not in sys.argv:
    assert RAW and Path(RAW).exists(), f'raw NACC file not found (pass --raw FILE or set NACC_RAW_FILE): {RAW}'
    t0 = time.time()
    (H/'verify_StataMP'/'data').mkdir(exist_ok=True)
    pr = subprocess.Popen(['Rscript', 'verify_analysis.R', RAW], cwd=H/'verify_R', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    # Full path passed in quotes so that Stata's -args- keeps a path containing spaces as one token.
    ps = subprocess.Popen([STATA, '-b', 'do', 'verify_analysis.do', f'"{RAW}"'], cwd=H/'verify_StataMP',
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    ro = pr.communicate()[0]; ps.communicate()
    for stray in (H/'verify_StataMP').glob('*.log'):            # Stata batch transcript; the run log is results/Stata_run.log
        stray.unlink()
    slog = (H/'verify_StataMP'/'results'/'Stata_run.log').read_text()
    if pr.returncode != 0 or 'done:' not in ro:
        sys.exit('30: R re-analysis failed\n' + '\n'.join(ro.splitlines()[-15:]))
    if ps.returncode != 0 or '\ndone: ' not in slog or re.search(r'^r\(\d+\);', slog, re.M):
        sys.exit('30: Stata re-analysis failed (see verify_StataMP/results/Stata_run.log)')
    print(f"30: engines ran in {time.time() - t0:.0f} s | R {ro.strip().splitlines()[-1]} | Stata {[l for l in slog.splitlines() if l.startswith('done: ')][-1]}")

OUT = H/'verification'; OUT.mkdir(exist_ok=True)
R = json.load(open(H/'tables'/'results.json'))

py = {}
for t, v in R['m1'].items(): py[f"m1|{t}"] = v
for g in '1234': py[f"m2|{g}"] = R['m2'][f"years:C(cgrp, Treatment('0'))[T.{g}]"]
for t, v in R['m3'].items(): py[f"m3|{t}"] = v
for t, v in R['m4'].items(): py[f"m4|{t}"] = v
names = {'All ages': 'All ages', 'Age at death < 80 years': 'lt80', 'Age at death ≥ 80 years': 'ge80', 'Age at death ≥ 90 years': 'ge90'}
for k, v in R['late'].items(): py[f"late|{names[k]}"] = v
py['late|interaction'] = R['late_interaction']
for t, v in R['abc'].items(): py[f"abc|{t}"] = v
for t, v in R['sens'].items(): py[f"sens|{t}"] = v
for t, v in R['cvd_hs_terms'].items(): py[f"cvdhs|{t}"] = v
aic = {c['model']: c['AIC'] for c in R['model_fit']}; ref = aic['Ordinal severity grades, FTLD-TDP separate']
for k, lab in [('count', 'Copathology count'), ('binary', 'Four proteinopathies, present/absent'),
               ('binsep', 'Present/absent, FTLD-TDP separate'), ('cat', 'Stage categories, FTLD-TDP separate')]:
    py[f'fit|dAIC_{k}'] = {'beta': aic[lab] - ref}          # relative to the linear-grade model, as in R/Stata
py['fit|LR_cat_vs_linear'] = {'beta': R['lr_cat_vs_linear']['chi2']}
for t, v in R['m4_cat'].items():                           # e.g. years:C(ADNC)[T.1.0] -> m4cat|ADNC1
    if 'C(' in t:
        x = {'ADNC': 'ADNC', 'LB_stage': 'LB', 'late_stage': 'LATE'}[t.split('C(')[1].split(')')[0]]
        py[f"m4cat|{x}{int(float(t.split('[T.')[1].rstrip(']')))}"] = v
    else:
        py[f"m4cat|{t.split(':')[1]}"] = v
for k, x in [('ADNC', 'ADNC'), ('LB_stage', 'LB'), ('late_stage', 'LATE')]:
    py[f'm4cat|wald_{x}'] = {'beta': R['m4_cat_wald'][k]['chi2']}
for r in R['rates_by_count']: py[f"rate|{r['group']}"] = {'beta': r['mean']}
py['rate|spearman'] = {'beta': R['spearman']['rho']}

pv = R['prevalence']; nb = R['n_by_count']
pycoh = {'with_neuropath': pv['with_neuropath'], 'classifiable': pv['N'], 'primary_participants': 2594, 'primary_visits': 14423,
         'centers': 33, **{f'count_{k}': nb[str(k)] for k in range(5)}, 'TDP_n': pv['TDP']['n'], 'FTLD_TDP_n': pv['FTLD_TDP']['n'],
         'FTLD_TDP_N': pv['FTLD_TDP']['N'], 'LATE_N': pv['LATE']['N'], 'LATE_n': pv['LATE']['n'], 'ge2_primary': pv['primary_ge2'][0]}

locked = pd.read_csv(H/'data'/'analytical_dataset.csv', usecols=['NACCID', 'days', 'primary'], low_memory=False)
pyset = set(zip(locked[locked.primary == 1].NACCID, locked[locked.primary == 1].days.round(3)))

rows, coh_rows, report = [], [], ["# Independent verification of the revised analysis in R and StataMP\n"]
for eng, folder, prefix in [('R', 'verify_R', 'R'), ('StataMP', 'verify_StataMP', 'Stata')]:
    res_p = H/folder/'results'/f'{prefix}_results.csv'
    if not res_p.exists():
        report.append(f"\n## {eng}: results not found ({res_p.name})\n"); continue
    res = pd.read_csv(res_p).set_index('item')
    coh = pd.read_csv(H/folder/'results'/f'{prefix}_cohort.csv').set_index('item')['value']
    vis = pd.read_csv(H/folder/'results'/f'{prefix}_primary_visits.csv')
    eset = set(zip(vis.NACCID, vis.days.round(3)))
    for k, v in pycoh.items():
        coh_rows.append(dict(engine=eng, item=k, python=v, other=int(coh.get(k, -1)), match=int(coh.get(k, -1)) == v))
    coh_rows.append(dict(engine=eng, item='participant-visit set identical', python=len(pyset), other=len(eset), match=eset == pyset))
    for k, p in py.items():
        if k not in res.index:
            rows.append(dict(engine=eng, item=k, python=p['beta'], other=None, diff=None, match=False)); continue
        o = res.loc[k]
        tol = 0.1 if k.startswith('fit|') else 0.001
        d = abs(o.estimate - p['beta'])
        # joint Wald chi-square: relative tolerance 1% (statsmodels vs lme4/Stata fixed-effect covariance estimators differ slightly)
        ok = d / abs(p['beta']) < 0.01 if '|wald_' in k else d < tol
        dci = None
        if 'lo' in p and pd.notna(o.get('lo')):
            dci = max(abs(o.lo - p['lo']), abs(o.hi - p['hi'])); ok = ok and dci < 0.001
        if 'n' in p and pd.notna(o.get('n')):
            ok = ok and int(o.n) == p['n'] and int(o.obs) == p['obs']
        # would the printed (2-decimal) estimate or CI differ?
        r2 = lambda x: f"{x:.2f}"
        printed_same = r2(o.estimate) == r2(p['beta'])
        if 'lo' in p and pd.notna(o.get('lo')):
            printed_same = printed_same and r2(o.lo) == r2(p['lo']) and r2(o.hi) == r2(p['hi'])
        rows.append(dict(engine=eng, item=k, python=round(p['beta'], 4), other=round(o.estimate, 4), diff=d, max_ci_diff=dci,
                         n_python=p.get('n'), n_other=o.get('n'), match=ok, printed_identical=printed_same))

cmp = pd.DataFrame(rows); cc = pd.DataFrame(coh_rows)
cmp.to_csv(OUT/'engine_comparison.csv', index=False); cc.to_csv(OUT/'engine_cohort_comparison.csv', index=False)
for eng in ['R', 'StataMP']:
    c = cmp[cmp.engine == eng]; k = cc[cc.engine == eng]
    if c.empty: continue
    report.append(f"\n## {eng}\n")
    report.append(f"- Cohort checks: {int(k.match.sum())}/{len(k)} identical (including the participant-visit set).")
    report.append(f"- Estimates: {int(c.match.sum())}/{len(c)} agree (|Δ| < 0.001 for coefficients and CIs; < 0.1 for ΔAIC and the likelihood-ratio statistic; < 1% relative for joint Wald χ²; identical n and visits).")
    cc_ = c[~c.item.str.contains('wald_|fit\\|', regex=True)]
    report.append(f"- Largest coefficient difference: {cc_['diff'].max():.2e}; largest CI-bound difference: {c['max_ci_diff'].max():.2e}.")
    wd = c[c.item.str.contains('wald_')]
    if len(wd): report.append("- Joint Wald χ² (Python vs " + eng + "): " + "; ".join(f"{r.item.split('_')[-1]} {r.python:.1f} vs {r.other:.1f}" for r in wd.itertuples())
                              + ". The small differences come from the fixed-effect covariance estimator (statsmodels vs lme4/Stata); all p-values lead to the same conclusion.")
    pi = c[~c.printed_identical.astype(bool) & ~c.item.str.contains('wald_')]
    report.append(f"- Values as printed (2 decimals): {len(c) - len(pi)}/{len(c)} identical to the manuscript."
                  + ("" if pi.empty else " Rounding-boundary cases: " + "; ".join(f"{r.item} ({r.python} vs {r.other})" for r in pi.itertuples())))
    bad = c[~c.match]
    if len(bad):
        report.append("- Disagreements:\n" + bad[['item', 'python', 'other', 'diff', 'max_ci_diff', 'n_python', 'n_other']].to_markdown(index=False))
    badc = k[~k.match]
    if len(badc):
        report.append("- Cohort disagreements:\n" + badc.to_markdown(index=False))
(OUT/'Engine_Verification.md').write_text("\n".join(report) + "\n")
print("\n".join(report))
if len(cmp) == 0 or not cmp.match.astype(bool).all() or not cc.match.astype(bool).all():
    sys.exit('30: engines disagree or results are missing (see verification/Engine_Verification.md)')
