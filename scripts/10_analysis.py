#!/usr/bin/env python3
"""10_analysis.py — all models and summaries reported in the manuscript and supplement.

Re-runs the original Models 1-3, individual annual rates (original definition: (last - first) CDR-SB
/ follow-up years, follow-up > 0.5 y), the pairwise-rate heatmap and the original sensitivity analyses
on the corrected analytic cohort, and adds only the reviewer-requested analyses:
  * Model 4: ordinal severity grades with FTLD-TDP separate from LATE-NC stage   (R2.2, R2.3, R1.4)
  * model-fit comparison binary vs graded, ML on a common sample                  (R2.3, R2.6)
  * LATE-NC by age at death, post hoc                                             (R1.9)
  * ABC components, each in its own model                                         (R2.2 amyloid)
  * restriction to participants with all four TDP-43 regions assessed             (R2.5)
  * original thresholds (Braak >= III; Lewy any) and unassessed-as-absent cohort   (R2.2, R2.5)
  * visits per annual interval                                                    (R1.7)
  * Model 4 with stage categories (dummy-coded vs stage 0) + omnibus Wald tests
  * last visit within 24 months of death (NACCINT <= 24)

Cohort (locked, data/analytical_dataset.csv, primary == 1): 2,594 participants, 14,423 visits; all four
proteinopathies classifiable; >= 2 visits with valid CDR-SB; complete covariates.
Covariates (age, sex, education, APOE e4) enter as main effects and interactions with time.
REML, L-BFGS (BFGS fallback), random intercept + slope for time. Seed 42.
"""
import json, os, random, warnings
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
warnings.filterwarnings('ignore')
np.random.seed(42); random.seed(42); os.environ["PYTHONHASHSEED"] = "42"

HERE = Path(__file__).resolve().parents[1]
DATA, OUT = HERE/'data', HERE/'tables'
df = pd.read_csv(DATA/'analytical_dataset.csv', low_memory=False)
npm = pd.read_csv(DATA/'neuropath_master.csv', low_memory=False)
P = df[df.primary == 1].copy()
assert P.NACCID.nunique() == 2594 and len(P) == 14423
COV = "age_c + female + educ_c + apoe4"
BASE = ['CDRSUM', 'years', 'age_c', 'female', 'educ_c', 'apoe4']
for d_ in (P, df):
    d_['late_stage'] = d_.LATE_stage.where(d_.FTLD_TDP.eq(0), np.nan)
    d_.loc[d_.FTLD_TDP.eq(1), 'late_stage'] = 0          # structural zero; FTLD-TDP entered separately
    d_['cgrp'] = d_.copath_count.fillna(-1).astype(int).astype(str)
R = {}


def fit(spec, data, extra, reml=True):
    d = data.dropna(subset=BASE + extra)
    vc = d.groupby('NACCID').years.nunique(); d = d[d.NACCID.isin(vc[vc >= 2].index)]
    for method in ('lbfgs', 'bfgs'):
        m = smf.mixedlm(spec, d, groups=d.NACCID, re_formula="~years").fit(reml=reml, method=method, maxiter=1000, disp=False)
        if m.converged and np.isfinite(m.bse_fe).all():
            break
    assert m.converged and np.isfinite(m.bse_fe).all(), spec
    assert np.linalg.eigvalsh(m.cov_re).min() > 0, spec
    return m, d


def terms(m, d, names):
    ci = m.conf_int()
    return {t: dict(beta=float(m.params[t]), lo=float(ci.loc[t, 0]), hi=float(ci.loc[t, 1]),
                    p=float(m.pvalues[t]), n=int(d.NACCID.nunique()), obs=int(len(d))) for t in names}

# ── Table 1 (original rows) ─────────────────────────────────────────────
b = P.sort_values('NACCVNUM').drop_duplicates('NACCID').copy()
b['nvis'] = b.NACCID.map(P.groupby('NACCID').size())
b['g'] = b.copath_count.astype(int).astype(str)
rows = []
for lab, v, kind in [('Age, years, mean (SD)', 'age_base', 'ms'), ('Female sex, n (%)', 'female', 'pct'),
                     ('Education, years, mean (SD)', 'educ_yrs', 'ms'), ('APOE ε4 carrier, n (%)', 'apoe4', 'pct'),
                     ('Baseline CDR-SB, mean (SD)', 'CDRSUM', 'ms'), ('Baseline MMSE, mean (SD)', 'NACCMMSE', 'ms'),
                     ('Visits, median [IQR]', 'nvis', 'medi')]:
    r = {'Characteristic': lab}; vals = []
    for g in '01234':
        x = b.loc[b.g == g, v].dropna(); vals.append(x)
        r[g] = (f"{x.mean():.1f} ({x.std():.1f})" if kind == 'ms' else
                f"{int(x.sum()):,} ({100*x.mean():.1f})" if kind == 'pct' else
                f"{x.median():.0f} [{x.quantile(.25):.0f}–{x.quantile(.75):.0f}]")
    p = stats.chi2_contingency(pd.crosstab(b.g, b[v]))[1] if kind == 'pct' else stats.f_oneway(*vals).pvalue
    r['P'] = "<0.001" if p < 0.001 else f"{p:.3f}"
    rows.append(r)
pd.DataFrame(rows).to_csv(OUT/'table1.csv', index=False)
R['n_by_count'] = {g: int((b.g == g).sum()) for g in '01234'}
R['cdrsb_median_iqr'] = [float(b.CDRSUM.median()), float(b.CDRSUM.quantile(.25)), float(b.CDRSUM.quantile(.75))]
R['n_cdr0'] = int((b.CDRGLOB == 0).sum())
x = b.NACCDAGE_x.dropna()                                  # age at death (Discussion: "82 years here")
R['age_at_death'] = dict(mean=float(x.mean()), sd=float(x.std()), median=float(x.median()), n=int(len(x)))

# ── Prevalence: all 3,391 with four classifiable proteinopathies ───────
comp = npm[npm.n_path_assessed == 4].copy()
pv = {'N': len(comp), 'with_neuropath': len(npm)}
for lab, col in [('ADNC', 'p_ADNC'), ('LBD', 'p_LBD'), ('TDP', 'p_TDP'), ('FTLDtau', 'p_FTLDtau'), ('FTLD_TDP', 'FTLD_TDP')]:
    x = comp[col].dropna(); pv[lab] = dict(n=int(x.sum()), N=len(x), pct=100*x.mean())
st = comp[comp.FTLD_TDP.eq(0) & comp.LATE_stage.notna()]
pv['LATE'] = dict(N=len(st), n=int((st.LATE_stage >= 1).sum()), pct=100*(st.LATE_stage >= 1).mean(),
                  stages={int(k): int(v) for k, v in st.LATE_stage.value_counts().sort_index().items()})
pv['count'] = {int(k): int(v) for k, v in comp.copath_count.value_counts().sort_index().items()}
pv['ge2_pct'] = 100*(comp.copath_count >= 2).mean()
cb = comp.path_combo.value_counts()
pv['combos'] = [[k, int(v), 100*v/len(comp)] for k, v in cb.items()]
pv['primary_ge2'] = [int((b.copath_count >= 2).sum()), 100*(b.copath_count >= 2).mean()]
R['prevalence'] = pv

# ── Models 1-3 (original) ───────────────────────────────────────────────
m1, d1 = fit(f"CDRSUM ~ years * (copath_count + {COV})", P, ['copath_count'])
R['m1'] = terms(m1, d1, ['Intercept', 'years', 'copath_count', 'years:copath_count', 'age_c', 'female', 'educ_c', 'apoe4'])
m2, d2 = fit(f"CDRSUM ~ years * (C(cgrp, Treatment('0')) + {COV})", P, ['cgrp'])
R['m2'] = terms(m2, d2, [f"years:C(cgrp, Treatment('0'))[T.{g}]" for g in '1234'])
PB = ['p_ADNC', 'p_TDP', 'p_LBD', 'p_FTLDtau']
m3, d3 = fit(f"CDRSUM ~ years * ({' + '.join(PB)} + {COV})", P, PB)
R['m3'] = terms(m3, d3, PB + [f'years:{t}' for t in PB])

# ── Model 4: ordinal severity, FTLD-TDP separate (R2.2/R2.3/R1.4) ───────
G = ['ADNC', 'LB_stage', 'late_stage', 'FTLD_TDP', 'p_FTLDtau']
m4, d4 = fit(f"CDRSUM ~ years * ({' + '.join(G)} + {COV})", P, G)
R['m4'] = terms(m4, d4, [f'years:{t}' for t in G])

# Model 4 with stage categories: each stage vs stage 0 (no linearity assumption across stages)
GC = ['C(ADNC)', 'C(LB_stage)', 'C(late_stage)', 'FTLD_TDP', 'p_FTLDtau']
m4c, d4c = fit(f"CDRSUM ~ years * ({' + '.join(GC)} + {COV})", P, G)
tn = [n for n in m4c.fe_params.index if n.startswith('years:') and any(k in n for k in G)]
R['m4_cat'] = terms(m4c, d4c, tn)
def wald(m, idx):
    b_ = m.fe_params[idx].values; V = m.cov_params().loc[idx, idx].values
    x2 = float(b_ @ np.linalg.solve(V, b_)); return dict(chi2=x2, df=len(idx), p=float(stats.chi2.sf(x2, len(idx))))
R['m4_cat_wald'] = {k: wald(m4c, [n for n in tn if f'({k})' in n]) for k in ['ADNC', 'LB_stage', 'late_stage']}

# ── ML comparison on a common sample ────────────────────────────────────
C = P.dropna(subset=BASE + PB + G).copy()
cmp = []
C['late_bin'] = (C.late_stage >= 1).astype(float)
llf = {}
for lab, rhs in [('Copathology count', 'copath_count'), ('Four proteinopathies, present/absent', ' + '.join(PB)),
                 ('Present/absent, FTLD-TDP separate', 'p_ADNC + p_LBD + late_bin + FTLD_TDP + p_FTLDtau'),
                 ('Ordinal severity grades, FTLD-TDP separate', ' + '.join(G)),
                 ('Stage categories, FTLD-TDP separate', ' + '.join(GC))]:
    m, d = fit(f"CDRSUM ~ years * ({rhs} + {COV})", C, [], reml=False)
    cmp.append(dict(model=lab, k=int(len(m.params)), AIC=float(m.aic), BIC=float(m.bic), n=int(d.NACCID.nunique()), obs=len(d)))
    llf[lab] = (float(m.llf), len(m.fe_params))
best = min(c['AIC'] for c in cmp)
for c in cmp: c['dAIC'] = c['AIC'] - best
R['model_fit'] = cmp
# likelihood-ratio test: stage categories vs linear grades (nested; departure from linearity across stages)
(l1, k1), (l0, k0) = llf['Stage categories, FTLD-TDP separate'], llf['Ordinal severity grades, FTLD-TDP separate']
R['lr_cat_vs_linear'] = dict(chi2=2*(l1 - l0), df=k1 - k0, p=float(stats.chi2.sf(2*(l1 - l0), k1 - k0)))

# ── LATE-NC by age at death (post hoc; R1.9) ───────────────────────────
NOF = P[P.FTLD_TDP == 0].copy()
LS = ['LATE_stage', 'ADNC', 'LB_stage', 'p_FTLDtau']
R['late'] = {}
for lab, sub in [('All ages', NOF), ('Age at death < 80 years', NOF[NOF.NACCDAGE_x < 80]),
                 ('Age at death ≥ 80 years', NOF[NOF.NACCDAGE_x >= 80]), ('Age at death ≥ 90 years', NOF[NOF.NACCDAGE_x >= 90])]:
    m, d = fit(f"CDRSUM ~ years * ({' + '.join(LS)} + {COV})", sub, LS)
    R['late'][lab] = terms(m, d, ['years:LATE_stage'])['years:LATE_stage']
NOF['old'] = (NOF.NACCDAGE_x >= 80).astype(float).where(NOF.NACCDAGE_x.notna())
m, d = fit(f"CDRSUM ~ years*LATE_stage*old + years*(ADNC + LB_stage + p_FTLDtau + {COV})", NOF, LS + ['old'])
R['late_interaction'] = terms(m, d, ['years:LATE_stage:old'])['years:LATE_stage:old']

# ── ABC components, each in its own model (R2.2 amyloid) ───────────────
ADJ = ['LB_stage', 'late_stage', 'FTLD_TDP', 'p_FTLDtau']
R['abc'] = {}
for v in ['A_score', 'B_score', 'C_score', 'ADNC']:
    m, d = fit(f"CDRSUM ~ years * ({v} + {' + '.join(ADJ)} + {COV})", P, [v] + ADJ)
    R['abc'][v] = terms(m, d, [f'years:{v}'])[f'years:{v}']

# ── Individual annual rates (original definition) ──────────────────────
def rate(g):
    g = g.sort_values('years'); dt = g.years.iloc[-1] - g.years.iloc[0]
    return (g.CDRSUM.iloc[-1] - g.CDRSUM.iloc[0]) / dt if dt > 0.5 else np.nan
rt = P.groupby('NACCID').apply(rate).rename('rate')
b = b.set_index('NACCID').join(rt)
rs = []
for g in '01234':
    x = b.loc[b.g == g, 'rate'].dropna()
    rs.append(dict(group=int(g), n=len(x), mean=float(x.mean()), sd=float(x.std()), median=float(x.median())))
R['rates_by_count'] = rs
v = b.dropna(subset=['rate'])
rho, prho = stats.spearmanr(v.copath_count, v.rate)
R['spearman'] = dict(rho=float(rho), p=float(prho), n=len(v))
b[['g', 'rate', 'path_combo', 'p_ADNC', 'p_LBD', 'p_TDP', 'p_FTLDtau']].to_csv(OUT/'individual_rates.csv')

# heatmap: diagonal = isolated pathology; off-diagonal = carries both (original definition)
H = []
labs = [('ADNC', 'p_ADNC'), ('LBD', 'p_LBD'), ('TDP-43', 'p_TDP'), ('FTLD-tau', 'p_FTLDtau')]
for li, ci_ in labs:
    for lj, cj in labs:
        if li == lj:
            s = v[(v[ci_] == 1) & (v[[c for _, c in labs if c != ci_]].sum(axis=1) == 0)]
        else:
            s = v[(v[ci_] == 1) & (v[cj] == 1)]
        H.append(dict(row=li, col=lj, n=len(s), mean=float(s.rate.mean())))
R['heatmap'] = H

# ── Sensitivity analyses ────────────────────────────────────────────────
fu = P.groupby('NACCID').years.max(); nv = P.groupby('NACCID').size()
S = {}
def cnt(data, var='copath_count', extra=None):
    m, d = fit(f"CDRSUM ~ years * ({var} + {COV})", data, [var])
    return terms(m, d, [f'years:{var}'])[f'years:{var}']
S['ge3visits'] = cnt(P[P.NACCID.isin(nv[nv >= 3].index)])
S['ge5years'] = cnt(P[P.NACCID.isin(fu[fu >= 5].index)])
S['all_tdp_regions'] = cnt(P[P.tdp_n_assessed.eq(4)])
S['death_24m'] = cnt(P[P.NACCINT <= 24])                  # last visit within 24 months of death
R['naccint'] = dict(n_le24=int((b.NACCINT <= 24).sum()), n_gt24=int((b.NACCINT > 24).sum()),
                    pct_gt24=100*float((b.NACCINT > 24).mean()), median=float(b.NACCINT.median()),
                    q1=float(b.NACCINT.quantile(.25)), q3=float(b.NACCINT.quantile(.75)))
P['count_orig_thr'] = (P.braak_stage >= 3).astype(float).where(P.braak_stage.notna()) + (P.LB_stage >= 1).astype(float) + P.p_TDP + P.p_FTLDtau
S['original_thresholds'] = cnt(P, 'count_orig_thr')
m, d = fit(f"CDRSUM ~ years * (copath_count + CVD_any + HS + {COV})", P, ['copath_count', 'CVD_any', 'HS'])
S['cvd_hs'] = terms(m, d, ['years:copath_count'])['years:copath_count']
R['cvd_hs_terms'] = terms(m, d, ['years:CVD_any', 'years:HS'])
E = df.copy(); E['count_abs'] = E[['p_ADNC', 'p_LBD', 'p_TDP', 'p_FTLDtau']].fillna(0).sum(axis=1)
S['unassessed_absent'] = cnt(E, 'count_abs')
R['sens'] = S

# ── visits per annual interval; fixed-composition subsample (R1.7) ─────
P['ybin'] = P.years.round().astype(int)
R['visits_by_year'] = {int(k): int(v) for k, v in P.groupby('ybin').NACCID.nunique().items()}

json.dump(R, open(OUT/'results.json', 'w'), indent=1, default=float)
print(json.dumps({k: R[k] for k in ['m1', 'm2', 'm3', 'm4', 'm4_cat', 'm4_cat_wald', 'lr_cat_vs_linear', 'naccint', 'model_fit', 'late', 'late_interaction', 'abc', 'rates_by_count', 'spearman', 'sens', 'cvd_hs_terms']}, indent=1, default=lambda x: round(float(x), 4)))
