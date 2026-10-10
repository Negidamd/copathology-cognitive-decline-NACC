* verify_analysis.do — independent re-analysis in StataMP 18 for the npj Dementia R1 revision.
* 1. Rebuilds the analytic cohort directly from the raw NACC investigator file (no Python or R output used).
* 2. Refits every model in the manuscript and supplement with -mixed- (REML; ML for AIC comparison),
*    random intercept and slope for time, unstructured covariance.
* 3. Writes results/Stata_results.csv and results/Stata_cohort.csv for ../scripts/30_verify_engines.py.
* Usage: stata-mp -b do verify_analysis.do "/path/to/investigator_ftldlbd_nacc72.csv"
version 18
clear all
set more off
set seed 42
set maxvar 5000
args RAW
capture mkdir results
capture mkdir data
log using results/Stata_run.log, replace text

* ── import (only the items needed) ───────────────────────────────────────────────────────────────────
local keep NACCID NACCADC PACKET NACCVNUM NACCFDYS NACCAGE NACCDAGE SEX EDUC NACCNE4S CDRSUM CDRGLOB NACCMMSE ///
    NPTHAL NACCBRAA NACCNEUR NPADNC NACCLEWY NPTDPB NPTDPC NPTDPD NPTDPE NPFTDTDP NPFTDTAU NPHIPSCL NPINF NACCARTE NACCINT
import delimited using "`RAW'", varnames(1) case(preserve) clear
keep `keep'
foreach v of local keep {
    if !inlist("`v'", "NACCID", "PACKET") capture destring `v', replace force
}
compress
save data/raw_subset.dta, replace

* valid-range helper: gen new = old if inrange(old, lo, hi)
capture program drop vr
program define vr
    args new old lo hi
    gen double `new' = `old' if inrange(`old', `lo', `hi') & !missing(`old')
end

* ── neuropathology: last record with any core item ────────────────────────────────────────────────────
vr braak NACCBRAA 0 6
vr cerad NACCNEUR 0 3
vr thal NPTHAL 0 5
vr adrep NPADNC 0 3
vr lewy NACCLEWY 0 4
vr ftd NPFTDTDP 0 1
vr tau NPFTDTAU 0 1
vr tB NPTDPB 0 1
vr tC NPTDPC 0 1
vr tD NPTDPD 0 1
vr tE NPTDPE 0 1
vr hs NPHIPSCL 0 3
vr inf NPINF 0 1
vr art NACCARTE 0 3
vr age_death NACCDAGE 18 120
vr naccint NACCINT 0 600
keep if !missing(braak) | !missing(cerad) | !missing(lewy) | !missing(adrep)
sort NACCID NACCVNUM
by NACCID: keep if _n == _N
recode thal (0=0) (1 2=1) (3=2) (4 5=3), gen(A)
recode braak (0=0) (1 2=1) (3 4=2) (5 6=3), gen(B)
gen double C = cerad
gen double abc = .
local ok !missing(A) & !missing(B) & !missing(C)
replace abc = 0 if `ok' & A == 0 & C == 0
replace abc = 1 if `ok' & A >= 1 & B <= 1
replace abc = 1 if `ok' & A == 1 & B >= 2 & C <= 1
replace abc = 2 if `ok' & A == 1 & B >= 2 & C >= 2
replace abc = 2 if `ok' & A >= 2 & B >= 2
replace abc = 3 if `ok' & A == 3 & B == 3 & C >= 2
gen double ADNC = cond(!missing(adrep), adrep, abc)
gen double LB_stage = lewy if lewy <= 3 & !missing(lewy)
gen tdp_n = !missing(tB) + !missing(tC) + !missing(tD) + !missing(tE)
gen anypos = tB == 1 | tC == 1 | tD == 1 | tE == 1
gen allneg = tB == 0 & tC == 0 & tD == 0 & tE == 0
gen double p_TDP = .
replace p_TDP = 0 if allneg & ftd == 0
replace p_TDP = 1 if anypos | ftd == 1
gen double LATE_stage = .
replace LATE_stage = cond(tE == 1, 3, cond(tC == 1 | tD == 1, 2, cond(tB == 1, 1, 0))) if tdp_n == 4 & ftd == 0
gen double FTLD_TDP = ftd
gen double p_FTLDtau = tau
gen double p_ADNC = (ADNC >= 2) if !missing(ADNC)
gen double p_LBD = (LB_stage >= 2) if !missing(LB_stage)
gen n_assessed = !missing(p_ADNC) + !missing(p_LBD) + !missing(p_TDP) + !missing(p_FTLDtau)
gen double copath_count = p_ADNC + p_LBD + p_TDP + p_FTLDtau if n_assessed == 4
gen double HS = (hs > 0) if !missing(hs)
gen double CVD_any = .
replace CVD_any = 0 if inf == 0 & art < 2 & !missing(art)
replace CVD_any = 1 if inf == 1 | (art >= 2 & !missing(art))
rename (A B C braak) (A_score B_score C_score braak_stage)
keep NACCID ADNC LB_stage LATE_stage FTLD_TDP p_ADNC p_LBD p_TDP p_FTLDtau n_assessed copath_count HS CVD_any ///
     age_death naccint A_score B_score C_score braak_stage tdp_n
save data/npm.dta, replace

* cohort checks on the neuropathology file
count
local with_np = r(N)
count if n_assessed == 4
local classifiable = r(N)
count if n_assessed == 4 & p_TDP == 1
local tdp_n = r(N)
count if n_assessed == 4 & FTLD_TDP == 1
local ftdp_n = r(N)
count if n_assessed == 4 & !missing(FTLD_TDP)
local ftdp_N = r(N)
count if n_assessed == 4 & FTLD_TDP == 0 & !missing(LATE_stage)
local late_N = r(N)
count if n_assessed == 4 & FTLD_TDP == 0 & !missing(LATE_stage) & LATE_stage >= 1
local late_n = r(N)

* ── visits ─────────────────────────────────────────────────────────────────────────────────────────────
use data/raw_subset.dta, clear
merge m:1 NACCID using data/npm.dta, keep(match) nogenerate keepusing(NACCID)
keep if inlist(PACKET, "I", "F", "IT")
vr cdr CDRSUM 0 18
gen double days = NACCFDYS
keep if !missing(cdr) & !missing(days) & days >= -0.05*365.25 & days <= 30*365.25
sort NACCID days NACCVNUM
by NACCID days: keep if _n == _N
by NACCID: gen nvis_days = _N
keep if nvis_days >= 2
sort NACCID days
by NACCID: gen double d0 = days[1]
foreach v in NACCAGE EDUC SEX NACCNE4S {
    by NACCID: gen double b_`v' = `v'[1]
}
gen double age_base = b_NACCAGE if inrange(b_NACCAGE, 18, 120)
gen double educ_yrs = b_EDUC if inrange(b_EDUC, 0, 36)
gen double female = cond(b_SEX == 1, 0, cond(b_SEX == 2, 1, .))
gen double apoe4 = (b_NACCNE4S > 0) if inrange(b_NACCNE4S, 0, 2)
gen double years = (days - d0) / 365.25
merge m:1 NACCID using data/npm.dta, keep(match) nogenerate
gen primary = n_assessed == 4 & !missing(age_base, educ_yrs, female, apoe4)
bysort NACCID (days): gen first = _n == 1
summarize age_base if primary & first, meanonly
gen double age_c = age_base - r(mean)
summarize educ_yrs if primary & first, meanonly
gen double educ_c = educ_yrs - r(mean)
gen double late_stage = cond(FTLD_TDP == 1, 0, cond(FTLD_TDP == 0, LATE_stage, .))
gen double count_abs = cond(missing(p_ADNC), 0, p_ADNC) + cond(missing(p_LBD), 0, p_LBD) + ///
    cond(missing(p_TDP), 0, p_TDP) + cond(missing(p_FTLDtau), 0, p_FTLDtau)
gen double count_orig = (braak_stage >= 3) + (LB_stage >= 1) + p_TDP + p_FTLDtau if !missing(braak_stage, LB_stage)
egen long pid = group(NACCID)
save data/visits_all.dta, replace

* primary-cohort checks
count if primary
local p_visits = r(N)
count if primary & first
local p_n = r(N)
quietly levelsof NACCADC if primary, local(ctr)
local centers : word count `ctr'
forvalues k = 0/4 {
    count if primary & first & copath_count == `k'
    local c`k' = r(N)
}
count if primary & first & copath_count >= 2
local ge2 = r(N)
tempname ch
postfile `ch' str30 item double value using results/Stata_cohort.dta, replace
foreach k in with_np classifiable p_n p_visits centers c0 c1 c2 c3 c4 tdp_n ftdp_n ftdp_N late_N late_n ge2 {
    local lab `k'
    if "`k'" == "with_np" local lab with_neuropath
    if "`k'" == "p_n" local lab primary_participants
    if "`k'" == "p_visits" local lab primary_visits
    if inlist("`k'", "c0", "c1", "c2", "c3", "c4") local lab = "count_" + substr("`k'", 2, 1)
    if "`k'" == "tdp_n" local lab TDP_n
    if "`k'" == "ftdp_n" local lab FTLD_TDP_n
    if "`k'" == "ftdp_N" local lab FTLD_TDP_N
    if "`k'" == "late_N" local lab LATE_N
    if "`k'" == "late_n" local lab LATE_n
    if "`k'" == "ge2" local lab ge2_primary
    post `ch' ("`lab'") (``k'')
}
postclose `ch'
preserve
use results/Stata_cohort.dta, clear
export delimited using results/Stata_cohort.csv, replace
restore
preserve
keep if primary
keep NACCID days
export delimited using results/Stata_primary_visits.csv, replace
restore

* ── models ─────────────────────────────────────────────────────────────────────────────────────────────
tempname post
postfile `post' str40 item double(estimate lo hi p n obs) using results/Stata_results.dta, replace
local COV c.age_c c.female c.educ_c c.apoe4

* fit: args = sample condition, rhs (factor-variable list), required vars, reml|ml
capture program drop fitm
program define fitm
    args cond rhs need method
    capture drop _use _nuniq
    gen _use = (`cond') & !missing(cdr, years, age_c, female, educ_c, apoe4 `need')
    bysort pid _use (years): gen _nuniq = sum(years != years[_n-1]) if _use
    bysort pid _use (_nuniq): replace _use = 0 if _nuniq[_N] < 2
    quietly mixed cdr c.years##(`rhs' c.age_c c.female c.educ_c c.apoe4) if _use || pid: years, cov(unstructured) `method' iterate(300)
end
capture program drop grab
program define grab
    args post key term
    local b = _b[`term']
    local s = _se[`term']
    post `post' ("`key'") (`b') (`b' - invnormal(0.975)*`s') (`b' + invnormal(0.975)*`s') (2*normal(-abs(`b'/`s'))) (el(e(N_g),1,1)) (e(N))
end

use data/visits_all.dta, clear
keep if primary
save data/primary.dta, replace

fitm 1 "c.copath_count" ", copath_count" reml
grab `post' "m1|Intercept" _cons
grab `post' "m1|years" years
grab `post' "m1|copath_count" copath_count
grab `post' "m1|years:copath_count" c.years#c.copath_count
foreach v in age_c female educ_c apoe4 {
    grab `post' "m1|`v'" `v'
}
fitm 1 "ib0.copath_count" ", copath_count" reml
forvalues k = 1/4 {
    grab `post' "m2|`k'" c.years#`k'.copath_count
}
fitm 1 "c.p_ADNC c.p_TDP c.p_LBD c.p_FTLDtau" ", p_ADNC, p_TDP, p_LBD, p_FTLDtau" reml
foreach v in p_ADNC p_TDP p_LBD p_FTLDtau {
    grab `post' "m3|`v'" `v'
    grab `post' "m3|years:`v'" c.years#c.`v'
}
local G c.ADNC c.LB_stage c.late_stage c.FTLD_TDP c.p_FTLDtau
fitm 1 "`G'" ", ADNC, LB_stage, late_stage, FTLD_TDP, p_FTLDtau" reml
foreach v in ADNC LB_stage late_stage FTLD_TDP p_FTLDtau {
    grab `post' "m4|years:`v'" c.years#c.`v'
}
* Model 4 with stage categories (reference = stage 0) and joint Wald tests of each pathology's stage x time terms
local GC i.ADNC i.LB_stage i.late_stage c.FTLD_TDP c.p_FTLDtau
fitm 1 "`GC'" ", ADNC, LB_stage, late_stage, FTLD_TDP, p_FTLDtau" reml
forvalues k = 1/3 {
    grab `post' "m4cat|ADNC`k'" c.years#`k'.ADNC
    grab `post' "m4cat|LB`k'" c.years#`k'.LB_stage
    grab `post' "m4cat|LATE`k'" c.years#`k'.late_stage
}
grab `post' "m4cat|FTLD_TDP" c.years#c.FTLD_TDP
grab `post' "m4cat|p_FTLDtau" c.years#c.p_FTLDtau
foreach v in ADNC LB_stage late_stage {
    quietly test c.years#1.`v' c.years#2.`v' c.years#3.`v'
    local lab = cond("`v'" == "ADNC", "ADNC", cond("`v'" == "LB_stage", "LB", "LATE"))
    post `post' ("m4cat|wald_`lab'") (r(chi2)) (.) (.) (r(p)) (.) (.)
}

* ML comparison on a common sample
gen _common = !missing(p_ADNC, p_TDP, p_LBD, p_FTLDtau, ADNC, LB_stage, late_stage, FTLD_TDP)
gen double late_bin = (late_stage >= 1) if !missing(late_stage)
foreach m in count binary binsep graded cat {
    if "`m'" == "count"  local r c.copath_count
    if "`m'" == "binary" local r c.p_ADNC c.p_TDP c.p_LBD c.p_FTLDtau
    if "`m'" == "binsep" local r c.p_ADNC c.p_LBD c.late_bin c.FTLD_TDP c.p_FTLDtau
    if "`m'" == "graded" local r `G'
    if "`m'" == "cat"    local r `GC'
    fitm "_common" "`r'" "" ml
    local ll_`m' = e(ll)
    quietly estat ic
    matrix S = r(S)
    local aic_`m' = S[1,5]
}
foreach m in count binary binsep cat {
    post `post' ("fit|dAIC_`m'") (`aic_`m'' - `aic_graded') (.) (.) (.) (.) (.)
}
post `post' ("fit|LR_cat_vs_linear") (2*(`ll_cat' - `ll_graded')) (.) (.) (.) (.) (.)

* LATE-NC by age at death (post hoc)
local LS c.LATE_stage c.ADNC c.LB_stage c.p_FTLDtau
fitm "FTLD_TDP == 0" "`LS'" ", LATE_stage, ADNC, LB_stage, p_FTLDtau" reml
grab `post' "late|All ages" c.years#c.LATE_stage
fitm "FTLD_TDP == 0 & age_death < 80" "`LS'" ", LATE_stage, ADNC, LB_stage, p_FTLDtau" reml
grab `post' "late|lt80" c.years#c.LATE_stage
fitm "FTLD_TDP == 0 & age_death >= 80 & !missing(age_death)" "`LS'" ", LATE_stage, ADNC, LB_stage, p_FTLDtau" reml
grab `post' "late|ge80" c.years#c.LATE_stage
fitm "FTLD_TDP == 0 & age_death >= 90 & !missing(age_death)" "`LS'" ", LATE_stage, ADNC, LB_stage, p_FTLDtau" reml
grab `post' "late|ge90" c.years#c.LATE_stage
gen double old = (age_death >= 80) if !missing(age_death)
capture drop _use _nuniq
gen _use = FTLD_TDP == 0 & !missing(cdr, years, age_c, female, educ_c, apoe4, LATE_stage, ADNC, LB_stage, p_FTLDtau, old)
bysort pid _use (years): gen _nuniq = sum(years != years[_n-1]) if _use
bysort pid _use (_nuniq): replace _use = 0 if _nuniq[_N] < 2
quietly mixed cdr c.years##c.LATE_stage##c.old c.years##(c.ADNC c.LB_stage c.p_FTLDtau `COV') if _use || pid: years, cov(unstructured) reml
grab `post' "late|interaction" c.years#c.LATE_stage#c.old

* ABC components
local ADJ c.LB_stage c.late_stage c.FTLD_TDP c.p_FTLDtau
foreach v in A_score B_score C_score ADNC {
    fitm 1 "c.`v' `ADJ'" ", `v', LB_stage, late_stage, FTLD_TDP, p_FTLDtau" reml
    grab `post' "abc|`v'" c.years#c.`v'
}

* sensitivity analyses
bysort pid: gen _nv = _N
bysort pid (years): gen double _fu = years[_N]
fitm "_nv >= 3" "c.copath_count" ", copath_count" reml
grab `post' "sens|ge3visits" c.years#c.copath_count
fitm "_fu >= 5" "c.copath_count" ", copath_count" reml
grab `post' "sens|ge5years" c.years#c.copath_count
fitm "tdp_n == 4" "c.copath_count" ", copath_count" reml
grab `post' "sens|all_tdp_regions" c.years#c.copath_count
fitm "naccint <= 24" "c.copath_count" ", copath_count" reml
grab `post' "sens|death_24m" c.years#c.copath_count
fitm 1 "c.count_orig" ", count_orig" reml
grab `post' "sens|original_thresholds" c.years#c.count_orig
fitm 1 "c.copath_count c.CVD_any c.HS" ", copath_count, CVD_any, HS" reml
grab `post' "sens|cvd_hs" c.years#c.copath_count
grab `post' "cvdhs|years:CVD_any" c.years#c.CVD_any
grab `post' "cvdhs|years:HS" c.years#c.HS
use data/visits_all.dta, clear
fitm 1 "c.count_abs" ", count_abs" reml
grab `post' "sens|unassessed_absent" c.years#c.count_abs

* individual rates (original definition) and Spearman
use data/primary.dta, clear
sort pid years
by pid: gen double _dt = years[_N] - years[1]
by pid: gen double rate = (cdr[_N] - cdr[1]) / _dt if _dt > 0.5
by pid: keep if _n == 1
forvalues k = 0/4 {
    summarize rate if copath_count == `k', meanonly
    post `post' ("rate|`k'") (r(mean)) (.) (.) (.) (.) (.)
}
spearman copath_count rate
post `post' ("rate|spearman") (r(rho)) (.) (.) (.) (.) (.)
postclose `post'

use results/Stata_results.dta, clear
export delimited using results/Stata_results.csv, replace
display "done: " _N " estimates; Stata " c(stata_version) " " c(edition_real)
log close
