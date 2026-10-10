# verify_analysis.R — independent re-analysis in R for the npj Dementia R1 revision.
#
# 1. Rebuilds the analytic cohort directly from the raw NACC investigator file (no Python output is used
#    to construct it), following the definitions in the revised Methods.
# 2. Refits every model reported in the manuscript and supplement with lme4 (REML; ML for AIC comparison).
# 3. Writes results/R_results.csv (item, estimate, lo, hi, p, n, obs) and results/R_cohort.csv.
# Comparison with the Python results is done by ../scripts/30_verify_engines.py.
#
# Usage: Rscript verify_analysis.R /path/to/investigator_ftldlbd_nacc72.csv
suppressPackageStartupMessages({library(data.table); library(lme4)})
args <- commandArgs(trailingOnly = TRUE); RAW <- args[1]
stopifnot(file.exists(RAW))
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
dir.create(file.path(here, "results"), showWarnings = FALSE)
set.seed(42)

cols <- c("NACCID","NACCADC","PACKET","NACCVNUM","NACCFDYS","NACCAGE","NACCDAGE","SEX","EDUC","NACCNE4S","CDRSUM",
          "CDRGLOB","NACCMMSE","NPTHAL","NACCBRAA","NACCNEUR","NPADNC","NACCLEWY","NPTDPB","NPTDPC","NPTDPD","NPTDPE",
          "NPFTDTDP","NPFTDTAU","NPHIPSCL","NPINF","NACCARTE","NACCINT")
d <- fread(RAW, select = cols, showProgress = FALSE)
val <- function(x, lo, hi) { x <- suppressWarnings(as.numeric(x)); fifelse(!is.na(x) & x >= lo & x <= hi, x, NA_real_) }
d[, `:=`(braak = val(NACCBRAA,0,6), cerad = val(NACCNEUR,0,3), thal = val(NPTHAL,0,5), adrep = val(NPADNC,0,3),
         lewy = val(NACCLEWY,0,4), ftd = val(NPFTDTDP,0,1), tau = val(NPFTDTAU,0,1), tB = val(NPTDPB,0,1),
         tC = val(NPTDPC,0,1), tD = val(NPTDPD,0,1), tE = val(NPTDPE,0,1))]
setorder(d, NACCID, NACCVNUM)

# ── neuropathology: last record with any core neuropathology item ─────────────────────────────────────────
npm <- d[!is.na(braak) | !is.na(cerad) | !is.na(lewy) | !is.na(adrep)][, .SD[.N], by = NACCID]
npm[, A := c(0,1,1,2,3,3)[thal + 1]][, B := c(0,1,1,2,2,3,3)[braak + 1]][, C := cerad]
npm[, abc := NA_real_]                                         # NIA-AA 2012 Table 3
npm[!is.na(A) & !is.na(B) & !is.na(C) & A == 0 & C == 0, abc := 0]
npm[!is.na(A) & !is.na(B) & !is.na(C) & A >= 1 & B <= 1, abc := 1]
npm[!is.na(A) & !is.na(B) & !is.na(C) & A == 1 & B >= 2 & C <= 1, abc := 1]
npm[!is.na(A) & !is.na(B) & !is.na(C) & A == 1 & B >= 2 & C >= 2, abc := 2]
npm[!is.na(A) & !is.na(B) & !is.na(C) & A >= 2 & B >= 2, abc := 2]
npm[!is.na(A) & !is.na(B) & !is.na(C) & A == 3 & B == 3 & C >= 2, abc := 3]
npm[, ADNC := fifelse(!is.na(adrep), adrep, abc)]
npm[, LB_stage := fifelse(!is.na(lewy) & lewy <= 3, lewy, NA_real_)]
npm[, tdp_n := (!is.na(tB)) + (!is.na(tC)) + (!is.na(tD)) + (!is.na(tE))]
npm[, anypos := (tB %in% 1) | (tC %in% 1) | (tD %in% 1) | (tE %in% 1)]
npm[, allneg := (tB %in% 0) & (tC %in% 0) & (tD %in% 0) & (tE %in% 0)]
npm[, p_TDP := fifelse(anypos | ftd %in% 1, 1, fifelse(allneg & ftd %in% 0, 0, NA_real_))]
npm[, LATE_stage := NA_real_]
npm[tdp_n == 4 & ftd %in% 0, LATE_stage := fifelse(tE == 1, 3, fifelse(tC == 1 | tD == 1, 2, fifelse(tB == 1, 1, 0)))]
npm[, `:=`(FTLD_TDP = ftd, p_FTLDtau = tau, p_ADNC = as.numeric(ADNC >= 2), p_LBD = as.numeric(LB_stage >= 2))]
npm[, n_assessed := (!is.na(p_ADNC)) + (!is.na(p_LBD)) + (!is.na(p_TDP)) + (!is.na(p_FTLDtau))]
npm[, copath_count := fifelse(n_assessed == 4, p_ADNC + p_LBD + p_TDP + p_FTLDtau, NA_real_)]
hs <- npm[, val(NPHIPSCL,0,3)]; npm[, HS := fifelse(is.na(hs), NA_real_, as.numeric(hs > 0))]
inf <- npm[, val(NPINF,0,1)]; art <- npm[, val(NACCARTE,0,3)]
npm[, CVD_any := fifelse((inf %in% 1) | (!is.na(art) & art >= 2), 1, fifelse((inf %in% 0) & !is.na(art) & art < 2, 0, NA_real_))]
npm[, age_death := val(NACCDAGE,18,120)]
npm[, naccint := val(NACCINT,0,600)]      # months from last UDS visit to death
npm[, `:=`(A_score = A, B_score = B, C_score = C, braak_stage = braak)]

# ── visits ───────────────────────────────────────────────────────────────────────────────────────────────
v <- d[NACCID %in% npm$NACCID & PACKET %in% c("I","F","IT")]
v[, `:=`(cdr = val(CDRSUM,0,18), mmse = val(NACCMMSE,0,30), cdrglob = val(CDRGLOB,0,3), days = suppressWarnings(as.numeric(NACCFDYS)))]
v <- v[!is.na(cdr) & !is.na(days) & days >= -0.05*365.25 & days <= 30*365.25]
setorder(v, NACCID, days, NACCVNUM)
v <- v[, .SD[.N], by = .(NACCID, days)]
setorder(v, NACCID, days, NACCVNUM)
v <- v[NACCID %in% v[, .(n = uniqueN(days)), by = NACCID][n >= 2, NACCID]]
first <- v[, .SD[1], by = NACCID][, .(NACCID, d0 = days, age_base = val(NACCAGE,18,120), educ_yrs = val(EDUC,0,36),
          female = fifelse(SEX == 1, 0, fifelse(SEX == 2, 1, NA_real_)), apoe_raw = val(NACCNE4S,0,2), cdr0 = cdr, mmse0 = mmse, glob0 = cdrglob)]
first[, apoe4 := fifelse(is.na(apoe_raw), NA_real_, as.numeric(apoe_raw > 0))]
v <- merge(v[, .(NACCID, NACCADC, days, cdr)], first, by = "NACCID")
v[, years := (days - d0) / 365.25]
v <- merge(v, npm[, .(NACCID, ADNC, LB_stage, LATE_stage, FTLD_TDP, p_ADNC, p_LBD, p_TDP, p_FTLDtau, n_assessed, copath_count,
                       HS, CVD_any, age_death, naccint, A_score, B_score, C_score, braak_stage, tdp_n)], by = "NACCID")
v[, primary := n_assessed == 4 & !is.na(age_base) & !is.na(educ_yrs) & !is.na(female) & !is.na(apoe4)]
b0 <- v[primary == TRUE, .SD[1], by = NACCID]
v[, `:=`(age_c = age_base - mean(b0$age_base), educ_c = educ_yrs - mean(b0$educ_yrs))]
v[, late_stage := fifelse(FTLD_TDP %in% 1, 0, fifelse(FTLD_TDP %in% 0, LATE_stage, NA_real_))]
P <- v[primary == TRUE]
fwrite(P[, .(NACCID, days)], file.path(here, "results", "R_primary_visits.csv"))

# ── cohort and descriptive checks ────────────────────────────────────────────────────────────────────────
comp <- npm[n_assessed == 4]
st <- comp[FTLD_TDP %in% 0 & !is.na(LATE_stage)]
coh <- data.table(item = c("with_neuropath","classifiable","primary_participants","primary_visits","centers",
                           paste0("count_", 0:4), "TDP_n","FTLD_TDP_n","FTLD_TDP_N","LATE_N","LATE_n","ge2_primary"),
  value = c(nrow(npm), nrow(comp), uniqueN(P$NACCID), nrow(P), uniqueN(P$NACCADC),
            sapply(0:4, function(k) P[, .SD[1], by = NACCID][copath_count == k, .N]),
            sum(comp$p_TDP), sum(comp$FTLD_TDP, na.rm = TRUE), sum(!is.na(comp$FTLD_TDP)), nrow(st), sum(st$LATE_stage >= 1),
            P[, .SD[1], by = NACCID][copath_count >= 2, .N]))
fwrite(coh, file.path(here, "results", "R_cohort.csv"))
print(coh)

# ── models ───────────────────────────────────────────────────────────────────────────────────────────────
COV <- "age_c + female + educ_c + apoe4"
out <- list()
fitm <- function(rhs, data, need, reml = TRUE) {
  dd <- data[complete.cases(data[, c("cdr","years","age_c","female","educ_c","apoe4", need), with = FALSE])]
  dd <- dd[NACCID %in% dd[, .(n = uniqueN(years)), by = NACCID][n >= 2, NACCID]]
  f <- as.formula(paste("cdr ~ years * (", rhs, "+", COV, ") + (1 + years | NACCID)"))
  m <- lmer(f, data = dd, REML = reml, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  list(m = m, n = uniqueN(dd$NACCID), obs = nrow(dd))
}
grab <- function(key, f, term) {
  cf <- summary(f$m)$coefficients; est <- cf[term, 1]; se <- cf[term, 2]
  out[[key]] <<- data.table(item = key, estimate = est, lo = est - 1.959964*se, hi = est + 1.959964*se,
                            p = 2*pnorm(-abs(est/se)), n = f$n, obs = f$obs)
}
f1 <- fitm("copath_count", P, "copath_count")
for (t in c("(Intercept)","years","copath_count","years:copath_count","age_c","female","educ_c","apoe4"))
  grab(paste0("m1|", sub("\\(Intercept\\)", "Intercept", t)), f1, t)
P[, cgrp := factor(copath_count, levels = 0:4)]
f2 <- fitm("cgrp", P, "cgrp")
for (k in 1:4) grab(paste0("m2|", k), f2, paste0("years:cgrp", k))
PB <- c("p_ADNC","p_TDP","p_LBD","p_FTLDtau")
f3 <- fitm(paste(PB, collapse = " + "), P, PB)
for (t in PB) { grab(paste0("m3|", t), f3, t); grab(paste0("m3|years:", t), f3, paste0("years:", t)) }
G <- c("ADNC","LB_stage","late_stage","FTLD_TDP","p_FTLDtau")
f4 <- fitm(paste(G, collapse = " + "), P, G)
for (t in G) grab(paste0("m4|years:", t), f4, paste0("years:", t))
# Model 4 with stage categories (reference = stage 0) and joint Wald tests of each pathology's stage x time terms
P[, `:=`(ADNCf = factor(ADNC), LBf = factor(LB_stage), LATEf = factor(late_stage))]
f4c <- fitm("ADNCf + LBf + LATEf + FTLD_TDP + p_FTLDtau", P, G)
for (k in 1:3) for (x in c("ADNC","LB","LATE")) grab(paste0("m4cat|", x, k), f4c, paste0("years:", x, "f", k))
grab("m4cat|FTLD_TDP", f4c, "years:FTLD_TDP"); grab("m4cat|p_FTLDtau", f4c, "years:p_FTLDtau")
for (x in c("ADNC","LB","LATE")) {
  b <- fixef(f4c$m); i <- grep(paste0("^years:", x, "f"), names(b)); V <- as.matrix(vcov(f4c$m))[i, i]
  x2 <- as.numeric(b[i] %*% solve(V, b[i]))
  out[[paste0("m4cat|wald_", x)]] <- data.table(item = paste0("m4cat|wald_", x), estimate = x2, p = pchisq(x2, length(i), lower.tail = FALSE))
}

# ML model comparison on a common sample
Cs <- P[complete.cases(P[, c(PB, G), with = FALSE])]
Cs[, late_bin := as.numeric(late_stage >= 1)]
mls <- lapply(list(count = "copath_count", binary = paste(PB, collapse = " + "),
                   binsep = "p_ADNC + p_LBD + late_bin + FTLD_TDP + p_FTLDtau", graded = paste(G, collapse = " + "),
                   cat = "ADNCf + LBf + LATEf + FTLD_TDP + p_FTLDtau"), function(r) fitm(r, Cs, character(0), reml = FALSE)$m)
aic <- sapply(mls, AIC)
for (k in c("count","binary","binsep","cat"))
  out[[paste0("fit|dAIC_", k)]] <- data.table(item = paste0("fit|dAIC_", k), estimate = aic[k] - aic["graded"])
out[["fit|LR_cat_vs_linear"]] <- data.table(item = "fit|LR_cat_vs_linear", estimate = 2 * as.numeric(logLik(mls$cat) - logLik(mls$graded)))

# LATE-NC by age at death (post hoc)
NOF <- P[FTLD_TDP %in% 0]; LS <- c("LATE_stage","ADNC","LB_stage","p_FTLDtau")
for (s in list(list("All ages", NOF), list("lt80", NOF[age_death < 80]), list("ge80", NOF[age_death >= 80]), list("ge90", NOF[age_death >= 90])))
  grab(paste0("late|", s[[1]]), fitm(paste(LS, collapse = " + "), s[[2]], LS), "years:LATE_stage")
NOF[, old := fifelse(is.na(age_death), NA_real_, as.numeric(age_death >= 80))]
dd <- NOF[complete.cases(NOF[, c("cdr","years","age_c","female","educ_c","apoe4", LS, "old"), with = FALSE])]
dd <- dd[NACCID %in% dd[, .(n = uniqueN(years)), by = NACCID][n >= 2, NACCID]]
mi <- lmer(as.formula(paste("cdr ~ years*LATE_stage*old + years*(ADNC + LB_stage + p_FTLDtau +", COV, ") + (1 + years | NACCID)")),
           data = dd, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
grab("late|interaction", list(m = mi, n = uniqueN(dd$NACCID), obs = nrow(dd)), "years:LATE_stage:old")

# ABC components
ADJ <- c("LB_stage","late_stage","FTLD_TDP","p_FTLDtau")
for (t in c("A_score","B_score","C_score","ADNC"))
  grab(paste0("abc|", t), fitm(paste(c(t, ADJ), collapse = " + "), P, c(t, ADJ)), paste0("years:", t))

# sensitivity analyses
nv <- P[, .N, by = NACCID]; fu <- P[, .(fu = max(years)), by = NACCID]
cnt <- function(key, data, var = "copath_count") grab(paste0("sens|", key), fitm(var, data, var), paste0("years:", var))
cnt("ge3visits", P[NACCID %in% nv[N >= 3, NACCID]])
cnt("ge5years", P[NACCID %in% fu[fu >= 5, NACCID]])
cnt("all_tdp_regions", P[tdp_n == 4])
cnt("death_24m", P[naccint <= 24])                 # last visit within 24 months of death
P[, count_orig := as.numeric(braak_stage >= 3) + as.numeric(LB_stage >= 1) + p_TDP + p_FTLDtau]
cnt("original_thresholds", P, "count_orig")
fc <- fitm("copath_count + CVD_any + HS", P, c("copath_count","CVD_any","HS"))
grab("sens|cvd_hs", fc, "years:copath_count"); grab("cvdhs|years:CVD_any", fc, "years:CVD_any"); grab("cvdhs|years:HS", fc, "years:HS")
v[, count_abs := fcoalesce(p_ADNC, 0) + fcoalesce(p_LBD, 0) + fcoalesce(p_TDP, 0) + fcoalesce(p_FTLDtau, 0)]
cnt("unassessed_absent", v, "count_abs")

# individual rates and Spearman (original definition)
rt <- P[order(NACCID, years), .(rate = { dt <- years[.N] - years[1]; if (dt > 0.5) (cdr[.N] - cdr[1]) / dt else NA_real_ },
                                 g = copath_count[1]), by = NACCID]
for (k in 0:4) out[[paste0("rate|", k)]] <- data.table(item = paste0("rate|", k), estimate = rt[g == k, mean(rate, na.rm = TRUE)])
out[["rate|spearman"]] <- data.table(item = "rate|spearman", estimate = rt[!is.na(rate), cor(g, rate, method = "spearman")])

res <- rbindlist(out, fill = TRUE)
fwrite(res, file.path(here, "results", "R_results.csv"))
writeLines(c(R.version.string, paste("lme4", packageVersion("lme4")), paste("data.table", packageVersion("data.table"))),
           file.path(here, "results", "R_session.txt"))
cat("done:", nrow(res), "estimates\n")
