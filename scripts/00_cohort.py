#!/usr/bin/env python3
"""00_cohort.py — the analytic cohort from the raw NACC investigator file (UDS + FTLD and LBD modules).
Writes data/analytical_dataset.csv, data/neuropath_master.csv and data/cohort_flow.json. An existing data/ is not
overwritten unless --overwrite-locked is given; use --outdir to rebuild elsewhere and compare.
Usage: python scripts/00_cohort.py --raw FILE [--outdir DIR]"""
from pathlib import Path
import numpy as np, pandas as pd, json, hashlib, os, argparse
BASE=Path(__file__).resolve().parent.parent
parser=argparse.ArgumentParser(description='Rebuild the corrected cohort from an authorized NACC investigator CSV.')
parser.add_argument('--raw',default=os.environ.get('NACC_RAW_FILE'),help='Raw investigator CSV path; alternatively set NACC_RAW_FILE.')
parser.add_argument('--outdir',default=None,help='Write data/ here instead of the project folder (to compare a rebuild with an existing dataset).')
parser.add_argument('--overwrite-locked',action='store_true',help='Replace the LOCKED project data/ (normally never needed).')
args=parser.parse_args()
OUT=Path(args.outdir).expanduser().resolve() if args.outdir else BASE
if OUT==BASE and (BASE/'data'/'analytical_dataset.csv').exists() and not args.overwrite_locked:
    parser.error('data/analytical_dataset.csv is LOCKED; use --outdir DIR to rebuild elsewhere and compare, or --overwrite-locked.')
if not args.raw:parser.error('--raw or NACC_RAW_FILE is required; individual-level NACC data are not distributed with code.')
RAW=Path(args.raw).expanduser().resolve()
COLS=['NACCID','NACCADC','PACKET','NACCVNUM','NACCFDYS','NACCAGE','NACCDAGE','SEX','EDUC','NACCNE4S','CDRSUM','CDRGLOB','NACCMMSE','NPFORMVER','NPTHAL','NACCBRAA','NACCNEUR','NPADNC','NACCLEWY','NPLBOD','NPTDPA','NPTDPB','NPTDPC','NPTDPD','NPTDPE','NPFTDTDP','NPFTDTAU','NPHIPSCL','NPINF','NACCARTE','NACCPICK','NACCCBD','NACCPROG','NPFTDT2','NPFTDT5','NPFTDT6','NPFTDT7','NPFTDT8','NPFTDT9','NPFTDT10','NACCINT']
(OUT/'data').mkdir(parents=True,exist_ok=True)
def valid(s,lo,hi):
    s=pd.to_numeric(s,errors='coerce');return s.where(s.between(lo,hi))
def old_abc(a,b,c):
    return np.where(a.isna()|b.isna()|c.isna(),np.nan,np.where((a==0)&(b==0)&(c==0),0,np.where(b<=1,np.where(a>=1,1,0),np.where(b==2,np.where(c>=1,2,1),np.where(c>=2,3,2)))))
def correct_abc(a,b,c):
    # NIA-AA 2012 Table 3. A0/C>0 combinations are inconsistent, not scored.
    ans=pd.Series(np.nan,index=a.index)
    ok=a.notna()&b.notna()&c.notna()
    ans[ok&(a==0)&(c==0)]=0
    ans[ok&(a>=1)&(b<=1)]=1
    ans[ok&(a==1)&(b>=2)&(c<=1)]=1
    ans[ok&(a==1)&(b>=2)&(c>=2)]=2
    ans[ok&(a>=2)&(b>=2)]=2
    ans[ok&(a==3)&(b==3)&(c>=2)]=3
    return ans
print('Reading raw data',flush=True)
d=pd.read_csv(RAW,usecols=COLS,low_memory=False)
print('Read',d.shape,flush=True)
for new,old,lo,hi in [('braak_stage','NACCBRAA',0,6),('cerad','NACCNEUR',0,3),('thal','NPTHAL',0,5),('reported_adnc','NPADNC',0,3),('lewy','NACCLEWY',0,4),('FTLD_TDP','NPFTDTDP',0,1),('p_FTLDtau','NPFTDTAU',0,1),('tB','NPTDPB',0,1),('tC','NPTDPC',0,1),('tD','NPTDPD',0,1),('tE','NPTDPE',0,1)]: d[new]=valid(d[old],lo,hi)
d=d.sort_values(['NACCID','NACCVNUM'])
npm=d[d[['braak_stage','cerad','lewy','reported_adnc']].notna().any(axis=1)].drop_duplicates('NACCID',keep='last').copy()
npm['A_score']=npm.thal.map({0:0,1:1,2:1,3:2,4:3,5:3});npm['B_score']=npm.braak_stage.map({0:0,1:1,2:1,3:2,4:2,5:3,6:3});npm['C_score']=npm.cerad
npm['abc_old']=old_abc(npm.A_score,npm.B_score,npm.C_score);npm['abc_corrected']=correct_abc(npm.A_score,npm.B_score,npm.C_score)
npm['ADNC_old']=npm.reported_adnc.fillna(npm.abc_old);npm['ADNC']=npm.reported_adnc.fillna(npm.abc_corrected)
npm['LB_stage']=npm.lewy.where(npm.lewy.between(0,3))
reg=['tB','tC','tD','tE'];npm['tdp_n_assessed']=npm[reg].notna().sum(axis=1);npm['tdp_any_positive']=npm[reg].eq(1).any(axis=1)
npm['old_tdp']=np.where((npm.tdp_n_assessed==0)&npm.FTLD_TDP.isna(),np.nan,(npm.tdp_any_positive|npm.FTLD_TDP.eq(1)).astype(float))
# Positive evidence determines presence; absence requires all four regions negative and no FTLD-TDP.
npm['p_TDP']=np.nan
npm.loc[npm.tdp_any_positive|npm.FTLD_TDP.eq(1),'p_TDP']=1
npm.loc[npm[reg].eq(0).all(axis=1)&npm.FTLD_TDP.eq(0),'p_TDP']=0
npm['LATE_stage']=np.nan
# Require all four measured regions and a recorded negative FTLD-TDP diagnosis.
complete_reg=npm[reg].notna().all(axis=1)&npm.FTLD_TDP.eq(0)
npm.loc[complete_reg,'LATE_stage']=np.select([npm.tE.eq(1),npm.tC.eq(1)|npm.tD.eq(1),npm.tB.eq(1)],[3,2,1],default=0)[complete_reg]
npm['LATE_NC']=npm.LATE_stage.ge(1).astype(float).where(npm.LATE_stage.notna())
npm['p_ADNC']=npm.ADNC.ge(2).astype(float).where(npm.ADNC.notna());npm['p_LBD']=npm.LB_stage.ge(2).astype(float).where(npm.LB_stage.notna())
path=['p_ADNC','p_LBD','p_TDP','p_FTLDtau'];npm['n_path_assessed']=npm[path].notna().sum(axis=1);npm['copath_count']=npm[path].sum(axis=1,min_count=4)
oldp=pd.DataFrame({'a':npm.ADNC_old.ge(2).astype(float).where(npm.ADNC_old.notna()),'b':npm.p_LBD,'c':npm.old_tdp,'d':npm.p_FTLDtau});npm['old_complete']=oldp.notna().all(axis=1);npm['old_count']=oldp.sum(axis=1,min_count=4)
# CVD/HS only sensitivity, no missing codes interpreted as absence.
hs=valid(npm.NPHIPSCL,0,3);npm['HS']=hs.gt(0).astype(float).where(hs.notna())
inf=valid(npm.NPINF,0,1);art=valid(npm.NACCARTE,0,3);npm['CVD_any']=np.where(inf.eq(1)|art.ge(2),1,np.where(inf.eq(0)&art.lt(2),0,np.nan))
npm['path_combo']=npm.apply(lambda r:' + '.join(lab for lab,col in zip(['ADNC','LBD','TDP-43','Other tauopathy'],path) if r[col]==1) if r.n_path_assessed==4 else 'Incomplete assessment',axis=1)
npm.loc[npm.copath_count.eq(0),'path_combo']='None above threshold'
vis=d[d.NACCID.isin(npm.NACCID)&d.PACKET.isin(['I','F','IT'])].copy()
vis['CDRSUM']=valid(vis.CDRSUM,0,18);vis['NACCMMSE']=valid(vis.NACCMMSE,0,30);vis['CDRGLOB']=valid(vis.CDRGLOB,0,3)
vis['days']=pd.to_numeric(vis.NACCFDYS,errors='coerce')
vis=vis[vis.CDRSUM.notna()&vis.days.between(-.05*365.25,30*365.25)].sort_values(['NACCID','days','NACCVNUM'])
dups=vis.duplicated(['NACCID','days']).sum();vis=vis.drop_duplicates(['NACCID','days'],keep='last')
vc=vis.groupby('NACCID').days.nunique();vis=vis[vis.NACCID.isin(vc[vc>=2].index)].copy()
# Landmark first eligible CDR-SB visit: consistent baseline age, education, sex, genotype and time.
first=vis.drop_duplicates('NACCID').copy();first['age_base']=valid(first.NACCAGE,18,120);first['educ_yrs']=valid(first.EDUC,0,36);first['female']=first.SEX.map({1:0.,2:1.});first['apoe4']=valid(first.NACCNE4S,0,2).gt(0).astype(float).where(valid(first.NACCNE4S,0,2).notna())
first=first.set_index('NACCID')
vis['years']=(vis.days-vis.NACCID.map(first.days))/365.25
for c in ['age_base','educ_yrs','female','apoe4']: vis[c]=vis.NACCID.map(first[c])
npcols=['NACCID','ADNC','ADNC_old','LB_stage','FTLD_TDP','LATE_stage','LATE_NC','p_ADNC','p_LBD','p_TDP','p_FTLDtau','n_path_assessed','copath_count','path_combo','HS','CVD_any','A_score','B_score','C_score','braak_stage','NACCDAGE','NPFORMVER','tdp_n_assessed','old_complete','old_count']
vis=vis.drop(columns=[c for c in npcols if c!='NACCID' and c in vis],errors='ignore').merge(npm[npcols],on='NACCID',validate='many_to_one')
vis['NACCDAGE_x']=valid(vis.NACCDAGE,18,120)
vis['path_complete']=vis.n_path_assessed.eq(4).astype(int)
vis['cov_complete']=vis[['age_base','educ_yrs','female','apoe4']].notna().all(axis=1).astype(int)
vis['primary']=(vis.path_complete.eq(1)&vis.cov_complete.eq(1)).astype(int)
b=vis[vis.primary.eq(1)].drop_duplicates('NACCID');means={'age_base':float(b.age_base.mean()),'educ_yrs':float(b.educ_yrs.mean())}
vis['age_c']=vis.age_base-means['age_base'];vis['educ_c']=vis.educ_yrs-means['educ_yrs']
flow={'raw_rows':len(d),'raw_participants':d.NACCID.nunique(),'with_neuropath':len(npm),'old_complete_assessment':int(npm.old_complete.sum()),'corrected_complete_assessment':int(npm.n_path_assessed.eq(4).sum()),'expanded_participants':vis.NACCID.nunique(),'expanded_visits':len(vis),'path_complete_longitudinal_participants':vis[vis.path_complete.eq(1)].NACCID.nunique(),'path_complete_longitudinal_visits':len(vis[vis.path_complete.eq(1)]),'primary_participants':len(b),'primary_visits':len(vis[vis.primary.eq(1)]),'duplicate_subject_days_removed':int(dups),'old_complete_tdp_false_negative':int((npm.old_complete&npm.old_tdp.eq(0)&npm.p_TDP.isna()).sum()),'abc_fallback_available':int(npm.reported_adnc.isna().sum()),'abc_fallback_value_changed':int((npm.reported_adnc.isna()&npm.abc_old.notna()&npm.abc_corrected.notna()&npm.abc_old.ne(npm.abc_corrected)).sum()),'abc_fallback_threshold_changed':int((npm.reported_adnc.isna()&npm.abc_old.notna()&npm.abc_corrected.notna()&npm.abc_old.ge(2).ne(npm.abc_corrected.ge(2))).sum()),'old_complete_regional_tdp_incomplete':int((npm.old_complete&npm.tdp_n_assessed.lt(4)).sum()),'center_count':b.NACCADC.nunique(),'centering':means,'group_counts':{str(int(k)):int(v) for k,v in b.copath_count.value_counts().sort_index().items()}}
flow['raw_sha256']=hashlib.file_digest(RAW.open('rb'),'sha256').hexdigest()
# provenance of the raw file (not derivable from the data)
flow.update(data_freeze='December 2025',data_distributed='March 2026')
vis.to_csv(OUT/'data/analytical_dataset.csv',index=False);npm.to_csv(OUT/'data/neuropath_master.csv',index=False)
json.dump(flow,(OUT/'data/cohort_flow.json').open('w'),indent=2,default=int)
print(json.dumps(flow,indent=2,default=int),flush=True)
print('NPFORMVER counts',npm.NPFORMVER.value_counts().to_dict())
print('NPFTDTAU component positive counts',npm[npm.p_FTLDtau.eq(1)][['NACCPICK','NACCCBD','NACCPROG','NPFTDT2','NPFTDT5','NPFTDT6','NPFTDT7','NPFTDT8','NPFTDT9','NPFTDT10']].eq(1).sum().to_dict())
