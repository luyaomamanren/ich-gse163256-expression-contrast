from pathlib import Path
import argparse, json, sys
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels
import statsmodels.formula.api as smf

p=argparse.ArgumentParser(); p.add_argument('--dataset',required=True); p.add_argument('--outdir',required=True)
a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(a.dataset); d.patient_id=d.patient_id.astype(str); d.day=pd.to_numeric(d.day,errors='coerce')
d=d[d.compartment.isin(['blood','hematoma']) & d.day.notna()].copy()

def mixed(formula,data):
    model=smf.mixedlm(formula,data=data,groups=data.patient_id); last=None
    for method in ['lbfgs','powell','cg']:
        try:
            r=model.fit(reml=False,method=method,maxiter=500,disp=False)
            if np.isfinite(r.params.get('Index',np.nan)): return r,method
        except Exception as e: last=e
    raise last

rows=[]
for outcome in ['Inflamm','Hypoxia']:
    r,method=mixed(f'{outcome} ~ Index + C(compartment) + day',d); ci=r.conf_int().loc['Index']
    rows.append(dict(model='overall_mixed',outcome=outcome,compartment='all',n_obs=len(d),n_patients=d.patient_id.nunique(),method=method,beta=r.params.Index,se=r.bse.Index,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=r.pvalues.Index,converged=bool(r.converged)))
    for comp in ['blood','hematoma']:
        x=d[d.compartment==comp].copy(); rr,meth=mixed(f'{outcome} ~ Index + day',x); ci=rr.conf_int().loc['Index']
        rows.append(dict(model='within_compartment_mixed',outcome=outcome,compartment=comp,n_obs=len(x),n_patients=x.patient_id.nunique(),method=meth,beta=rr.params.Index,se=rr.bse.Index,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=rr.pvalues.Index,converged=bool(rr.converged)))
    ols=smf.ols(f'{outcome} ~ Index + C(compartment) + day',d).fit(cov_type='cluster',cov_kwds={'groups':d.patient_id}); ci=ols.conf_int().loc['Index']
    rows.append(dict(model='overall_ols_cluster_robust',outcome=outcome,compartment='all',n_obs=len(d),n_patients=d.patient_id.nunique(),method='OLS clustered by patient',beta=ols.params.Index,se=ols.bse.Index,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=ols.pvalues.Index,converged=True))
pd.DataFrame(rows).to_csv(out/'gse163256_reanalysis_models.csv',index=False)

cor=[]
for outcome in ['Inflamm','Hypoxia']:
    for comp,x in d.groupby('compartment'):
        rho,pv=stats.spearmanr(x.Index,x[outcome]); cor.append(dict(outcome=outcome,compartment=comp,n_obs=len(x),n_patients=x.patient_id.nunique(),spearman_rho=rho,p=pv))
pd.DataFrame(cor).to_csv(out/'gse163256_unadjusted_correlations.csv',index=False)
d.groupby(['patient_id','compartment']).size().rename('n_observations').reset_index().to_csv(out/'gse163256_patient_compartment_manifest.csv',index=False)
(out/'gse163256_reanalysis_summary.json').write_text(json.dumps({'python':sys.version,'statsmodels':statsmodels.__version__,'n_observations':len(d),'n_patients':int(d.patient_id.nunique()),'n_blood':int((d.compartment=='blood').sum()),'n_hematoma':int((d.compartment=='hematoma').sum()),'days':sorted(map(int,d.day.unique()))},indent=2),encoding='utf-8')

