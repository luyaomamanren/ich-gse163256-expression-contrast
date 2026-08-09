from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

p=argparse.ArgumentParser()
p.add_argument('--expr',required=True); p.add_argument('--gmt',required=True); p.add_argument('--outdir',required=True)
a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

sets={}
for line in Path(a.gmt).read_text(encoding='utf-8').splitlines():
    f=line.rstrip().split('\t')
    if len(f)>=3: sets[f[0]]=set(f[2:])
names={'gly':'HALLMARK_GLYCOLYSIS','ox':'HALLMARK_OXIDATIVE_PHOSPHORYLATION','hyp':'HALLMARK_HYPOXIA','inf':'HALLMARK_INFLAMMATORY_RESPONSE'}
missing=[v for v in names.values() if v not in sets]
if missing: raise ValueError(f'Missing Hallmark sets: {missing}')
gs={k:sets[v] for k,v in names.items()}

expr=pd.read_csv(a.expr,index_col=0)
all_cols=set(expr.columns); base=[]
for c in expr.columns:
    root=c.rsplit('.',1)[0]
    base.append(root if c.rsplit('.',1)[-1].isdigit() and root in all_cols else c)
expr.columns=base
expr=expr.T.groupby(level=0).mean().T
available=set(expr.index.astype(str)); gs={k:v & available for k,v in gs.items()}

def score(g): return expr.loc[sorted(g)].mean(axis=0)
rows=[]
for s in expr.columns:
    f=s.split('_')
    if len(f)==3 and f[1] in {'blood','hematoma'}:
        try: rows.append((s,str(f[0]),f[1],int(float(f[2]))))
        except ValueError: pass
meta=pd.DataFrame(rows,columns=['sample','patient_id','compartment','day']).set_index('sample')

def dataset(use):
    d=meta.copy()
    for key,label in [('gly','Gly'),('ox','Ox'),('hyp','Hypoxia'),('inf','Inflamm')]: d[label]=score(use[key]).reindex(d.index)
    d['Index']=d.Gly-d.Ox
    return d.reset_index().dropna()

def fit_cluster(d,outcome,pred='Index'):
    m=smf.ols(f'{outcome} ~ {pred} + C(compartment) + day',d).fit(cov_type='cluster',cov_kwds={'groups':d.patient_id})
    ci=m.conf_int().loc[pred]
    return dict(beta=m.params[pred],se=m.bse[pred],ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=m.pvalues[pred],n=len(d),patients=d.patient_id.nunique())

clean={k:gs[k]-set().union(*(gs[j] for j in gs if j!=k)) for k in gs}
overlap=[]
keys=list(gs)
for i,x in enumerate(keys):
    for y in keys[i+1:]:
        inter=sorted(gs[x]&gs[y]); overlap.append({'set_a':x,'set_b':y,'n_overlap':len(inter),'genes':';'.join(inter)})
pd.DataFrame(overlap).to_csv(out/'signature_overlap_2025_1.csv',index=False)

res=[]
for label,use in [('original',gs),('all_overlap_pruned',clean)]:
    d=dataset(use)
    if label=='original': d.to_csv(out/'gse163256_locked_analysis_dataset.csv',index=False)
    for outcome in ['Inflamm','Hypoxia']:
        z=fit_cluster(d,outcome); z.update(signature=label,outcome=outcome,predictor='Gly-minus-Ox'); res.append(z)
        for pred in ['Gly','Ox']:
            z=fit_cluster(d,outcome,pred)
            if pred=='Ox':
                z['beta']=-z['beta']; lo,hi=-z['ci_high'],-z['ci_low']; z['ci_low'],z['ci_high']=lo,hi
            z.update(signature=label,outcome=outcome,predictor='Gly-only' if pred=='Gly' else 'negative-Ox-only'); res.append(z)
pd.DataFrame(res).to_csv(out/'signature_sensitivity_models.csv',index=False)

loo=[]
for predset in ['gly','ox']:
    for gene in sorted(clean[predset]):
        use={k:set(v) for k,v in clean.items()}; use[predset].remove(gene); d=dataset(use)
        for outcome in ['Inflamm','Hypoxia']:
            loo.append({'removed_from':predset,'removed_gene':gene,'outcome':outcome,**fit_cluster(d,outcome)})
pd.DataFrame(loo).to_csv(out/'signature_leave_one_gene_out.csv',index=False)

d=dataset(clean); cv=[]
for outcome in ['Inflamm','Hypoxia']:
    actual=[]; pred=[]
    for pid in sorted(d.patient_id.unique()):
        tr=d[d.patient_id!=pid]; te=d[d.patient_id==pid]
        m=smf.ols(f'{outcome} ~ Index + C(compartment) + day',tr).fit()
        actual.extend(te[outcome]); pred.extend(m.predict(te))
    actual=np.asarray(actual); pred=np.asarray(pred); sp=stats.spearmanr(actual,pred)
    cv.append({'outcome':outcome,'n':len(actual),'patients':d.patient_id.nunique(),'rmse':float(np.sqrt(np.mean((actual-pred)**2))),'r2_predictive':float(1-np.sum((actual-pred)**2)/np.sum((actual-actual.mean())**2)),'spearman_rho':float(sp.statistic),'spearman_p':float(sp.pvalue)})
pd.DataFrame(cv).to_csv(out/'signature_leave_one_patient_out_cv.csv',index=False)
(out/'signature_sensitivity_summary.json').write_text(json.dumps({'msigdb':'2025.1.Hs Hallmark symbols','gene_counts_original':{k:len(v) for k,v in gs.items()},'gene_counts_pruned':{k:len(v) for k,v in clean.items()},'overlap_pairs':overlap,'n_samples':len(d),'n_patients':int(d.patient_id.nunique())},indent=2),encoding='utf-8')

