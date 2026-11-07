"""Independent result consistency checks. Does not rerun fitted models or bootstrap."""
from pathlib import Path
import json,hashlib
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy import stats
from factor_analyzer import FactorAnalyzer
from factor_correction import fit
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;OUT=HERE/'results'
def load(name):return pd.read_csv(OUT/name)
checks=[]
def check(name,condition):
    checks.append(dict(check=name,passed=bool(condition)))
    if not condition:raise AssertionError(name)

v=json.loads((OUT/'VALIDATION_REPORT.json').read_text(encoding='utf-8'))
check('Pipeline status',v['status']=='PASS')
for rel,h in v['protected_hashes'].items():check('Immutable '+rel,hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==h)
for rel,h in v['output_sha256'].items():check('Output hash '+rel,hashlib.sha256((OUT/rel).read_bytes()).hexdigest()==h)
df=pd.read_csv(ROOT/'05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv',keep_default_na=False)
items=[f'Q{i:02d}_score' for i in range(1,15)];x=df[items].to_numpy(float)
levels={'NEVER':0,'RARELY':25,'SOMETIMES':50,'OFTEN':75,'ALWAYS':100}
for i in range(1,15):check('Numeric coding Q'+str(i),np.array_equal(df[f'Q{i:02d}_level'].map(levels),df[f'Q{i:02d}_score']))
sc=load('EXP03_Scoring_Comparison.csv')
check('No respondent private fields',set(sc.columns)=={'case','M1','M2','M3','M1_rank','M2_rank','M3_rank'})
check('M1 values retained',np.allclose(sc.M1,df.BDSF_score))
check('M2 formula',np.allclose(sc.M2,x.mean(1)))
coef=load('EXP03_M3_Coefficients.csv').regression_coefficient.to_numpy()
check('M3 independently reconstructed',np.allclose(sc.M3,((x-x.mean(0))/x.std(0))@coef))
for m in ['M1','M2','M3']:
    # Stored rounded CSV scores can turn tiny floating differences into ties; tolerate half a rank only if needed.
    check(m+' rank sum',abs(sc[m+'_rank'].sum()-159*160/2)<1e-7)
tr=load('EXP03_Category_Transitions.csv');t=tr[tr.to_method=='M2'].pivot(index='from_category',columns='to_category',values='n').reindex(index=['LOW','MEDIUM','HIGH'],columns=['LOW','MEDIUM','HIGH']).to_numpy()
check('Transition total',t.sum()==159)
changed=int(t.sum()-np.trace(t));cs=load('EXP03_Category_Summary.csv')
check('Category changed count',int(cs.n_changed.iloc[0])==changed)
w=(np.arange(3)[:,None]-np.arange(3)[None,:])**2
expected=np.outer(t.sum(1),t.sum(0))/t.sum();k=1-(w*t).sum()/(w*expected).sum()
check('Weighted kappa independent formula',abs(k-cs.quadratic_weighted_kappa.iloc[0])<1e-10)
metrics=load('EXP01_Factor_Comparison.csv');ld=load('EXP01_Factor_Loadings.csv');phis=load('EXP01_Factor_Correlations.csv');res=load('EXP01_Residual_Correlations.csv')
for _,r in metrics.iterrows():
    mask=(ld['sample']==r['sample'])&(ld.correlation==r.correlation)&(ld.factors==r.factors)
    a=ld[mask];l=a[[f'loading_F{i+1}' for i in range(int(r.factors))]].to_numpy()
    p=phis[(phis['sample']==r['sample'])&(phis.correlation==r.correlation)&(phis.factors==r.factors)].correlation_value.to_numpy().reshape(int(r.factors),int(r.factors))
    check('Oblique communality '+str(tuple(r[['sample','correlation','factors']])),np.allclose(np.diag(l@p@l.T),a.communality))
    rr=res[(res['sample']==r['sample'])&(res.correlation==r.correlation)&(res.factors==r.factors)].residual.to_numpy()
    check('RMSR '+str(tuple(r[['sample','correlation','factors']])),abs(np.sqrt(np.mean(rr**2))-r.offdiag_rmsr)<1e-10)
lang=load('EXP04_Language_Comparison.csv');check('10 language reliability rows',len(lang)==10)
rr=np.corrcoef(x,rowvar=False)
unrotated=FactorAnalyzer(n_factors=4,rotation=None,method='minres',is_corr_matrix=True);unrotated.fit(rr)
fixed=fit(rr,4)
check('Corrected oblique covariance equals unrotated covariance',np.allclose(fixed['L']@fixed['Phi']@fixed['L'].T,unrotated.loadings_@unrotated.loadings_.T,atol=1e-7))
check('Bootstrap counts',((lang.bootstrap_valid+lang.bootstrap_failed)==1000).all())
check('Valid omega bounds',lang.omega_raw_sum.between(0,1).all())
for p in OUT.glob('*.csv'):
    fields=set(pd.read_csv(p,nrows=0).columns)
    check('No direct participant fields '+p.name,not(fields&{'respondent_id','timestamp_iso','timestamp_raw','gender_code','age_group_code','email','name'}))
report={'status':'PASS','checks_count':len(checks),'checks':checks}
(OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(f'PASS: {len(checks)} independent checks')
