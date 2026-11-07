"""Continue the existing project: cache missing analyses and preserve all previous results."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import argparse, hashlib, importlib.util, json, sys, traceback, warnings, shutil, time
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import cohen_kappa_score
from methods import *
from factor_correction import fit

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
OUT=HERE/'results'; CACHE=HERE/'cache'; WORK=HERE/'work'
for d in [OUT,CACHE,WORK]:d.mkdir(exist_ok=True)
P=argparse.ArgumentParser();P.add_argument('--bootstrap',type=int,default=1000);P.add_argument('--poly-permutations',type=int,default=200);P.add_argument('--skip-article',action='store_true');args=P.parse_args()
assert args.bootstrap>=1000 and args.poly_permutations>=200
SEED=20261107; checks=[]; reuse=[]; warnings_log=[]; metrics=[]; loads=[]; phis=[]; residuals=[]; omega_rows=[]; poly_notes=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def writejson(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)),encoding='utf-8')
def csv(name,data):
    frame=data if isinstance(data,pd.DataFrame) else pd.DataFrame(data)
    frame.to_csv(OUT/name,index=False,float_format='%.12g',encoding='utf-8-sig');return frame
def check(name,ok,detail=''):
    checks.append(dict(check=name,status='PASS' if bool(ok) else 'FAIL',detail=str(detail)))
    if not ok:raise AssertionError(name+': '+str(detail))
def read(path):return pd.read_csv(ROOT/path,keep_default_na=False)
def stage(n,text):print(f'STAGE {n}: {text}',flush=True)
def protected():
    paths=[]
    for folder in ROOT.iterdir():
        if folder.is_dir() and folder.name[:2].isdigit() and int(folder.name[:2])<=15:
            paths += [p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    return {p.relative_to(ROOT).as_posix():sha(p) for p in paths}
before=protected()
fingerprint=hashlib.sha256(json.dumps({'inputs':before,'methods':sha(HERE/'methods.py'),'seed':SEED,'bootstrap':args.bootstrap,'poly_permutations':args.poly_permutations},sort_keys=True).encode()).hexdigest()
def cached(name,fn):
    p=CACHE/(name+'.npz');meta=CACHE/(name+'.json')
    if p.exists() and meta.exists() and json.loads(meta.read_text())['fingerprint']==fingerprint:
        reuse.append(dict(task=name,action='cached_missing_analysis'));return dict(np.load(p,allow_pickle=False))
    result=fn();np.savez_compressed(p,**result);writejson(meta,dict(fingerprint=fingerprint));return result
def export_model(label,kind,m,r,origin):
    row=dict(sample=label,correlation=kind,source=origin,**fit_metrics(r,m));metrics.append(row)
    l=m['L'];k=l.shape[1];h=np.diag(l@m['Phi']@l.T)
    for i,item in enumerate(ITEMS):
        intended=next(d for d,ix in DIMS.items() if i in ix)
        order=np.argsort(np.abs(l[i]))[::-1]
        loads.append(dict(sample=label,correlation=kind,factors=k,item=item,theoretical_dimension=intended,**{f'loading_F{j+1}':l[i,j] for j in range(k)},communality=h[i],uniqueness=1-h[i],largest_loading_factor=f'F{order[0]+1}',second_abs_loading=float(abs(l[i,order[1]])) if k>1 else 0))
    for i in range(k):
        for j in range(k):phis.append(dict(sample=label,correlation=kind,factors=k,factor1=f'F{i+1}',factor2=f'F{j+1}',correlation_value=m['Phi'][i,j]))
    e=r-l@m['Phi']@l.T
    for i in range(14):
        for j in range(i+1,14):residuals.append(dict(sample=label,correlation=kind,factors=k,item1=ITEMS[i],item2=ITEMS[j],residual=e[i,j]))
    warnings_log.extend(dict(sample=label,correlation=kind,factors=k,message=w) for w in m.get('warnings',[]))
    return row
def summarize(x,scores):
    sd=scores.std(ddof=1);delta=stats.t.ppf(.975,len(x)-1)*sd/np.sqrt(len(x))
    return dict(n=len(x),BDSF_mean=scores.mean(),BDSF_sd=sd,BDSF_ci_low=scores.mean()-delta,BDSF_ci_high=scores.mean()+delta,alpha=alpha(x),**diagnostics(np.corrcoef(x,rowvar=False),len(x)))
def parallel(x,label):
    def calc():
        rng=np.random.default_rng(SEED);null=np.array([np.linalg.eigvalsh(np.corrcoef(rng.normal(size=x.shape),rowvar=False))[::-1] for _ in range(500)])
        return dict(observed=np.linalg.eigvalsh(np.corrcoef(x,rowvar=False))[::-1],null95=np.quantile(null,.95,axis=0))
    a=cached('pa_'+label,calc);csv('PA_'+label+'.csv',[dict(component=i+1,observed=a['observed'][i],null95=a['null95'][i],retain=a['observed'][i]>a['null95'][i]) for i in range(14)])
    return int((a['observed']>a['null95']).sum())
def getpoly(x,label):
    def calc():
        raw,r,detail,note=polychoric(x)
        return dict(raw=raw,r=r,detail=json.dumps(detail),note=json.dumps(note))
    d=cached('poly_'+label,calc);note=json.loads(str(d['note']));poly_notes.append(dict(sample=label,**note))
    csv('POLY_'+label+'_pairs.csv',json.loads(str(d['detail'])))
    for suffix,mat in [('raw',d['raw']),('used',d['r'])]:csv('POLY_'+label+'_'+suffix+'.csv',pd.DataFrame(mat,columns=ITEMS).assign(item=ITEMS))
    return d['r']
def omega_boot(x,lang,scale,point):
    def calc():
        rng=np.random.default_rng(SEED+(0 if lang=='KZ' else 10000)+list({'ALL14':None,**DIMS}).index(scale));values=[];failures=[]
        for i in range(args.bootstrap):
            sample=x[rng.integers(0,len(x),len(x))]
            try:
                if np.any(sample.std(axis=0)==0):raise ValueError('constant item')
                model=fit(np.corrcoef(sample,rowvar=False),1)
                if any('converge' in w.lower() for w in model['warnings']):raise ValueError('nonconvergence')
                v=omega(model,sample.std(axis=0,ddof=1))
                if not np.isfinite(v):raise ValueError('inadmissible model')
                values.append(v)
            except (ValueError,np.linalg.LinAlgError) as e:failures.append(str(e))
        return dict(values=np.array(values),failures=json.dumps(failures))
    d=cached('omega_boot_'+lang+'_'+scale,calc);v=d['values'];valid=len(v)
    lo,hi=np.quantile(v,[.025,.975]) if valid>=.9*args.bootstrap else (np.nan,np.nan)
    return dict(ci95_low=lo,ci95_high=hi,bootstrap_requested=args.bootstrap,bootstrap_valid=valid,bootstrap_failed=args.bootstrap-valid,ci_status='OK' if valid>=.9*args.bootstrap else 'UNSTABLE_NOT_REPORTED')

try:
    writejson(OUT/'VALIDATION_REPORT.json',dict(status='RUNNING',fingerprint=fingerprint))
    stage(1,'Verify existing CSVs, raw-source mapping, immutable scoring and inclusion')
    numerical_tests=method_tests()
    df=read('05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv');raw=read('03_harmonization/03_01_BDSF_Kazakhstan_Harmonized_Dataset.csv')
    x=df[[q+'_score' for q in ITEMS]].to_numpy(float);scores=df.BDSF_score.to_numpy(float)
    check('159 adults and 180 original rows',len(df)==159 and len(raw)==180)
    check('Exact adult partition',set(df.respondent_id)==set(raw.loc[raw.age_group_code.ne('LT18'),'respondent_id']))
    check('No duplicated IDs',df.respondent_id.nunique()==159)
    check('All 14 values complete and correctly coded',np.isfinite(x).all() and np.isin(x,[0,25,50,75,100]).all())
    check('Language counts',df.source_language.value_counts().to_dict()=={'KZ':81,'RU':78})
    weights=np.array([1/(4*len(ix)) for i in range(14) for ix in DIMS.values() if i in ix])
    check('Original M1 unchanged',np.allclose(x@weights,scores,atol=1e-9))
    for d,ix in DIMS.items():check(d+' items unchanged',np.allclose(x[:,ix].mean(1),df[d+'_score']))
    # Verify raw answers against stored harmonized/scored data; no statistics pipeline rerun.
    for lang,word in [('KZ','Kazakh'),('RU','Russian')]:
        src=next((ROOT/'00_sources').glob('*_'+word+'_Responses.csv'));r=pd.read_csv(src,dtype=str,keep_default_na=False)
        h=raw[raw.source_language==lang].sort_values('source_row')
        check(lang+' original response texts retained',np.array_equal(r.iloc[:,4:].to_numpy(),h[[q+'_raw' for q in ITEMS]].to_numpy()))
    r=np.corrcoef(x,rowvar=False);summary=summarize(x,scores)
    oldvalid=read('08_construct_validity/08_04_BDSF_Kazakhstan_Validity_Indices.csv').set_index('metric').value.astype(float)
    oldrel=read('07_reliability_analysis/07_02_BDSF_Kazakhstan_Scale_Reliability.csv').set_index('scale')
    check('Baseline mean and alpha',abs(summary['BDSF_mean']-60.3347746331237)<1e-9 and abs(summary['alpha']-float(oldrel.loc['14_items','cronbach_alpha']))<1e-9)
    for key in ['KMO','Bartlett_chi2']:check('Baseline '+key,abs(summary[key]-oldvalid[key])<1e-7)
    csv('STAGE01_Data_Audit.csv',checks)
    stage(2,'Compute missing one-factor MINRES; reuse pooled Pearson parallel analysis')
    one=fit(r,1);export_model('N159','Pearson',one,r,'new')
    pa=read('08_construct_validity/08_05_BDSF_Kazakhstan_Parallel_Analysis.csv')
    check('Stored parallel observed eigenvalues',np.allclose(pa.observed_eigenvalue,np.linalg.eigvalsh(r)[::-1]))
    reuse.append(dict(task='N159 Pearson parallel 500 simulations',action='reused_verified_existing'))
    summary['parallel_retained']=int(oldvalid['Parallel_retained_factors']);summary['omega_raw_sum']=omega(one,x.std(0,ddof=1));summary['omega_M1_weighted']=omega(one,x.std(0,ddof=1),weights)
    stage(3,'Reuse existing four-factor solution and calculate new residual diagnostics')
    oldloads=read('08_construct_validity/08_03_BDSF_Kazakhstan_Factor_Loadings.csv')
    four=dict(L=oldloads[[f'loading_{d}' for d in DIMS]].to_numpy(float),Phi=read('08_construct_validity/08_06_BDSF_Kazakhstan_Factor_Correlations.csv')[list(DIMS)].to_numpy(float),warnings=[])
    check('Existing oblique communalities',np.allclose(np.diag(four['L']@four['Phi']@four['L'].T),oldloads.communality.astype(float)))
    fixed=align(fit(r,4),four['L'])
    check('Corrected model preserves original loadings',np.allclose(fixed['L'],four['L'],atol=1e-7))
    legacy_audit=[dict(sample='N159',old_rmsr=fit_metrics(r,four)['offdiag_rmsr'],corrected_rmsr=fit_metrics(r,fixed)['offdiag_rmsr'],max_communality_change=float(np.max(np.abs(np.diag(fixed['L']@fixed['Phi']@fixed['L'].T)-oldloads.communality.to_numpy(float)))),reason='factor_analyzer sorts loadings and structure but not phi; corrected via structure=L@Phi')]
    four=dict(L=four['L'],Phi=fixed['Phi'],warnings=[])
    export_model('N159','Pearson',four,r,'existing_loadings_Phi_order_corrected');reuse.append(dict(task='N159 four-factor EFA',action='loadings preserved; Phi order bug corrected in new outputs only'))
    rp=getpoly(x,'N159');polyone=fit(rp,1);polyfour=align(fit(rp,4),four['L'])
    export_model('N159','Polychoric',polyone,rp,'new');export_model('N159','Polychoric',polyfour,rp,'new')
    stage(4,'Compare models and run ordinal marginal-preserving permutation parallel analysis')
    def poly_pa():
        rng=np.random.default_rng(SEED);ev=[];adjustments=0
        for i in range(args.poly_permutations):
            null=np.column_stack([rng.permutation(x[:,j]) for j in range(14)])
            _,rn,_,note=polychoric(null);ev.append(np.linalg.eigvalsh(rn)[::-1]);adjustments+=note['smoothed']
            if (i+1)%25==0:print(f'  Polychoric permutation {i+1}/{args.poly_permutations}',flush=True)
        return dict(null=np.array(ev),smoothed=adjustments)
    pp=cached('poly_pa_N159',poly_pa);obs=np.linalg.eigvalsh(rp)[::-1];threshold=np.quantile(pp['null'],.95,axis=0)
    csv('EXP01_Polychoric_Parallel.csv',[dict(component=i+1,observed=obs[i],null95=threshold[i],retain=bool(obs[i]>threshold[i])) for i in range(14)])
    polyretained=int((obs>threshold).sum())
    stage(5,'Verify quality flags on all 180 records, including overlap')
    levels=[q+'_level' for q in ITEMS];dup=raw.duplicated(levels,keep=False);ts=raw.timestamp_iso.duplicated(keep=False);straight=raw[levels].nunique(axis=1).eq(1)
    flaglist=[';'.join((["STRAIGHTLINE_ALL_14"] if straight.iloc[i] else [])+(['DUPLICATE_RESPONSE_PATTERN'] if dup.iloc[i] else [])+(['DUPLICATE_TIMESTAMP'] if ts.iloc[i] else [])) for i in range(len(raw))]
    flagmap=dict(zip(raw.respondent_id,flaglist));verified=df.respondent_id.map(flagmap)
    check('All stored flags reproduce on original 180',verified.equals(df.quality_flags))
    flagged=verified.ne('').to_numpy();check('Exactly 20 flagged adults',int(flagged.sum())==20)
    csv('EXP02_Flag_Overlap.csv',pd.DataFrame({'flag_combination':verified.replace('','NONE')}).value_counts().reset_index(name='n'))
    csv('EXP02_Flag_Counts.csv',[dict(flag=f,n=int(verified.str.contains(f).sum()),scope='adults; algorithm evaluated on 180') for f in ['STRAIGHTLINE_ALL_14','DUPLICATE_RESPONSE_PATTERN','DUPLICATE_TIMESTAMP']])
    stage(6,'Reuse and extend N159 summary')
    summary['sample']='N159';summary['role']='PRIMARY'
    stage(7,'Compute N139 sensitivity only')
    x139=x[~flagged];r139=np.corrcoef(x139,rowvar=False);s139=summarize(x139,scores[~flagged]);check('Sensitivity n139',len(x139)==139)
    one139=fit(r139,1);four139=align(fit(r139,4),four['L'])
    export_model('N139','Pearson',one139,r139,'new_sensitivity');export_model('N139','Pearson',four139,r139,'new_sensitivity')
    s139.update(sample='N139',role='SENSITIVITY_ONLY',parallel_retained=parallel(x139,'N139'),omega_raw_sum=omega(one139,x139.std(0,ddof=1)),omega_M1_weighted=omega(one139,x139.std(0,ddof=1),weights))
    stage(8,'Compare nested samples without an independent-samples significance test')
    csv('EXP02_Quality_Sensitivity.csv',[summary,s139])
    csv('EXP02_Absolute_Relative_Changes.csv',[dict(metric=k,N159=summary[k],N139=s139[k],absolute_change=s139[k]-summary[k],relative_change_percent=100*(s139[k]-summary[k])/summary[k] if summary[k]!=0 else np.nan) for k in ['BDSF_mean','BDSF_sd','alpha','omega_raw_sum','omega_M1_weighted','KMO','parallel_retained']])
    csv('EXP02_Loading_Changes.csv',[dict(item=ITEMS[i],one_factor_N159=one['L'][i,0],one_factor_N139=one139['L'][i,0],**{f'F{j+1}_difference':four139['L'][i,j]-four['L'][i,j] for j in range(4)}) for i in range(14)])
    stage(9,'M1 original, M2 equal items, M3 one-factor regression scores')
    z=(x-x.mean(0))/x.std(0,ddof=0);coef=np.linalg.solve(r,one['L'][:,0]);m3=z@coef
    check('M3 direction anchored to overall score',stats.spearmanr(m3,scores).statistic>0)
    methods={'M1':scores,'M2':x.mean(1),'M3':m3};ranks={k:stats.rankdata(-v,method='average') for k,v in methods.items()}
    # Output-only anonymous row labels; never export respondent IDs, timestamps or demographics.
    comp=pd.DataFrame({'case':[f'P{i+1:03}' for i in range(159)],**methods,**{k+'_rank':v for k,v in ranks.items()}})
    csv('EXP03_Scoring_Comparison.csv',comp)
    csv('EXP03_M3_Coefficients.csv',[dict(item=q,regression_coefficient=coef[i],training_mean=x[:,i].mean(),training_sd_population=x[:,i].std(ddof=0)) for i,q in enumerate(ITEMS)])
    stage(10,'Rank stability including ties')
    rankstats=[]
    for a,b in [('M1','M2'),('M1','M3'),('M2','M3')]:
        delta=ranks[b]-ranks[a]
        rankstats.append(dict(method1=a,method2=b,spearman=stats.spearmanr(methods[a],methods[b]).statistic,kendall_tau_b=stats.kendalltau(methods[a],methods[b]).statistic,mean_abs_rank_difference=np.abs(delta).mean(),max_abs_rank_difference=np.abs(delta).max(),n_rank_changed=int((delta!=0).sum())))
    csv('EXP03_Rank_Stability.csv',rankstats)
    changed=comp.assign(M2_minus_M1_rank=ranks['M2']-ranks['M1'],M3_minus_M1_rank=ranks['M3']-ranks['M1'])
    csv('EXP03_Largest_Rank_Changes.csv',changed.assign(max_abs_change=lambda t:t[['M2_minus_M1_rank','M3_minus_M1_rank']].abs().max(axis=1)).sort_values('max_abs_change',ascending=False).head(15))
    stage(11,'Category transitions M1 to M2; no unjustified M3 score transformation')
    def cat(v):return np.where(np.round(v,2)<=33.32,'LOW',np.where(np.round(v,2)<=66.66,'MEDIUM','HIGH'))
    c1=cat(scores);c2=cat(methods['M2']);check('Original category rule',np.array_equal(c1,df.BDSF_level.to_numpy()))
    cats=['LOW','MEDIUM','HIGH'];transition=[dict(from_method='M1',to_method='M2',from_category=a,to_category=b,n=int(((c1==a)&(c2==b)).sum()),status='COMPUTED') for a in cats for b in cats]
    transition.append(dict(from_method='M1',to_method='M3',from_category='NOT_APPLICABLE',to_category='NOT_APPLICABLE',n=None,status='M3 is an uncalibrated latent score; no 0-100 thresholds'))
    csv('EXP03_Category_Transitions.csv',transition)
    kappa=cohen_kappa_score([cats.index(v) for v in c1],[cats.index(v) for v in c2],labels=[0,1,2],weights='quadratic')
    category_summary=dict(n_changed=int((c1!=c2).sum()),percent_changed=100*float((c1!=c2).mean()),quadratic_weighted_kappa=float(kappa))
    csv('EXP03_Category_Summary.csv',[dict(method=m,category=c,n=int((arr==c).sum()),**category_summary) for m,arr in [('M1',c1),('M2',c2)] for c in cats])
    stage(12,'Weighting sensitivity documented; M1 remains primary')
    stage(13,'Verify and reuse previous language estimates and hashes')
    oldlang=read('15_language_psychometrics/15_01_Reliability.csv');oldit=read('15_language_psychometrics/15_02_Item_Total.csv');oldld=read('15_language_psychometrics/15_03_Factor_Loadings.csv');olddist=read('15_language_psychometrics/15_04_Response_Distributions.csv')
    val=json.loads((ROOT/'15_language_psychometrics/15_05_Validation.json').read_text(encoding='utf-8'))
    check('Language inputs SHA256',val['source_sha256']==sha(ROOT/'05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv'))
    for name,hashvalue in val['outputs'].items():check('Stored language hash '+name,sha(ROOT/'15_language_psychometrics'/name)==hashvalue)
    language_models={};language_poly={};langrows=[]
    for lang in ['KZ','RU']:
        a=x[df.source_language.eq(lang).to_numpy()];rr=np.corrcoef(a,rowvar=False)
        stored=oldld[oldld.language==lang][[f'loading_{d}' for d in DIMS]].to_numpy(float)
        verified4=align(fit(rr,4),stored)
        check(lang+' previous loadings validation only',np.allclose(verified4['L'],stored,atol=1e-7))
        oldh=oldld[oldld.language==lang].communality.to_numpy(float)
        legacy_audit.append(dict(sample=lang,old_rmsr=None,corrected_rmsr=fit_metrics(rr,verified4)['offdiag_rmsr'],max_communality_change=float(np.max(np.abs(np.diag(verified4['L']@verified4['Phi']@verified4['L'].T)-oldh))),reason='preserved loadings; recovered consistent Phi from library structure matrix'))
        language_models[lang]=verified4
        for scale,ix in {'ALL14':list(range(14)),**DIMS}.items():
            b=a[:,ix];previous=oldlang[(oldlang.language==lang)&(oldlang.scale==scale)].iloc[0]
            check(lang+' alpha '+scale,abs(alpha(b)-float(previous.alpha))<1e-8)
            if scale!='ALL14':
                for j,idx in enumerate(ix):
                    cit=np.corrcoef(b[:,j],np.delete(b,j,axis=1).sum(1))[0,1]
                    st=oldit[(oldit.language==lang)&(oldit.item==ITEMS[idx])].iloc[0]
                    check(lang+' ITC '+ITEMS[idx],abs(cit-float(st.corrected_item_total_within_dimension))<1e-8)
        for i,q in enumerate(ITEMS):
            st=olddist[(olddist.language==lang)&(olddist.item==q)]
            check(lang+' distribution '+q,np.array_equal(st['count'].to_numpy(int),[(a[:,i]==v).sum() for v in [0,25,50,75,100]]))
        reuse.append(dict(task=lang+' alpha, item-total, four-factor, response frequencies',action='verified_and_reused'))
    stored_cong=read('15_language_psychometrics/15_03_Factor_Congruence.csv').tucker_congruence.to_numpy(float)
    lk=language_models['KZ']['L'];lr=language_models['RU']['L']
    check('Existing language Tucker coefficients',np.allclose(stored_cong,np.sum(lk*lr,0)/(np.linalg.norm(lk,axis=0)*np.linalg.norm(lr,axis=0)),atol=1e-8))
    stage(14,'Add language omega bootstrap intervals and polychoric sensitivity')
    for lang in ['KZ','RU']:
        a=x[df.source_language.eq(lang).to_numpy()];rr=np.corrcoef(a,rowvar=False);rp_lang=getpoly(a,lang)
        lp1=fit(rp_lang,1);lp4=align(fit(rp_lang,4),language_models[lang]['L']);language_poly[lang]=lp4
        export_model(lang,'Pearson',language_models[lang],rr,'verified_existing')
        export_model(lang,'Polychoric',lp1,rp_lang,'new');export_model(lang,'Polychoric',lp4,rp_lang,'new')
        for scale,ix in {'ALL14':list(range(14)),**DIMS}.items():
            print('  Omega bootstrap '+lang+' '+scale,flush=True)
            b=a[:,ix];model=fit(np.corrcoef(b,rowvar=False),1);w=omega(model,b.std(0,ddof=1));ci=omega_boot(b,lang,scale,w)
            poly_model=fit(rp_lang[np.ix_(ix,ix)],1);wp=omega(poly_model)
            previous=oldlang[(oldlang.language==lang)&(oldlang.scale==scale)].iloc[0]
            langrows.append(dict(language=lang,scale=scale,n=len(a),alpha_reused=float(previous.alpha),alpha_ci_low_reused=float(previous.ci95_low),alpha_ci_high_reused=float(previous.ci95_high),omega_raw_sum=w,**ci,omega_polychoric_latent=wp,polychoric_omega_ci='NOT_ESTIMATED',polychoric_onefactor_rmsr=fit_metrics(rp_lang[np.ix_(ix,ix)],poly_model)['offdiag_rmsr']))
    csv('EXP04_Language_Comparison.csv',langrows)
    stage(15,'Language structural stability and limits')
    cong=[]
    for label,a,b in [('KZ Pearson vs polychoric',language_models['KZ']['L'],language_poly['KZ']['L']),('RU Pearson vs polychoric',language_models['RU']['L'],language_poly['RU']['L']),('KZ vs RU polychoric',language_poly['KZ']['L'],align(language_poly['RU'],language_poly['KZ']['L'])['L'])]:
        for j in range(4):cong.append(dict(comparison=label,factor=f'F{j+1}',tucker_congruence=float(a[:,j]@b[:,j]/(np.linalg.norm(a[:,j])*np.linalg.norm(b[:,j])))))
    csv('EXP04_Factor_Stability.csv',cong)
    stage(16,'Integrate evidence, document columns, prepare revised article copy')
    csv('EXP01_Factor_Comparison.csv',metrics);csv('EXP01_Factor_Loadings.csv',loads);csv('EXP01_Factor_Correlations.csv',phis);csv('EXP01_Residual_Correlations.csv',residuals);csv('POLY_Matrix_Diagnostics.csv',poly_notes);csv('REUSE_LEDGER.csv',reuse)
    csv('EXP01_Legacy_Phi_Audit.csv',legacy_audit)
    check('All previous source and result files unchanged',before==protected())
    csv('STAGE_Completion.csv',[dict(stage=i,status='COMPLETE',scope='existing results verified and reused; missing analyses added') for i in range(1,17)])
    evidence=dict(summary159=summary,summary139=s139,poly_parallel_retained=polyretained,poly_parallel_permutations=args.poly_permutations,poly_null_matrices_smoothed=int(pp['smoothed']),rank=rankstats,categories=category_summary,language=langrows,model_comparison=metrics,poly_diagnostics=poly_notes,congruence=cong,bootstrap=args.bootstrap,legacy_phi_audit=legacy_audit)
    # Nonfinite omega/CI values are explicitly null, never presented as valid estimates.
    def finite(obj):
        if isinstance(obj,dict):return {k:finite(v) for k,v in obj.items()}
        if isinstance(obj,list):return [finite(v) for v in obj]
        if isinstance(obj,(float,np.floating)):return float(obj) if np.isfinite(obj) else None
        return obj
    writejson(OUT/'EVIDENCE.json',finite(evidence))
    from reporting import make_report, update_article
    make_report(HERE,ROOT,evidence,checks)
    if not args.skip_article:update_article(HERE,ROOT,evidence)
    check('Protected files unchanged after reporting',before==protected())
    csv('VALIDATION_CHECKS.csv',checks)
    writejson(OUT/'VALIDATION_REPORT.json',dict(status='PASS',stages_completed=16,checks=checks,numerical_tests=numerical_tests,input_fingerprint=fingerprint,protected_file_count=len(before),protected_hashes=before,model_warnings=warnings_log,bootstrap=args.bootstrap,poly_permutations=args.poly_permutations,article_created=not args.skip_article,output_sha256={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name not in ['VALIDATION_REPORT.json','INDEPENDENT_VERIFICATION.json']}))
    print('PASS: 16 stages complete; original results preserved',flush=True)
except Exception as e:
    writejson(OUT/'VALIDATION_REPORT.json',dict(status='FAIL',error=str(e),checks=checks,traceback=traceback.format_exc()))
    traceback.print_exc();sys.exit(1)
