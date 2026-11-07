"""Reproduce three supplementary experiments from the two original CSVs.

All writes stay in section 17. Prior project files, methods and primary M1 are preserved.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys,importlib.util,hashlib,json,subprocess,math,itertools,argparse,platform
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(HERE));sys.path.insert(0,str(ROOT/'16_extended_experiments'))
import numpy as np
import pandas as pd
from scipy import stats
from methods import alpha,diagnostics,fit,omega,ITEMS,DIMS
SEED=20261107
EXPS=['EXP01_Ordinal_CFA','EXP02_TopK_Stability','EXP03_Quality_Sensitivity']
WORK=HERE/'work';VALID=HERE/'validation'
for p in [WORK,VALID,*[HERE/v for v in EXPS]]:p.mkdir(parents=True,exist_ok=True)
CHECKS=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    CHECKS.append(dict(check=name,status='PASS' if bool(ok) else 'FAIL',detail=str(detail)))
    if not ok:raise AssertionError(name+': '+str(detail))
def save(p,rows):
    f=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
    f.to_csv(p,index=False,float_format='%.12g',encoding='utf-8-sig');return f
def jsonsave(p,x):
    p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False,default=lambda o:o.item() if isinstance(o,np.generic) else str(o)),encoding='utf-8')
def snapshot():
    files=[]
    for d in ROOT.iterdir():
        if d.is_dir() and d.name[:2].isdigit() and int(d.name[:2])<=16:
            files.extend(p for p in d.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    return {p.relative_to(ROOT).as_posix():sha(p) for p in files}
def raw_data():
    # Reuse definitions, not precomputed results. Existing harmonizer writes only to supplied section-17 work directory.
    spec=importlib.util.spec_from_file_location('bdsf_original_pipeline',ROOT/'13_reproducibility_package/13_pipeline.py')
    pipe=importlib.util.module_from_spec(spec);spec.loader.exec_module(pipe)
    qmap=pd.read_csv(ROOT/'13_reproducibility_package/inputs/question_map.csv',keep_default_na=False)
    pipe.check_reference(qmap)
    df,frames=pipe.harmonize(ROOT,qmap,WORK/'reconstructed_inputs')
    check('Raw 180 records',len(frames['03_harmonization/03_01_BDSF_Kazakhstan_Harmonized_Dataset.csv'])==180)
    check('159 adults',len(df)==159)
    check('No minors in analytical sample',not df.age_group_code.eq('LT18').any())
    check('Unique original technical IDs',df.respondent_id.nunique()==159)
    check('Languages 81 and 78',df.source_language.value_counts().to_dict()=={'KZ':81,'RU':78})
    x=df[[q+'_score' for q in ITEMS]].to_numpy(float)
    check('Valid ordinal encoding',np.isfinite(x).all() and np.isin(x,[0,25,50,75,100]).all())
    m1=np.column_stack([x[:,ix].mean(1) for ix in DIMS.values()]).mean(1)
    check('Original M1 formula',np.allclose(m1,df.BDSF_score.to_numpy(),atol=1e-12))
    old=pd.read_csv(ROOT/'05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv',keep_default_na=False).set_index('respondent_id').loc[df.respondent_id]
    check('Original stored scores match raw reconstruction',np.allclose(old.BDSF_score,m1,atol=1e-9))
    check('Original quality flags match raw reconstruction',np.array_equal(old.quality_flags,df.quality_flags))
    check('Exactly 20 flagged adults',int(df.quality_flags.ne('').sum())==20)
    check('Baseline mean',abs(m1.mean()-60.3347746331237)<1e-9)
    check('Baseline alpha',abs(alpha(x)-.928379270943)<1e-9)
    save(VALID/'raw_source_hashes.csv',[dict(file=p.name,sha256=sha(p)) for p in (ROOT/'00_sources').glob('*_Responses.csv')])
    return df,x,m1

def run_cfa(x):
    print('EXP01: ordinal CFA, lavaan WLSMV',flush=True)
    out=HERE/EXPS[0];ordinal=WORK/'ordinal_items.csv'
    save(ordinal,pd.DataFrame((x/25).astype(int),columns=ITEMS))
    exe=HERE/'runtime/R/bin/Rscript.exe';lib=HERE/'runtime/R-library'
    if not exe.exists():raise FileNotFoundError('Section 17 R runtime missing; run setup_runtime.ps1')
    env=os.environ.copy();env['R_LIBS_USER']=str(lib);env['R_HOME']=str(HERE/'runtime/R')
    env['LC_ALL']='C'
    env['R_LIBS']='';env['R_LIBS_SITE']=''
    proc=subprocess.run([str(exe),'--vanilla',str(out/'models.R'),str(ordinal),str(out),str(lib)],env=env,text=True,encoding='utf-8',errors='replace',capture_output=True)
    (out/'execution.log').write_text(proc.stdout+'\n'+proc.stderr,encoding='utf-8')
    check('R WLSMV execution',proc.returncode==0,proc.stdout+proc.stderr)
    result=pd.read_csv(out/'EXP01_CFA_RESULTS.csv')
    check('Two CFA models evaluated',set(result.model)=={'OneFactor','FourFactor'})
    check('WLSMV explicitly used',result.estimator.eq('WLSMV').all() and result.test.str.split(';').map(lambda v:'scaled.shifted' in v).all() and result.standard_errors.eq('robust.sem').all())
    check('CFA convergence',result.converged.all())
    required=['CFI_scaled','TLI_scaled','RMSEA_scaled','RMSEA90_low_scaled','RMSEA90_high_scaled','SRMR']
    check('Required CFA fit measures finite',np.isfinite(result[required].to_numpy(float)).all())
    for _,r in result.iterrows():
        check(r.model+' admissibility assessed',not pd.isna(r.admissible))
        # Inadmissible models remain a scientific result, not silently accepted or deleted.
        check(r.model+' RMSEA CI ordered',r.RMSEA90_low_scaled<=r.RMSEA90_high_scaled)
    save(HERE/'EXP01_CFA_RESULTS.csv',result)
    return result

def run_topk(df,x,m1):
    print('EXP02: Top-K low-score selection stability',flush=True)
    out=HERE/EXPS[1];r=np.corrcoef(x,rowvar=False);one=fit(r,1)
    coefficients=np.linalg.solve(r,one['L'][:,0]);m3=((x-x.mean(0))/x.std(0,ddof=0))@coefficients
    check('M3 oriented positively',stats.spearmanr(m1,m3).statistic>0)
    old=pd.read_csv(ROOT/'16_extended_experiments/results/EXP03_Scoring_Comparison.csv')
    check('M3 agrees with section 16',np.allclose(m3,old.M3,atol=1e-8))
    oldcoef=pd.read_csv(ROOT/'16_extended_experiments/results/EXP03_M3_Coefficients.csv').regression_coefficient.to_numpy()
    check('M3 coefficients preserved',np.allclose(coefficients,oldcoef,atol=1e-8))
    methods={'M1':m1,'M2':x.mean(1),'M3':m3}
    # 10 decimals removes machine-precision differences only; actual scores remain unrounded in export.
    # The minimum nonzero attainable M1/M2 score gap is much larger than this tolerance.
    keys={m:np.round(a,10) for m,a in methods.items()}
    ranks={m:stats.rankdata(a,method='average') for m,a in keys.items()}
    ids=df.respondent_id.astype(str).to_numpy();anon={v:f'P{i+1:03}' for i,v in enumerate(sorted(ids))}
    orders={m:sorted(range(159),key=lambda i:(keys[m][i],ids[i])) for m in methods}
    result=[];selection=[];ties=[]
    for p in [.1,.2,.3]:
        k=math.ceil(p*159);sets={m:set(order[:k]) for m,order in orders.items()}
        for m,indexes in sets.items():
            check(f'{m} Top{k} unique IDs',len(indexes)==k and len({ids[i] for i in indexes})==k)
            for pos,i in enumerate(orders[m][:k],1):selection.append(dict(budget_percent=int(p*100),K=k,method=m,selection_position=pos,case=anon[ids[i]],score=methods[m][i]))
            threshold=keys[m][orders[m][k-1]];equal=np.where(keys[m]==threshold)[0]
            ties.append(dict(method=m,K=k,cutoff_score=threshold,tied_at_cutoff=len(equal),selected_from_tie=len(set(equal)&indexes),tie_split=len(set(equal)&indexes)<len(equal)))
        for a,b in [('M1','M2'),('M1','M3'),('M2','M3')]:
            intersection=len(sets[a]&sets[b]);union=len(sets[a]|sets[b]);j=intersection/union;delta=np.abs(ranks[a]-ranks[b])
            check(f'Jaccard bounds {a} {b} {k}',0<=j<=1)
            result.append(dict(method1=a,method2=b,budget_percent=int(p*100),K=k,n=159,spearman=stats.spearmanr(keys[a],keys[b]).statistic,kendall_tau_b=stats.kendalltau(keys[a],keys[b]).statistic,jaccard=j,overlap_n=intersection,overlap_percent=100*intersection/k,mean_abs_rank_change=delta.mean(),max_abs_rank_change=delta.max()))
    scoretable=pd.DataFrame({'case':[anon[v] for v in ids],**methods,**{m+'_ascending_rank':a for m,a in ranks.items()}})
    save(out/'scores.csv',scoretable);save(out/'topk_selections.csv',selection);save(out/'cutoff_ties.csv',ties)
    save(out/'M3_coefficients.csv',[dict(item=q,coefficient=coefficients[i],mean=x[:,i].mean(),sd_population=x[:,i].std(ddof=0)) for i,q in enumerate(ITEMS)])
    result=save(out/'results.csv',result);save(HERE/'EXP02_TOPK_RESULTS.csv',result)
    check('Top-K result finite',np.isfinite(result.select_dtypes('number').to_numpy()).all())
    # Independent set arithmetic against the exported selection lists.
    ss=pd.DataFrame(selection)
    for _,row in result.iterrows():
        a=set(ss[(ss.method==row.method1)&(ss.K==row.K)]['case']);b=set(ss[(ss.method==row.method2)&(ss.K==row.K)]['case'])
        check(f'Exported sets agree {row.method1} {row.method2} {row.K}',len(a&b)==row.overlap_n and abs(len(a&b)/len(a|b)-row.jaccard)<1e-12)
    return result

def run_quality(df,x,m1):
    print('EXP03: individual and combined quality-flag sensitivity',flush=True)
    out=HERE/EXPS[2]
    flags={f:df.quality_flags.str.contains(f).to_numpy() for f in ['DUPLICATE_RESPONSE_PATTERN','DUPLICATE_TIMESTAMP','STRAIGHTLINE_ALL_14']}
    masks={'S0':np.ones(159,bool),'S1':~flags['DUPLICATE_RESPONSE_PATTERN'],'S2':~flags['DUPLICATE_TIMESTAMP'],'S3':~flags['STRAIGHTLINE_ALL_14'],'S4':df.quality_flags.eq('').to_numpy()}
    check('S4 has 139',int(masks['S4'].sum())==139)
    overlap=[]
    for a,b in itertools.combinations(flags,2):overlap.append(dict(flag1=a,flag2=b,flag1_n=int(flags[a].sum()),flag2_n=int(flags[b].sum()),intersection_n=int((flags[a]&flags[b]).sum()),identical_sets=bool(np.array_equal(flags[a],flags[b]))))
    save(out/'flag_overlaps.csv',overlap)
    save(out/'flag_combinations.csv',df.quality_flags.replace('','NONE').value_counts().rename_axis('flags').reset_index(name='n'))
    cache={};rows=[];pa_rows=[]
    for scenario,mask in masks.items():
        check(scenario+' is subset of S0',mask.shape==(159,) and set(df.loc[mask,'respondent_id'])<=set(df.respondent_id))
        key=tuple(np.flatnonzero(mask))
        if key in cache:
            result=dict(cache[key]);result['same_data_as']='S1';rows.append(dict(scenario=scenario,**result));continue
        a=x[mask];scores=m1[mask];n=len(a)
        if n<=15:raise ValueError(scenario+' sample too small; cannot silently compute statistics')
        r=np.corrcoef(a,rowvar=False)
        check(scenario+' finite positive definite Pearson matrix',np.isfinite(r).all() and np.linalg.eigvalsh(r).min()>0)
        f=fit(r,1);w=omega(f,a.std(0,ddof=1));check(scenario+' admissible one-factor omega',np.isfinite(w) and 0<=w<=1)
        # Same Gaussian null PA method and seed as original project, recomputed from original CSV-derived sample.
        rng=np.random.default_rng(20261030);obs=np.linalg.eigvalsh(r)[::-1]
        null=np.array([np.linalg.eigvalsh(np.corrcoef(rng.normal(size=a.shape),rowvar=False))[::-1] for _ in range(500)])
        threshold=np.quantile(null,.95,axis=0);retained=int((obs>threshold).sum())
        for i in range(14):pa_rows.append(dict(scenario=scenario,component=i+1,observed_eigenvalue=obs[i],random95=threshold[i],retain=bool(obs[i]>threshold[i])))
        result=dict(n=n,excluded_n=159-n,BDSF_mean=float(scores.mean()),BDSF_median=float(np.median(scores)),BDSF_sd=float(scores.std(ddof=1)),alpha=float(alpha(a)),omega_raw_sum=float(w),**diagnostics(r,n),parallel_retained_factors=retained,same_data_as='')
        cache[key]=result;rows.append(dict(scenario=scenario,**result))
    t=pd.DataFrame(rows);base=t.iloc[0]
    for col in ['BDSF_mean','BDSF_median','BDSF_sd','alpha','omega_raw_sum','KMO']:
        t[col+'_delta_vs_S0']=t[col]-base[col];t[col+'_relative_percent_vs_S0']=100*(t[col]-base[col])/base[col]
    check('All quality results finite',np.isfinite(t.select_dtypes('number').to_numpy()).all())
    check('S0 baseline alpha',abs(t.iloc[0].alpha-.928379270943)<1e-9)
    old=pd.read_csv(ROOT/'16_extended_experiments/results/EXP02_Quality_Sensitivity.csv').set_index('n')
    for n in [159,139]:
        now=t[t.n==n].iloc[0]
        check(f'{n} agrees with previous mean alpha omega',all(abs(float(now[c])-float(old.loc[n,c]))<1e-8 for c in ['BDSF_mean','alpha','omega_raw_sum']))
    save(out/'results.csv',t);save(out/'parallel_analysis.csv',pa_rows);save(HERE/'EXP03_QUALITY_RESULTS.csv',t)
    return t

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--only',type=int,choices=[1,2,3]);args=parser.parse_args()
    protected=snapshot();jsonsave(VALID/'STATUS.json',dict(status='RUNNING'))
    try:
        df,x,m1=raw_data();results={}
        if args.only in [None,1]:results[1]=run_cfa(x)
        if args.only in [None,2]:results[2]=run_topk(df,x,m1)
        if args.only in [None,3]:results[3]=run_quality(df,x,m1)
        check('Previous files unchanged',snapshot()==protected)
        from reports import build_reports
        build_reports(HERE,results)
        if args.only is None:
            test=subprocess.run([sys.executable,str(VALID/'tests/test_results.py')],capture_output=True,text=True,encoding='utf-8',errors='replace')
            (VALID/'independent_tests.log').write_text(test.stdout+'\n'+test.stderr,encoding='utf-8')
            check('Seven independent test groups',test.returncode==0,test.stdout+test.stderr)
        versions={}
        from importlib.metadata import version
        for package in ['numpy','pandas','scipy','factor-analyzer','scikit-learn','matplotlib','openpyxl']:versions[package]=version(package)
        jsonsave(VALID/'environment.json',dict(python=sys.version,platform=platform.platform(),packages=versions,seed=SEED,parallel_analysis_seed=20261030))
        (HERE/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in versions.items())+'\n',encoding='utf-8')
        save(VALID/'checks.csv',CHECKS)
        lines=['# Автоматтандырылған тексерулер','',f'Тексеру саны: {len(CHECKS)}. Бәрі PASS.','',*['- '+c['status']+' — '+c['check'] for c in CHECKS]]
        (HERE/'ALL_VALIDATION_CHECKS.md').write_text('\n'.join(lines),encoding='utf-8')
        jsonsave(VALID/'STATUS.json',dict(status='PASS' if args.only is None else 'PARTIAL_PASS',experiments_completed=list(results),checks=len(CHECKS),previous_files_count=len(protected),previous_hashes=protected,output_hashes={p.relative_to(HERE).as_posix():sha(p) for p in HERE.rglob('*.csv') if 'runtime' not in p.parts and 'work' not in p.parts}))
        print(f'PASS: experiments {list(results)}; {len(CHECKS)} checks; {len(protected)} prior files unchanged',flush=True)
    except Exception as e:
        import traceback
        jsonsave(VALID/'STATUS.json',dict(status='FAIL',error=str(e),checks=CHECKS,traceback=traceback.format_exc()));raise
if __name__=='__main__':main()
