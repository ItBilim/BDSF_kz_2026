"""Additional reviewer experiments; independently decode raw CSVs, preserve stages 00-18."""
from pathlib import Path
import os,sys,csv,json,hashlib,subprocess,collections
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;OUT=HERE/'results'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,x):pd.DataFrame(x).to_csv(OUT/name,index=False,float_format='%.12g',encoding='utf-8-sig')
def main():
    OUT.mkdir(exist_ok=True)
    (OUT/'STATUS.json').write_text('{"status":"RUNNING"}')
    inputs=list((ROOT/'00_sources').glob('*Responses.csv'))
    before={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    mapping={}
    for words in [('Ешқашан','Сирек','Кейде','Жиі','Әрқашан'),('Никогда','Редко','Иногда','Часто','Всегда'),('Мүлдем көз жеткізбеймін','Сирек көз жеткіземін','Кейде көз жеткіземін','Жиі көз жеткіземін','Әрқашан толық көз жеткіземін'),('Никогда не проверяю','Редко проверяю','Иногда проверяю','Часто проверяю','Всегда тщательно проверяю')]:mapping.update(dict(zip(words,range(5))))
    qmap=pd.read_csv(ROOT/'13_reproducibility_package/inputs/question_map.csv')
    raw=[]
    for lang,label,col in [('KZ','Kazakh','kazakh_question'),('RU','Russian','russian_question')]:
        p=next((ROOT/'00_sources').glob('*_'+label+'_Responses.csv'))
        with p.open(encoding='utf-8-sig',newline='') as f:
            rr=csv.reader(f);header=next(rr)
            assert [s.strip() for s in header[4:]]==qmap[col].tolist()
            for i,r in enumerate(rr,1):
                assert len(r)==18 and all(r)
                raw.append(dict(id=f'{lang}{i:03}',language=lang,minor=r[3] in ['18 жасқа дейін','Младше 18 лет'],timestamp=r[0],x=tuple(mapping[v] for v in r[4:])))
    assert len(raw)==180 and sum(r['minor'] for r in raw)==21
    patterns=collections.Counter(r['x'] for r in raw);times=collections.Counter(r['timestamp'] for r in raw)
    adults=[r for r in raw if not r['minor']];assert len(adults)==159
    x=np.array([r['x'] for r in adults],float)*25
    dims=[slice(0,4),slice(4,7),slice(7,11),slice(11,14)]
    score=np.mean([x[:,s].mean(1) for s in dims],axis=0)
    old=pd.read_csv(ROOT/'05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv').set_index('respondent_id').loc[[r['id'] for r in adults]]
    assert np.allclose(score,old.BDSF_score,atol=1e-9,rtol=0)
    assert np.array_equal(x,old[[f'Q{i:02}_score' for i in range(1,15)]].to_numpy())
    f1=np.array([patterns[r['x']]>1 for r in adults]);f2=np.array([times[r['timestamp']]>1 for r in adults]);f3=np.array([len(set(r['x']))==1 for r in adults])
    masks={'S0':np.ones(159,bool),'S1':~f1,'S2':~f2,'S3':~f3,'S4':~(f1|f2|f3)}
    assert [int(m.sum()) for m in masks.values()]==[159,141,157,141,139]
    assert np.array_equal(f1,f3)
    reasons=['None','Duplicate response pattern','Duplicate timestamp','Straightline all 14','Union of all three flags']
    previous=pd.read_csv(ROOT/'18_final_sensitivity_topk/RQ2_SENSITIVITY.csv').set_index('scenario')
    rows=[]
    for (s,m),why in zip(masks.items(),reasons):
        a=x[m];alpha=14/13*(1-a.var(0,ddof=1).sum()/a.sum(1).var(ddof=1));scores=score[m]
        assert abs(alpha-previous.loc[s,'Cronbach_alpha'])<1e-9
        assert abs(scores.mean()-previous.loc[s,'M'])<1e-8
        z=previous.loc[s].to_dict();z.update(scenario=s,filter=why,excluded_n=int((~m).sum()),median=float(np.median(scores)),low=int((np.round(scores,2)<=33.32).sum()),medium=int(((np.round(scores,2)>33.32)&(np.round(scores,2)<=66.66)).sum()),high=int((np.round(scores,2)>66.66).sum()),KZ=sum(r['language']=='KZ' for r,keep in zip(adults,m) if keep),RU=sum(r['language']=='RU' for r,keep in zip(adults,m) if keep))
        rows.append(z)
    save('SENSITIVITY_COMPLETE.csv',rows)
    ordinal=pd.DataFrame((x/25).astype(int),columns=[f'Q{i:02}' for i in range(1,15)]);ordinal['language']=[r['language'] for r in adults]
    for s,m in masks.items():ordinal[s]=m.astype(int)
    save('ordinal_input.csv',ordinal)
    counts=[]
    for lang in ['KZ','RU']:
        for q in ordinal.columns[:14]:
            for value in range(5):counts.append(dict(language=lang,item=q,category=value,n=int(((ordinal.language==lang)&(ordinal[q]==value)).sum())))
    save('CATEGORY_COUNTS.csv',counts)
    empty=[r for r in counts if r['n']==0]
    assert empty==[dict(language='RU',item='Q01',category=0,n=0)],empty
    (OUT/'INVARIANCE_SCOPE.json').write_text(json.dumps({'original_5_category_invariance':'NOT_ESTIMABLE: RU Q01 category 0 absent','diagnostic_only':'Merge Q01 categories 0 and 1 in BOTH languages for MG models only; primary scoring and sensitivity unchanged','empty_cells':empty},indent=2))
    env=os.environ.copy();env.update(R_HOME=str(ROOT/'17_q3_additional_experiments/runtime/R'),R_LIBS='',R_LIBS_SITE='',R_LIBS_USER=str(HERE/'runtime/R-library'),LC_ALL='C')
    # R needs a space-free temp string; relative temp remains inside project.
    # UTF-8 locale is required for Kazakh/Cyrillic paths on Windows UCRT R.
    temp=HERE/'work/rtemp';temp.mkdir(parents=True,exist_ok=True)
    env.update(TMPDIR='rtemp',TEMP='rtemp',TMP='rtemp',LC_ALL='English_United States.utf8')
    proc=subprocess.run([str(ROOT/'17_q3_additional_experiments/runtime/R/bin/Rscript.exe'),'--vanilla',str(HERE/'models.R'),str(ROOT),str(OUT)],cwd=temp.parent,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (OUT/'execution.log').write_text(proc.stdout+'\n'+proc.stderr,encoding='utf-8')
    if proc.returncode:raise RuntimeError('R analysis failed; see execution.log')
    cfa=pd.read_csv(OUT/'CFA_ALL_MODELS.csv');assert len(cfa)==18
    oldcfa=pd.read_csv(ROOT/'17_q3_additional_experiments/EXP01_CFA_RESULTS.csv').set_index('model')
    checks=[]
    for model in ['OneFactor','FourFactor']:
        now=cfa.set_index('key').loc['S0_'+model]
        for a,b in [('CFI','CFI_scaled'),('RMSEA','RMSEA_scaled'),('SRMR','SRMR')]:
            delta=abs(now[a]-oldcfa.loc[model,b]);checks.append(dict(model=model,metric=a,absolute_difference=delta,tolerance=1e-5,passed=bool(delta<1e-5)))
    assert all(r['passed'] for r in checks)
    save('CFA_REPRODUCTION_CHECKS.csv',checks)
    assert before=={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    status=dict(status='PASS',raw_n=180,excluded_minors=21,adults=159,independent_raw_decoder=True,old_scores_reproduced=True,models_attempted=18,models_converged=int(cfa.converged.sum()),source_sha256=before,output_sha256={p.relative_to(OUT).as_posix():sha(p) for p in OUT.rglob('*.csv')})
    (OUT/'STATUS.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in status.items() if not k.endswith('sha256')},indent=2),flush=True)
if __name__=='__main__':
    try:main()
    except Exception as e:
        OUT.mkdir(exist_ok=True);(OUT/'STATUS.json').write_text(json.dumps({'status':'FAILED','error':str(e)}));raise
