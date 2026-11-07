"""Finish RQ2/RQ3 only, preserving all earlier stages and reusing verified results."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys,json,hashlib,argparse,itertools,math,traceback
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'16_extended_experiments'))
import numpy as np
import pandas as pd
from scipy import stats
from importlib.metadata import version
from methods import ITEMS,polychoric,alpha,diagnostics

SEED=20261107;PERMUTATIONS=200
DETAIL=HERE/'diagnostics';CACHE=HERE/'cache'
for p in [HERE,DETAIL,CACHE]:p.mkdir(parents=True,exist_ok=True)
CHECKS=[];NOTES=[];SOURCES={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(name,ok,details=''):
    CHECKS.append(dict(check=name,status='PASS' if bool(ok) else 'FAIL',details=str(details)))
    if not ok:raise AssertionError(name+': '+str(details))
def save(p,data):
    t=data if isinstance(data,pd.DataFrame) else pd.DataFrame(data)
    t.to_csv(p,index=False,encoding='utf-8-sig',float_format='%.12g');return t
def dump(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)),encoding='utf-8')
def read(rel):
    p=ROOT/rel;SOURCES[rel]=sha(p);return pd.read_csv(p,keep_default_na=False)
def protected():
    return {p.relative_to(ROOT).as_posix():sha(p) for d in ROOT.iterdir() if d.is_dir() and d.name[:2].isdigit() and int(d.name[:2])<=17 for p in d.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
def ledger_verify(stage,validation):
    v=json.loads((ROOT/validation).read_text(encoding='utf-8'))
    check(stage+' validation status',v['status']=='PASS')
    entries=v.get('output_hashes',v.get('output_sha256',{}))
    base=ROOT/stage
    if stage.startswith('16'):base=base/'results'
    for rel,h in entries.items():check(stage+' verified output '+rel,sha(base/rel)==h)
    return v

def new_poly_pa(x,scenario,force=False):
    key=hashlib.sha256(x.astype('<f8').tobytes()+f'{SEED}:{PERMUTATIONS}'.encode()+sha(ROOT/'16_extended_experiments/methods.py').encode()).hexdigest()
    cache=CACHE/(scenario+'.npz');meta=CACHE/(scenario+'.json')
    reusable=cache.exists() and meta.exists() and json.loads(meta.read_text())['fingerprint']==key
    if reusable and not force:
        c=dict(np.load(cache,allow_pickle=False));print(scenario+': verified checkpoint reused',flush=True)
    else:
        raw,r,pairs,note=polychoric(x)
        rng=np.random.default_rng(SEED);sim=[];smoothed=0
        for i in range(PERMUTATIONS):
            null=np.column_stack([rng.permutation(x[:,j]) for j in range(14)])
            _,matrix,_,nd=polychoric(null);sim.append(np.linalg.eigvalsh(matrix)[::-1]);smoothed+=int(nd['smoothed'])
            if (i+1)%25==0:print(f'{scenario}: polychoric permutation {i+1}/{PERMUTATIONS}',flush=True)
        c=dict(raw=raw,used=r,null=np.array(sim),pairs=json.dumps(pairs),note=json.dumps(note),null_smoothed=smoothed)
        np.savez_compressed(cache,**c);dump(meta,dict(fingerprint=key,input_n=len(x),seed=SEED,permutations=PERMUTATIONS))
    check(scenario+' finite polychoric and null matrices',np.isfinite(c['raw']).all() and np.isfinite(c['used']).all() and np.isfinite(c['null']).all())
    check(scenario+' valid used polychoric matrix',np.linalg.eigvalsh(c['used']).min()>0 and np.allclose(c['used'],c['used'].T) and np.allclose(np.diag(c['used']),1))
    # Independently regenerate two deterministic permutations, rather than merely trust cached output.
    rng=np.random.default_rng(SEED)
    for i in range(2):
        null=np.column_stack([rng.permutation(x[:,j]) for j in range(14)])
        _,matrix,_,_=polychoric(null)
        check(scenario+f' repeat permutation {i+1}',np.allclose(np.linalg.eigvalsh(matrix)[::-1],c['null'][i],atol=1e-10,rtol=0))
    observed=np.linalg.eigvalsh(c['used'])[::-1];threshold=np.quantile(c['null'],.95,axis=0)
    retained=observed>threshold
    pa=pd.DataFrame(dict(component=range(1,15),observed_eigenvalue=observed,null95=threshold,retain=retained))
    save(DETAIL/(scenario+'_POLY_PARALLEL.csv'),pa)
    for label in ['raw','used']:save(DETAIL/(scenario+'_POLY_'+label+'.csv'),pd.DataFrame(c[label],columns=ITEMS).assign(item=ITEMS))
    save(DETAIL/(scenario+'_POLY_PAIRS.csv'),json.loads(str(c['pairs'])))
    note=json.loads(str(c['note']));note.update(scenario=scenario,n=len(x),null_matrices_smoothed=int(c['null_smoothed']))
    return int(retained.sum()),note

def run(force=False,verify_repeat=False):
    before=protected();old_numerical={p.relative_to(HERE).as_posix():sha(p) for p in HERE.rglob('*.csv') if 'cache' not in p.parts}
    dump(HERE/'STATUS.json',dict(status='RUNNING'))
    try:
        v17=ledger_verify('17_q3_additional_experiments','17_q3_additional_experiments/validation/STATUS.json')
        v16=ledger_verify('16_extended_experiments','16_extended_experiments/results/VALIDATION_REPORT.json')
        # Link existing estimates to unchanged raw sources and to the actual analytical dataset.
        for rel,h in v17['previous_hashes'].items():
            if rel.startswith(('00_sources/','03_harmonization/','05_scoring/','13_reproducibility_package/inputs/')):
                check('Prior validated source unchanged '+rel,sha(ROOT/rel)==h)
        df=read('05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv');raw=read('03_harmonization/03_01_BDSF_Kazakhstan_Harmonized_Dataset.csv')
        ids=df.respondent_id.astype(str).to_numpy();x=df[[q+'_score' for q in ITEMS]].to_numpy(float)
        check('S0=159 and raw=180',len(df)==159 and len(raw)==180)
        check('All 21 minors excluded',raw.age_group_code.eq('LT18').sum()==21 and set(ids)==set(raw.loc[raw.age_group_code.ne('LT18'),'respondent_id']))
        check('14 valid complete items',np.isfinite(x).all() and np.isin(x,[0,25,50,75,100]).all())
        check('159 unique IDs',len(set(ids))==159)
        m1=(x[:,:4].mean(1)+x[:,4:7].mean(1)+x[:,7:11].mean(1)+x[:,11:].mean(1))/4
        check('M1 unchanged',np.allclose(m1,df.BDSF_score,atol=1e-9))
        # Recompute only quality flags for verification; existing statistics are reused.
        levels=[q+'_level' for q in ITEMS]
        flagraw=pd.DataFrame({'DUPLICATE_RESPONSE_PATTERN':raw.duplicated(levels,keep=False),'DUPLICATE_TIMESTAMP':raw.timestamp_iso.duplicated(keep=False),'STRAIGHTLINE_ALL_14':raw[levels].nunique(axis=1).eq(1)},index=raw.index)
        flagraw.index=raw.respondent_id;flags=flagraw.loc[ids]
        for name in flags:check('Flag definition '+name,np.array_equal(flags[name],df.quality_flags.str.contains(name)))
        union=flags.any(axis=1).to_numpy();masks={'S0':np.ones(159,bool),'S1':~flags.DUPLICATE_RESPONSE_PATTERN.to_numpy(),'S2':~flags.DUPLICATE_TIMESTAMP.to_numpy(),'S3':~flags.STRAIGHTLINE_ALL_14.to_numpy(),'S4':~union}
        check('S4=139 and union counted once',masks['S4'].sum()==139 and union.sum()==20)
        overlap=[]
        for a,b in itertools.combinations(flags.columns,2):
            aa=flags[a].to_numpy();bb=flags[b].to_numpy();overlap.append(dict(flag1=a,flag2=b,flag1_n=int(aa.sum()),flag2_n=int(bb.sum()),intersection_n=int((aa&bb).sum()),union_n=int((aa|bb).sum()),same_ids=bool(np.array_equal(aa,bb))))
        save(HERE/'RQ2_FILTER_OVERLAP.csv',overlap)
        old=read('17_q3_additional_experiments/EXP03_QUALITY_RESULTS.csv').set_index('scenario')
        oldpa=read('17_q3_additional_experiments/EXP03_Quality_Sensitivity/parallel_analysis.csv')
        scores=read('17_q3_additional_experiments/EXP02_TopK_Stability/scores.csv')
        coef=read('17_q3_additional_experiments/EXP02_TopK_Stability/M3_coefficients.csv')
        cases={v:f'P{i+1:03}' for i,v in enumerate(sorted(ids))};mapback={v:k for k,v in cases.items()}
        check('Same 159 participants for all methods',set(scores['case'])==set(mapback) and scores['case'].nunique()==159)
        scores=scores.set_index('case').loc[[cases[v] for v in ids]]
        check('M1 previous result',np.allclose(scores.M1,m1,atol=1e-9))
        check('M2 previous result',np.allclose(scores.M2,x.mean(1),atol=1e-9))
        check('M3 exact standardization means and population SD',np.allclose(coef['mean'],x.mean(0)) and np.allclose(coef.sd_population,x.std(0,ddof=0)))
        m3=((x-x.mean(0))/x.std(0,ddof=0))@coef.coefficient.to_numpy(float)
        check('M3 previous regression coefficients reproduce scores',np.allclose(m3,scores.M3,atol=1e-8))
        check('M3 positive security direction',stats.spearmanr(m3,m1).statistic>0)
        # No new factor fitting, scoring model or category for M3.
        previous=read('17_q3_additional_experiments/EXP02_TOPK_RESULTS.csv')
        selection17=read('17_q3_additional_experiments/EXP02_TopK_Stability/topk_selections.csv')
        selected=selection17[['budget_percent','K','method','selection_position','case']].copy()
        selected['technical_id']=selected['case'].map(mapback);selected=selected.drop(columns='case')
        check('No missing technical IDs',selected.technical_id.notna().all())
        # Validate the existing deterministic selection using full precision M1/M2 and original M3 formula.
        methods={'M1':m1,'M2':x.mean(1),'M3':m3}
        for (method,k),sub in selected.groupby(['method','K']):
            order=sorted(range(159),key=lambda i:(round(float(methods[method][i]),10),ids[i]))
            expect=[ids[i] for i in order[:k]]
            check(f'{method} Top{k} exact unique deterministic IDs',len(sub)==k and sub.technical_id.nunique()==k and sub.sort_values('selection_position').technical_id.tolist()==expect)
        results=[]
        for k in [16,32,48]:
            for a,b in [('M1','M2'),('M1','M3'),('M2','M3')]:
                aa=set(selected[(selected.K==k)&(selected.method==a)].technical_id);bb=set(selected[(selected.K==k)&(selected.method==b)].technical_id)
                intersection=len(aa&bb);union_size=len(aa|bb);j=intersection/union_size;ov=intersection/k;changed=k-intersection
                check(f'Top{k} {a}-{b} valid overlap bounds',0<=j<=1 and 0<=ov<=1 and len(aa-bb)==changed and len(bb-aa)==changed)
                prior=previous[(previous.K==k)&(previous.method1==a)&(previous.method2==b)].iloc[0]
                check(f'Top{k} {a}-{b} previous result',abs(j-float(prior.jaccard))<1e-10 and intersection==int(prior.overlap_n))
                results.append(dict(K=k,comparison=a+'–'+b,Jaccard=j,Overlap_at_K=ov,N_changed=changed,intersection_n=intersection,union_n=union_size,Spearman_rho=float(prior.spearman),Kendall_tau_b=float(prior.kendall_tau_b)))
        top=save(HERE/'RQ3_TOPK_RESULTS.csv',results);save(HERE/'RQ3_TOPK_SELECTED.csv',selected)
        check('Nine Top-K comparisons',len(top)==9)
        # Reuse S0 polychoric parallel analysis from section 16 after input/output provenance checks.
        s0pa=read('16_extended_experiments/results/EXP01_Polychoric_Parallel.csv')
        poly_counts={'S0':int((s0pa.observed>s0pa.null95).sum())};poly_diagnostics=[]
        s0note=read('16_extended_experiments/results/POLY_Matrix_Diagnostics.csv');s0note=s0note[s0note['sample']=='N159'].iloc[0]
        poly_diagnostics.append(dict(scenario='S0',n=159,min_eigenvalue_raw=float(s0note.min_eigenvalue_raw),smoothed=bool(s0note.smoothed),max_abs_adjustment=float(s0note.max_abs_adjustment),null_matrices_smoothed=int(np.load(ROOT/'16_extended_experiments/cache/poly_pa_N159.npz')['smoothed'])))
        save(DETAIL/'S0_POLY_PARALLEL.csv',s0pa.rename(columns={'observed':'observed_eigenvalue'}))
        completed={tuple(np.flatnonzero(masks['S0'])):'S0'}
        quality=[]
        for scenario,mask in masks.items():
            a=x[mask];row=old.loc[scenario]
            check(scenario+' n matches previous and subset of adults',len(a)==int(row.n) and set(ids[mask])<=set(ids))
            check(scenario+' mean SD alpha verified',abs(m1[mask].mean()-float(row.BDSF_mean))<1e-9 and abs(m1[mask].std(ddof=1)-float(row.BDSF_sd))<1e-9 and abs(alpha(a)-float(row.alpha))<1e-9)
            diag=diagnostics(np.corrcoef(a,rowvar=False),len(a));check(scenario+' KMO Bartlett verified',abs(diag['KMO']-float(row.KMO))<1e-9 and abs(diag['Bartlett_chi2']-float(row.Bartlett_chi2))<1e-7)
            pa_scenario='S1' if scenario=='S3' and np.array_equal(masks['S1'],mask) else scenario
            pearson=oldpa[oldpa.scenario==pa_scenario]
            check(scenario+' prior Pearson eigenvalues',len(pearson)==14 and np.allclose(pearson.observed_eigenvalue,np.linalg.eigvalsh(np.corrcoef(a,rowvar=False))[::-1],atol=1e-9))
            key=tuple(np.flatnonzero(mask));origin=completed.get(key)
            if scenario!='S0':
                if origin:
                    poly_counts[scenario]=poly_counts[origin]
                    reference=next(v for v in poly_diagnostics if v['scenario']==origin);poly_diagnostics.append(dict(reference,scenario=scenario))
                    prior_pa=pd.read_csv(DETAIL/(origin+'_POLY_PARALLEL.csv'));save(DETAIL/(scenario+'_POLY_PARALLEL.csv'),prior_pa)
                    print(scenario+': identical participant set, reused '+origin,flush=True)
                else:
                    count,note=new_poly_pa(a,scenario,force);poly_counts[scenario]=count;poly_diagnostics.append(note);completed[key]=scenario
            quality.append(dict(scenario=scenario,n=len(a),M=float(row.BDSF_mean),SD=float(row.BDSF_sd),Cronbach_alpha=float(row.alpha),McDonald_omega=float(row.omega_raw_sum),KMO=float(row.KMO),Bartlett_chi2=float(row.Bartlett_chi2),Bartlett_df=int(row.Bartlett_df),Bartlett_p=float(row.Bartlett_p),Pearson_PA_factors=int(row.parallel_retained_factors),Polychoric_PA_factors=poly_counts[scenario],delta_M=float(row.BDSF_mean)-float(old.loc['S0','BDSF_mean']),delta_alpha=float(row.alpha)-float(old.loc['S0','alpha']),delta_omega=float(row.omega_raw_sum)-float(old.loc['S0','omega_raw_sum']),statistics_source='17_q3_additional_experiments',polychoric_source='16 existing' if scenario=='S0' else ('18 '+origin+' identical sample' if origin else '18 new permutation PA')))
        q=save(HERE/'RQ2_SENSITIVITY.csv',quality);save(DETAIL/'POLYCHORIC_DIAGNOSTICS.csv',poly_diagnostics)
        check('All required estimates finite',np.isfinite(q.select_dtypes('number').to_numpy()).all() and np.isfinite(top.select_dtypes('number').to_numpy()).all())
        check('Previous S0 controls',abs(q.iloc[0].M-60.33)<.01 and abs(q.iloc[0].Cronbach_alpha-.928)<.001 and abs(q.iloc[0].McDonald_omega-.929)<.001)
        check('Previous S4 controls',abs(q.iloc[-1].M-60.79)<.01 and abs(q.iloc[-1].Cronbach_alpha-.896)<.001 and abs(q.iloc[-1].McDonald_omega-.897)<.001)
        if verify_repeat:
            check('Previous numerical output exists for repeat check',bool(old_numerical))
            for rel,h in old_numerical.items():check('Repeated numerical CSV '+rel,sha(HERE/rel)==h)
        check('All previous project files unchanged',before==protected())
        write_report(q,top,pd.DataFrame(overlap))
        requirements={p:version(p) for p in ['numpy','pandas','scipy','factor-analyzer','scikit-learn']}
        (HERE/'requirements.txt').write_text('\n'.join(f'{p}=={v}' for p,v in requirements.items())+'\n',encoding='utf-8')
        status=dict(status='PASS',checks=len(CHECKS),protected_files=len(before),seed=SEED,polychoric_permutations=PERMUTATIONS,repeat_checked=verify_repeat,dependencies=requirements,sources=SOURCES,original_files_unchanged=True,output_sha256={p.relative_to(HERE).as_posix():sha(p) for p in HERE.rglob('*.csv')})
        dump(HERE/'STATUS.json',status)
        lines=['# Автоматты тексеру','',f'PASS. Тексеру саны: {len(CHECKS)}. Бұрынғы {len(before)} файл өзгермеді.','',
          'Бұрынғы S0–S4 және Top-K нәтижелерімен сандық ауытқу анықталмады. Бұрынғы есептеулер қайта жазылмады.',
          'S1 және S3 бірдей 18 адамды белгілейді; олар екі тәуелсіз әсер емес. S4 union арқылы тек 20 бірегей адамды шығарады.',
          'Екі қосымша permutation әр жаңа polychoric сценарийде қайта есептеліп, сақталған нәтижемен салыстырылды.',
          'Бірдей кодты қайта іске қосқандағы барлық сандық CSV хештері тексерілді.' if verify_repeat else 'Барлық CSV хештері сақталды; толық қайталау тексеруі --verify-repeat арқылы жасалады.',
          'Polychoric PA: 200 маргиналдық үлестірімді сақтайтын permutation, 95-процентиль шегі. S0 — бұрынғы 16 бөлімнен; S1/S2/S4 — жаңа; S3=S1. Eigenvalue clipping қажет болса диагностикада көрсетіледі. Бұл толық measurement invariance не тәуелсіз валидация емес.',
          'McDonald ω — бұрынғы бір факторлы congeneric модельде 14 сұрақтың бақыланатын тең қосындысына арналған баға; негізгі M1 сенімділігімен автоматты тең емес.',
          'M3=Z R⁻¹λ; Z үшін сол 159 адамның item mean және population SD(ddof=0) пайдаланылды. Бұрынғы коэффициенттер өзгермеді. Іріктеу кілті бұрынғыдай 10 ондыққа дөңгелектеліп, тең мәндер техникалық ID арқылы реттеледі.',
          'N_changed=K−intersection — бір жақтағы ауысу; екі жиынның симметриялық айырмасы 2×N_changed. Overlap@K үлес ретінде 0–1, пайыз емес.',
          'RQ3_TOPK_SELECTED.csv тек техникалық ID, әдіс, бюджет, K және таңдау орнын қамтиды. Шикі жеке жауап, уақыт, жас немесе жыныс жаңа нәтижелерге көшірілген жоқ.',
          '',*['- '+c['status']+' — '+c['check'] for c in CHECKS]]
        (HERE/'VALIDATION_REPORT.md').write_text('\n\n'.join(lines),encoding='utf-8')
        print(f'PASS: {len(CHECKS)} checks; {len(before)} previous files unchanged',flush=True)
    except Exception as e:
        dump(HERE/'STATUS.json',dict(status='FAIL',error=str(e),checks=CHECKS,traceback=traceback.format_exc()))
        (HERE/'VALIDATION_REPORT.md').write_text('# FAIL\n\n'+str(e)+'\n\n'+traceback.format_exc(),encoding='utf-8');raise

def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def write_report(q,t,overlap):
    tq=table(['Сценарий','n','M','SD','α','ω','KMO','Bartlett χ²','p','PA Pearson','PA polychoric','ΔM','Δα','Δω'],[[r.scenario,int(r.n),f'{r.M:.2f}',f'{r.SD:.2f}',f'{r.Cronbach_alpha:.3f}',f'{r.McDonald_omega:.3f}',f'{r.KMO:.3f}',f'{r.Bartlett_chi2:.2f}',f'{r.Bartlett_p:.2e}',int(r.Pearson_PA_factors),int(r.Polychoric_PA_factors),f'{r.delta_M:.3f}',f'{r.delta_alpha:.4f}',f'{r.delta_omega:.4f}'] for _,r in q.iterrows()])
    tt=table(['K','Салыстыру','Jaccard','Overlap@K','Өзгерген қатысушылар'],[[int(r.K),r.comparison,f'{r.Jaccard:.4f}',f'{r.Overlap_at_K:.4f}',int(r.N_changed)] for _,r in t.iterrows()])
    s=q.set_index('scenario');stable=t.groupby('comparison').agg(mean_J=('Jaccard','mean'),mean_O=('Overlap_at_K','mean'),max_changed=('N_changed','max'))
    best=stable[stable.mean_J==stable.mean_J.max()].index.tolist()
    lines=['# RQ2','',tq,'',
      f"Жеке сүзгілер ішінде S1 және S3 көбірек әсер етті: екеуі бірдей 18 жауапты алып тастайды (n={int(s.loc['S1','n'])}). ΔM={s.loc['S1','delta_M']:.3f}, Δα={s.loc['S1','delta_alpha']:.4f}, Δω={s.loc['S1','delta_omega']:.4f}. Олардың әсерін бір-бірінен ажыратуға болмайды. S2 тек екі жауапқа қатысты: ΔM={s.loc['S2','delta_M']:.3f}, Δα={s.loc['S2','delta_alpha']:.4f}, Δω={s.loc['S2','delta_omega']:.4f}.",
      f"Біріктірілген S4: n=139, M={s.loc['S4','M']:.2f}, α={s.loc['S4','Cronbach_alpha']:.3f}, ω={s.loc['S4','McDonald_omega']:.3f}. Бастапқы S0: n=159, M={s.loc['S0','M']:.2f}, α={s.loc['S0','Cronbach_alpha']:.3f}, ω={s.loc['S0','McDonald_omega']:.3f}. Орташа мән аз өзгергенімен, S1/S3/S4 ішкі келісімділікті төмендетеді.",
      'Pearson фактор саны: '+', '.join(f'{r.scenario}={int(r.Pearson_PA_factors)}' for _,r in q.iterrows())+'. Polychoric фактор саны: '+', '.join(f'{r.scenario}={int(r.Polychoric_PA_factors)}' for _,r in q.iterrows())+'.',
      'Негізгі орташа, α, ω және Pearson нәтижелері бұрынғы S0–S4 есебімен сәйкес. Polychoric нәтижесі — сол іріктемеге қосымша сезімталдық тексеруі. Жалауша жалған жауапты дәлелдемейді; уақыт белгісінің қайталануы бір адамның қайталап қатысуының дәлелі емес. S1–S4 — S0 ішкі жиындары, сондықтан бұл тәуелсіз топтар салыстыруы емес.',
      '', '# RQ3','',tt,'',
      'Үш бюджет бойынша орташа Jaccard ең жоғары жұп: '+', '.join(best)+f' ({stable.mean_J.max():.4f}). Бұл — сипаттамалық қорытынды, статистикалық мәнділік тесті емес.']
    for k,percent in [(16,10),(32,20),(48,30)]:
        part=t[t.K==k];lines.append(f'Top-{percent}% (K={k}): '+', '.join(f'{r.comparison}: {int(r.N_changed)} ауысу, overlap={r.Overlap_at_K*100:.2f}%' for _,r in part.iterrows())+'.')
    lines+=['Толық рейтинг корреляциясы жоғары: Spearman ρ='+f'{t.Spearman_rho.min():.4f}–{t.Spearman_rho.max():.4f}'+'. Бұл Top-K толық сәйкестігін кепілдемейді: бюджет шегіндегі аз ранг айырмасы құрамды өзгерте алады. Тең балл техникалық ID бойынша реттелді; оның өзі тәуекел айырмашылығын дәлелдемейді.',
      'Өзгерген қатысушылар саны K−қиылысу ретінде бір жақтағы ауысуды білдіреді; екі жиынның симметриялық айырмасы одан екі есе үлкен. M3 үшін жаңа қауіпсіздік санаты енгізілмеді. Бұл нақты киберқауіпсіздік оқиғаларын болжауды емес, бағалау алгоритмінің тұрақтылығын сипаттайды.']
    (HERE/'FINAL_RESULTS.md').write_text('\n\n'.join(lines)+'\n',encoding='utf-8')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--recompute-poly',action='store_true');parser.add_argument('--verify-repeat',action='store_true');args=parser.parse_args()
    run(args.recompute_poly,args.verify_repeat)
