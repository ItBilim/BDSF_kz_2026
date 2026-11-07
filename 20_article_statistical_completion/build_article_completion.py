"""Complete article tables from raw CSV and immutable archived experiments."""
from pathlib import Path
import os,sys,csv,json,hashlib,collections,argparse
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import pandas as pd
H=Path(__file__).resolve().parent;R=H.parent
sys.path.insert(0,str(R/'16_extended_experiments'))
from methods import alpha,omega,fit,diagnostics,polychoric
SOURCES={};DIFF=[];CHECK=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def source(rel):
    p=R/rel;SOURCES[rel]=sha(p);return p
def read(rel):return pd.read_csv(source(rel),float_precision='round_trip')
def check(name,ok):
    CHECK.append(dict(check=name,passed=bool(ok)))
    if not ok:raise AssertionError(name)
def save(name,df):pd.DataFrame(df).to_csv(H/name,index=False,encoding='utf-8-sig',float_format='%.17g')
def comparison(scope,metric,new,old,tol):
    diff=float(new)-float(old)
    DIFF.append(dict(scope=scope,metric=metric,new=float(new),archived=float(old),difference=diff,tolerance=tol,within_tolerance=abs(diff)<=tol))
    check(scope+'/'+metric,abs(diff)<=tol)
def table(df):
    def fmt(x,col):
        if pd.isna(x):return 'NA'
        if isinstance(x,(int,np.integer)):return str(x)
        if isinstance(x,(float,np.floating)):
            if x!=0 and abs(x)<.001 and (str(col).endswith('_p') or col=='p'):return f'{x:.3e}'
            return f'{x:.3f}'
        return str(x)
    return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(x,c) for x,c in zip(row,df.columns))+' |' for row in df.itertuples(index=False,name=None))
def main():
    (H/'STATUS.json').write_text('{"status":"RUNNING"}')
    # Explicit independent raw response decoder; no previous harmonizer is called.
    mapping={}
    for words in [('Ешқашан','Сирек','Кейде','Жиі','Әрқашан'),('Никогда','Редко','Иногда','Часто','Всегда'),('Мүлдем көз жеткізбеймін','Сирек көз жеткіземін','Кейде көз жеткіземін','Жиі көз жеткіземін','Әрқашан толық көз жеткіземін'),('Никогда не проверяю','Редко проверяю','Иногда проверяю','Часто проверяю','Всегда тщательно проверяю')]:mapping.update(dict(zip(words,[0,25,50,75,100])))
    qmap=read('13_reproducibility_package/inputs/question_map.csv');raw=[]
    for lang,label,col in [('KZ','Kazakh','kazakh_question'),('RU','Russian','russian_question')]:
        p=next((R/'00_sources').glob('*_'+label+'_Responses.csv'));source(p.relative_to(R).as_posix())
        with p.open(encoding='utf-8-sig',newline='') as f:
            rr=csv.reader(f);header=next(rr);check(lang+' question order',[v.strip() for v in header[4:]]==qmap[col].tolist())
            for i,row in enumerate(rr,1):
                assert len(row)==18 and all(row)
                raw.append(dict(id=f'{lang}{i:03}',lang=lang,minor=row[3] in ['18 жасқа дейін','Младше 18 лет'],timestamp=row[0],x=tuple(mapping[v] for v in row[4:])))
    check('180 raw / 21 minors',len(raw)==180 and sum(r['minor'] for r in raw)==21)
    pat=collections.Counter(r['x'] for r in raw);times=collections.Counter(r['timestamp'] for r in raw)
    adults=[r for r in raw if not r['minor']];check('159 retained',len(adults)==159)
    x=np.array([r['x'] for r in adults],float);scores=np.column_stack([x[:,:4].mean(1),x[:,4:7].mean(1),x[:,7:11].mean(1),x[:,11:].mean(1)]).mean(1)
    stored=read('05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv').set_index('respondent_id').loc[[r['id'] for r in adults]]
    check('all item scores exact',np.array_equal(x,stored[[f'Q{i:02}_score' for i in range(1,15)]].to_numpy()))
    check('all M1 scores',np.allclose(scores,stored.BDSF_score,rtol=0,atol=1e-9))
    dup=np.array([pat[r['x']]>1 for r in adults]);time=np.array([times[r['timestamp']]>1 for r in adults]);straight=np.array([len(set(r['x']))==1 for r in adults])
    check('S1/S3 identical 18 exclusions',np.array_equal(dup,straight) and dup.sum()==18)
    save('BDSF_FILTER_AUDIT.csv',[dict(respondent_id=r['id'],S1_excluded=bool(dup[i]),S2_excluded=bool(time[i]),S3_excluded=bool(straight[i])) for i,r in enumerate(adults)])
    masks={'S0':np.ones(159,bool),'S1':~dup,'S2':~time,'S3':~straight,'S4':~(dup|time|straight)}
    old=read('18_final_sensitivity_topk/RQ2_SENSITIVITY.csv').set_index('scenario')
    rows=[];cache={};pa_rows=[]
    source('16_extended_experiments/methods.py')
    for s,mask in masks.items():
        key=tuple(np.flatnonzero(mask))
        if key in cache:row=dict(cache[key]);row['scenario']=s;rows.append(row);continue
        print('Recalculate '+s,flush=True)
        a=x[mask];r=np.corrcoef(a,rowvar=False);n=len(a)
        row=dict(scenario=s,n=n,M=float(scores[mask].mean()),SD=float(scores[mask].std(ddof=1)),Cronbach_alpha=float(alpha(a)),McDonald_omega=float(omega(fit(r,1),a.std(0,ddof=1))),**diagnostics(r,n))
        rng=np.random.default_rng(20261030)
        null=np.array([np.linalg.eigvalsh(np.corrcoef(rng.normal(size=a.shape),rowvar=False))[::-1] for _ in range(500)])
        eigen=np.linalg.eigvalsh(r)[::-1];limit=np.quantile(null,.95,axis=0)
        row['Pearson_PA_factors']=int((eigen>limit).sum())
        _,pr,_,_=polychoric(a);obs=np.linalg.eigvalsh(pr)[::-1]
        if s=='S0':p=source('16_extended_experiments/cache/poly_pa_N159.npz')
        else:
            p=source('18_final_sensitivity_topk/cache/'+s+'.npz');meta=json.loads(source('18_final_sensitivity_topk/cache/'+s+'.json').read_text())
            fingerprint=hashlib.sha256(a.astype('<f8').tobytes()+b'20261107:200'+sha(R/'16_extended_experiments/methods.py').encode()).hexdigest()
            check(s+' poly cache fingerprint',meta['fingerprint']==fingerprint)
        checkpoint=dict(np.load(p,allow_pickle=False));check(s+' poly null shape',checkpoint['null'].shape==(200,14))
        rng=np.random.default_rng(20261107)
        for i in range(2):
            perm=np.column_stack([rng.permutation(a[:,j]) for j in range(14)])
            _,m,_,_=polychoric(perm)
            check(s+f' poly permutation {i+1}',np.allclose(np.linalg.eigvalsh(m)[::-1],checkpoint['null'][i],atol=1e-10,rtol=0))
        threshold=np.quantile(checkpoint['null'],.95,axis=0);row['Polychoric_PA_factors']=int((obs>threshold).sum())
        for j in range(14):pa_rows.append(dict(scenario=s,component=j+1,Pearson_observed=eigen[j],Pearson_null95=limit[j],Polychoric_observed=obs[j],Polychoric_null95=threshold[j]))
        for metric in ['M','SD','Cronbach_alpha','McDonald_omega','KMO','Bartlett_chi2','Bartlett_df','Pearson_PA_factors','Polychoric_PA_factors']:
            comparison(s,metric,row[metric],old.loc[s,metric],1e-8)
        comparison(s,'Bartlett_log_p',np.log(row['Bartlett_p']),np.log(old.loc[s,'Bartlett_p']),1e-8)
        cache[key]=dict(row);rows.append(row)
    full=pd.DataFrame(rows)
    for short,col in [('M','M'),('alpha','Cronbach_alpha'),('omega','McDonald_omega')]:full['delta_'+short]=full[col]-full.iloc[0][col]
    save('BDSF_S0_S4_FULL_RESULTS.csv',full);save('BDSF_PA_CHECKS.csv',pa_rows)
    # Extract published CFA estimates; do not refit, round, or replace originals.
    loads=[];cors=[]
    for model in ['OneFactor','FourFactor']:
        st=read(f'19_reviewer_experiments/results/parameters/S0_{model}_standardized.csv')
        l=st[st.op=='=~'][['lhs','rhs','est.std','se','ci.lower','ci.upper']].copy();l.columns=['factor','item','loading','SE','CI95_lower','CI95_upper'];l.insert(0,'model',model);loads.append(l)
        check(model+' 14 loadings',len(l)==14)
        c=st[(st.op=='~~')&(st.lhs!=st.rhs)&st.lhs.isin(['D1','D2','D3','D4'])][['lhs','rhs','est.std','se','ci.lower','ci.upper']].copy();c.columns=['factor1','factor2','correlation','SE','CI95_lower','CI95_upper'];c.insert(0,'model',model);cors.append(c)
        save('BDSF_CFA_'+model.upper()+'_LOADINGS.csv',l)
    loads=pd.concat(loads,ignore_index=True);cors=pd.concat(cors,ignore_index=True);check('six factor correlations',len(cors)==6)
    save('BDSF_CFA_LOADINGS.csv',loads);save('BDSF_FACTOR_CORRELATIONS.csv',cors)
    cfa=read('19_reviewer_experiments/results/CFA_ALL_MODELS.csv');save('BDSF_CFA_MODEL_FIT.csv',cfa[~cfa.key.str.startswith('MG')])
    old17=read('17_q3_additional_experiments/EXP01_CFA_RESULTS.csv').set_index('model')
    for model in ['OneFactor','FourFactor']:
        row=cfa.set_index('key').loc['S0_'+model]
        for a,b in [('CFI','CFI_scaled'),('TLI','TLI_scaled'),('RMSEA','RMSEA_scaled'),('SRMR','SRMR'),('CFI_robust','CFI_robust'),('TLI_robust','TLI_robust'),('RMSEA_robust','RMSEA_robust')]:comparison(model,a,row[a],old17.loc[model,b],1e-5)
    mg=cfa[cfa.key.str.startswith('MG')].copy();mg['scope']='DIAGNOSTIC ONLY: Q01 categories 0/1 merged in BOTH languages; original scoring unchanged'
    mg['reason']=np.where(mg.key.str.contains('FourFactor'),'INVALID: KZ latent covariance matrix not positive definite','Poor configural fit prevents confirmation; exact LRT is exploratory')
    mg.loc[mg.key=='MG_OneFactor_configural','reason']+='; near-singular parameter vcov eigenvalue 1.14e-16 warning'
    mg['supports_invariance']=False;mg['step_LRT_p']=np.nan
    for ix,row in mg.iterrows():
        if not row['key'].endswith('configural'):
            lrt=read('19_reviewer_experiments/results/'+row['key']+'_LRT.csv');mg.loc[ix,'step_LRT_p']=lrt.iloc[-1]['Pr(>Chisq)']
    targets=[('MG_OneFactor_configural','RMSEA',.1261,.00005),('MG_OneFactor_configural','CFI_robust',.7751,.00005),('MG_OneFactor_metric','step_LRT_p',.03071,.000005),('MG_OneFactor_scalar','step_LRT_p',.000922,.0000005),('MG_FourFactor_configural','min_latent_eigen',-.0522,.00005)]
    for model,col,value,tol in targets:comparison(model,col,mg.set_index('key').loc[model,col],value,tol)
    save('BDSF_MG_CFA_DIAGNOSTICS.csv',mg)
    # Recheck both current files and the original relocation ledger; absent remote
    # work directory is explicitly recorded, never represented as a new rerun.
    ledger=read('19_reviewer_experiments/results/RELOCATION_HASHES.csv');check('48 scientific ledger entries',len(ledger)==48)
    relocation=json.loads(source('19_reviewer_experiments/results/RELOCATION_STATUS.json').read_text(encoding='utf-8'))
    relocated=Path(relocation['relocated_to'])/'19_reviewer_experiments/results';records=[]
    for row in ledger.to_dict('records'):
        p=source('19_reviewer_experiments/results/'+row['file']);q=relocated/row['file']
        current=sha(p);check('archived original '+row['file'],current==row['original_sha256'])
        check('ledger pair '+row['file'],row['original_sha256']==row['relocated_sha256'])
        other=sha(q) if q.exists() else ''
        if q.exists():check('current relocated '+row['file'],other==current)
        records.append(dict(file=row['file'],current_sha256=current,recorded_sha256=row['original_sha256'],relocated_available=q.exists(),relocated_current_sha256=other,identical=q.exists() and other==current))
    save('BDSF_48_FILE_CHECKS.csv',records)
    save('BDSF_NUMERICAL_DIFFERENCES.csv',DIFF);save('BDSF_VALIDATION_CHECKS.csv',CHECK)
    text='# Мақалаға арналған толық кестелер\n\n159 ересек өзгермеді; 180 бастапқы жауаптан тек 21 кәмелетке толмаған шығарылды. CSV файлдары 17 мәнді цифрмен сақталған; төмендегі кестелер үш ондық таңбамен берілген. Өте кіші p мәндері нөлге айналдырылмай ғылыми нотациямен көрсетілген. SD — ddof=1.\n\n## 9-кесте. S0–S4\n\n'+table(full[['scenario','n','M','SD','Cronbach_alpha','McDonald_omega','delta_M','delta_alpha','delta_omega']])
    text+='\n\n## 9-кестенің жалғасы. Факторлану және PA\n\n'+table(full[['scenario','KMO','Bartlett_chi2','Bartlett_df','Bartlett_p','Pearson_PA_factors','Polychoric_PA_factors']])
    text+='\n\nS0 — барлық 159 ересек; S1 — қайталанған жауап үлгілерінсіз; S2 — қайталанған уақыт белгілерінсіз; S3 — 14 сұраққа бірдей жауап бергендерсіз; S4 — үш жалаушаның кемінде бірі бар қатысушыларсыз. S1 және S3 бірдей 18 қатысушыны шығарады. ω — бұрынғы бір факторлы MINRES моделіне негізделген raw-item-sum omega; M1 өлшем салмақтарына арналған жаңа omega емес.\n'
    for model in ['OneFactor','FourFactor']:text+='\n## '+model+' CFA жүктемелері\n\n'+table(loads[loads.model==model].drop(columns='model'))+'\n'
    text+='\n## Төрт фактордың барлық корреляциялары\n\n'+table(cors.drop(columns='model'))
    text+='\n\n## Негізгі CFA сәйкестігі\n\n'+table(cfa[cfa.key.str.startswith('S0')][['key','converged','admissible','CFI','TLI','RMSEA','RMSEA_low','RMSEA_high','SRMR','CFI_robust','TLI_robust','RMSEA_robust']])
    text+='\n\nRMSEA CI — 90%; жүктемелер мен корреляциялардың CI — 95%. Төрт факторлы модель салыстырмалы жақсырақ, бірақ жоғары факторлық корреляциялар және robust/RMSEA индекстері төрт дербес өлшемді толық растауға мүмкіндік бермейді.\n\n## Диагностикалық MG ordinal CFA\n\n'+table(mg[['key','converged','admissible','RMSEA','CFI_robust','min_latent_eigen','step_LRT_p']])
    text+='\n\nRU Q01=0 санаты бос. Тек MG диагностикасында екі тілде де Q01=0/1 біріктірілген. Бастапқы бес санатты нұсқаның толық инварианттылығы тексерілген/расталған деп көрсетілмейді. Төрт факторлы MG модельдерінің барлығында KZ латенттік ковариациясы оң анықталған емес; олардың p мәндері растау дәлелі емес. Бір факторлы configural модельде нашар fit және параметрлік vcov тұрақтылығы ескертуі бар. NA — есептелмеген/қолжетімсіз индекс, нөл емес.\n'
    (H/'BDSF_ARTICLE_TABLES.md').write_text(text,encoding='utf-8')
    matched=sum(r['identical'] for r in records)
    report=f'''# Қайта өндіру есебі

Скрипт: `build_article_completion.py`; іске қосу: жоба түбірінен `20_article_statistical_completion\\RUN_TABLES.cmd`.

Бастапқы екі CSV-ден 180 жауап оқылды, 21 кәмелетке толмаған шығарылды, 159 ересектің 14 жауабы және M1 балдары сақталған нәтижелерге сәйкес келді. Сапа жалаушалары бастапқы 180 жауапта анықталды. S1 және S3 бірдей 18 ересекті шығарады; қалған n=141.

SD қайта есептелді: S1={full.loc[1,'SD']:.17g}; S2={full.loc[2,'SD']:.17g}; S3={full.loc[3,'SD']:.17g}. Мақалада: 19.090; 21.884; 19.090.

M, SD, α, ω, KMO, Bartlett және Pearson PA бастапқы CSV-ден қайта есептелді. Pearson PA: 500 Gaussian simulations, seed=20261030, 95% шек, бұрынғы әдіс. Polychoric: бастапқы жауаптардан матрица жаңадан есептелді; сақталған 200 permutation нөлдік спектр қолданылды (seed=20261107). Әр бірегей сценарийде алғашқы екі permutation жаңадан есептеліп checkpoint-пен салыстырылды; барлық 200 permutation қайта орындалған жоқ. S3=S1. Бұл есеп бұрынғы есептеу әдісін өзгертпейді.

CFA жаңа бағаланған жоқ: 19-бөлімнің толық дәлдіктегі сақталған CSV бағалары шығарылды. Scaled және robust көрсеткіштер сол күйінде сақталды. 17-бөліммен өте шағын сандық айырмашылықтар бар: бұрын delta, 19-бөлімде theta параметрлеуі қолданылған. Қайта оңтайландыру айырмашылықтары 1e-5 төзімділік ішінде; ешбір бұрынғы нәтиже алмастырылмады. Әр нақты айырмашылық BDSF_NUMERICAL_DIFFERENCES.csv ішінде; SD/басқа мәндердегі ұсақ айырмашылықтардың қосымша себебі — бұрынғы CSV-лердің 12 мәнді цифрмен сақталуы.

MG бастапқы есептің бес мәні өз көрсетілген дөңгелектеу дәлдігінде расталды. Төрт факторлы барлық MG модельдері жарамсыз. Конвергенция модельдің ғылыми жарамдылығына тең емес. Бір факторлы MG configural fit нашар және vcov жақын-сингуляр ескертуі бар; инварианттылық расталған жоқ.

48 ғылыми CSV: бастапқы файлдардың барлығы бұрынғы хешке сәйкес. Осы компьютердегі сақталған басқа жолдағы көшірмеден дәл сәйкес файлдар: {matched}/48. Бұл жаңа толық CFA rerun емес, бұрын қайта есептелген екі файл жиынын бүгінгі қайта тексеру. Жұмыс көшірмесі басқа компьютерде жоқ болса, скрипт оны NOT AVAILABLE ретінде тіркейді; архивтік сәйкестік тек тарихи дәлел ретінде қалады.

Тәуелсіз зерттеуші/бөлек компьютер/таза VM тексеруі орындалған жоқ. Желіні ОЖ деңгейінде бұғаттау болмаған. Есеп осы шектеуді жойылған деп көрсетпейді. Бұрынғы деректер мен нәтижелер тек оқылды. Қате болса STATUS=FAILED және nonzero exit code. PASS — техникалық тексеру мәртебесі, ғылыми гипотезаларды растау емес.

Кіріс файлдардың SHA-256 тізімі: BDSF_INPUT_MANIFEST.csv. Толық бастапқы бағалар алты сұралған негізгі файлда және қосымша диагностикалық CSV-лерде. Бастапқы мақала DOCX өзгертілген жоқ; мақалаға арналған кестелер BDSF_ARTICLE_TABLES.md ішінде.
'''
    (H/'BDSF_REPRODUCIBILITY_REPORT.md').write_text(report,encoding='utf-8')
    check('all sources unchanged',all(sha(R/k)==v for k,v in SOURCES.items()))
    save('BDSF_INPUT_MANIFEST.csv',[dict(path=k,sha256=v) for k,v in sorted(SOURCES.items())])
    save('BDSF_VALIDATION_CHECKS.csv',CHECK)
    (H/'STATUS.json').write_text(json.dumps(dict(status='PASS',n=159,checks=len(CHECK),relocated_csv_matches=matched,inputs_unchanged=True),indent=2))
    print(full[['scenario','n','M','SD','Cronbach_alpha','McDonald_omega']].to_string(index=False),flush=True)
    print('PASS',len(CHECK),'checks; relocation matches',matched,flush=True)
if __name__=='__main__':
    try:main()
    except Exception as e:
        (H/'STATUS.json').write_text(json.dumps(dict(status='FAILED',error=str(e))))
        save('BDSF_NUMERICAL_DIFFERENCES.csv',DIFF);raise
