from pathlib import Path
import hashlib, json, sys, warnings
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.utils.validation import check_array
import factor_analyzer.factor_analyzer as fam
from factor_analyzer import FactorAnalyzer

def compatible(*a, **k):
    if 'force_all_finite' in k: k['ensure_all_finite'] = k.pop('force_all_finite')
    return check_array(*a, **k)
fam.check_array = compatible
ROOT = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
OUT = ROOT/'15_language_psychometrics'
OUT.mkdir(exist_ok=True)
source = ROOT/'05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv'
df = pd.read_csv(source)
excluded = pd.read_csv(ROOT/'04_cleaning/04_03_BDSF_Kazakhstan_Exclusion_Log.csv')
assert len(df)==159 and len(excluded)==21 and df.respondent_id.is_unique
assert not set(df.respondent_id)&set(excluded.respondent_id)
assert not (df.age_group_code=='LT18').any()
assert df.source_language.value_counts().to_dict()=={'KZ':81,'RU':78}
items=[f'Q{i:02}' for i in range(1,15)]
groups={'D1':items[:4],'D2':items[4:7],'D3':items[7:11],'D4':items[11:]}
x=df[[i+'_score' for i in items]].to_numpy(float)
assert np.isfinite(x).all() and np.isin(x,[0,25,50,75,100]).all()
assert np.allclose(np.column_stack([df[[i+'_score' for i in v]].mean(axis=1) for v in groups.values()]).mean(axis=1),df.BDSF_score)
rng=np.random.default_rng(20261107)
def alpha(a):
    v=a.sum(axis=1).var(ddof=1)
    return a.shape[1]/(a.shape[1]-1)*(1-a.var(axis=0,ddof=1).sum()/v) if v>0 else np.nan
def save(name, rows):
    t=pd.DataFrame(rows); t.to_csv(OUT/name,index=False,encoding='utf-8-sig',float_format='%.10g'); return t
reliability=[]; diagnostics=[]; distributions=[]; loadings=[]; fits={}; fitnotes=[]
for lang in ['KZ','RU']:
    sub=df[df.source_language==lang]; a=sub[[i+'_score' for i in items]].to_numpy(float)
    for scale, qs in {**groups,'ALL14':items}.items():
        b=sub[[i+'_score' for i in qs]].to_numpy(float)
        boot=np.array([alpha(b[rng.integers(0,len(b),len(b))]) for _ in range(2000)])
        assert np.isfinite(boot).all()
        lo,hi=np.quantile(boot,[.025,.975])
        reliability.append(dict(language=lang,scale=scale,n=len(b),items=len(qs),alpha=alpha(b),ci95_low=lo,ci95_high=hi))
    for dim,qs in groups.items():
        for q in qs:
            val=sub[q+'_score']; rest=sub[[j+'_score' for j in qs if j!=q]].sum(axis=1)
            diagnostics.append(dict(language=lang,item=q,dimension=dim,n=len(sub),corrected_item_total_within_dimension=val.corr(rest)))
    for q in items:
        for score in [0,25,50,75,100]:
            count=int((sub[q+'_score']==score).sum())
            distributions.append(dict(language=lang,item=q,score=score,count=count,percent=100*count/len(sub),n=len(sub)))
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        fa=FactorAnalyzer(n_factors=4,method='minres',rotation='oblimin'); fa.fit(a)
    fits[lang]=fa
    fitnotes.append(dict(language=lang,warnings=[str(w.message) for w in ws]))
# Fix factor identities using intended dimensions in KZ, then match RU to KZ by absolute cosine.
L=fits['KZ'].loadings_
aff=np.array([np.abs(L[[items.index(i) for i in qs]]).mean(axis=0) for qs in groups.values()])
_,order=linear_sum_assignment(-aff)
ref=L[:,order]
sign=np.array([1 if ref[[items.index(i) for i in qs],j].sum()>=0 else -1 for j,qs in enumerate(groups.values())])
ref=ref*sign
raw=fits['RU'].loadings_
cos=ref.T@raw/(np.linalg.norm(ref,axis=0)[:,None]*np.linalg.norm(raw,axis=0)[None,:])
_,ro=linear_sum_assignment(-np.abs(cos)); rs=np.sign(cos[np.arange(4),ro]); ru=raw[:,ro]*rs
aligned={'KZ':ref,'RU':ru}
congruence=[]
for j,dim in enumerate(groups):
    congruence.append(dict(factor=dim,tucker_congruence=float(ref[:,j]@ru[:,j]/(np.linalg.norm(ref[:,j])*np.linalg.norm(ru[:,j])))))
for lang in ['KZ','RU']:
    fa=fits[lang]; comm=np.diag(fa.loadings_@fa.phi_@fa.loadings_.T)
    for i,q in enumerate(items):
        loadings.append(dict(language=lang,item=q,**{f'loading_{d}':aligned[lang][i,j] for j,d in enumerate(groups)},communality=comm[i],uniqueness=1-comm[i]))
rel=save('15_01_Reliability.csv',reliability)
it=save('15_02_Item_Total.csv',diagnostics)
save('15_03_Factor_Loadings.csv',loadings)
con=save('15_03_Factor_Congruence.csv',congruence)
dist=save('15_04_Response_Distributions.csv',distributions)
assert (dist.groupby(['language','item'])['count'].sum().loc['KZ']==81).all()
assert (dist.groupby(['language','item'])['count'].sum().loc['RU']==78).all()
diff=[]
for q in items:
    kz=dist[(dist.language=='KZ')&(dist.item==q)].percent.to_numpy()
    ru_p=dist[(dist.language=='RU')&(dist.item==q)].percent.to_numpy()
    diff.append(dict(item=q,total_variation_percent=np.abs(kz-ru_p).sum()/2,max_category_difference_pp=np.abs(kz-ru_p).max()))
dd=save('15_04_Distribution_Differences.csv',diff).sort_values('total_variation_percent',ascending=False)
lines=['# Қазақша және орысша нұсқалардың қосымша психометриялық талдауы','',
'Іріктеме: 159 ересек (KZ = 81; RU = 78). 21 кәмелетке толмаған қатысушы енгізілмеді. Сапа жалаушалары бойынша қосымша шығару болған жоқ. Негізгі баллдар, салмақтар және жіктеу өзгерген жоқ.',
'','## 1. Ішкі келісімділік — Cronbach α','',
'Бірінші көрсеткіштің атауы бастапқы сұрауда көрсетілмегендіктен, Cronbach α қабылданды. Шикі α; әр тілде қатысушыларды 2000 рет bootstrap қайта іріктеу, percentile 95% CI, seed 20261107. ALL14 — барлық сұрақтың диагностикалық көрсеткіші; ол жаңа BDSF баллы емес.','',
'| Шкала | KZ α [95% CI] | RU α [95% CI] |','|---|---|---|']
for s in [*groups,'ALL14']:
    vals=[]
    for l in ['KZ','RU']:
        r=rel[(rel.scale==s)&(rel.language==l)].iloc[0]; vals.append(f'{r.alpha:.3f} [{r.ci95_low:.3f}; {r.ci95_high:.3f}]')
    lines.append('| '+s+' | '+' | '.join(vals)+' |')
lines+=['','## 2. Түзетілген item–total correlation','',
'Әр сұрақ өз өлшеміндегі қалған сұрақтар қосындысымен Pearson корреляциясы арқылы салыстырылды. Сұрақтың өзі қосындыдан шығарылды.']
for l in ['KZ','RU']:
    t=it[it.language==l].corrected_item_total_within_dimension
    lines.append(f'{l}: ең төменгі r = {t.min():.3f}; ең жоғарғы r = {t.max():.3f}.')
lines+=['','## 3. Факторлық жүктемелер','',
'Әр тілге жеке төрт факторлы MINRES EFA, Pearson матрицасы, oblimin айналдыруы есептелді. Төрт фактор алдын ала салыстыру үшін бекітілді; бұл төрт фактордың деректермен дәлелденгенін білдірмейді. KZ факторлары күтілетін өлшемдердің орташа абсолюттік жүктемесіне қарай бірегей сәйкестендірілді. RU факторлары KZ бағандарына абсолюттік косинусты максимизациялау арқылы сәйкестендіріліп, таңбалары түзетілді. D1–D4 белгілері салыстыру белгілері, өлшем валидтілігінің кепілі емес.','']
for r in congruence: lines.append(f"- {r['factor']}: Tucker сәйкестік коэффициенті = {r['tucker_congruence']:.3f}.")
lines+=['','## 4. Сұрақтардың жауап таралуы','',
'14 сұрақтың әрқайсысы үшін бес жауап санатының саны мен тіл ішіндегі пайызы берілді. Total variation = абсолюттік пайыз айырмалары қосындысының жартысы; 0 — бірдей, 100 — толық бөлінген таралулар. Бұл сипаттамалық салыстыру, статистикалық мәнділік тесті емес.','']
for _,r in dd.head(3).iterrows(): lines.append(f"- {r['item']}: total variation = {r.total_variation_percent:.2f}%; ең үлкен санат айырмасы = {r.max_category_difference_pp:.2f} пайыздық тармақ.")
lines+=['','## Шектеулер және ғылыми мәні','',
'Төрт талдау осы іріктемеде тілдер арасындағы ішкі келісімділік, сұрақ–өлшем байланысы, зерттеушілік фактор құрылымы және жауап таралуының ұқсастығы мен айырмашылығын сипаттайды. Тіл топтары тәуелсіз, рандомизацияланбаған: айырмашылықтарды аударма әсеріне ғана жатқызуға болмайды. 81/78 қатысушыдан алынған төрт факторлы шешімдер тұрақсыз болуы мүмкін. Likert жауаптарына Pearson тәсілі қолданылды; ordinal CFA жүргізілген жоқ. Bootstrap интервалдары жеке көрсеткіштерге арналған, көптік салыстыруға түзетілмеген. Толық measurement invariance, configural/metric/scalar invariance немесе DIF дәлелденген жоқ. Негізгі зерттеу қорытындысы өзгертілмейді.',
'','Әдіс құжаттамасы: https://factor-analyzer.readthedocs.io/en/latest/factor_analyzer.html','',
'Қайта іске қосу (жоба түбірінен): `runtime\\python\\python.exe -I -X utf8 15_language_psychometrics\\run_analysis.py`.',
'Бұл қосымша модуль дайын 05_scoring деректерін оқиды. Негізгі RUN_REVIEW.cmd және мақала автоматты түрде өзгертілген жоқ.']
(OUT/'15_00_Report_KK.md').write_text('\n'.join(lines),encoding='utf-8')
audit={'status':'PASS','adults':159,'excluded':21,'KZ':81,'RU':78,'seed':20261107,'bootstrap_resamples':2000,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'factor_fit_notes':fitnotes,'factor_matching':{'KZ_order':order.tolist(),'RU_order':ro.tolist()},'outputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.glob('*.csv'))}}
(OUT/'15_05_Validation.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print((OUT/'15_00_Report_KK.md').read_text(encoding='utf-8'))
