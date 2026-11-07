"""Kazakh research report and a revised COPY of the existing manuscript."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

def num(v,n=3):
    return 'анықталмады' if v is None or not np.isfinite(float(v)) else f'{float(v):.{n}f}'

def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])

def narrative(e):
    a=e['summary159'];b=e['summary139'];cmp=e['model_comparison']
    models={(r['correlation'],r['factors']):r for r in cmp if r['sample']=='N159'}
    r12=e['rank'][0];r13=e['rank'][1];cat=e['categories']
    lk=next(r for r in e['language'] if r['language']=='KZ' and r['scale']=='ALL14');lr=next(r for r in e['language'] if r['language']=='RU' and r['scale']=='ALL14')
    return dict(
      rq1=f"Pearson параллель талдауы {int(a['parallel_retained'])}, polychoric permutation талдауы {e['poly_parallel_retained']} факторды қолдады. Бір факторлы Pearson шешімінің RMSR көрсеткіші {num(models['Pearson',1]['offdiag_rmsr'])}, төрт факторлы шешімде {num(models['Pearson',4]['offdiag_rmsr'])}; polychoric шешімдерде тиісінше {num(models['Polychoric',1]['offdiag_rmsr'])} және {num(models['Polychoric',4]['offdiag_rmsr'])}. Төрт факторлы модельдің қалдықтары азаюы оның теориялық D1–D4 құрылымын растауға жеткіліксіз: модель еркін параметрлері көбірек. Бір басым факторға дәлел бар, бірақ бұл бір факторлы модельдің мінсіз сәйкестігі деген сөз емес.",
      rq2=f"Негізгі n=159 іріктемеде BDSF M={num(a['BDSF_mean'],2)}, SD={num(a['BDSF_sd'],2)}, α={num(a['alpha'])}, ω={num(a['omega_raw_sum'])}. 20 жалаушалы жауапты тек сезімталдық талдауынан шығарғанда n=139, M={num(b['BDSF_mean'],2)}, SD={num(b['BDSF_sd'],2)}, α={num(b['alpha'])}, ω={num(b['omega_raw_sum'])}. Орташа балл өзгерісі {num(b['BDSF_mean']-a['BDSF_mean'],2)} балл ({num(100*(b['BDSF_mean']-a['BDSF_mean'])/a['BDSF_mean'],2)}%). Параллель талдаудағы фактор саны {int(a['parallel_retained'])} → {int(b['parallel_retained'])}. Бұл іріктемелер өзара тәуелді; тәуелсіз топтарға арналған мәнділік тесті қолданылған жоқ. Жалауша жауаптың жарамсыздығын дәлелдемейді.",
      rq3=f"M1–M2 Spearman ρ={num(r12['spearman'],4)}, Kendall τb={num(r12['kendall_tau_b'],4)}; орташа абсолюттік ранг өзгерісі {num(r12['mean_abs_rank_difference'],2)}, ең үлкені {num(r12['max_abs_rank_difference'],1)} орын. M1–M3 үшін ρ={num(r13['spearman'],4)}, орташа өзгеріс {num(r13['mean_abs_rank_difference'],2)}, ең үлкені {num(r13['max_abs_rank_difference'],1)} орын. M2 қолданғанда {cat['n_changed']} қатысушының санаты өзгерді ({num(cat['percent_changed'],2)}%); квадраттық салмақталған κ={num(cat['quadratic_weighted_kappa'],4)}. Жоғары корреляция жеке қатысушының рангы мен шек маңындағы санаты өзгермейді дегенді білдірмейді. M3 үшін сыртқы калибрлеу болмағандықтан санаттау жасалған жоқ. M1 негізгі әдіс ретінде сақталды; оның оңтайлылығы дәлелденген жоқ.",
      rq4=f"Бұрынғы тілдік нәтижелер сақталды: KZ n=81, α={num(lk['alpha_reused'])}; RU n=78, α={num(lr['alpha_reused'])}. Бір факторлы raw-sum ω: KZ {num(lk['omega_raw_sum'])} [{num(lk['ci95_low'])}; {num(lk['ci95_high'])}], RU {num(lr['omega_raw_sum'])} [{num(lr['ci95_low'])}; {num(lr['ci95_high'])}]. Жауап таралуындағы ең үлкен сипаттамалық айырмашылықтар бұрынғыдай Q13, Q03, Q07 сұрақтарында. Pearson төрт факторлы шешімдерінің Tucker коэффициенттері 0.449–0.621. Тіл топтарында ішкі келісімділік жоғары болғанымен, факторлық құрылымның теңдігі дәлелденбеді. Тіл мен қатысушылар құрамы әсерін бұл дизайн ажыратпайды."
    )

def scientific_tables(out,e):
    c=pd.read_csv(out/'EXP01_Factor_Comparison.csv');c=c[c['sample']=='N159']
    modelrows=[[r.correlation,int(r.factors),num(r.offdiag_rmsr),num(r.max_abs_residual),int(r.free_covariance_parameters),int(r.heywood_items)] for _,r in c.iterrows()]
    l=pd.read_csv(out/'EXP01_Factor_Loadings.csv');l=l[(l['sample']=='N159')&(l.correlation=='Pearson')&(l.factors==4)]
    loadrows=[[r['item'],r.theoretical_dimension]+[num(r[f'loading_F{i}']) for i in range(1,5)] for _,r in l.iterrows()]
    qualrows=[[int(s['n']),num(s['BDSF_mean'],2),num(s['BDSF_sd'],2),num(s['alpha']),num(s['omega_raw_sum']),int(s['parallel_retained'])] for s in [e['summary159'],e['summary139']]]
    rankrows=[[r['method1']+'–'+r['method2'],num(r['spearman'],4),num(r['kendall_tau_b'],4),num(r['mean_abs_rank_difference'],2),num(r['max_abs_rank_difference'],1)] for r in e['rank']]
    t=pd.read_csv(out/'EXP03_Category_Transitions.csv');t=t[t.to_method=='M2']
    transrows=[[a]+[int(t[(t.from_category==a)&(t.to_category==b)].n.iloc[0]) for b in ['LOW','MEDIUM','HIGH']] for a in ['LOW','MEDIUM','HIGH']]
    langrows=[[r['language'],r['scale'],num(r['alpha_reused']),num(r['omega_raw_sum']),f"{num(r['ci95_low'])}–{num(r['ci95_high'])}",num(r['omega_polychoric_latent'])] for r in e['language']]
    return [
      ('Бір және төрт факторлы шешімдер',['Корреляция','Фактор','RMSR','max |e|','Параметр','Heywood'],modelrows),
      ('Негізгі төрт факторлы Pearson жүктемелері',['Сұрақ','Өлшем','F1','F2','F3','F4'],loadrows),
      ('Сапа жалаушаларына сезімталдық',['n','M','SD','α','ω','PA'],qualrows),
      ('Рейтинг тұрақтылығы',['Әдістер','Spearman','Kendall','Орташа |Δr|','max |Δr|'],rankrows),
      ('M1 санатынан M2 санатына ауысу',['M1 → M2','LOW','MEDIUM','HIGH'],transrows),
      ('Тілдер бойынша ішкі келісімділік',['Тіл','Шкала','α','ω raw','95% CI','ω latent'],langrows)]

def make_report(here,root,e,checks):
    out=here/'results';n=narrative(e);tabs=scientific_tables(out,e)
    lines=['# BDSF жобасының төрт қосымша эксперименті',
      '## 1 Зерттеудің мақсаты',
      'Қолданыстағы нәтижелерді сақтап, BDSF құрылымын, сапа жалаушаларына, салмақтау әдісіне және сауалнама тіліне сезімталдығын бағалау. Бұл сол деректердің кеңейтілген қайта талдауы, тәуелсіз валидация емес.',
      '## 2 Зерттеу сұрақтары',
      'RQ1 — бір немесе төрт факторға эмпирикалық қолдау; RQ2 — жалаушалы жауаптардың әсері; RQ3 — рейтинг пен санаттардың есептеу әдісіне тұрақтылығы; RQ4 — екі тілдегі психометриялық қасиеттер.',
      '## 3 Деректер мен әдістер',
      '180 бастапқы жауап өзгертілмеді. Негізгі іріктеме — 159 ересек: қазақша 81, орысша 78. Кәмелетке толмаған 21 жазба негізгі және барлық қосымша талдауларға кірмейді. Сапа жалаушасы бар 20 ересек негізгі талдауда қалады; тек RQ2 сезімталдық іріктемесі n=139. Жалаушалар бастапқы 180 жазбада бағаланған қолданыстағы алгоритммен тексерілді. Бос жауаптарды толтыру не қосымша сүзу қолданылған жоқ.',
      'D1=Q01–Q04; D2=Q05–Q07; D3=Q08–Q11; D4=Q12–Q14. Жауап балдары 0/25/50/75/100. Негізгі M1 — төрт өлшемнің тең орташа мәні. M2 — 14 сұрақтың орташа мәні. M3 — бір факторлы Pearson MINRES шешімінің regression баллы: Z R⁻¹λ; Z әр сұрақтың осы 159 қатысушыдағы орташа мәні мен population SD арқылы стандартталған. Фактор таңбасы жалпы баллмен оң бағытта бекітілді. M3 іріктемеге тәуелді, 0–100 шкаласы емес. Ранг 1 — ең жоғары нәтиже, тең баллдар орташа рангпен беріледі. M1/M2 санаты үшін бастапқы екі ондыққа дөңгелектеу және 33.32/66.66 шектері сақталды; M3 санаттауы негізсіз болғандықтан қолданылмады.',
      'EFA — MINRES, төрт фактор үшін oblimin. Негізгі Pearson төрт факторлы шешім мен оның Φ матрицасы бұрынғы CSV-ден алынды. Қалдық rᵢⱼ−(ΛΦΛᵀ)ᵢⱼ; RMSR тек диагональдан тыс 91 жұп бойынша есептелді. Коммуналдылық diag(ΛΦΛᵀ). EFA еркін ковариациялық параметрлері p(k+1)−k(k−1)/2; df=p(p+1)/2−параметрлер. AIC, CFI, TLI, RMSEA не CFA χ² MINRES нәтижесіне жасанды қосылған жоқ. Төрт фактор қосымша салыстыру моделі; PA таңдаған шешім ретінде көрсетілмейді. Факторларды сәйкестендіру тек баған реті мен таңбасын өзгертеді; сұрақ құрамы, баллдар өзгермейді. Модельдер ішіндегі корреляциялар сақталды.',
      'Polychoric — екі сатылы бағалау: шектер эмпирикалық маргиналдық жиіліктерден алынады; әр жұп үшін бивариатты стандартты қалыпты үлестірімдегі ρ шартты ықтималдықты максимизациялау арқылы бағаланады (−0.995…0.995). Бос санаттар жұптағы бақыланған категориялардан алынып тасталады, жауап жазбалары шығарылмайды; бос ұяшықтарға псевдосанау қосылмайды. CDF детерминдік 96 түйінді Gauss–Legendre интегралымен есептеліп, тәуелсіз adaptive quadrature және белгілі ρ=0.6 симуляциясымен тексерілді. Оң анықталмаған матрицада меншікті мәндер 10⁻⁶-дан төмен болса көтеріліп, диагональ 1-ге нормаланады. Бастапқы және қолданылған матрица, түзету шамасы мен жұп диагностикасы сақталған. Бұл түзету нәтижені шектеуі мүмкін [2,3].',
      f"Негізгі Pearson параллель талдауы бұрынғы 500 қалыпты кездейсоқ матрицадан қайта пайдаланылды; n=139 үшін осы тәсіл жаңа есептелді. Polychoric PA: n=159 үшін әр баған жеке араластырылған {e['poly_parallel_permutations']} permutation, маргиналдық үлестірім сақталады, әр қайталауда polychoric қайта бағаланады. 95-процентильдік eigenvalue шегі қолданылады. Seed=20261107. Бұл ordinal permutation null бұрынғы Gaussian null-ға дәл тең емес.",
      f"McDonald's ω бір факторлы congeneric MINRES модельге шартталған. Raw-sum ω = (Σσᵢλᵢ)² / [(Σσᵢλᵢ)² + Σσᵢ²(1−λᵢ²)]. Бұл барлық сұрақтардың тең қосындысына қатысты, негізгі M1 сенімділігімен автоматты тең емес: M1 үшін бөлек салмақталған ω сақталды. Polychoric ω латентті стандартталған жауаптарға қатысты; оны бақыланатын жауаптардың сенімділігіне тікелей теңеуге болмайды [2]. Әр тіл мен шкалаға {e['bootstrap']} person-level bootstrap, raw-sum ω үшін percentile 95% CI; шектер тіл ішіндегі топқа шартталған, көптік салыстыруға түзетілмеген. Конвергенциясы жоқ, тұрақты сұрағы бар немесе Heywood шешімдері жарамсыз деп саналды; жарамдылар 90%-дан аз болса CI берілмейді. Жарамсыз қайталаулар саны кестеде ашық көрсетілген. Polychoric ω CI бағаланған жоқ; ол қосымша нүктелік диагностика.",
      '## 4 Эксперимент 1 Факторлық құрылым',
      'Техникалық аудитте factor_analyzer 0.5.1 кітапханасы loadings_ және structure_ бағандарын қайта реттеп, phi_ матрицасын сол ретпен өзгертпейтіні анықталды. Сондықтан бұрынғы сақталған Φ мен коммуналдылықтарда сәйкессіздік болды. Жаңа қосымшада Φ structure_=ΛΦ теңдігінен қалпына келтіріліп, симметриясы, диагоналі және жүктемелердің өзгермегені тексерілді. Бұрынғы CSV файлдары өзгертілмеді; түзетілген Φ, коммуналдылық және қалдықтар жаңа файлдарда ғана. Түзетудің шамалары EXP01_Legacy_Phi_Audit.csv кестесінде. Негізгі BDSF, α, бір факторлы ω, жүктемелердің өзі және бұрынғы Tucker коэффициенттері өзгермейді. Бастапқы конвейерде бұл ақау әлі бар, сондықтан оның төрт факторлы Φ/коммуналдылық файлдарын түзетілген қосымшамен шатастырмау керек. Кітапхананың ресми коды: https://factor-analyzer.readthedocs.io/en/latest/_modules/factor_analyzer/factor_analyzer.html',
      n['rq1'],table(tabs[0][1],tabs[0][2]),table(tabs[1][1],tabs[1][2]),
      'F1–F4 — бұрынғы негізгі модельдегі D1–D4 бағандарына сәйкес салыстыру белгілері; мазмұндық жарамдылықтың дәлелі емес. Бір фактордың меншікті мәні мен жалпы ортақ дисперсиясы бір ұғым емес; кестедегі common_variance_share_of_total — коммуналдылықтардың орташа мәні. Қиғаш факторлардың дисперсия үлестерін тәуелсіз бөліктер сияқты қосып түсіндіруге болмайды.',
      '## 5 Эксперимент 2 Деректер сапасы',n['rq2'],table(tabs[2][1],tabs[2][2]),
      'STRAIGHTLINE_ALL_14, DUPLICATE_RESPONSE_PATTERN және DUPLICATE_TIMESTAMP жалаушалары бір-бірімен қабаттасуы мүмкін; себептер саны қатысушылар санына тең болмауы мүмкін. n=139 — сапасы дәлелденген «дұрыс» іріктеме емес, алдын ала анықталған қосымша сүзудің сезімталдық сценарийі.',
      '## 6 Эксперимент 3 Салмақтау',n['rq3'],table(tabs[3][1],tabs[3][2]),table(tabs[4][1],tabs[4][2]),
      'M1→M3 санат ауысуы — NOT_APPLICABLE. M3-ті min–max арқылы созу немесе percentile-ді қауіпсіздік баллы деп атау негізделмегендіктен қолданылмады. Ең көп ранг өзгерткен жазбалар тек P001… техникалық талдау кодымен берілген; ат, бастапқы ID, уақыт белгісі, жас, жыныс экспортталмайды.',
      '## 7 Эксперимент 4 Екі тіл',n['rq4'],table(tabs[5][1],tabs[5][2]),
      'Жеке өлшемдердің α, ω және CI кестеде берілді; интервалдардың қабаттасуы формалды теңдік тесті емес. Pearson–polychoric және KZ–RU Tucker коэффициенттері EXP04_Factor_Stability.csv файлында. Фактор бағандары алдын ала белгіленген reference бойынша сәйкестендірілді; бұл процедура олардың айырмашылығын жоққа шығармайды.',
      '## 8 Біріктірілген талқылау',
      'Нәтижелер жалпы индекс үшін ішкі келісімділік дәлелін береді, бірақ жоғары α/ω төрт дербес құрылымның валидтілігін дәлелдемейді. RQ1–RQ4 осы іріктемеге қатысты сипаттамалық және сезімталдық дәлелдерімен жауапталды. Салмақтау мен сапа сүзуі туралы шешімдер негізгі есептеуді үнсіз алмастырмайды. Қандай да бір әдісті жақсырақ нәтиже көрсеткені үшін негізгі әдіске ауыстыру жасалған жоқ.',
      '## 9 Шектеулер',
      'Іріктеме ықтималдық емес, тіл топтары шағын және рандомизацияланбаған. Қайта талдау тәуелсіз валидация емес. Polychoric латентті қалыптылықты болжайды; сирек ұяшықтар мен матрицаны түзету тұрақтылықты төмендетуі мүмкін. Төрт факторлы шешімде Heywood жағдайлары болса, олар модельдің шектеуі ретінде сақталады. ω бір факторлық модель мен қателердің корреляциясыздығына шартталған; жоғары мәні модель fit дәлелі емес. Ағымдағы bootstrap ω белгісіздігін бағалайды, фактор санын және 14×4 жүктемелердің толық bootstrap тұрақтылығын бағаламайды. Жалпыға қолданылатын қауіпсіздік шектерінің сыртқы валидациясы жоқ. Жеке келісім рәсімдері туралы жаңа мәлімдеме жасалмады.',
      '## 10 Ғылыми қорытынды',
      'Қолдау: сұрақтарда жалпы ортақ өзгергіштік және ішкі келісімділік бар. Қолдау жеткіліксіз: төрт дербес өлшемнің толық бөлінуі, тілдік measurement invariance, 25% салмақтардың эмпирикалық оңтайлылығы. Анық емес: тіл айырмашылығының аудармадан не қатысушы құрамынан болуы және қауіпсіздік санаттарының сыртқы критериймен сәйкестігі. Тәуелсіз үлкен іріктемеде ordinal CFA, тілдер бойынша configural/metric/scalar invariance, DIF, test–retest және объективті қауіпсіздік көрсеткіштерімен сыртқы валидация қажет.',
      '## 11 Мақаладағы өзгерістер',
      'Бөлек ARTICLE_Updated_Four_Experiments.docx көшірмесінде әдістер толықтырылды, нәтижелерге 7–12 кестелер және төрт эксперимент енгізілді, үш тілдегі аңдатпаға негізгі жаңа қорытындылар қосылды, талқылау, шектеулер және қорытынды жаңартылды. Бастапқы мақала, 00–15 кезеңдерінің файлдары және негізгі конвейер өзгертілмеді. Бұл жаңа каталог GitHub-қа автоматты жарияланбайды.',
      '## Әдіс дереккөздері',
      '[1] FactorAnalyzer MINRES және oblimin ресми құжаттамасы: https://factor-analyzer.readthedocs.io/en/latest/factor_analyzer.html',
      '[2] Revelle psych polychoric құжаттамасы, ordinal reliability шектеулері: https://www.personality-project.org/r/psych/help/tetrachor.html',
      '[3] Revelle psych correlation matrix smoothing: https://www.personality-project.org/r/psych/help/cor.smooth.html',
      '[4] Revelle omega tutorial: https://www.personality-project.org/r/tutorials/HowTo/omega.tutorial/omega.pdf']
    (out/'FINAL_Research_Report_KK.md').write_text('\n\n'.join(lines)+'\n',encoding='utf-8')
    schema=['# Нәтижелерді оқу және қайта іске қосу','Жоба түбірінен: `16_extended_experiments\\RUN_EXPERIMENTS.cmd`. Бұрынғы RUN_REVIEW.cmd өзгермейді. Жаңа талдау алдында бастапқы кезеңдердің дайын CSV нәтижелері болуы керек.','Есептеудің қайталаулары cache ішінде кіріс SHA256, methods.py, seed және параметрлер бойынша тексеріліп қайта пайдаланылады. Негізгі нәтижелер қайта жасалмайды. Негізгі кодтағы логика өзгерсе cache жойып емес, жаңа fingerprint-пен жаңартылуы тиіс. Қате болса nonzero exit және VALIDATION_REPORT.json=FAIL; PASS болмаған нәтижелер аяқталған пакет емес.','ω — raw-score sum үшін; ω latent — polychoric модель үшін. ci95 — raw ω bootstrap percentile; бос ұяшық — бағалау жарамсыз не қолданылмайды. PA — параллель талдауда қалған компоненттер саны; RMSR — off-diagonal residual root mean square. n — қатысушылар/ұяшық саны, контекстке сәйкес. Барлық айырма N139−N159 немесе method2−method1, кесте атауына сәйкес.','M3 санаттары әдейі берілмейді; EXP03_Category_Transitions.csv файлындағы NOT_APPLICABLE жолы есептеу қатесі емес. Негізгі BDSF M1 ретінде сақталған.','## Файлдар және бағандар']
    for p in sorted(out.glob('*.csv')):
        t=pd.read_csv(p);schema.append('### '+p.name+'\n\n'+', '.join('`'+c+'`' for c in t.columns)+f'\n\nЖол саны: {len(t)}.')
    (here/'README_KK.md').write_text('\n\n'.join(schema),encoding='utf-8')
    definitions={
      'sample':'Іріктеме: N159 негізгі, N139 сапаға сезімталдық, KZ/RU тіл топтары.',
      'correlation':'Қолданылған матрица: Pearson немесе Polychoric.',
      'source':'Есептің шыққан жері: бұрынғы нәтижеден алынған не жаңа есептеу.',
      'factors':'Модельде бағаланған факторлар саны; фактор санын таңдау қорытындысы емес.',
      'offdiag_rmsr':'91 диагональдан тыс корреляция қалдығы квадраттарының орташа мәнінің түбірі.',
      'max_abs_residual':'Диагональдан тыс ең үлкен абсолюттік қалдық корреляция.',
      'residual_pairs_gt_005':'Абсолюттік қалдық корреляциясы 0.05-тен үлкен жұп саны.',
      'mean_communality':'Коммуналдылықтардың орташа мәні.',
      'sum_common_variance':'Коммуналдылықтар қосындысы; стандартталған сұрақтар бойынша ортақ дисперсия.',
      'common_variance_share_of_total':'Ортақ дисперсия / барлық стандартталған сұрақ дисперсиясы; mean_communality-мен бірдей.',
      'heywood_items':'Коммуналдылығы 0-ден төмен не 1-ден кем емес сұрақтар саны.',
      'efa_df':'EFA ковариациялық модель еркіндік дәрежесі; χ² fit тесті есептелмеген.',
      'free_covariance_parameters':'EFA ковариациялық параметр саны, ротациялық анықталмаушылық ескерілген.',
      'max_abs_factor_correlation':'Факторлар арасындағы ең үлкен абсолюттік Φ корреляциясы.',
      'item':'Q01–Q14 сұрағы; item1/item2 — сұрақ жұбы.',
      'theoretical_dimension':'Алдын ала бекітілген сұрақтың D1–D4 құрамы.',
      'communality':'diag(ΛΦΛᵀ); қиғаш айналдыруда жай квадраттар қосындысы емес.',
      'uniqueness':'1 − коммуналдылық.',
      'largest_loading_factor':'Абсолюттік жүктемесі ең үлкен баған; сұрақтың жаңа өлшемі емес.',
      'second_abs_loading':'Екінші ең үлкен абсолюттік жүктеме.',
      'factor1':'Φ жолының факторы.','factor2':'Φ бағанының факторы.','correlation_value':'Факторлар корреляциясының Φ элементі.',
      'residual':'Бақыланатын корреляция − модельдің ортақ корреляциясы.',
      'n':'Жол контекстіндегі қатысушылар саны; санат/жалауша кестесінде сол топтың саны.',
      'BDSF_mean':'Бастапқы M1 орташа мәні.','BDSF_sd':'M1 таңдама стандарттық ауытқуы, ddof=1.',
      'BDSF_ci_low':'M1 орташа мәнінің t-based 95% CI төменгі шегі.','BDSF_ci_high':'M1 орташа мәнінің t-based 95% CI жоғарғы шегі.',
      'alpha':'Шикі Cronbach α.','KMO':'Pearson матрицасына негізделген жалпы KMO.',
      'Bartlett_chi2':'Pearson Bartlett сфералық тест статистикасы.','Bartlett_df':'Bartlett тестінің еркіндік дәрежесі.','Bartlett_p':'Bartlett тестінің түзетілмеген p мәні.',
      'parallel_retained':'Параллель талдауда 95% кездейсоқ шектен асатын компонент саны.',
      'omega_raw_sum':'Бақыланатын сұрақтардың тең қосындысының бір факторлы модельге шартталған ω.',
      'omega_M1_weighted':'Бастапқы 25% өлшемдік M1 салмақтарына арналған бір факторлы ω.',
      'role':'PRIMARY не SENSITIVITY_ONLY.',
      'case':'P001… дербестендірілмеген талдау коды; бастапқы ID экспортталмайды.',
      'M1':'Бастапқы өлшемдердің тең орташа мәні, 0–100.','M2':'14 сұрақтың тең орташа мәні, 0–100.','M3':'Бір факторлы regression score, Z R⁻¹λ; 0–100 емес.',
      'method1':'Салыстырудың бірінші әдісі.','method2':'Салыстырудың екінші әдісі.',
      'spearman':'Тең рангтарды ескеретін Spearman корреляциясы.','kendall_tau_b':'Тең мәндерді ескеретін Kendall τb.',
      'mean_abs_rank_difference':'159 адам бойынша |rank2−rank1| орташа мәні.','max_abs_rank_difference':'Ең үлкен |rank2−rank1|.','n_rank_changed':'Рангы өзгерген адамдар саны, тең рангтарды қоса.',
      'from_method':'Бастапқы әдіс.','to_method':'Салыстырылатын әдіс.','from_category':'Бастапқы санат.','to_category':'Соңғы санат.',
      'n_changed':'M1→M2 санаты өзгерген адамдар саны.','percent_changed':'n_changed / 159 × 100.',
      'quadratic_weighted_kappa':'LOW/MEDIUM/HIGH ретімен квадраттық салмақталған Cohen κ.',
      'language':'KZ — қазақша, RU — орысша.','scale':'ALL14 — 14 сұрақ жиынтығы; D1–D4 — жеке мазмұндық бағыттар.',
      'alpha_reused':'15-кезеңнен тексеріліп қайта алынған α.','alpha_ci_low_reused':'Бұрынғы 2000 bootstrap α интервалының төменгі шегі.','alpha_ci_high_reused':'Бұрынғы 2000 bootstrap α интервалының жоғарғы шегі.',
      'ci95_low':'Raw-sum ω percentile 95% CI төменгі шегі.','ci95_high':'Raw-sum ω percentile 95% CI жоғарғы шегі.',
      'bootstrap_requested':'Жоспарланған қайта іріктеу саны.','bootstrap_valid':'Жарамды модель алынған қайта іріктеу саны.','bootstrap_failed':'Тұрақты сұрақ, конвергенция не Heywood себебімен жарамсыз қайта іріктеу саны.',
      'ci_status':'OK: кемінде 90% жарамды; әйтпесе CI берілмейді.',
      'omega_polychoric_latent':'Латентті стандартталған жауаптарға арналған polychoric ω, raw-score сенімділігі емес.',
      'polychoric_omega_ci':'NOT_ESTIMATED: polychoric ω bootstrap интервалы бағаланбаған.',
      'polychoric_onefactor_rmsr':'Тиісті тіл/шкала polychoric бір факторлы моделінің off-diagonal RMSR.',
      'comparison':'Қай матрица не тілдер салыстырылғаны.','factor':'Сәйкестендірілген фактор бағаны.','tucker_congruence':'Баған жүктемелерінің косинусы; ротация реті/таңбасы сәйкестендірілген.',
      'min_eigenvalue_raw':'Түзетуге дейінгі ең кіші eigenvalue.','smoothed':'Матрицаға eigenvalue clipping қолданылған-қолданылмағаны.','max_abs_adjustment':'Матрицаны түзетудің ең үлкен абсолюттік элемент айырмасы.',
      'metric':'Салыстырылатын көрсеткіш.','N159':'Негізгі мән.','N139':'Сезімталдық мән.','absolute_change':'N139−N159.','relative_change_percent':'100×(N139−N159)/N159.',
      'flag':'Алгоритмдік сапа белгісі.','flag_combination':'Бір адамдағы белгілер жиынтығы; NONE — жоқ.','scope':'Белгіні есептеу популяциясы.',
      'component':'Кему ретімен eigenvalue нөмірі.','observed':'Бақыланған eigenvalue.','null95':'Кездейсоқ null үлестірімінің 95-процентилі.','retain':'observed > null95.',
      'rho':'Екі сатылы polychoric корреляция бағасы.','zero_cells':'Жұп кестедегі нөлдік жиілік ұяшықтары.','cells':'Жұп кестедегі ұяшық саны.','boundary':'Баға ±0.995 шекарасына жақын ба.','negative_loglik':'Шартты теріс log likelihood, шектер бекітілген.',
      'regression_coefficient':'R⁻¹λ құрамындағы M3 коэффициенті.','training_mean':'M3 стандарттауға пайдаланған сұрақ орташа мәні.','training_sd_population':'M3 стандарттау SD, ddof=0.',
      'one_factor_N159':'Негізгі бір факторлы Pearson жүктеме.','one_factor_N139':'Сезімталдық бір факторлы Pearson жүктеме.',
      'task':'Қайта пайдаланылған есептеу.','action':'Бұрынғы нәтиже не fingerprint-пен тексерілген cache.',
      'stage':'1–16 кезең нөмірі.','status':'Есептің не тексерудің күйі.','check':'Тексерілген талап.','detail':'Тексеру туралы қосымша мәлімет.'}
    dictionary=['# CSV бағандарының сөздігі']
    for p in sorted(out.glob('*.csv')):
        dictionary.append('## '+p.name)
        for c in pd.read_csv(p,nrows=0).columns:
            description=definitions.get(c)
            if description is None:
                if c.startswith('loading_F'):description='MINRES pattern жүктемесі; баған сәйкестендірілген фактор.'
                elif c.endswith('_rank'):description='Ранг не екі ранг айырмасы; ең жоғары баллдың рангы 1, теңдер орташа рангпен.'
                elif c.endswith('_difference'):description='Тиісті фактордың N139−N159 жүктеме айырмасы, реті/таңбасы сәйкестендірілген.'
                elif c in ['item1','item2']:description='Корреляциялық жұптағы Q01–Q14 сұрағы.'
                elif c=='category':description='LOW/MEDIUM/HIGH, бастапқы дөңгелектеу мен шектер бойынша.'
                elif c=='max_abs_change':description='M2−M1 және M3−M1 ранг айырмаларының үлкен абсолюттік мәні.'
                else:description='Тексеру не әдістің атауы; мән контексті README_KK.md және ғылыми есепте берілген.'
            dictionary.append(f'- `{c}`: {description}')
    (here/'DATA_DICTIONARY_KK.md').write_text('\n\n'.join(dictionary),encoding='utf-8')

def update_article(here,root,e):
    from docx import Document
    from docx.shared import Pt,RGBColor,Cm
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    source=root/'12_article_results/manuscript_20261106/BDSF_Kazakhstan_New_Article_2026.docx'
    original_hash=hashlib.sha256(source.read_bytes()).hexdigest();doc=Document(source);n=narrative(e);out=here/'results'
    def find(text):return next(p for p in doc.paragraphs if p.text==text)
    def before(anchor,text,style=None):
        p=anchor.insert_paragraph_before(text,style=style);p.paragraph_format.widow_control=True;return p
    methods=find('Нәтижелер')
    before(methods,'Қосымша төрт эксперименттің әдістері','Heading 2')
    before(methods,'Техникалық аудит factor_analyzer 0.5.1 нұсқасында ротациядан кейін жүктемелер мен Φ бағандарының реті сәйкес келмейтінін көрсетті [12]. Қосымша талдауда Φ structure_=ΛΦ теңдігінен қалпына келтірілді; жүктемелер, негізгі баллдар және бастапқы файлдар өзгертілмеді. Осы мақаладағы жаңа қалдықтар мен коммуналдылықтар түзетілген Φ матрицасына негізделеді. Айырмалар EXP01_Legacy_Phi_Audit.csv файлында берілген.')
    before(methods,'Негізгі 159 ересек өзгеріссіз сақталды. Қосымша талдау төрт сұраққа бағытталды: бір және төрт факторлы құрылым; сапа жалаушаларын сүзуге сезімталдық; салмақтауға рейтинг пен санаттардың тұрақтылығы; екі тілдегі психометриялық қасиеттер. 20 жалаушалы жауапты алып тастаған n=139 жиынтығы тек сезімталдық сценарийі болды. Талдау тәуелсіз валидация емес.')
    before(methods,'Pearson MINRES бір факторлы шешімі төрт факторлы oblimin шешімімен диагональдан тыс RMSR, қалдық корреляциялар, коммуналдылық және модель күрделілігі бойынша салыстырылды. Қосымша polychoric корреляциялар эмпирикалық маргиналдық шектер мен бивариатты қалыпты likelihood арқылы екі сатылы бағаланды [12,13]. Қажет болса оң анықталмаған матрицаларға eigenvalue clipping қарастырылды; негізгі және екі тілдік матрицада мұндай түзету қажет болған жоқ. Бастапқы және қолданылған матрицалар сақталды. Polychoric параллель талдауы маргиналдық жиіліктерді сақтайтын 200 permutation негізінде орындалды.')
    before(methods,'M1=0.25(D1+D2+D3+D4); M2=14 сұрақтың орташа мәні; M3=Z R⁻¹λ, бір факторлы MINRES regression score. M3 сол іріктеменің стандартталған жауаптарынан есептелді және сыртқы калибрлеу болмағандықтан 0–100 қауіпсіздік санаттарына түрлендірілмеді. Рейтинг үшін Spearman ρ, Kendall τb, орташа және ең үлкен абсолюттік ранг айырмасы есептелді. Тең баллдарға орташа ранг берілді. M1–M2 санат ауысуы бастапқы шектермен, екі ондыққа дөңгелектеп, квадраттық салмақталған Cohen κ арқылы бағаланды.')
    before(methods,'McDonald ω бір факторлы congeneric модель бойынша есептелді [14]. Бақыланатын сұрақтардың тең қосындысы үшін ω=(Σσᵢλᵢ)²/[(Σσᵢλᵢ)²+Σσᵢ²(1−λᵢ²)]. Әр тілдегі ALL14 және D1–D4 үшін 1000 person-level bootstrap percentile 95% CI есептелді. Жарамсыз және Heywood қайталаулары саналып, 90%-дан аз жарамды қайталау болса интервал берілмеді. Polychoric ω — латентті стандартталған жауаптардың қосымша нүктелік бағасы; ол бақыланатын балл сенімділігімен теңестірілмейді. Негізгі 25% салмақтау өзгертілмеді.')
    results=find('Талқылау')
    for i,key in enumerate(['rq1','rq2','rq3','rq4'],1):
        before(results,f'Қосымша эксперимент {i}','Heading 2');before(results,n[key])
    def insert_table(title,headers,rows,index):
        before(results,f'{index} кесте {title}')
        t=doc.add_table(rows=1,cols=len(headers));t.style='Table Grid';t.autofit=False
        for j,h in enumerate(headers):t.rows[0].cells[j].text=h
        for row in rows:
            cells=t.add_row().cells
            for j,v in enumerate(row):cells[j].text=str(v)
        for rowidx,row in enumerate(t.rows):
            for cell in row.cells:
                cell.width=Cm(17/len(headers))
                for p in cell.paragraphs:
                    p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.space_after=Pt(3);p.paragraph_format.space_before=Pt(3)
                    p.paragraph_format.keep_with_next=(rowidx<2 or len(t.rows)-3<=rowidx<len(t.rows)-1)
                    for run in p.runs:run.font.name='Times New Roman';run.font.size=Pt(9);run.font.bold=rowidx==0;run.font.color.rgb=RGBColor(0,0,0)
                if rowidx==0:
                    shading=OxmlElement('w:shd');shading.set(qn('w:fill'),'E8EEF4');cell._tc.get_or_add_tcPr().append(shading)
            trpr=row._tr.get_or_add_trPr();no=OxmlElement('w:cantSplit');trpr.append(no)
        repeat=OxmlElement('w:tblHeader');t.rows[0]._tr.get_or_add_trPr().append(repeat)
        borders=OxmlElement('w:tblBorders')
        for edge in ['top','left','bottom','right','insideH','insideV']:
            el=OxmlElement('w:'+edge);el.set(qn('w:val'),'single');el.set(qn('w:sz'),'4');el.set(qn('w:color'),'D9D9D9');borders.append(el)
        t._tbl.tblPr.append(borders);results._p.addprevious(t._tbl)
        before(results,'')
    for i,(title,heads,rows) in enumerate(scientific_tables(out,e),7):insert_table(title,heads,rows,i)
    before(results,'Ескерту: F1–F4 бағандары бұрынғы D1–D4 сәйкестендіруіне негізделген; фактор атауы теориялық валидтіліктің дәлелі емес. ω raw — бақыланатын жауаптардың тең қосындысы, ω latent — polychoric модельдегі латентті стандартталған жауаптар. CI тек ω raw үшін. PA — параллель талдауда қолдау алған компоненттер саны. M3 санаттау жасалмағандықтан 11-кестеде тек M1→M2 көрсетілген.')
    for p in doc.paragraphs:
        if 'Барлық талдау бір ғана 159 ересек іріктемеде жүргізілді' in p.text:
            p.text=p.text.replace('Барлық талдау бір ғана 159 ересек іріктемеде жүргізілді','Негізгі талдаулар 159 ересек іріктемеде жүргізілді')
        if 'барлық 14 жауапты бірден орташалау қолданылмайды' in p.text:
            p.text=p.text.replace('барлық 14 жауапты бірден орташалау қолданылмайды','негізгі M1 формуласында барлық 14 жауапты бірден орташалау қолданылмайды')
        if 'мұның бағалауға ықпалы бөлек іріктеме арқылы тексерілген жоқ' in p.text:
            p.text=p.text.replace('мұның бағалауға ықпалы бөлек іріктеме арқылы тексерілген жоқ','мұның бағалауға ықпалы n=139 қосымша сезімталдық сценарийінде тексерілді, бірақ жалаушаның өзі жауаптың жарамсыздығын дәлелдемейді')
        if p.text.startswith('Талдаудың жалғыз шығару шарты'):
            p.text=p.text.replace('Талдаудың жалғыз шығару шарты','Негізгі талдаудың жалғыз шығару шарты')+' Қосымша сапаға сезімталдық талдауында ғана жалаушасыз 139 ересек бөлек қарастырылды.'
        elif p.text.startswith('Тең салмақтау әр мазмұндық бағыттың'):
            p.text='Тең салмақтау әр мазмұндық бағыттың үлесін бірдей ұстайды; ол хаттамада бекітілген авторлық шешім [7], эмпирикалық оңтайлылық кепілі емес. '+n['rq3']
        elif p.text.startswith('Бес деңгейлі жауаптар Pearson корреляцияларымен'):
            p.text='Бес деңгейлі жауаптар Pearson және polychoric корреляцияларымен талданды. Polychoric бағалар сирек ұяшықтарға, латентті қалыптылық болжамына және қажет болған матрица түзетуіне сезімтал. ω бір факторлы модельге шартталған; жоғары сенімділік құрылымдық fit дәлелі емес. Тіл топтарындағы n=81/78 төрт факторлы шешімді тұрақсыз етуі мүмкін. Тәуелсіз ordinal CFA, тілдік өлшеу инварианттылығы, DIF, test–retest және объективті критерийлермен валидтілік тексерілген жоқ.'
        elif p.text.startswith('Аңдатпа.'):
            p.text+=' Қосымша төрт эксперимент бір/төрт факторлы шешімдерді, n=159/139 сезімталдығын, M1–M3 рейтингтерін және тілдердегі ω мен polychoric құрылымын салыстырды. '+n['rq2'].split('Бұл іріктемелер')[0]+f" M1–M2 бойынша {e['categories']['n_changed']} қатысушы санатын өзгертті."
        elif p.text.startswith('Abstract.'):
            p.text+=f" Supplementary analyses compared one- and four-factor Pearson/polychoric solutions, the primary sample with 139 unflagged adults, three scoring methods, and language-specific omega. In the sensitivity sample, the mean was {num(e['summary139']['BDSF_mean'],2)}, alpha {num(e['summary139']['alpha'])}, and raw-sum omega {num(e['summary139']['omega_raw_sum'])}. Equal item weighting changed the categories of {e['categories']['n_changed']} participants. The primary sample and scoring method were retained; these analyses do not establish measurement invariance."
        elif p.text.startswith('Аннотация.'):
            p.text+=f" Дополнительные анализы сопоставили одно- и четырёхфакторные решения Pearson/polychoric, выборки 159/139, методы M1–M3 и ω по языкам. В выборке чувствительности средний балл равен {num(e['summary139']['BDSF_mean'],2)}, α={num(e['summary139']['alpha'])}, ω={num(e['summary139']['omega_raw_sum'])}. Равные веса вопросов изменили категорию у {e['categories']['n_changed']} участников. Основная выборка и формула сохранены; инвариантность измерения не доказана."
    limits=find('Шектеулер');before(limits,n['rq2']);before(limits,'Polychoric нәтижелері мен ω қосымша модельдік дәлел береді. Оларды төрт өлшемді немесе тілдік теңдікті растау үшін іріктеп түсіндіруге болмайды. Сезімталдық талдауы зерттеу дизайнындағы шектеулерді жоймайды.')
    availability=find('Деректер мен кодтың қолжетімділігі');before(availability,'Қосымша төрт эксперимент негізгі индексті өзгертпей оның құрылымдық және алгоритмдік шектеулерін сандық түрде сипаттады. '+n['rq1'])
    refs=find('Әдебиеттер');before(refs,'Қосымша код, кестелер және 16 кезеңдік тексеру 16_extended_experiments каталогында сақталған. Оның жергілікті нұсқасы негізгі жарияланған репозиторийден бөлек; осы қосымшаның GitHub-та жарияланғаны туралы мәлімдеме жасалмайды.')
    last=next(p for p in doc.paragraphs if p.text.startswith('[11]'))
    # Add references after [11], without moving or renumbering the existing sources.
    for text in reversed(['[12] FactorAnalyzer. MINRES and oblimin source and API documentation. https://factor-analyzer.readthedocs.io/en/latest/_modules/factor_analyzer/factor_analyzer.html','[13] Revelle W. psych documentation: polychoric correlations and smoothing. https://www.personality-project.org/r/psych/help/tetrachor.html','[14] Revelle W. Using R and the psych package to find omega. https://www.personality-project.org/r/tutorials/HowTo/omega.tutorial/omega.pdf']):
        p=doc.add_paragraph(text);last._p.addnext(p._p)
    target=out/'ARTICLE_Updated_Four_Experiments.docx';doc.save(target)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==original_hash
    assert len(doc.tables)==12
    (out/'ARTICLE_CHANGE_MAP.md').write_text('Бастапқы мақала өзгертілмеді.\n\nӘдістер: 4 қосымша әдіс абзацы.\nНәтижелер: 4 эксперимент және 7–12 кестелер.\nАңдатпа: қазақша, орысша, ағылшынша толықтырылды.\nТалқылау: сапа және салмақтау сезімталдығы.\nШектеулер: polychoric/omega және measurement invariance.\nҚорытынды: бір/төрт фактор дәлелдері.\nӘдебиеттер: 12–14 әдіс дереккөздері.\n',encoding='utf-8')
