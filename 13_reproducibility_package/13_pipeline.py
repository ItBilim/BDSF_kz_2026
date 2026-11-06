"""Rebuild BDSF results from two raw CSV exports. Fail before publish on any check."""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json,os,platform,re,shutil,subprocess,sys,tempfile,uuid
from datetime import datetime,timezone,timedelta
import numpy as np
import pandas as pd
import openpyxl

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('statistics_bdsf',HERE/'13_statistics.py');st=importlib.util.module_from_spec(sp);sp.loader.exec_module(st)
sp=importlib.util.spec_from_file_location('verify_bdsf',HERE/'13_01_reproduce_analysis.py');core=importlib.util.module_from_spec(sp);sp.loader.exec_module(core)
LEVELS={'NEVER':0,'RARELY':25,'SOMETIMES':50,'OFTEN':75,'ALWAYS':100}
DEMOS={
'KZ':{'headers':['Отметка времени','Жынысыңыз','Сіз студентсіз бе?','Жас аралығыңыз'],'gender':{'Ер':'MALE','Әйел':'FEMALE'},'student':{'Иә, студентпін':'STUDENT'},'age':{'18 жасқа дейін':'LT18','18–20 жас':'18_20','21–23 жас':'21_23','27 жас және одан жоғары':'GE27'}},
'RU':{'headers':['Отметка времени','1. Ваш пол','Являетесь ли Вы студентом?','3. Укажите Вашу возрастную категорию'],'gender':{'Мужской':'MALE','Женский':'FEMALE'},'student':{'Да, я студент(ка)':'STUDENT'},'age':{'Младше 18 лет':'LT18','18–20 лет':'18_20','21–23 года':'21_23','27 лет и старше':'GE27'}}}
ANSWERS={}
for words in [('Ешқашан','Сирек','Кейде','Жиі','Әрқашан'),('Никогда','Редко','Иногда','Часто','Всегда'),('Мүлдем көз жеткізбеймін','Сирек көз жеткіземін','Кейде көз жеткіземін','Жиі көз жеткіземін','Әрқашан толық көз жеткіземін'),('Никогда не проверяю','Редко проверяю','Иногда проверяю','Часто проверяю','Всегда тщательно проверяю')]: ANSWERS.update(dict(zip(words,LEVELS)))

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest().upper()
def write(path,text):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
def jsonwrite(path,obj):write(path,json.dumps(obj,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
def run(cmd):subprocess.run([str(x) for x in cmd],check=True)
def checked_map(value,mapping,where):
    key=str(value).strip()
    if key not in mapping:raise ValueError(f'Unknown category at {where}: {value!r}')
    return mapping[key]
def longstring(values):
    best=now=1
    for a,b in zip(values,values[1:]):now=now+1 if a==b else 1;best=max(best,now)
    return best

def harmonize(project,qmap,build):
    records=[];audits=[];cats=[];missing=[];inventory=[]
    for lang,name in [('KZ','Kazakh'),('RU','Russian')]:
        path=next((project/'00_sources').glob(f'*_{name}_Responses.csv'))
        df=pd.read_csv(path,dtype=str,keep_default_na=False);cfg=DEMOS[lang]
        question_col='kazakh_question' if lang=='KZ' else 'russian_question'
        expected=cfg['headers']+qmap[question_col].tolist()
        if [str(x).strip() for x in df.columns]!=expected:raise ValueError(f'{lang}: missing, reordered, duplicated or unknown input columns')
        if df.eq('').any().any():raise ValueError(f'{lang}: missing response; no silent row exclusion allowed')
        inventory.append(dict(file_name=path.name,sha256=sha(path),rows=len(df),columns=len(df.columns)))
        for pos,col in enumerate(df.columns):
            audits.append(dict(language=lang,column_position=pos+1,column_name=col,dtype='string',rows=len(df),unique_values=df[col].nunique()))
            missing.append(dict(language=lang,column_name=col,missing_count=int(df[col].str.strip().eq('').sum())))
        for j in range(14):
            for answer,n in df.iloc[:,j+4].value_counts().items():cats.append(dict(language=lang,item_code=f'Q{j+1:02d}',response_category=answer,count=n))
        for idx,row in df.iterrows():
            v=row.tolist();rid=f'{lang}{idx+1:03d}'
            if not v[0].endswith(' GMT+5'):raise ValueError(f'{rid}: unrecognized timestamp timezone')
            dt=datetime.strptime(v[0].removesuffix(' GMT+5'),'%Y/%m/%d %I:%M:%S %p').replace(tzinfo=timezone(timedelta(hours=5)))
            r=dict(respondent_id=rid,source_language=lang,source_file=path.name,source_row=idx+2,timestamp_raw=v[0],timestamp_iso=dt.isoformat(),gender_raw=v[1],gender_code=checked_map(v[1],cfg['gender'],rid),student_status_raw=v[2],student_status_code=checked_map(v[2],cfg['student'],rid),age_raw=v[3],age_group_code=checked_map(v[3],cfg['age'],rid))
            for j in range(14):r[f'Q{j+1:02d}_raw']=v[j+4];r[f'Q{j+1:02d}_level']=checked_map(v[j+4],ANSWERS,f'{rid}/Q{j+1:02d}')
            r.update(harmonization_status='PASS',harmonization_issues='');records.append(r)
    raw=pd.DataFrame(records)
    if len(raw)!=180 or raw.respondent_id.nunique()!=180:raise ValueError('This registered study requires exactly 180 unique input records')
    frame={
    '01_audit/01_01_file_inventory.csv':pd.DataFrame(inventory),
    '01_audit/01_02_dataset_structure.csv':pd.DataFrame(audits),
    '01_audit/01_04_missing_values_report.csv':pd.DataFrame(missing),
    '01_audit/01_05_kz_ru_question_alignment.csv':qmap.copy(),
    '01_audit/01_06_response_category_alignment.csv':pd.DataFrame(cats),
    '03_harmonization/03_01_BDSF_Kazakhstan_Harmonized_Dataset.csv':raw,
    '03_harmonization/03_03_BDSF_Kazakhstan_Bilingual_Question_Map.csv':qmap.copy()}
    levels=[f'Q{i:02d}_level' for i in range(1,15)];clean=raw.copy();dup=clean.duplicated(levels,keep=False);ts=clean.timestamp_iso.duplicated(keep=False)
    clean['analysis_inclusion']=np.where(clean.age_group_code.eq('LT18'),'EXCLUDE','INCLUDE')
    clean['exclusion_reasons']=np.where(clean.age_group_code.eq('LT18'),'UNDER_18_NO_PARENTAL_PERMISSION','')
    clean['unique_response_levels']=clean[levels].nunique(axis=1);clean['longstring_max']=clean[levels].apply(lambda r:longstring(r.tolist()),axis=1)
    clean['quality_flags']=[';'.join((["STRAIGHTLINE_ALL_14"] if r.unique_response_levels==1 else [])+(['DUPLICATE_RESPONSE_PATTERN'] if dup.iloc[i] else [])+(['DUPLICATE_TIMESTAMP'] if ts.iloc[i] else [])) for i,r in clean.iterrows()]
    # Flags are evaluated on the complete collected data and never exclude an adult.
    excluded=clean.loc[clean.age_group_code.eq('LT18'),['respondent_id','source_language','source_row','exclusion_reasons']]
    adults=clean.loc[clean.analysis_inclusion.eq('INCLUDE')].copy().reset_index(drop=True)
    if len(adults)!=159 or len(excluded)!=21:raise ValueError('Registered population must be 159 adults plus 21 LT18 records')
    frame['01_audit/01_03_duplicate_report.csv']=clean.loc[dup|ts,['respondent_id','source_language','quality_flags']]
    frame['04_cleaning/04_01_BDSF_Kazakhstan_Clean_Dataset.csv']=adults.copy()
    frame['04_cleaning/04_03_BDSF_Kazakhstan_Exclusion_Log.csv']=excluded
    for q in range(1,15):adults[f'Q{q:02d}_score']=adults[f'Q{q:02d}_level'].map(LEVELS)
    for d,cols in st.DIMS.items():adults[d+'_score']=adults[cols].mean(axis=1)
    adults['BDSF_score']=adults[[d+'_score' for d in st.DIMS]].mean(axis=1)
    adults['BDSF_level']=core.classify(adults.BDSF_score).astype(str)
    frame['05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv']=adults
    qc=adults[['respondent_id','BDSF_score']].copy();qc['recalc_BDSF']=sum(adults[cols].sum(axis=1)/len(cols)/4 for cols in st.DIMS.values());qc['max_abs_difference']=abs(qc.BDSF_score-qc.recalc_BDSF);qc['qc_status']=np.where(qc.max_abs_difference<1e-10,'PASS','FAIL')
    if not qc.qc_status.eq('PASS').all():raise ValueError('Independent scoring verification failed')
    frame['05_scoring/05_03_BDSF_Kazakhstan_Scoring_QC.csv']=qc
    for relative,data in frame.items():
        dest=build/relative;dest.parent.mkdir(parents=True,exist_ok=True);data.to_csv(dest,index=False,float_format='%.12g')
    return adults,frame

def check_reference(qmap):
    wb=openpyxl.load_workbook(HERE/'inputs/reference_scoring_matrix.xlsx',read_only=True,data_only=True)
    rows=list(wb['Scoring'].values); wb.close()
    if len(rows)!=15:raise ValueError('Reference matrix must have 14 items')
    for (_,q),r in zip(qmap.iterrows(),rows[1:]):
        if q.dimension_code!=r[0] or q.source_question_code!=r[4] or [r[i] for i in [7,9,11,13,15]]!=[0,25,50,75,100]:raise ValueError('Reference question/dimension/answer score mismatch')

def reports(build,s,adult):
    texts={
    '01_audit/01_07_protocol_compliance_check.md':'# Хаттама\n\nАвтор мөрлі хаттама үзіндісін кейін енгізетінін хабарлады. Оның мазмұны мен ересектердің келісім рәсімі автоматты түрде расталды деп есептелмейді.',
    '01_audit/01_08_audit_report.md':'# Бастапқы деректер аудиті\n\nЕкі CSV: KZ=94, RU=86; барлығы 180. Бағандар мен барлық жауап санаттары бекітілген картаға сәйкес тексерілді. Танылмаған мән болса конвейер тоқтайды. Қайталану белгілері ақпараттық сипатта.',
    '03_harmonization/03_04_BDSF_Kazakhstan_Harmonization_Rules.md':'# Үйлестіру\n\nСұрақтар inputs/question_map.csv картасымен, жауаптар кодтағы анық сәйкестікпен тексеріледі. Бастапқы жауап мәтіндері сақталады. Уақыт белдеуі +05:00. Бос не белгісіз мәндер үнсіз алмастырылмайды.',
    '03_harmonization/03_05_BDSF_Kazakhstan_Harmonization_Report.md':'# Үйлестіру нәтижесі\n\n180 жазба; 180 бірегей ID; 14 сұрақ. Бұл бастапқы толық жиынтық, аналитикалық іріктеме емес.',
    '04_cleaning/04_04_BDSF_Kazakhstan_Cleaning_Log.md':'# Шығару журналы\n\n180 − 21 LT18 = 159 ересек. Басқа шығару жоқ. Шығарылған ID және себептер 04_03 файлында.',
    '04_cleaning/04_05_BDSF_Kazakhstan_Cleaning_Report.md':f'# Тазалау\n\n159 ересек толық сақталды. Сапа жалаушасы бар {adult.quality_flags.ne("").sum()} ересек те талдауға енгізілді.',
    '05_scoring/05_04_BDSF_Kazakhstan_Scoring_Log.md':'# Бағалау\n\n0/25/50/75/100 жауап балдары; D1=Q01–Q04, D2=Q05–Q07, D3=Q08–Q11, D4=Q12–Q14. BDSF=(D1+D2+D3+D4)/4. Деңгей үшін ғана екі ондыққа дөңгелектенеді. Тең өлшемдік салмақ қолданыстағы спецификациядан сақталған; екінші Excel-де жеке өлшемдік салмақ бағаны жоқ.',
    '06_descriptive_analysis/06_06_BDSF_Kazakhstan_Descriptive_Report.md':f'# Сипаттамалық нәтиже\n\nn=159; BDSF={s["BDSF"]["mean"]:.6f}; SD={s["BDSF"]["sd"]:.6f}; 95% CI [{s["BDSF"]["ci95_low"]:.6f}, {s["BDSF"]["ci95_high"]:.6f}].',
    '07_reliability_analysis/07_04_BDSF_Kazakhstan_Reliability_Report.md':'# Сенімділік\n\nБір ғана n=159 іріктеме.\n\n'+'\n'.join(f'- {r["scale"]}: alpha={r["cronbach_alpha"]:.6f}' for r in s['reliability']),
    '08_construct_validity/08_07_BDSF_Kazakhstan_Validity_Report.md':f'# Конструктивтік валидтілік\n\nKMO={s["validity"]["KMO"]:.6f}; Bartlett={s["validity"]["Bartlett_chi2"]:.6f}, df=91. Parallel analysis: {s["validity"]["parallel_factors"]} фактор. Төрт факторлық MINRES/oblimin зерттеушілік шешім. REVIEW: {", ".join(s["validity"]["review_items"])}.\n\nҚиғаш айналдырудағы коммуналдылық diag(L Phi Lᵀ) арқылы есептелді; бұрынғы квадраттар қосындысы түзетілді. Бағалау формуласы мен негізгі нәтижелер өзгермейді. Төрт факторлық құрылым толық расталған жоқ.',
    '09_group_comparison/09_07_BDSF_Kazakhstan_Group_Comparison_Report.md':'# Топтық салыстыру\n\nТіл және жыныс бойынша әр салыстыру барлық 159 ересекті қамтиды. Welch t және әр салыстыру түріндегі бес нәтиже үшін Holm түзетуі қолданылды. Жас тек сипаттамалық түрде көрсетілді. Мәнді емес нәтиже эквиваленттілікті дәлелдемейді.',
    '11_figures_and_tables/11_09_Figure_and_Table_Guide.md':'# Жарияланым материалдары\n\nБарлық кесте мен сурет бір ғана 159 ересек іріктемеден жасалған. 11_02 өлшемдер, 11_03 индикаторлар, 11_04 сенімділік, 11_05 EFA. Бір деректер құрылымы CSV, XLSX, графиктер және мақала үшін пайдаланылады.'}
    for path,text in texts.items():write(build/path,text+'\n')

def workbook_specs(build,frames,qmap):
    books={
    '03_harmonization/03_02_BDSF_Kazakhstan_Codebook.xlsx':[('Questions',qmap),('Levels',pd.DataFrame(list(LEVELS.items()),columns=['level','score']))],
    '04_cleaning/04_02_BDSF_Kazakhstan_Validation_Report.xlsx':[('Population',pd.DataFrame([['Collected',180],['LT18 excluded',21],['Adults included',159]],columns=['population','n']))],
    '05_scoring/05_02_BDSF_Kazakhstan_Scoring_Report.xlsx':[],
    '06_descriptive_analysis/06_01_BDSF_Kazakhstan_Descriptive_Statistics.xlsx':[],
    '07_reliability_analysis/07_01_BDSF_Kazakhstan_Reliability_Analysis.xlsx':[],
    '08_construct_validity/08_01_BDSF_Kazakhstan_Construct_Validity.xlsx':[],
    '09_group_comparison/09_01_BDSF_Kazakhstan_Group_Comparisons.xlsx':[],
    '11_figures_and_tables/11_01_BDSF_Kazakhstan_Publication_Tables.xlsx':[]}
    for book,sheets in books.items():
        stage=book.split('/')[0]
        for path,df in frames.items():
            if path.split('/')[0]!=stage:continue
            prefix=Path(path).name[:5]
            if stage[:2] in ['03'] or prefix in ['04_01','05_01']:continue
            sheets.append((prefix,df))
    payload={}
    for book,sheets in books.items():
        if not sheets:raise ValueError('Workbook cannot be empty: '+book)
        payload[book]=[dict(name=name,rows=[df.columns.tolist()]+df.astype(object).where(df.notna(),None).values.tolist()) for name,df in sheets]
    jsonwrite(build/'workbooks.json',payload)

def verify_outputs(build,frames,summary):
    checks=[]
    def record(name,ok,detail=''):checks.append(dict(check=name,status='PASS' if ok else 'FAIL',detail=detail))
    for rel,df in frames.items():
        actual=pd.read_csv(build/rel,keep_default_na=False)
        record('CSV '+rel,len(actual)==len(df) and list(actual.columns)==list(df.columns))
        same=True
        for col in df.columns:
            if pd.api.types.is_numeric_dtype(df[col]) and not pd.api.types.is_bool_dtype(df[col]):
                same &= np.allclose(pd.to_numeric(actual[col],errors='coerce'),pd.to_numeric(df[col],errors='coerce'),rtol=1e-10,atol=1e-10,equal_nan=True)
            else:same &= actual[col].fillna('').astype(str).tolist()==df[col].fillna('').astype(str).tolist()
        record('CSV cell values '+rel,bool(same))
    raw=frames['03_harmonization/03_01_BDSF_Kazakhstan_Harmonized_Dataset.csv'];adult=frames['05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv'];clean=frames['04_cleaning/04_01_BDSF_Kazakhstan_Clean_Dataset.csv'];excluded=frames['04_cleaning/04_03_BDSF_Kazakhstan_Exclusion_Log.csv']
    record('Exact adult ID partition',set(adult.respondent_id)==set(raw.loc[raw.age_group_code.ne('LT18'),'respondent_id']) and set(excluded.respondent_id)==set(raw.loc[raw.age_group_code.eq('LT18'),'respondent_id']))
    record('159 clean IDs equal scored IDs',len(clean)==159 and set(clean.respondent_id)==set(adult.respondent_id))
    record('No duplicate IDs',adult.respondent_id.nunique()==159)
    record('Flags do not exclude adults',adult.quality_flags.ne('').sum()==20)
    record('Score range',bool(adult[st.ITEMS+['BDSF_score']].ge(0).all().all() and adult[st.ITEMS+['BDSF_score']].le(100).all().all()))
    specs=json.loads((build/'workbooks.json').read_text(encoding='utf-8'))
    for book,sheets in specs.items():
        wb=openpyxl.load_workbook(build/book,read_only=True,data_only=True)
        record('Workbook sheets '+book,len(sheets)>0 and wb.sheetnames==[x['name'] for x in sheets])
        for spec in sheets:
            rows=list(wb[spec['name']].values);expected=spec['rows'];ok=len(rows)==len(expected)
            for actual,wanted in zip(rows,expected):
                for a,b in zip(actual,wanted):
                    if b is None or b=='':ok &= a is None or a==''
                    elif isinstance(b,(int,float)) and not isinstance(b,bool):ok &= isinstance(a,(int,float)) and np.isclose(a,b,rtol=1e-10,atol=1e-10)
                    else:ok &= str(a)==str(b)
            record('XLSX '+book+'/'+spec['name'],ok)
        wb.close()
    for p in (build/'11_figures_and_tables').glob('*.png'):record('Figure '+p.name,p.stat().st_size>10000)
    # Full manuscript is also checked against the original author and reference text.
    from docx import Document
    full=Document(build/'12_article_results/12_01_BDSF_Kazakhstan_Full_Updated_Article.docx');original=Document(HERE/'inputs/original_article.docx');text='\n'.join(p.text for p in full.paragraphs)+'\n'+'\n'.join(c.text for t in full.tables for row in t.rows for c in row.cells)
    for i in list(range(2,7))+list(range(66,83))+list(range(84,89))+list(range(96,101)):record('Preserved original paragraph '+str(i),original.paragraphs[i].text in text)
    for token in ['176 студент','177 респондент','63.20','63,20','0.4966','0.5416','174.4425','Varimax','HC3','no-flags']:record('No stale manuscript token '+token,token not in text)
    record('Article mean',f'{summary["BDSF"]["mean"]:.2f}'.replace('.',',') in text)
    record('Article N','159' in text)
    record('Figures embedded',len(full.inline_shapes)==3)
    record('Only adult population',summary['n']==159)
    if any(r['status']=='FAIL' for r in checks):raise ValueError('Output validation failed: '+str([r for r in checks if r['status']=='FAIL']))
    return checks

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--project-root',type=Path,default=HERE.parent);parser.add_argument('--build-only',action='store_true');args=parser.parse_args()
    project=args.project_root.resolve();build=project/'work'/('build_'+uuid.uuid4().hex[:12]);build.mkdir(parents=True);print('BUILD_DIRECTORY='+str(build),flush=True)
    sources=list((project/'00_sources').glob('*'))+list((HERE/'inputs').glob('*'));before={p.relative_to(project).as_posix():sha(p) for p in sources if p.is_file()}
    try:
        qmap=pd.read_csv(HERE/'inputs/question_map.csv');check_reference(qmap)
        adult,frames=harmonize(project,qmap,build);summary,analytic=st.analyse(adult,build);frames.update(analytic)
        actual=core.calculate(adult);expected=json.loads((HERE/'13_03_expected_results.json').read_text(encoding='utf-8'));baseline=core.verify(actual,expected)
        if any(r['status']!='PASS' for r in baseline):raise ValueError('Baseline disagreement: '+str([r for r in baseline if r['status']!='PASS']))
        reports(build,summary,adult);jsonwrite(build/'results.json',summary);workbook_specs(build,frames,qmap)
        run([sys.executable,'-I','-X','utf8',HERE/'13_export_workbooks.py',build])
        run([sys.executable,'-I','-X','utf8',HERE/'13_article.py','--build',build,'--original',HERE/'inputs/original_article.docx'])
        checks=verify_outputs(build,frames,summary)
        if any(sha(project/p)!=h for p,h in before.items()):raise ValueError('Immutable input modified during build')
        out=build/'14_final_quality_control';out.mkdir(exist_ok=True)
        pd.DataFrame(checks).to_csv(out/'14_02_BDSF_Kazakhstan_Final_QC_Checks.csv',index=False)
        pd.DataFrame(baseline).to_csv(out/'14_08_Numerical_Baseline_Checks.csv',index=False)
        versions={}
        from importlib.metadata import version
        for package in ['numpy','pandas','scipy','scikit-learn','factor-analyzer','matplotlib','python-docx','openpyxl']:versions[package]=version(package)
        status=dict(technical_status='PASS',checks_passed=len(checks),baseline_passed=len(baseline),n=159,excluded_under18=21,seed=st.SEED,python=platform.python_version(),dependencies=versions,ethics_status='AUTHOR_PROTOCOL_PENDING_INSERTION',consent_evidence='NOT_VERIFIED',original_archive_available=False,scientific_validity='FOUR_FACTOR_STRUCTURE_NOT_CONFIRMED',source_hashes=before)
        jsonwrite(out/'14_03_BDSF_Kazakhstan_Final_QC_Summary.json',status)
        write(out/'14_04_BDSF_Kazakhstan_Final_QC_Report.md',f'# Конвейер тексерімі\n\nТехникалық тексеру: {len(checks)}/{len(checks)} PASS. Негізгі эталон: {len(baseline)}/{len(baseline)} PASS. Бастапқы CSV-ден 180 жауап оқылып, тек 21 LT18 шығарылды; 159 ересек сақталды.\n\nМөрлі хаттаманы автор кейін енгізеді; келісім рәсімдері автоматты түрде расталған жоқ. Ескі ZIP қолжетімділігі: {status["original_archive_available"]}. Төрт факторлық құрылым толық расталмаған. Бастапқы әдебиеттер сақталды; олардың библиографиялық дұрыстығына жаңа сыртқы аудит жүргізілген жоқ.\n\nЕсептер толық дәлдікпен есептелді. Қиғаш EFA коммуналдылығы diag(L Phi Lᵀ) формуласына түзетілді.\n')
        write(out/'14_05_BDSF_Kazakhstan_Submission_Checklist.md','# Авторлық тексеру\n\n- Мөрлі хаттама үзіндісін енгізу және оның қолданылу аясын растау.\n- Ересектер келісімі мен құпиялық рәсімдерін нақты құжатпен көрсету.\n- Әдебиеттердің библиографиялық дәлдігін және журнал форматын тексеру.\n- Жаңа конвейер репозиторийге жарияланды деген тұжырымды жариялау орындалмайынша қоспау.\n')
        write(out/'14_06_BDSF_Kazakhstan_Correction_Log.md','# Өзгерістер\n\nДайын балдарды оқитын тексерім толық бастапқы CSV конвейерімен толықтырылды. Барлық нәтижелер бір 159 ересек іріктемеден шығарылады. Толық мақала бастапқы авторлар, кіріспе және әдебиеттер сақталып жаңартылды. Ескі 10 индикаторлық стандарт сәйкестігі кестесі расталмаған жаңа сәйкестіктермен алмастырылған жоқ; орнына 14 сұрақтың нақты құрам кестесі берілді. Қиғаш EFA коммуналдылығы түзетілді.\n')
        jsonwrite(build/'13_reproducibility_package/reproduced_results/13_07_reproduced_results.json',dict(provenance=dict(input='00_sources/*.csv',source_hashes=before),results=actual))
        dest=build/'13_reproducibility_package/reproduced_results/13_08_verification_report.csv';pd.DataFrame(baseline).to_csv(dest,index=False)
        # Manifest uses relative paths; portable previews are published with the results.
        manifest=[]
        for p in sorted(build.glob('*/*')):
            if p.is_file() and p.suffix in ['.csv','.md','.json','.png','.docx','.xlsx'] and not p.name.endswith('.inspect.ndjson'):manifest.append(dict(path=p.relative_to(build).as_posix(),sha256=sha(p),size_bytes=p.stat().st_size))
        pd.DataFrame(manifest).to_csv(out/'14_01_BDSF_Kazakhstan_Final_File_Manifest.csv',index=False)
        jsonwrite(build/'13_reproducibility_package/reproduced_results/13_09_full_analysis_results.json',summary)
        qa=build/'review_materials';qa.rename(build/'qa_stage');qa.mkdir();shutil.move(str(build/'qa_stage'),str(qa/'generated'))
        if args.build_only:print('BUILD_ONLY_PASS='+str(build));return 0
        publish(project,build)
        print(f'PIPELINE PASS: raw=180; adult=159; excluded=21; baseline={len(baseline)}/{len(baseline)}; outputs={len(checks)}/{len(checks)}')
        print('FULL_ARTICLE='+str(project/'12_article_results/12_01_BDSF_Kazakhstan_Full_Updated_Article.docx'))
        return 0
    except Exception as exc:
        jsonwrite(build/'FAILED.json',dict(status='FAILED_NOT_PUBLISHED',error=str(exc)));raise

def publish(project,build):
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f');archive=project/'work/recovery'/('before_run_'+stamp);archive.mkdir(parents=True)
    stages=['01_audit','03_harmonization','04_cleaning','05_scoring','06_descriptive_analysis','07_reliability_analysis','08_construct_validity','09_group_comparison','11_figures_and_tables','12_article_results','14_final_quality_control','13_reproducibility_package/reproduced_results','review_materials/generated']
    installed=[];saved=[]
    try:
        for name in stages:
            target=(project/name).resolve();source=build/name
            if not target.is_relative_to(project) or not source.is_dir():raise ValueError('Unsafe or absent stage: '+name)
            old=archive/name;old.parent.mkdir(parents=True,exist_ok=True)
            if target.exists():shutil.move(str(target),str(old));saved.append((target,old))
            shutil.move(str(source),str(target));installed.append(target)
    except Exception:
        for target in reversed(installed):
            recovery=build/'rollback'/target.relative_to(project);recovery.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(target),str(recovery))
        for target,old in reversed(saved):shutil.move(str(old),str(target))
        raise
    print('ARCHIVE='+str(archive))

if __name__=='__main__':raise SystemExit(main())
