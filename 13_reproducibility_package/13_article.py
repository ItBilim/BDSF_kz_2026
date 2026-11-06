"""Merge new results into the original article's scientific and author context."""
from pathlib import Path
import argparse,json
import pandas as pd
from docx import Document
from docx.shared import Inches,Pt,RGBColor,Cm
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

def build_article(build,original):
    s=json.loads((build/'results.json').read_text(encoding='utf-8'));old=Document(original);d=Document()
    section=d.sections[0];section.page_width=Cm(21);section.page_height=Cm(29.7);section.top_margin=Cm(2);section.bottom_margin=Cm(2);section.left_margin=Cm(2.2);section.right_margin=Cm(2.2)
    normal=d.styles['Normal'];normal.font.name='Times New Roman';normal.font.size=Pt(11);normal.paragraph_format.space_after=Pt(6);normal.paragraph_format.line_spacing=1.08
    for sty in ['Title','Heading 1','Heading 2']:
        d.styles[sty].font.name='Times New Roman';d.styles[sty].font.color.rgb=RGBColor(0,0,0);d.styles[sty].font.size=Pt(14 if sty=='Title' else 12);d.styles[sty].font.bold=True
    for style in d.styles:
        for border in list(style.element.iter(qn('w:pBdr'))):border.getparent().remove(border)
    def label(code):return {'Language':'Сауалнама тілі','Gender':'Жыныс','Age':'Жас','KZ':'Қазақша','RU':'Орысша','MALE':'Ер','FEMALE':'Әйел','18_20':'18–20','21_23':'21–23','GE27':'27 және одан жоғары','14_items':'14 сұрақ','4_dimensions':'4 өлшем'}.get(code,code)
    def p(text,style=None):return d.add_paragraph(text,style)
    def heading(text):return d.add_heading(text,level=1)
    def fmt(value,n=2):return f'{value:.{n}f}'.replace('.',',')
    def table(title,headers,rows,widths=None):
        cap=p(title);cap.paragraph_format.keep_with_next=True;cap.runs[0].bold=True
        t=d.add_table(rows=1,cols=len(headers));t.autofit=False
        if widths:
            for col,w in zip(t.columns,widths):col.width=Cm(w)
        for cell,text in zip(t.rows[0].cells,headers):cell.text=str(text)
        for row in rows:
            for cell,text in zip(t.add_row().cells,row):cell.text=str(text)
        for ri,row in enumerate(t.rows):
            for ci,cell in enumerate(row.cells):
                if widths:cell.width=Cm(widths[ci])
                tcpr=cell._tc.get_or_add_tcPr();sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'17365D' if ri==0 else ('EAF0F5' if ri%2==0 else 'FFFFFF'));tcpr.append(sh)
                borders=OxmlElement('w:tcBorders')
                for edge in ['top','left','bottom','right']:
                    b=OxmlElement('w:'+edge);b.set(qn('w:val'),'single');b.set(qn('w:sz'),'4');b.set(qn('w:color'),'D9D9D9');borders.append(b)
                tcpr.append(borders)
                for para in cell.paragraphs:
                    para.paragraph_format.space_after=Pt(3);para.paragraph_format.space_before=Pt(3)
                    if ri==0:para.paragraph_format.keep_with_next=True
                    for run in para.runs:
                        run.font.size=Pt(10);run.bold=ri==0;run.font.color.rgb=RGBColor(255,255,255) if ri==0 else RGBColor(0,0,0)
            trpr=row._tr.get_or_add_trPr();cant=OxmlElement('w:cantSplit');trpr.append(cant)
        repeat=OxmlElement('w:tblHeader');t.rows[0]._tr.get_or_add_trPr().append(repeat)
        p('')
        return t
    def figure(name,caption,width=15):
        para=d.add_paragraph();para.paragraph_format.keep_with_next=True;para.add_run().add_picture(str(build/'11_figures_and_tables'/name),width=Cm(width))
        cap=p(caption);cap.paragraph_format.keep_with_next=False;cap.alignment=1
        for r in cap.runs:r.font.size=Pt(10)
    b=s['BDSF'];alpha={r['scale']:r['cronbach_alpha'] for r in s['reliability']};v=s['validity']
    kzabstract=f'Аңдатпа. Зерттеудің мақсаты халықаралық стандарттар мен ғылыми әдебиеттердің тұжырымдамалық негізінде студенттердің цифрлық қауіпсіздік мінез-құлқын бағалауға арналған BDSF моделін апробациялау болды. Жаңартылған құрал төрт өлшемге біріктірілген 14 сұрақтан тұрады. Екітілді сауалнамадағы 180 жауаптың ішінен ата-ана немесе заңды өкіл рұқсаты жоқ 18 жасқа толмаған 21 жауап осы қайта талдауда шығарылып, 159 ересек қатысушы сақталды. Әр өлшем жалпы көрсеткішке 25% үлес қосады. BDSF орташа мәні {fmt(b["mean"])} (SD={fmt(b["sd"])}; 95% CI {fmt(b["ci95_low"])}–{fmt(b["ci95_high"])}) болды. 14 сұрақ үшін Cronbach alpha={fmt(alpha["14_items"],3)}. KMO={fmt(v["KMO"],3)}; Bartlett χ²(91)={fmt(v["Bartlett_chi2"])} (p<0,001). Parallel analysis {v["parallel_factors"]} басым факторды қолдады; теориялық төрт өлшемнің толық ажырауы расталмады. Тіл және жыныс бойынша Holm түзетуінен кейін мәнді айырмашылық анықталмады. Ғылыми үлес бағалаудың ашық алгоритмін және оның сенімділігі мен құрылымдық шектеулерін бір деректер жиынында көрсетуге негізделеді.'
    enabstract=f'Abstract. This study evaluated an updated Behavioral Digital Security Assessment Framework using the conceptual context of international standards and previous research. Fourteen items were grouped into four equally weighted dimensions. Of 180 bilingual survey responses, 21 from participants under 18 were excluded in this reanalysis because parental or legal-guardian permission was unavailable, leaving 159 adults. The mean BDSF score was {b["mean"]:.2f} (SD={b["sd"]:.2f}; 95% CI {b["ci95_low"]:.2f}–{b["ci95_high"]:.2f}). Cronbach alpha for the 14 items was {alpha["14_items"]:.3f}. KMO was {v["KMO"]:.3f}; Bartlett χ²(91)={v["Bartlett_chi2"]:.2f} (p<0.001). Parallel analysis supported {v["parallel_factors"]} dominant factor; the proposed four-dimensional structure was not fully confirmed. No language or gender differences remained significant after Holm correction. The contribution is a reproducible scoring workflow with explicit evidence on internal consistency and structural limitations. Independent confirmatory validation remains necessary.'
    ruabstract=f'Аннотация. Проведена апробация обновлённой модели BDSF в концептуальном контексте международных стандартов и научной литературы. Четырнадцать вопросов объединены в четыре измерения с равным весом 25%. Из 180 ответов при повторном анализе исключён 21 ответ участника младше 18 лет ввиду отсутствия разрешения родителя или законного представителя; сохранены все 159 взрослых участников. Средний BDSF составил {fmt(b["mean"])} (SD={fmt(b["sd"])}; 95% ДИ {fmt(b["ci95_low"])}–{fmt(b["ci95_high"])}). Cronbach alpha для 14 вопросов — {fmt(alpha["14_items"],3)}. KMO={fmt(v["KMO"],3)}; критерий Бартлетта χ²(91)={fmt(v["Bartlett_chi2"])} (p<0,001). Параллельный анализ поддержал один доминирующий фактор; полное разделение четырёх теоретических измерений не подтверждено. После поправки Holm значимых различий по языку и полу не обнаружено. Научный вклад состоит в воспроизводимом алгоритме оценки и совместном анализе надёжности и структурных ограничений инструмента. Требуется независимая подтверждающая валидация.'
    p(old.paragraphs[0].text)
    for i in range(2,7):p(old.paragraphs[i].text)
    p(old.paragraphs[8].text,'Title');p(kzabstract);p(old.paragraphs[12].text)
    heading('Кіріспе')
    for i in range(15,21):
        text=old.paragraphs[i].text
        if i==16:text=text.replace('10 мінез-құлық','14 мінез-құлық')
        p(text)
    p('Халықаралық киберқауіпсіздік стандарттары (ХКС) [12–17] осы жұмыста тұжырымдамалық негіз ретінде қарастырылады. Олар студенттердің мінез-құлқын өлшейтін дайын психометриялық шкала ретінде қабылданбайды. Жаңартылған 14 сұрақ төрт мазмұндық өлшемге бекітілді. Ескі он индикаторлық кестенің нормативтік сәйкестіктері жаңа сұрақтарға автоматты түрде көшірілген жоқ; төменде тек тексерілген сұрақ құрамы берілген.')
    qmap=pd.read_csv(build/'03_harmonization/03_03_BDSF_Kazakhstan_Bilingual_Question_Map.csv')
    labels=['Құпиялық параметрлерін қайта қарау','Қолданбалар рұқсаттарын басқару','Жеке ақпаратты жариялау ережелері','Қатаң құпиялық параметрлері','Қауіпсіз қосылымды тексеру','Қауіпсіз интернет желісін таңдау','Қосылатын қолданбаларды тексеру','Осал қолданбаларды жаңарту','Жаңа киберқауіптер туралы ақпарат іздеу','Зиянды бағдарламалардан қорғау','Құрылғының негізгі қауіпсіздігі','Аккаунт белсенділігін бақылау','Құпиясөз құпиялығын сақтау','Аккаунтты қалпына келтіруді дайындау']
    table('1 кесте Жаңартылған құралдың құрамы',['Өлшем','Сұрақ','Индикатор','Дереккөз коды'],[[r.dimension_code,r.item_code,labels[i],r.source_question_code] for i,r in enumerate(qmap.itertuples())],[1.7,1.7,10.8,2.4])
    heading('Материалдар мен тәсілдер')
    p('Көлденең қималық сауалнама Google Forms арқылы қазақ және орыс тілдерінде жүргізілді. Уақыт белгілеріне сәйкес жауаптар 2026 жылғы 20–24 қазанда жиналған. Бастапқы CSV файлдарында 180 жауап болды. Осы қайта талдауда 18 жасқа толмаған 21 жауап ата-ана немесе заңды өкіл рұқсаты болмағандықтан шығарылды. Қалған 159 ересектің барлығы сақталды; сапа жалаушалары қосымша шығару өлшемі болған жоқ. Интервенция, рандомизация немесе бақылау тобы қолданылған жоқ.')
    p('Мөрлі хаттама үзіндісі авторда бар және кейін енгізіледі. Оның нақты нөмірі, күні және осы жаңартылған сауалнамаға қолданылу аясы түпнұсқа құжат енгізілгеннен кейін көрсетіледі. Ересектердің ақпараттандырылған келісімі мен құпиялық рәсімдері туралы расталмаған мәлімет бұл мәтінде орындалған рәсім ретінде берілмейді.')
    table('2 кесте Аналитикалық іріктеме',['Сипаттама','Санат','n','%'],[[label(r.Characteristic),label(r.Category),r.n,fmt(r.Percent,1)] for r in pd.read_csv(build/'11_figures_and_tables/11_07_Table_1_Sample_Characteristics.csv').itertuples()],[4,7,2,3.6])
    p('Жалпы алгоритм бастапқы жауаптарды үйлестіру, жауап балдарын есептеу, өлшемдерді орташалау және жалпы индексті алу қадамдарын қамтиды. Жаңартылған сұрақтардың өлшемге тиесілілігі және жауап балдары екінші бағалау матрицасымен салыстырылды. Жауаптар NEVER, RARELY, SOMETIMES, OFTEN, ALWAYS деңгейлеріне және тиісінше 0, 25, 50, 75, 100 балға ауыстырылды. Q10 арнайы тұжырымдары сол бес мағыналық деңгейге сәйкестендірілді. Бос немесе танылмаған мән нөлге ауыстырылмайды: конвейер нақты қатемен тоқтайды.')
    p('D1=mean(Q01–Q04), D2=mean(Q05–Q07), D3=mean(Q08–Q11), D4=mean(Q12–Q14). Жалпы BDSF=(D1+D2+D3+D4)/4. Әр өлшемнің 25% салмағы қолданыстағы спецификациядан өзгеріссіз сақталды; жауап балдары берілген екінші Excel-де жеке өлшемдік салмақ бағаны жоқ. Деңгейді анықтау үшін BDSF екі ондыққа дөңгелектенеді: 0,00–33,32 төмен, 33,33–66,66 орташа, 66,67–100,00 жоғары. Бұл шектер сыртқы критериймен валидацияланған нормативтер емес.')
    p('Орташа, медиана, стандарттық ауытқу және Student t үлестіріміне негізделген 95% сенімділік аралықтары есептелді. Ішкі келісімділік Cronbach alpha, стандартталған alpha, corrected item–total correlation және alpha if item deleted арқылы бағаланды. Сұрақтар осы диагностикалар негізінде автоматты түрде алынып тасталған жоқ.')
    p(f'Факторлық талдауға жарамдылық KMO және Bartlett тестімен бағаланды. Pearson корреляциялары бойынша parallel analysis 500 қалыпты кездейсоқ деректер жиынының 95-процентильдік меншікті мәндерімен салыстырылды; seed={s["seed"]}. Теориялық құрылым үшін төрт факторлық MINRES және oblimin айналдыру қолданылды. Фактор бағандары күтілген өлшемдердегі абсолют жүктемелердің орташа шамасына қарай сәйкестендірілді; бұл қадам төрт өлшемді растау болып саналмайды. Қиғаш шешімнің коммуналдылығы diag(L Phi Lᵀ) бойынша есептелді. Күтілген фактор басым болып, жүктеме ≥0,50 және ең үлкен қосалқы жүктеме <0,30 болса STRONG; күтілген фактор басым және жүктеме ≥0,40 болса ACCEPTABLE; қалғандары REVIEW деп белгіленді.')
    p('Тіл және жыныс бойынша Welch t-test, орташа айырманың 95% CI және Hedges g есептелді. Әр салыстыру түріндегі бес нәтиже үшін Holm түзетуі қолданылды. Әр салыстыру барлық 159 ересекті қамтиды. Бір адамдық 27+ тобы болғандықтан, жас құрамы тек сипаттамалық түрде берілді. Есептеулер Python арқылы, Excel есептері бір есептелген деректер құрылымынан жасалды. Бір команда бастапқы CSV-ден есептер мен осы мақаланың сандық мазмұнын қайта қалыптастырады.')
    heading('Нәтижелер')
    p(f'Негізгі іріктеме n={s["n"]}: қазақша жауаптар {s["language"]["KZ"]}, орысша {s["language"]["RU"]}. Жалпы BDSF M={fmt(b["mean"])}; медиана={fmt(b["median"])}; SD={fmt(b["sd"])}; 95% CI {fmt(b["ci95_low"])}–{fmt(b["ci95_high"])}. Төмен деңгейде {s["levels"]["LOW"]}, орташа деңгейде {s["levels"]["MEDIUM"]}, жоғары деңгейде {s["levels"]["HIGH"]} қатысушы болды.')
    table('3 кесте Өлшемдер мен жалпы көрсеткіш',['Көрсеткіш','n','M','SD','95% CI'],[[name,z['n'],fmt(z['mean']),fmt(z['sd']),fmt(z['ci95_low'])+'–'+fmt(z['ci95_high'])] for name,z in [*s['dimensions'].items(),('BDSF',b)]],[3.2,1.5,2.6,2.6,6.7])
    figure('11_02_Figure_1_Dimension_Means.png','1 сурет Өлшемдер мен BDSF орташа мәндері')
    ind=pd.read_csv(build/'06_descriptive_analysis/06_02_BDSF_Kazakhstan_Indicator_Statistics.csv');high=ind.loc[ind['mean'].idxmax()];low=ind.loc[ind['mean'].idxmin()]
    p(f'Өлшемдер арасында D4 орташа мәні ең жоғары ({fmt(s["dimensions"]["D4"]["mean"])}), D3 ең төмен ({fmt(s["dimensions"]["D3"]["mean"])}). Индикаторлар арасында {high.variable} ең жоғары (M={fmt(high["mean"])}), {low.variable} ең төмен (M={fmt(low["mean"])}) көрсеткіш көрсетті.')
    table('4 кесте Ішкі келісімділік',['Шкала','Сұрақ/өлшем саны','n','Cronbach alpha'],[[label(r['scale']),r['items'],r['n'],fmt(r['cronbach_alpha'],3)] for r in s['reliability']],[4.2,5,2.4,5])
    figure('11_04_Figure_3_Reliability.png','2 сурет Негізгі іріктемедегі ішкі келісімділік')
    p(f'KMO={fmt(v["KMO"],3)}; Bartlett χ²(91)={fmt(v["Bartlett_chi2"])} (p<0,001). Parallel analysis {v["parallel_factors"]} басым факторды ұсынды. Төрт факторлық шешімде қайта қарауды қажет ететін сұрақтар: {", ".join(v["review_items"])}. Сондықтан жоғары alpha төрт тәуелсіз өлшемнің валидтілігін дәлелдемейді.')
    figure('11_05_Figure_4_Factor_Loadings.png','3 сурет Төрт факторлық зерттеушілік шешім',12)
    for r in s['comparisons']:
        if r['outcome']=='BDSF':p(f'{"Тіл" if r["comparison_type"]=="Language" else "Жыныс"} бойынша {r["group_1"]}: M={fmt(r["mean_1"])} (n={r["n_1"]}); {r["group_2"]}: M={fmt(r["mean_2"])} (n={r["n_2"]}). Айырма {fmt(r["difference_group2_minus_group1"])}; 95% CI {fmt(r["ci95_low"])}–{fmt(r["ci95_high"])}; Holm p={fmt(r["p_holm"],3)}; Hedges g={fmt(r["hedges_g"],3)}.')
    p('Holm түзетуінен кейін D1–D4 және жалпы BDSF үшін тіл немесе жыныс бойынша мәнді айырмашылық анықталмады. Бұл тұжырым топтардың эквиваленттілігі немесе өлшеу инварианттылығы дәлелденді дегенді білдірмейді.')
    heading('Талқылау')
    p('Жаңартылған деректер D4 және D3 мазмұндық бағыттары арасында айырмашылық бар екенін көрсетті. Бұл көрсеткіштер өзіндік есепке негізделеді; қауіпсіздік туралы білім немесе нақты қорғану әрекеті бөлек өлшенбегендіктен, олардың арасындағы алшақтық дәлелденді деп айтуға болмайды. Бастапқы әдебиетте қарастырылған хабардарлық, мінез-құлық және оқыту мәселелері [1–11] зерттеудің контекстін береді, бірақ осы көлденең сауалнама оқыту тиімділігін бағаламайды.')
    p(f'Жаңартылған 14 сұрақтық жиынтықта alpha={fmt(alpha["14_items"],3)} болғанымен, бір басым фактордың сақталуы төрт өлшемді толық тәуелсіз конструкциялар ретінде түсіндіруді шектейді. D1–D4 әзірге жалпы индекстің мазмұндық профильдері ретінде қарастырылады. REVIEW сұрақтары мазмұндық және екі тілдегі тұжырым тұрғысынан қайта қаралуы тиіс; оларды осы іріктемеде жақсы статистика алу үшін автоматты шығару қолданылған жоқ.')
    p('Зерттеудің ғылыми үлесі жаңартылған мінез-құлық сұрақтарын ашық есептеу тәртібіне біріктіруде және ішкі келісімділік пен құрылымдық дәлелдердің айырмасын көрсетуінде. Бір ғана ересек іріктемеге негізделген конвейер бастапқы жауаптан мақала нәтижесіне дейінгі байланысты тексеруге мүмкіндік береді. Тіл мен жыныс бойынша нәтижелер сипаттамалық ассоциациялар ретінде түсіндіріледі.')
    heading('Шектеулер')
    p('Іріктеме ықтималдық емес және демографиялық тұрғыдан теңгерілмеген. Өзіндік есеп нақты киберқауіпсіздік әрекеттерімен сыртқы валидациядан өткен жоқ. Бес деңгейлі жауаптар Pearson корреляциясымен талданды; polychoric корреляциялар мен тәуелсіз CFA келесі зерттеуде қарастырылуы тиіс. Қазақ және орыс нұсқаларының өлшеу инварианттылығы тексерілген жоқ. Жіктеу шектері сыртқы критериймен расталмаған. Нәтижелер Қазақстандағы барлық студенттерге жалпыланбайды. Хаттаманың жаңартылған сауалнамаға қолданылу аясы және келісім рәсімдері құжат енгізілгеннен кейін нақтыланады.')
    heading('Қорытынды')
    p(f'Жаңартылған BDSF 159 ересек қатысушыда орташа жалпы көрсеткіш көрсетті (M={fmt(b["mean"])}). Тек 18 жасқа толмаған 21 жауап шығарылды; барлық қалған қатысушылар сақталды. D4 ең жоғары, D3 ең төмен орташа балға ие болды. Ішкі келісімділік жоғары болғанымен, төрт факторлық құрылым толық расталмады. Қазіргі нәтиже ересек студенттердің цифрлық қауіпсіздік мінез-құлқын сипаттайтын алдын ала үздіксіз индекс ретінде қолданылуы мүмкін; тәуелсіз растау қажет.')
    p('Қайта өндіру коды мен нәтижелер жоба құрамында беріледі. Осы жаңартылған нұсқаның ашық репозиторийге жарияланғаны тексерілген жоқ.')
    heading('Әдебиеттер')
    for i in range(66,83):p(f'[{i-65}] '+old.paragraphs[i].text)
    d.add_page_break()
    for i in range(84,89):p(old.paragraphs[i].text)
    p(old.paragraphs[90].text,'Title');p(enabstract);p(old.paragraphs[94].text)
    d.add_page_break()
    for i in range(96,101):p(old.paragraphs[i].text)
    p(old.paragraphs[102].text,'Title');p(ruabstract);p(old.paragraphs[106].text)
    # Add page fields; no administrative identity block from the journal cover sheet.
    footer=section.footer.paragraphs[0];footer.alignment=2;field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');footer._p.append(field)
    dest=build/'12_article_results';dest.mkdir(exist_ok=True);d.save(dest/'12_01_BDSF_Kazakhstan_Full_Updated_Article.docx')
    (dest/'12_02_BDSF_Kazakhstan_Article_Replacement_Map.md').write_text('# Толық мақаланы жаңарту\n\nАвторлар үш тілде, кіріспенің ғылыми контексті және 17 әдебиет сақталды. Әдістер, нәтижелер, талқылау, шектеулер, қорытынды және үш аңдатпа results.json нәтижелерінен жасалды. Бастапқы файл өзгертілмеді. Ескі 10 индикаторлық нормативтік кесте мен ескі статистикалық кестелер ауыстырылды. Журналға арналған жеке әкімшілік байланыс/жеке сәйкестендіру блогы ғылыми мақалаға көшірілген жоқ және түпнұсқада сақталады.\n',encoding='utf-8')
    (dest/'12_03_BDSF_Kazakhstan_Statistical_Reporting_Checklist.md').write_text('# Есептілік\n\nБір іріктеме n=159; 21 LT18 шығарылды. Барлық сандық мәтіндер конвейер нәтижелерінен қалыптасады. 3 сурет, 4 кесте. Автор мөрлі хаттаманы кейін енгізеді; нөмірі мен күні ойдан қосылмайды. Техникалық PASS ғылыми валидтілік не этикалық мақұлдау болып саналмайды.\n',encoding='utf-8')
    print('Full article generated from original context and calculated results')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--build',type=Path,required=True);a.add_argument('--original',type=Path,required=True);args=a.parse_args();build_article(args.build,args.original)
