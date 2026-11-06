"""Single adult sample statistics. Input: freshly scored in-memory dataframe only."""
from pathlib import Path
import math,json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from factor_analyzer import FactorAnalyzer
import factor_analyzer.factor_analyzer as fam
from sklearn.utils.validation import check_array

def compatible_check(*args,**kwargs):
    if 'force_all_finite' in kwargs: kwargs['ensure_all_finite']=kwargs.pop('force_all_finite')
    return check_array(*args,**kwargs)
fam.check_array=compatible_check
ITEMS=[f'Q{i:02d}_score' for i in range(1,15)]
DIMS={'D1':ITEMS[:4],'D2':ITEMS[4:7],'D3':ITEMS[7:11],'D4':ITEMS[11:]}
SEED=20261030

def describe(s):
    s=pd.Series(s,dtype=float);n=len(s);sd=s.std(ddof=1);se=sd/np.sqrt(n);t=stats.t.ppf(.975,n-1) if n>1 else float('nan')
    return dict(n=n,mean=s.mean(),median=s.median(),sd=sd,min=s.min(),q1=s.quantile(.25),q3=s.quantile(.75),max=s.max(),se=se,ci95_low=s.mean()-t*se,ci95_high=s.mean()+t*se)

def alpha(x):
    x=pd.DataFrame(x);k=x.shape[1]
    return k/(k-1)*(1-x.var(ddof=1).sum()/x.sum(axis=1).var(ddof=1))

def holm(p):
    p=np.asarray(p); order=np.argsort(p);adj=np.empty(len(p));running=0
    for rank,i in enumerate(order): running=max(running,min(1,p[i]*(len(p)-rank)));adj[i]=running
    return adj

def welch(a,b):
    a=np.asarray(a);b=np.asarray(b);n1=len(a);n2=len(b);v1=a.var(ddof=1);v2=b.var(ddof=1);diff=b.mean()-a.mean();se=np.sqrt(v1/n1+v2/n2)
    df=(v1/n1+v2/n2)**2/((v1/n1)**2/(n1-1)+(v2/n2)**2/(n2-1));t=-diff/se;crit=stats.t.ppf(.975,df)
    pooled=np.sqrt(((n1-1)*v1+(n2-1)*v2)/(n1+n2-2));g=diff/pooled*(1-3/(4*(n1+n2)-9))
    return dict(n_1=n1,mean_1=a.mean(),sd_1=np.sqrt(v1),n_2=n2,mean_2=b.mean(),sd_2=np.sqrt(v2),difference_group2_minus_group1=diff,ci95_low=diff-crit*se,ci95_high=diff+crit*se,welch_t=t,df=df,p_raw=2*stats.t.sf(abs(t),df),hedges_g=g)

def analyse(adult,root):
    frames={}
    def save(stage,name,df):
        path=root/stage/name;path.parent.mkdir(parents=True,exist_ok=True);df.to_csv(path,index=False,float_format='%.12g');frames[path.relative_to(root).as_posix()]=df
    outcomes=[f'{d}_score' for d in DIMS]+['BDSF_score']
    ind=pd.DataFrame([dict(variable=q[:3],**describe(adult[q])) for q in ITEMS])
    dim=pd.DataFrame([dict(variable=q.replace('_score',''),**describe(adult[q])) for q in outcomes])
    save('06_descriptive_analysis','06_02_BDSF_Kazakhstan_Indicator_Statistics.csv',ind)
    save('06_descriptive_analysis','06_03_BDSF_Kazakhstan_Dimension_Statistics.csv',dim)
    groups=[];sample=[]
    for typ,col in [('Language','source_language'),('Gender','gender_code'),('Age','age_group_code')]:
        for group,data in adult.groupby(col):
            groups.append(dict(group_type=typ,group_code=group,**describe(data.BDSF_score)))
            sample.append(dict(Characteristic=typ,Category=group,n=len(data),Percent=len(data)/len(adult)*100))
    save('06_descriptive_analysis','06_04_BDSF_Kazakhstan_Group_Descriptives.csv',pd.DataFrame(groups))
    counts=adult.BDSF_level.value_counts().reindex(['LOW','MEDIUM','HIGH'],fill_value=0)
    save('06_descriptive_analysis','06_05_BDSF_Kazakhstan_Level_Distribution.csv',pd.DataFrame(dict(level=counts.index,n=counts.values,percent=counts.values/len(adult)*100)))
    rel=[];diagnostics=[]
    for scale,cols in {**DIMS,'14_items':ITEMS,'4_dimensions':outcomes[:4]}.items():
        x=adult[cols];corr=x.corr();rel.append(dict(sample='Adults_all',scale=scale,n=len(x),items=len(cols),cronbach_alpha=alpha(x),standardized_alpha=alpha(x.apply(stats.zscore)),mean_interitem_correlation=corr.where(~np.eye(len(cols),dtype=bool)).stack().mean()))
        if scale in DIMS:
            for item in cols:
                rest=x.drop(columns=item);cit=x[item].corr(rest.sum(axis=1))
                diagnostics.append(dict(sample='Adults_all',scale=scale,item=item[:3],n=len(x),item_mean=x[item].mean(),item_sd=x[item].std(ddof=1),corrected_item_total_correlation=cit,alpha_if_deleted=alpha(rest),review_flag='REVIEW' if cit<.3 else 'OK'))
    reliability=pd.DataFrame(rel)
    save('07_reliability_analysis','07_02_BDSF_Kazakhstan_Scale_Reliability.csv',reliability)
    save('07_reliability_analysis','07_03_BDSF_Kazakhstan_Item_Diagnostics.csv',pd.DataFrame(diagnostics))
    x=adult[ITEMS].to_numpy(float);corr=np.corrcoef(x,rowvar=False);inv=np.linalg.inv(corr);partial=-inv/np.sqrt(np.outer(np.diag(inv),np.diag(inv)));np.fill_diagonal(partial,0);off=corr.copy();np.fill_diagonal(off,0)
    kmo=(off**2).sum()/((off**2).sum()+(partial**2).sum());chi=-(len(x)-1-(2*14+5)/6)*np.log(np.linalg.det(corr));df=91;bp=stats.chi2.sf(chi,df)
    save('08_construct_validity','08_02_BDSF_Kazakhstan_Item_Correlation_Matrix.csv',pd.DataFrame(corr,index=[q[:3] for q in ITEMS],columns=[q[:3] for q in ITEMS]).rename_axis('item').reset_index())
    rng=np.random.default_rng(SEED);obs=np.linalg.eigvalsh(corr)[::-1];sim=np.array([np.linalg.eigvalsh(np.corrcoef(rng.normal(size=x.shape),rowvar=False))[::-1] for _ in range(500)]);r95=np.quantile(sim,.95,axis=0)
    pa=pd.DataFrame(dict(component=range(1,15),observed_eigenvalue=obs,random_95th_eigenvalue=r95,retain=obs>r95))
    save('08_construct_validity','08_05_BDSF_Kazakhstan_Parallel_Analysis.csv',pa)
    fa=FactorAnalyzer(n_factors=4,rotation='oblimin',method='minres');fa.fit(x);expected=['D1']*4+['D2']*3+['D3']*4+['D4']*3
    similarity=np.array([np.abs(fa.loadings_[np.array(expected)==d]).mean(axis=0) for d in DIMS]);rows,cols=linear_sum_assignment(-similarity)
    load=fa.loadings_[:,cols];phi=fa.phi_[np.ix_(cols,cols)]
    # Oblique common variance is diag(L Phi L'), not the orthogonal sum of squares.
    communality=np.diag(load@phi@load.T);loads=[]
    for i,item in enumerate(ITEMS):
        ed=list(DIMS).index(expected[i]);absrow=np.abs(load[i]);dom=absrow.argmax();cross=max(np.delete(absrow,ed));intended=absrow[ed]
        align='STRONG' if dom==ed and intended>=.5 and cross<.3 else ('ACCEPTABLE' if dom==ed and intended>=.4 else 'REVIEW')
        loads.append(dict(item=item[:3],expected_dimension=expected[i],**{f'loading_{d}':load[i,j] for j,d in enumerate(DIMS)},communality=communality[i],uniqueness=1-communality[i],dominant_factor=list(DIMS)[dom],intended_abs_loading=intended,max_cross_loading=cross,alignment=align))
    save('08_construct_validity','08_03_BDSF_Kazakhstan_Factor_Loadings.csv',pd.DataFrame(loads))
    save('08_construct_validity','08_06_BDSF_Kazakhstan_Factor_Correlations.csv',pd.DataFrame(phi,index=DIMS,columns=DIMS).rename_axis('factor').reset_index())
    save('08_construct_validity','08_04_BDSF_Kazakhstan_Validity_Indices.csv',pd.DataFrame([('N',len(x)),('KMO',kmo),('Bartlett_chi2',chi),('Bartlett_df',91),('Bartlett_p',bp),('Parallel_retained_factors',int(pa.retain.sum()))],columns=['metric','value']))
    comparisons=[]
    for typ,col,g1,g2 in [('Language','source_language','KZ','RU'),('Gender','gender_code','MALE','FEMALE')]:
        rows=[]
        for outcome in outcomes: rows.append(dict(comparison_type=typ,outcome=outcome.replace('_score',''),group_1=g1,group_2=g2,**welch(adult.loc[adult[col].eq(g1),outcome],adult.loc[adult[col].eq(g2),outcome])))
        for row,p in zip(rows,holm([r['p_raw'] for r in rows])): row.update(p_holm=p,significant_holm_0_05='YES' if p<.05 else 'NO')
        comparisons+=rows
    save('09_group_comparison','09_02_BDSF_Kazakhstan_Two_Group_Comparisons.csv',pd.DataFrame(comparisons))
    save('11_figures_and_tables','11_07_Table_1_Sample_Characteristics.csv',pd.DataFrame(sample))
    main=[dict(Section='Descriptive',Metric=r.variable,Estimate=r['mean'],Secondary=r.sd,CI_low=r.ci95_low,CI_high=r.ci95_high,Note='Mean; secondary=SD') for _,r in dim.iterrows()]
    main += [dict(Section='Reliability',Metric=r['scale'],Estimate=r['cronbach_alpha'],Secondary=r['standardized_alpha'],CI_low=None,CI_high=None,Note='Cronbach alpha; secondary=standardized alpha') for r in rel]
    save('11_figures_and_tables','11_08_Table_2_Main_Results.csv',pd.DataFrame(main))
    out=root/'11_figures_and_tables';plt.style.use('seaborn-v0_8-whitegrid')
    fig,ax=plt.subplots(figsize=(8,4));ax.bar(dim.variable,dim['mean'],yerr=[dim['mean']-dim.ci95_low,dim.ci95_high-dim['mean']],capsize=4,color=['#287db2']*4+['#e67e22']);ax.set_ylim(0,100);ax.set_ylabel('Mean score (95% CI)');ax.set_title('Adult sample n=159');fig.tight_layout();fig.savefig(out/'11_02_Figure_1_Dimension_Means.png',dpi=240);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4));ax.bar(ind.variable,ind['mean'],color='#287db2');ax.set_ylim(0,100);ax.set_ylabel('Mean score');ax.set_title('14 indicators; n=159');fig.tight_layout();fig.savefig(out/'11_03_Figure_2_Indicator_Means.png',dpi=240);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4));ax.bar(reliability.scale,reliability.cronbach_alpha,color='#287db2');ax.axhline(.7,ls='--',c='gray');ax.set_ylim(0,1);ax.set_ylabel('Cronbach alpha');ax.set_title('Adult sample n=159');fig.tight_layout();fig.savefig(out/'11_04_Figure_3_Reliability.png',dpi=240);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,6));ax.grid(False);im=ax.imshow(load,cmap='RdYlGn',vmin=-1,vmax=1,aspect='auto');ax.set_xticks(range(4),DIMS);ax.set_yticks(range(14),[q[:3] for q in ITEMS]);
    for i in range(14):
        for j in range(4): ax.text(j,i,f'{load[i,j]:.2f}',ha='center',va='center',fontsize=9)
    fig.colorbar(im,ax=ax,label='Oblimin pattern loading');ax.set_title('Four-factor exploratory solution');fig.tight_layout();fig.savefig(out/'11_05_Figure_4_Factor_Loadings.png',dpi=240);plt.close(fig)
    summary=dict(n=len(adult),excluded_under18=21,language=adult.source_language.value_counts().to_dict(),gender=adult.gender_code.value_counts().to_dict(),age=adult.age_group_code.value_counts().to_dict(),BDSF=describe(adult.BDSF_score),dimensions={d:describe(adult[d+'_score']) for d in DIMS},levels=counts.to_dict(),reliability=rel,validity=dict(KMO=kmo,Bartlett_chi2=chi,Bartlett_df=91,Bartlett_p=bp,parallel_factors=int(pa.retain.sum()),review_items=[r['item'] for r in loads if r['alignment']=='REVIEW']),comparisons=comparisons,seed=SEED)
    return summary,frames
