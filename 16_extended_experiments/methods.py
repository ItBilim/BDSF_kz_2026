"""Numerical methods for supplementary BDSF analyses; no source files are edited."""
import warnings
import numpy as np
from scipy import stats, optimize, special, integrate
from factor_analyzer import FactorAnalyzer
import factor_analyzer.factor_analyzer as fam
from sklearn.utils.validation import check_array

def compatible(*a, **kw):
    if 'force_all_finite' in kw: kw['ensure_all_finite']=kw.pop('force_all_finite')
    return check_array(*a, **kw)
fam.check_array=compatible
ITEMS=[f'Q{i:02d}' for i in range(1,15)]
DIMS={'D1':list(range(4)), 'D2':list(range(4,7)), 'D3':list(range(7,11)), 'D4':list(range(11,14))}

def alpha(x):
    return x.shape[1]/(x.shape[1]-1)*(1-x.var(axis=0,ddof=1).sum()/x.sum(axis=1).var(ddof=1))

def diagnostics(r,n):
    inv=np.linalg.inv(r); p=-inv/np.sqrt(np.outer(np.diag(inv),np.diag(inv)))
    np.fill_diagonal(p,0); off=r.copy();np.fill_diagonal(off,0)
    k=r.shape[0]; chi=-(n-1-(2*k+5)/6)*np.linalg.slogdet(r)[1]
    return dict(KMO=float(np.square(off).sum()/(np.square(off).sum()+np.square(p).sum())),Bartlett_chi2=float(chi),Bartlett_df=k*(k-1)//2,Bartlett_p=float(stats.chi2.sf(chi,k*(k-1)//2)))

def fit(r,k):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        fa=FactorAnalyzer(n_factors=k,rotation='oblimin' if k>1 else None,method='minres',is_corr_matrix=True)
        fa.fit(r)
    l=fa.loadings_.copy(); phi=fa.phi_ if k>1 else np.eye(1)
    if k==1 and l.sum()<0:l=-l
    return dict(L=l,Phi=phi,warnings=[str(v.message) for v in w])

def align(model, reference):
    l=model['L']; sim=reference.T@l/(np.linalg.norm(reference,axis=0)[:,None]*np.linalg.norm(l,axis=0)[None,:])
    _,order=optimize.linear_sum_assignment(-np.abs(sim)); sign=np.where(sim[np.arange(len(order)),order]>=0,1,-1)
    return dict(L=l[:,order]*sign,Phi=model['Phi'][np.ix_(order,order)]*np.outer(sign,sign),warnings=model.get('warnings',[]),order=order.tolist(),sign=sign.tolist())

def fit_metrics(r,m):
    l=m['L']; phi=m['Phi']; common=l@phi@l.T; h=np.diag(common); off=np.triu_indices(len(r),1)
    residual=r-common; vals=residual[off]; k=l.shape[1]; p=len(r)
    return dict(factors=k,offdiag_rmsr=float(np.sqrt(np.mean(vals**2))),max_abs_residual=float(np.max(np.abs(vals))),residual_pairs_gt_005=int((np.abs(vals)>.05).sum()),mean_communality=float(h.mean()),sum_common_variance=float(h.sum()),common_variance_share_of_total=float(h.mean()),heywood_items=int(((h<0)|(h>=1)).sum()),efa_df=int(((p-k)**2-p-k)/2),free_covariance_parameters=int(p*(k+1)-k*(k-1)/2),max_abs_factor_correlation=float(np.max(np.abs(phi-np.eye(k)))))

def omega(m,sd=None,weights=None):
    # Model-implied omega for a sum: a' Lambda Phi Lambda' a / a' Sigma_model a.
    # One-factor model here; observed raw-score omega uses a = item SD.
    l=m['L'];h=np.diag(l@m['Phi']@l.T)
    if np.any(h>=1) or np.any(h<0) or not np.isfinite(h).all():return np.nan
    a=np.ones(len(l)) if sd is None else np.array(sd)
    if weights is not None:a=a*np.array(weights)
    common=float(a@l@m['Phi']@l.T@a); error=float(np.sum(a*a*(1-h)))
    return common/(common+error)

# Bivariate normal CDF via Plackett's integral over rho, deterministic Gauss-Legendre quadrature.
NODES,WEIGHTS=np.polynomial.legendre.leggauss(96)
def bvncdf_grid(tx,ty,rho):
    h=np.asarray(tx)[:,None];k=np.asarray(ty)[None,:]; out=special.ndtr(h)*special.ndtr(k)
    mask=np.isfinite(h)&np.isfinite(k)
    hh=np.broadcast_to(h,out.shape)[mask,None]; kk=np.broadcast_to(k,out.shape)[mask,None]
    t=rho*(NODES+1)/2; v=1-t*t
    density=np.exp(-(hh*hh-2*t*hh*kk+kk*kk)/(2*v))/(2*np.pi*np.sqrt(v))
    out[mask]+=rho/2*(density@WEIGHTS)
    return out

def poly_pair(x,y):
    _,ix=np.unique(x,return_inverse=True);_,iy=np.unique(y,return_inverse=True)
    table=np.zeros((ix.max()+1,iy.max()+1),int);np.add.at(table,(ix,iy),1)
    if min(table.shape)<2:raise ValueError('Constant item cannot have a polychoric correlation')
    tx=np.r_[-np.inf,stats.norm.ppf(np.cumsum(table.sum(1))[:-1]/table.sum()),np.inf]
    ty=np.r_[-np.inf,stats.norm.ppf(np.cumsum(table.sum(0))[:-1]/table.sum()),np.inf]
    def objective(r):
        c=bvncdf_grid(tx,ty,r);p=np.diff(np.diff(c,axis=0),axis=1)
        return -float(np.sum(table*np.log(np.maximum(p,1e-14))))
    opt=optimize.minimize_scalar(objective,bounds=(-.995,.995),method='bounded',options={'xatol':1e-7})
    if not opt.success:raise RuntimeError('Polychoric optimization failed')
    return float(opt.x),dict(zero_cells=int((table==0).sum()),cells=int(table.size),boundary=bool(abs(opt.x)>.994),negative_loglik=float(opt.fun))

def smooth(r):
    ev,v=np.linalg.eigh(r); minimum=float(ev.min())
    if minimum>=1e-8:return r.copy(),dict(min_eigenvalue_raw=minimum,smoothed=False,max_abs_adjustment=0.)
    repaired=(v*np.maximum(ev,1e-6))@v.T; d=np.sqrt(np.diag(repaired));repaired=repaired/np.outer(d,d)
    return repaired,dict(min_eigenvalue_raw=minimum,smoothed=True,max_abs_adjustment=float(np.max(np.abs(repaired-r))))

def polychoric(x):
    r=np.eye(x.shape[1]);detail=[]
    for i in range(x.shape[1]):
        for j in range(i+1,x.shape[1]):
            r[i,j],note=poly_pair(x[:,i],x[:,j]);r[j,i]=r[i,j]
            detail.append(dict(item1=ITEMS[i],item2=ITEMS[j],rho=r[i,j],**note))
    adjusted,note=smooth(r)
    return r,adjusted,detail,note

def method_tests():
    # Independent adaptive quadrature checks the numerical CDF, including high rho.
    errors=[]
    for h,k,r in [(0.,0.,.8),(.3,-.7,-.6),(1.2,.4,.99),(-1.1,-.5,.995)]:
        expected=integrate.quad(lambda z:stats.norm.pdf(z)*stats.norm.cdf((k-r*z)/np.sqrt(1-r*r)),-np.inf,h,epsabs=1e-11)[0]
        errors.append(abs(bvncdf_grid([h],[k],r)[0,0]-expected))
    assert max(errors)<1e-7,errors
    rng=np.random.default_rng(911);z=rng.multivariate_normal([0,0],[[1,.6],[.6,1]],size=12000)
    ordinal=np.digitize(z,[-1.,-.3,.3,1.]);estimated,_=poly_pair(ordinal[:,0],ordinal[:,1])
    assert abs(estimated-.6)<.035,estimated
    assert abs(omega(dict(L=np.full((4,1),.7),Phi=np.eye(1)))-7.84/9.88)<1e-12
    return dict(cdf_max_error=max(errors),synthetic_polychoric_true=.6,synthetic_polychoric_estimate=estimated,omega_identity='PASS')
