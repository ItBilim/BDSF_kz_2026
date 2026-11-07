"""Fix factor_analyzer 0.5.1 Phi ordering without altering its installation.

FactorAnalyzer.fit sorts loadings_ and structure_, but not phi_.
Since structure_ = loadings_ @ Phi in consistent factor coordinates,
recover Phi from that identity, verify it, and preserve loadings verbatim.
This changes only k>1 models; cached one-factor omega/polychoric matrices are unaffected.
"""
import warnings
import numpy as np
from factor_analyzer import FactorAnalyzer
from methods import fit as original_fit

def fit(r,k):
    if k==1:return original_fit(r,k)
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        fa=FactorAnalyzer(n_factors=k,rotation='oblimin',method='minres',is_corr_matrix=True);fa.fit(r)
    l=fa.loadings_.copy();phi=np.linalg.lstsq(l,fa.structure_,rcond=None)[0]
    assert np.allclose(l@phi,fa.structure_,atol=1e-9)
    assert np.allclose(phi,phi.T,atol=1e-9) and np.allclose(np.diag(phi),1,atol=1e-9)
    phi=(phi+phi.T)/2
    assert np.linalg.eigvalsh(phi).min()>0
    return dict(L=l,Phi=phi,warnings=[str(w.message) for w in ws],phi_order_correction=float(np.max(np.abs(phi-fa.phi_))))
