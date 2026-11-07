"""Independent consistency tests on saved research outputs; no refitting of models."""
from pathlib import Path
import unittest,math
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]

class ResearchResults(unittest.TestCase):
    def test_primary_formula(self):
        d=pd.read_csv(ROOT/'work/reconstructed_inputs/05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv')
        self.assertEqual(len(d),159)
        a=d[[f'Q{i:02d}_score' for i in range(1,15)]].to_numpy()
        expected=(a[:,:4].sum(1)/4+a[:,4:7].sum(1)/3+a[:,7:11].sum(1)/4+a[:,11:].sum(1)/3)/4
        np.testing.assert_allclose(expected,d.BDSF_score,rtol=0,atol=1e-9)

    def test_topk_selection_is_deterministic(self):
        scores=pd.read_csv(ROOT/'EXP02_TopK_Stability/scores.csv');selected=pd.read_csv(ROOT/'EXP02_TopK_Stability/topk_selections.csv')
        # P-codes preserve original lexicographic ID order, allowing an independent tie-break check.
        for (method,k),group in selected.groupby(['method','K']):
            expected=scores.assign(key=scores[method].round(10)).sort_values(['key','case']).head(k)['case'].tolist()
            self.assertEqual(group.sort_values('selection_position')['case'].tolist(),expected)
            self.assertEqual(len(set(group['case'])),k)
            shuffled=scores.sample(frac=1,random_state=42).assign(key=lambda t:t[method].round(10)).sort_values(['key','case']).head(k)['case'].tolist()
            self.assertEqual(expected,shuffled)

    def test_overlap_arithmetic(self):
        s=pd.read_csv(ROOT/'EXP02_TopK_Stability/topk_selections.csv');r=pd.read_csv(ROOT/'EXP02_TOPK_RESULTS.csv')
        self.assertEqual(len(r),9)
        for _,row in r.iterrows():
            self.assertEqual(row.K,math.ceil(row.budget_percent/100*159))
            a=set(s[(s.method==row.method1)&(s.K==row.K)]['case']);b=set(s[(s.method==row.method2)&(s.K==row.K)]['case'])
            self.assertEqual(row.overlap_n,len(a&b));self.assertAlmostEqual(row.jaccard,len(a&b)/len(a|b),places=10)
            self.assertAlmostEqual(row.overlap_percent,100*len(a&b)/row.K,places=8)

    def test_quality_subsets_and_counts(self):
        d=pd.read_csv(ROOT/'work/reconstructed_inputs/05_scoring/05_01_BDSF_Kazakhstan_Scored_Dataset.csv',keep_default_na=False)
        results=pd.read_csv(ROOT/'EXP03_QUALITY_RESULTS.csv').set_index('scenario')
        for s,flag in [('S1','DUPLICATE_RESPONSE_PATTERN'),('S2','DUPLICATE_TIMESTAMP'),('S3','STRAIGHTLINE_ALL_14')]:
            sub=d[~d.quality_flags.str.contains(flag)]
            self.assertEqual(len(sub),results.loc[s,'n'])
            self.assertAlmostEqual(sub.BDSF_score.mean(),results.loc[s,'BDSF_mean'],places=8)
        self.assertEqual(results.loc['S4','n'],139)
        np.testing.assert_allclose(results.loc['S1',['BDSF_mean','alpha','omega_raw_sum']].to_numpy(float),results.loc['S3',['BDSF_mean','alpha','omega_raw_sum']].to_numpy(float))

    def test_cfa_options_and_model_structure(self):
        t=pd.read_csv(ROOT/'EXP01_CFA_RESULTS.csv')
        self.assertTrue(t.converged.all());self.assertTrue(t.ordered_items.eq(14).all())
        self.assertTrue(t.parameter_estimator.eq('DWLS').all());self.assertTrue(t.test.str.split(';').map(lambda v:'scaled.shifted' in v).all())
        self.assertTrue(t.standard_errors.str.contains('robust').all())
        for model,nf in [('OneFactor',1),('FourFactor',4)]:
            ld=pd.read_csv(ROOT/f'EXP01_Ordinal_CFA/{model}_standardized_loadings.csv')
            self.assertEqual(len(ld),14);self.assertEqual(ld.lhs.nunique(),nf)
            self.assertEqual(ld.rhs.nunique(),14)
            if nf==4:
                expected={'D1':['Q01','Q02','Q03','Q04'],'D2':['Q05','Q06','Q07'],'D3':['Q08','Q09','Q10','Q11'],'D4':['Q12','Q13','Q14']}
                for factor,items in expected.items():self.assertEqual(set(ld[ld.lhs==factor].rhs),set(items))

    def test_cfa_admissibility_matches_exported_parameters(self):
        t=pd.read_csv(ROOT/'EXP01_CFA_RESULTS.csv').set_index('model')
        for model in t.index:
            phi=pd.read_csv(ROOT/f'EXP01_Ordinal_CFA/{model}_factor_correlations.csv',index_col=0).to_numpy()
            theta=pd.read_csv(ROOT/f'EXP01_Ordinal_CFA/{model}_residual_variances.csv').residual_variance.to_numpy()
            ld=pd.read_csv(ROOT/f'EXP01_Ordinal_CFA/{model}_standardized_loadings.csv')['est.std'].to_numpy()
            valid=np.linalg.eigvalsh(phi).min()>0 and np.all(theta>=-1e-8) and np.all(np.abs(ld)<=1+1e-8)
            self.assertEqual(bool(t.loc[model,'admissible']),bool(valid))
            np.testing.assert_allclose(phi,phi.T,atol=1e-10);np.testing.assert_allclose(np.diag(phi),1,atol=1e-10)

    def test_fit_values_and_intervals(self):
        t=pd.read_csv(ROOT/'EXP01_CFA_RESULTS.csv')
        columns=['CFI_scaled','TLI_scaled','RMSEA_scaled','SRMR','RMSEA90_low_scaled','RMSEA90_high_scaled']
        self.assertTrue(np.isfinite(t[columns]).all().all())
        self.assertTrue((t.RMSEA90_low_scaled<=t.RMSEA_scaled).all())
        self.assertTrue((t.RMSEA_scaled<=t.RMSEA90_high_scaled).all())

if __name__=='__main__':unittest.main(verbosity=2)
