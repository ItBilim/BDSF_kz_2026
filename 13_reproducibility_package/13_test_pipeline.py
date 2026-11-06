"""Boundary and provenance checks; all mutations occur in disposable test fixtures."""
from pathlib import Path
import importlib.util,shutil,unittest,uuid
import pandas as pd
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('pipeline',HERE/'13_pipeline.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.base=HERE.parent/'work/tests'/('BDSF_test_'+uuid.uuid4().hex[:10]);self.base.mkdir(parents=True);self.project=self.base/'project';self.project.mkdir();shutil.copytree(HERE.parent/'00_sources',self.project/'00_sources');self.build=self.base/'build';self.build.mkdir();self.q=pd.read_csv(HERE/'inputs/question_map.csv')
    def tearDown(self):
        # Fixtures remain inside work/tests for traceability; no deletion of user data.
        pass
    def mutate(self,col,value):
        path=next((self.project/'00_sources').glob('*Kazakh_Responses.csv'));df=pd.read_csv(path,dtype=str);df.iloc[0,col]=value;df.to_csv(path,index=False)
    def test_unknown_age_fails(self):
        self.mutate(3,'UNKNOWN AGE')
        with self.assertRaisesRegex(ValueError,'Unknown category'):m.harmonize(self.project,self.q,self.build)
    def test_unknown_answer_fails(self):
        self.mutate(4,'UNREGISTERED RESPONSE')
        with self.assertRaisesRegex(ValueError,'Unknown category'):m.harmonize(self.project,self.q,self.build)
    def test_blank_answer_fails(self):
        self.mutate(4,'')
        with self.assertRaisesRegex(ValueError,'missing response'):m.harmonize(self.project,self.q,self.build)
    def test_missing_column_fails(self):
        path=next((self.project/'00_sources').glob('*Kazakh_Responses.csv'));df=pd.read_csv(path);df.drop(columns=df.columns[-1]).to_csv(path,index=False)
        with self.assertRaisesRegex(ValueError,'input columns'):m.harmonize(self.project,self.q,self.build)
    def test_partition_and_raw_preservation(self):
        before={p.name:m.sha(p) for p in (self.project/'00_sources').iterdir()};adult,frames=m.harmonize(self.project,self.q,self.build)
        raw=frames['03_harmonization/03_01_BDSF_Kazakhstan_Harmonized_Dataset.csv'];ex=frames['04_cleaning/04_03_BDSF_Kazakhstan_Exclusion_Log.csv']
        self.assertEqual(set(adult.respondent_id),set(raw.loc[raw.age_group_code.ne('LT18'),'respondent_id']))
        self.assertEqual(set(ex.respondent_id),set(raw.loc[raw.age_group_code.eq('LT18'),'respondent_id']))
        self.assertEqual(len(adult),159);self.assertEqual(len(ex),21);self.assertEqual(adult.quality_flags.ne('').sum(),20)
        for p in (self.project/'00_sources').iterdir():self.assertEqual(m.sha(p),before[p.name])
        for lang,name in [('KZ','Kazakh'),('RU','Russian')]:
            source=pd.read_csv(next((self.project/'00_sources').glob(f'*{name}_Responses.csv')),dtype=str,keep_default_na=False)
            subset=raw.loc[raw.source_language.eq(lang)]
            for i in range(14):self.assertEqual(subset[f'Q{i+1:02d}_raw'].tolist(),source.iloc[:,i+4].tolist())
    def test_score_boundaries(self):
        scores=pd.Series([0,33.32,33.33,66.66,66.67,100])
        self.assertEqual(m.core.classify(scores).astype(str).tolist(),['LOW','LOW','MEDIUM','MEDIUM','HIGH','HIGH'])
    def test_failed_build_does_not_publish(self):
        self.mutate(3,'unknown');marker=self.project/'14_final_quality_control';marker.mkdir();(marker/'previous.txt').write_text('unchanged')
        with self.assertRaises(ValueError):m.harmonize(self.project,self.q,self.build)
        self.assertEqual((marker/'previous.txt').read_text(),'unchanged');self.assertFalse((self.project/'05_scoring').exists())
    def test_reference_dimensions_and_scores(self):m.check_reference(self.q)

if __name__=='__main__':unittest.main(verbosity=2)
