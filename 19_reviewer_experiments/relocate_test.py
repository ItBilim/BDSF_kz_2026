from pathlib import Path
import os,sys,shutil,subprocess,hashlib,json,csv
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
DST=ROOT/'work'/'Тәуелсіз тексеру 19'/'Көшірілген жоба 03'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert DST.resolve().is_relative_to((ROOT/'work').resolve())
    if '--verify-existing' in sys.argv:
        assert json.loads((DST/'19_reviewer_experiments/results/STATUS.json').read_text())['status']=='PASS'
        return verify()
    DST.mkdir(parents=True,exist_ok=True)
    folders=['00_sources','05_scoring','runtime','19_reviewer_experiments/runtime','17_q3_additional_experiments/runtime/R','17_q3_additional_experiments/runtime/R-library']
    for folder in folders:
        shutil.copytree(ROOT/folder,DST/folder,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    files=['13_reproducibility_package/inputs/question_map.csv','17_q3_additional_experiments/EXP01_CFA_RESULTS.csv','18_final_sensitivity_topk/RQ2_SENSITIVITY.csv','19_reviewer_experiments/run.py','19_reviewer_experiments/models.R','19_reviewer_experiments/RUN_REVIEWER.cmd']
    for rel in files:
        p=DST/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,p)
    # No old section-19 result or precomputed ordinal data is copied.
    out=DST/'19_reviewer_experiments/results'
    if out.exists():
        raise RuntimeError('Use a fresh relocation directory; existing results would contaminate this test')
    env=os.environ.copy()
    for k in ['PYTHONPATH','PYTHONHOME','R_HOME','R_LIBS','R_LIBS_USER','R_LIBS_SITE']:env.pop(k,None)
    env['PATH']=str(Path(os.environ['SystemRoot'])/'System32');env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    proc=subprocess.run([str(DST/'runtime/python/python.exe'),'-I','-X','utf8',str(DST/'19_reviewer_experiments/run.py')],cwd=DST,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (HERE/'results/RELOCATION_EXECUTION.log').write_text(proc.stdout+'\n'+proc.stderr,encoding='utf-8')
    assert proc.returncode==0,proc.stderr
    verify()
def verify():
    out=DST/'19_reviewer_experiments/results'
    original=HERE/'results';checks=[]
    for p in original.rglob('*.csv'):
        rel=p.relative_to(original)
        if rel.name in ['RELOCATION_HASHES.csv']:continue
        q=out/rel
        assert q.exists(),str(rel)
        if rel.name=='R_PACKAGES.csv':
            # Locale changes package-list row ordering, not package versions.
            def packages(path):
                with path.open(encoding='utf-8-sig',newline='') as f:return sorted(tuple(sorted(r.items())) for r in csv.DictReader(f))
            assert packages(p)==packages(q),'Package versions/license records differ'
            continue
        checks.append({'file':rel.as_posix(),'original_sha256':sha(p),'relocated_sha256':sha(q),'identical':sha(p)==sha(q)})
    assert all(c['identical'] for c in checks)
    with (original/'RELOCATION_HASHES.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(checks[0]));w.writeheader();w.writerows(checks)
    report=dict(status='PASS',csv_files_identical=len(checks),package_metadata_equal_after_row_sort=True,package_metadata_byte_difference='Locale-dependent row order only',fresh_raw_reconstruction=True,isolated_python=True,system_dependency_path_removed=True,bundled_R_and_Python=True,relocated_to=str(DST),network_os_blocked=False,separate_computer=False,independent_researcher=False,scope='Section 19 raw scoring, sensitivity and 18 CFA models; not the complete legacy manuscript pipeline')
    (original/'RELOCATION_STATUS.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
