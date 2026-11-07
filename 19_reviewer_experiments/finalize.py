from pathlib import Path
import hashlib,csv,json,subprocess,sys
H=Path(__file__).resolve().parent;R=H.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=R/'work/github_publish/PACKAGE_MANIFEST.csv'
rows=list(csv.DictReader(manifest.open(encoding='utf-8')))
checked=[];mismatch=[]
for row in rows:
    name=row['path'];part=name.split('/')[0]
    if not(part[:2].isdigit() and int(part[:2])<=18):continue
    if name=='13_reproducibility_package/13_05_REPRODUCIBILITY_README.md':continue
    p=R/name
    if not p.exists() or sha(p)!=row['sha256']:mismatch.append(name)
    checked.append(name)
assert not mismatch,mismatch
rel=json.loads((H/'results/RELOCATION_STATUS.json').read_text(encoding='utf-8'))
assert rel['status']=='PASS'
subprocess.run([sys.executable,str(H/'write_report.py')],check=True)
(H/'results/FINAL_VALIDATION.json').write_text(json.dumps(dict(status='PASS',prior_published_files_unchanged=len(checked),explicit_exception='Existing local 13_05 README differs from previously published distribution; not modified by section 19',relocation=rel),ensure_ascii=False,indent=2),encoding='utf-8')
with (H/'MANIFEST.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f);w.writerow(['path','bytes','sha256'])
    for p in sorted(H.rglob('*')):
        if p.is_file() and p.name!='MANIFEST.csv' and 'work' not in p.relative_to(H).parts and '__pycache__' not in p.parts:
            w.writerow([p.relative_to(H).as_posix(),p.stat().st_size,sha(p)])
print('FINAL PASS; prior files unchanged:',len(checked))
