from pathlib import Path
import subprocess,sys
raise SystemExit(subprocess.call([sys.executable,str(Path(__file__).resolve().parents[1]/'run_all.py'),'--only','3']))
