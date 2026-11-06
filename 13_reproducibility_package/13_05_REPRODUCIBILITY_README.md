# BDSF Kazakhstan reviewer package

Windows x64. Extract/copy the complete package to a writable directory, then run from that directory:

```powershell
.\RUN_REVIEW.cmd
```

No Python, Node, Office, pip installation or internet connection is needed for computation.
The bundled Python 3.14.7 uses isolated relative imports; versions and licenses are in runtime/.
00_sources contains immutable raw sources; 13_reproducibility_package/inputs contains the question map, reference scoring matrix and original article.
The single pipeline regenerates CSV, 8 XLSX workbooks, 4 scientific figures, stage reports and the full updated DOCX article.
It keeps all 159 adults and excludes only the 21 respondents coded LT18 from the 180 collected responses.
Scores, dimensions, weights and thresholds are unchanged. Quality flags do not exclude adults.
24 baseline checks and 134 output checks must pass before results replace the published stages.
Failure returns a nonzero exit code. Previous outputs remain available; a failed build is marked FAILED_NOT_PUBLISHED.
Temporary builds, test copies and recoverable backups are confined to work/. They are not research outputs and are excluded from the reviewer ZIP.
The author retains pre-portability recovery files locally under work/recovery; backups are not published to GitHub.
review_materials/portability contains actual relocation, isolation and failure-test logs.
The relocation test uses a path containing spaces and Kazakh/Russian characters, an empty dependency PATH, Python isolated mode, and a Python audit guard blocking networking and author paths.
This is a test on the same Windows computer, not an independent clean VM or a second computer. Native-library file access is not OS-sandboxed.
The supported target is Windows x64; Linux/macOS/ARM64 are not tested.
DOCX is regenerated without Office. PDF page rendering is not part of this pipeline. Fonts used to create scientific PNGs and workbook previews are bundled.
The author-updated protocol PDF is included; its accepted SHA-256 is recorded in review_materials/protocol_revision.json. Its ethical scope and consent evidence have not been independently verified. Some generated manuscript/report wording still refers to pending insertion and requires a separate content review.
BDSF3.zip was not available and is not a required input. No ethics approval, date, number or consent statement has been invented.
The four-factor structure is not fully confirmed; packaging does not change this scientific limitation.

GitHub distribution: clone this repository or download and extract its source ZIP, then run RUN_REVIEW.cmd. The bundled runtime is tracked directly. Build/recovery folders and the redundant reviewer ZIP are excluded. Third-party components retain their own licenses in runtime; the root LICENSE does not replace those licenses.
