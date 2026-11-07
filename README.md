# BDSF Kazakhstan research package

Updated 7 November 2026. Primary sample: 159 adults out of 180 responses; 21 under-18 responses excluded. The primary M1 score, question allocation and 25% dimension weights are unchanged. Quality flags remain in the primary sample; reduced samples are sensitivity analyses only.

## Latest results and code

| Section | Content | Report |
| --- | --- | --- |
| 15_language_psychometrics | Original bilingual alpha, item–total, exploratory loadings and response distributions | [Report](15_language_psychometrics/15_00_Report_KK.md) |
| 16_extended_experiments | Factor comparison, n=159/139 sensitivity, M1–M3 and language omega/polychoric analyses | [Research report](16_extended_experiments/results/FINAL_Research_Report_KK.md) |
| 17_q3_additional_experiments | Ordinal CFA using actual lavaan WLSMV, Top-K and individual quality-filter scenarios | [Research report](17_q3_additional_experiments/FINAL_Q3_REPORT.md) |
| 18_final_sensitivity_topk | Final RQ2/RQ3 tables, polychoric parallel analysis for S0–S4, overlap and changed selections | [Final results](18_final_sensitivity_topk/FINAL_RESULTS.md) |

Latest revised full manuscript: [ARTICLE_Updated_Four_Experiments.docx](16_extended_experiments/results/ARTICLE_Updated_Four_Experiments.docx). Section 17 and 18 findings are supplied as separate article-ready tables and reports; they have not been inserted into that DOCX.

## Run on Windows x64

Clone or download the complete repository to a writable directory. From its root:

```powershell
.\18_final_sensitivity_topk\RUN_FINAL.cmd
.\18_final_sensitivity_topk\RUN_FINAL.cmd --verify-repeat
```

The final RQ2/RQ3 module verifies and reuses existing results and checkpoints. To recompute its new polychoric permutations, use `--recompute-poly`.

Other sections can be run independently:

```powershell
.\16_extended_experiments\RUN_EXPERIMENTS.cmd
.\17_q3_additional_experiments\RUN_Q3.cmd
```

Section 17 reconstructs its analytical inputs from the original CSV files. Section 16 and 18 continue from validated existing stage outputs; they are not standalone packages. Do not run these commands concurrently. Each run writes its own output and validation records. File-hash checks intentionally fail when prerequisite files differ from the validated version; do not bypass those checks or copy old PASS reports over changed results.

The Python runtime is bundled under `runtime/`. R 4.5.2 and the lavaan dependencies are bundled under `17_q3_additional_experiments/runtime/`. Normal computation needs no system Python/R, Office or network. Package versions, source URLs and license notices are included. Linux, macOS, ARM64 and a separate clean computer have not been tested. Word was used only for manual manuscript layout QA, not by the computational pipelines.

## Preserved legacy pipeline and known correction

`RUN_REVIEW.cmd` preserves the original primary pipeline and its earlier stage outputs. It regenerates primary CSV/XLSX/figures/manuscript files and can replace the `12_article_results` directory, moving old files under local `work/` recovery. Preserve the current manuscript source before using that legacy regeneration path.

The legacy four-factor `factor_analyzer` Phi ordering defect remains in the unchanged original code. Corrected factor correlations, communalities and residuals are in section 16, with the explicit audit in [EXP01_Legacy_Phi_Audit.csv](16_extended_experiments/results/EXP01_Legacy_Phi_Audit.csv). Original loadings, primary BDSF scores and alpha were preserved. Section 17 CFA uses lavaan independently. Do not interpret an old legacy PASS report as validation of the corrected Phi ordering.

## Scientific and source scope

The updated signed departmental protocol PDF is included in `00_sources`. It documents the questionnaire/scoring methodology; it does not by itself establish individual informed consent or independent ethics-committee approval. No consent statement or approval record has been invented. The unavailable old BDSF3.zip is not a required input.

These are analyses of the same collected sample, not independent validation. A dominant factor, high reliability or high ranking correlations do not establish four separate constructs, measurement invariance or prediction of cybersecurity incidents. Follow each report's admissibility checks and limitations.

All required research code, sources, outputs, aggregate checkpoints and runtimes are tracked. Temporary `work/` directories, personal environments, build/recovery copies and installation logs are excluded. Section 18 selection lists contain technical IDs only; no new raw-response copies are published from working folders.

`PACKAGE_MANIFEST.csv` records SHA-256 and byte size for all tracked distribution files except itself. Historical section-specific manifests and validation snapshots are retained as records of their respective runs. Third-party components retain their own licenses; the root LICENSE does not replace them.

## Reviewer diagnostics and completed article tables

- [Section 19: full reviewer response](19_reviewer_experiments/REVIEWER_RESPONSE_KK.md): raw-data verification, complete CFA parameters, diagnostic ordinal language models, and a same-computer relocation check. Four-factor MG models were inadmissible; language equivalence was not established. This is not external independent replication.
- [Section 20: article-ready Table 9 and factor diagnostics](20_article_statistical_completion/BDSF_ARTICLE_TABLES.md): full S0–S4 statistics, missing SD values, loadings/SE/CI, correlations, diagnostics, and the 48-file audit.
- Run the table completion: `20_article_statistical_completion\RUN_TABLES.cmd`.

Sections 19 and 20 are separate supplements; the full manuscript DOCX has not been overwritten. The original n=159 scoring and previous results are unchanged. All new input hashes, numerical differences and audit limitations are included.
