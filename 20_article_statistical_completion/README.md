# BDSF article statistical completion

Primary sample remains 159 adults; 21 under-18 responses are excluded. All previous results are read-only.

Run from the repository root: `20_article_statistical_completion\RUN_TABLES.cmd`.
The bundled Python and existing section-16 numerical methods are used. No R/Office/network is needed for this table extraction and numerical verification.

Deliverables:
- [Full S0–S4 statistics](BDSF_S0_S4_FULL_RESULTS.csv)
- [CFA loadings, SE and 95% CI](BDSF_CFA_LOADINGS.csv)
- [All six factor correlations](BDSF_FACTOR_CORRELATIONS.csv)
- [MG ordinal CFA diagnostics](BDSF_MG_CFA_DIAGNOSTICS.csv)
- [Reproducibility report](BDSF_REPRODUCIBILITY_REPORT.md)
- [Article-ready tables](BDSF_ARTICLE_TABLES.md)

Full machine-precision calculations are written with 17 significant digits; extracted CFA values retain the precision available in the original R CSV. Article tables use three decimal places; very small p-values use scientific notation.

Pearson PA is recalculated. Polychoric observed matrices and two deterministic permutations per unique scenario are recalculated; the verified existing 200-permutation checkpoints supply the full null distribution. CFA is extracted, not refitted. All numerical discrepancies and tolerances are explicitly listed in BDSF_NUMERICAL_DIFFERENCES.csv.

The 48-file audit verifies archived hashes and, when available, the original same-computer relocation directory. That work directory is not distributed. On other computers its absence is reported explicitly; it is not treated as evidence of a new independent replication. The original 48-file ledger and previous scientific results are retained under section 19. Neither the old DOCX nor old experiment files are overwritten.

Section 19 contains diagnostic language models with Q01 adjacent-category merging; primary scoring is unchanged. The invalid four-factor MG models must not be presented as confirmation. See the detailed report for limitations.
