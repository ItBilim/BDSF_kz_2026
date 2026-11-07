# Нәтижелерді оқу және қайта іске қосу

Жоба түбірінен: `16_extended_experiments\RUN_EXPERIMENTS.cmd`. Бұрынғы RUN_REVIEW.cmd өзгермейді. Жаңа талдау алдында бастапқы кезеңдердің дайын CSV нәтижелері болуы керек.

Есептеудің қайталаулары cache ішінде кіріс SHA256, methods.py, seed және параметрлер бойынша тексеріліп қайта пайдаланылады. Негізгі нәтижелер қайта жасалмайды. Негізгі кодтағы логика өзгерсе cache жойып емес, жаңа fingerprint-пен жаңартылуы тиіс. Қате болса nonzero exit және VALIDATION_REPORT.json=FAIL; PASS болмаған нәтижелер аяқталған пакет емес.

ω — raw-score sum үшін; ω latent — polychoric модель үшін. ci95 — raw ω bootstrap percentile; бос ұяшық — бағалау жарамсыз не қолданылмайды. PA — параллель талдауда қалған компоненттер саны; RMSR — off-diagonal residual root mean square. n — қатысушылар/ұяшық саны, контекстке сәйкес. Барлық айырма N139−N159 немесе method2−method1, кесте атауына сәйкес.

M3 санаттары әдейі берілмейді; EXP03_Category_Transitions.csv файлындағы NOT_APPLICABLE жолы есептеу қатесі емес. Негізгі BDSF M1 ретінде сақталған.

## Файлдар және бағандар

### EXP01_Factor_Comparison.csv

`sample`, `correlation`, `source`, `factors`, `offdiag_rmsr`, `max_abs_residual`, `residual_pairs_gt_005`, `mean_communality`, `sum_common_variance`, `common_variance_share_of_total`, `heywood_items`, `efa_df`, `free_covariance_parameters`, `max_abs_factor_correlation`

Жол саны: 12.

### EXP01_Factor_Correlations.csv

`sample`, `correlation`, `factors`, `factor1`, `factor2`, `correlation_value`

Жол саны: 117.

### EXP01_Factor_Loadings.csv

`sample`, `correlation`, `factors`, `item`, `theoretical_dimension`, `loading_F1`, `communality`, `uniqueness`, `largest_loading_factor`, `second_abs_loading`, `loading_F2`, `loading_F3`, `loading_F4`

Жол саны: 168.

### EXP01_Legacy_Phi_Audit.csv

`sample`, `old_rmsr`, `corrected_rmsr`, `max_communality_change`, `reason`

Жол саны: 3.

### EXP01_Polychoric_Parallel.csv

`component`, `observed`, `null95`, `retain`

Жол саны: 14.

### EXP01_Residual_Correlations.csv

`sample`, `correlation`, `factors`, `item1`, `item2`, `residual`

Жол саны: 1092.

### EXP02_Absolute_Relative_Changes.csv

`metric`, `N159`, `N139`, `absolute_change`, `relative_change_percent`

Жол саны: 7.

### EXP02_Flag_Counts.csv

`flag`, `n`, `scope`

Жол саны: 3.

### EXP02_Flag_Overlap.csv

`flag_combination`, `n`

Жол саны: 3.

### EXP02_Loading_Changes.csv

`item`, `one_factor_N159`, `one_factor_N139`, `F1_difference`, `F2_difference`, `F3_difference`, `F4_difference`

Жол саны: 14.

### EXP02_Quality_Sensitivity.csv

`n`, `BDSF_mean`, `BDSF_sd`, `BDSF_ci_low`, `BDSF_ci_high`, `alpha`, `KMO`, `Bartlett_chi2`, `Bartlett_df`, `Bartlett_p`, `parallel_retained`, `omega_raw_sum`, `omega_M1_weighted`, `sample`, `role`

Жол саны: 2.

### EXP03_Category_Summary.csv

`method`, `category`, `n`, `n_changed`, `percent_changed`, `quadratic_weighted_kappa`

Жол саны: 6.

### EXP03_Category_Transitions.csv

`from_method`, `to_method`, `from_category`, `to_category`, `n`, `status`

Жол саны: 10.

### EXP03_Largest_Rank_Changes.csv

`case`, `M1`, `M2`, `M3`, `M1_rank`, `M2_rank`, `M3_rank`, `M2_minus_M1_rank`, `M3_minus_M1_rank`, `max_abs_change`

Жол саны: 15.

### EXP03_M3_Coefficients.csv

`item`, `regression_coefficient`, `training_mean`, `training_sd_population`

Жол саны: 14.

### EXP03_Rank_Stability.csv

`method1`, `method2`, `spearman`, `kendall_tau_b`, `mean_abs_rank_difference`, `max_abs_rank_difference`, `n_rank_changed`

Жол саны: 3.

### EXP03_Scoring_Comparison.csv

`case`, `M1`, `M2`, `M3`, `M1_rank`, `M2_rank`, `M3_rank`

Жол саны: 159.

### EXP04_Factor_Stability.csv

`comparison`, `factor`, `tucker_congruence`

Жол саны: 12.

### EXP04_Language_Comparison.csv

`language`, `scale`, `n`, `alpha_reused`, `alpha_ci_low_reused`, `alpha_ci_high_reused`, `omega_raw_sum`, `ci95_low`, `ci95_high`, `bootstrap_requested`, `bootstrap_valid`, `bootstrap_failed`, `ci_status`, `omega_polychoric_latent`, `polychoric_omega_ci`, `polychoric_onefactor_rmsr`

Жол саны: 10.

### PA_N139.csv

`component`, `observed`, `null95`, `retain`

Жол саны: 14.

### POLY_KZ_pairs.csv

`item1`, `item2`, `rho`, `zero_cells`, `cells`, `boundary`, `negative_loglik`

Жол саны: 91.

### POLY_KZ_raw.csv

`Q01`, `Q02`, `Q03`, `Q04`, `Q05`, `Q06`, `Q07`, `Q08`, `Q09`, `Q10`, `Q11`, `Q12`, `Q13`, `Q14`, `item`

Жол саны: 14.

### POLY_KZ_used.csv

`Q01`, `Q02`, `Q03`, `Q04`, `Q05`, `Q06`, `Q07`, `Q08`, `Q09`, `Q10`, `Q11`, `Q12`, `Q13`, `Q14`, `item`

Жол саны: 14.

### POLY_Matrix_Diagnostics.csv

`sample`, `min_eigenvalue_raw`, `smoothed`, `max_abs_adjustment`

Жол саны: 3.

### POLY_N159_pairs.csv

`item1`, `item2`, `rho`, `zero_cells`, `cells`, `boundary`, `negative_loglik`

Жол саны: 91.

### POLY_N159_raw.csv

`Q01`, `Q02`, `Q03`, `Q04`, `Q05`, `Q06`, `Q07`, `Q08`, `Q09`, `Q10`, `Q11`, `Q12`, `Q13`, `Q14`, `item`

Жол саны: 14.

### POLY_N159_used.csv

`Q01`, `Q02`, `Q03`, `Q04`, `Q05`, `Q06`, `Q07`, `Q08`, `Q09`, `Q10`, `Q11`, `Q12`, `Q13`, `Q14`, `item`

Жол саны: 14.

### POLY_RU_pairs.csv

`item1`, `item2`, `rho`, `zero_cells`, `cells`, `boundary`, `negative_loglik`

Жол саны: 91.

### POLY_RU_raw.csv

`Q01`, `Q02`, `Q03`, `Q04`, `Q05`, `Q06`, `Q07`, `Q08`, `Q09`, `Q10`, `Q11`, `Q12`, `Q13`, `Q14`, `item`

Жол саны: 14.

### POLY_RU_used.csv

`Q01`, `Q02`, `Q03`, `Q04`, `Q05`, `Q06`, `Q07`, `Q08`, `Q09`, `Q10`, `Q11`, `Q12`, `Q13`, `Q14`, `item`

Жол саны: 14.

### REUSE_LEDGER.csv

`task`, `action`

Жол саны: 19.

### STAGE01_Data_Audit.csv

`check`, `status`, `detail`

Жол саны: 15.

### STAGE_Completion.csv

`stage`, `status`, `scope`

Жол саны: 16.

### VALIDATION_CHECKS.csv

`check`, `status`, `detail`

Жол саны: 101.