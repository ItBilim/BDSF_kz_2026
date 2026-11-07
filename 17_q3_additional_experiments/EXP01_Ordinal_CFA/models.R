# Actual ordinal CFA in lavaan: WLSMV = DWLS estimates + robust SE + scaled/shifted test.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
input <- args[1]; output <- args[2]; lib <- args[3]
.libPaths(c(lib,.Library), include.site=FALSE)
library(lavaan,lib.loc=lib)
set.seed(20261107)
d <- read.csv(input, check.names=FALSE, fileEncoding='UTF-8-BOM')
items <- sprintf('Q%02d',1:14)
stopifnot(nrow(d)==159, identical(names(d),items), !anyNA(d))
stopifnot(all(as.matrix(d) %in% 0:4))
models <- list(OneFactor=paste('G =~',paste(items,collapse=' + ')),
 FourFactor='D1 =~ Q01 + Q02 + Q03 + Q04
 D2 =~ Q05 + Q06 + Q07
 D3 =~ Q08 + Q09 + Q10 + Q11
 D4 =~ Q12 + Q13 + Q14')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
results <- list(); log <- character()
for (name in names(models)) {
  warnings <- character()
  f <- withCallingHandlers(tryCatch(cfa(models[[name]],data=d,ordered=items,
       estimator='WLSMV',parameterization='delta',std.lv=TRUE,
       control=list(iter.max=10000)), error=function(e)e),
       warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  if (inherits(f,'error')) stop(paste(name,conditionMessage(f)))
  conv <- lavInspect(f,'converged')
  opts <- lavInspect(f,'options')
  post <- suppressWarnings(lavInspect(f,'post.check'))
  # Nonconvergence is reported as a failed analysis, not a valid fit.
  if (!conv) stop(paste(name,'did not converge'))
  fm <- suppressWarnings(fitMeasures(f,fit.measures='all',fm.args=list(rmsea.ci.level=.90)))
  get <- function(n) if(n %in% names(fm)) unname(fm[n]) else NA_real_
  theta <- lavInspect(f,'theta'); covlv <- lavInspect(f,'cov.lv'); corlv <- lavInspect(f,'cor.lv')
  st <- standardizedSolution(f); ld <- st[st$op=='=~',]
  negative_residuals <- sum(diag(theta)< -1e-8)
  out_of_range_loadings <- sum(abs(ld$est.std)>1+1e-8)
  min_latent_eigen <- min(eigen(covlv,symmetric=TRUE,only.values=TRUE)$values)
  excessive_correlations <- sum(abs(corlv[upper.tri(corlv)])>1+1e-8)
  admissible <- isTRUE(post) && negative_residuals==0 && out_of_range_loadings==0 && min_latent_eigen>0 && excessive_correlations==0
  results[[name]] <- data.frame(model=name,n=159,estimator='WLSMV',
    parameter_estimator=opts$estimator,test=paste(opts$test,collapse=';'),parameterization=opts$parameterization,
    standard_errors=opts$se,ordered_items=length(lavNames(f,'ov.ord')),
    converged=conv,admissible=admissible,post_check=post,
    chisq_scaled=get('chisq.scaled'),df_scaled=get('df.scaled'),p_scaled=get('pvalue.scaled'),
    CFI_scaled=get('cfi.scaled'),TLI_scaled=get('tli.scaled'),RMSEA_scaled=get('rmsea.scaled'),
    RMSEA90_low_scaled=get('rmsea.ci.lower.scaled'),RMSEA90_high_scaled=get('rmsea.ci.upper.scaled'),
    SRMR=get('srmr'),CFI_robust=get('cfi.robust'),TLI_robust=get('tli.robust'),
    RMSEA_robust=get('rmsea.robust'),RMSEA90_low_robust=get('rmsea.ci.lower.robust'),RMSEA90_high_robust=get('rmsea.ci.upper.robust'),
    negative_residual_variances=negative_residuals,standardized_loadings_gt1=out_of_range_loadings,
    min_latent_covariance_eigenvalue=min_latent_eigen,latent_correlations_gt1=excessive_correlations,
    max_abs_latent_correlation=if(nrow(corlv)>1)max(abs(corlv[upper.tri(corlv)])) else 0,
    npar=get('npar'),warnings=paste(unique(warnings),collapse=' | '))
  write.csv(ld,file.path(output,paste0(name,'_standardized_loadings.csv')),row.names=FALSE,na='')
  write.csv(corlv,file.path(output,paste0(name,'_factor_correlations.csv')),row.names=TRUE,na='')
  pe <- parameterEstimates(f,standardized=TRUE)
  write.csv(pe,file.path(output,paste0(name,'_parameters.csv')),row.names=FALSE,na='')
  write.csv(data.frame(item=items,residual_variance=diag(theta)),file.path(output,paste0(name,'_residual_variances.csv')),row.names=FALSE)
  samp <- lavInspect(f,'sampstat')
  write.csv(samp$cov,file.path(output,paste0(name,'_polychoric.csv')),row.names=TRUE)
  capture.output(summary(f,fit.measures=TRUE,standardized=TRUE),file=file.path(output,paste0(name,'_summary.txt')))
  log <- c(log,paste(name,unique(warnings),sep=': '))
}
result <- do.call(rbind,results);row.names(result)<-NULL
write.csv(result,file.path(output,'EXP01_CFA_RESULTS.csv'),row.names=FALSE,na='')
write.csv(result,file.path(output,'results.csv'),row.names=FALSE,na='')
writeLines(log,file.path(output,'warnings.txt'))
capture.output(sessionInfo(),file=file.path(output,'R_session_info.txt'))
writeLines(.libPaths(),file.path(output,'R_library_paths.txt'))
pk <- installed.packages(lib.loc=lib,fields='License')
write.csv(pk[,intersect(c('Package','Version','License'),colnames(pk)),drop=FALSE],file.path(output,'R_package_versions.csv'),row.names=FALSE)
cat('Ordinal CFA complete; WLSMV models:',nrow(result),'\n')
