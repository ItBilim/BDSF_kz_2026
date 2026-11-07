args <- commandArgs(TRUE)
root <- args[1]; out <- args[2]
.libPaths(c(file.path(root,'19_reviewer_experiments/runtime/R-library'),file.path(root,'17_q3_additional_experiments/runtime/R-library'),.Library),include.site=FALSE)
library(lavaan);library(semTools)
d <- read.csv(file.path(out,'ordinal_input.csv'),fileEncoding='UTF-8-BOM')
items <- sprintf('Q%02d',1:14)
models <- list(OneFactor=paste('G =~',paste(items,collapse=' + ')),FourFactor='D1 =~ Q01 + Q02 + Q03 + Q04
D2 =~ Q05 + Q06 + Q07
D3 =~ Q08 + Q09 + Q10 + Q11
D4 =~ Q12 + Q13 + Q14')
dir.create(file.path(out,'parameters'),showWarnings=FALSE)
rows <- list(); fits <- list(); warnings.all <- list(); comparisons <- list()
runfit <- function(key,model,data,group=NULL) {
  ww <- character()
  f <- withCallingHandlers(tryCatch(cfa(model,data=data,ordered=items,group=group,estimator='WLSMV',parameterization='theta',control=list(iter.max=10000)),error=function(e)e),warning=function(w){ww <<- c(ww,conditionMessage(w));invokeRestart('muffleWarning')})
  warnings.all[[key]] <<- ww
  if(inherits(f,'error')) {rows[[key]] <<- data.frame(key=key,converged=FALSE,error=conditionMessage(f));return(NULL)}
  conv <- lavInspect(f,'converged')
  if(!conv) {rows[[key]] <<- data.frame(key=key,converged=FALSE,error='Nonconvergence');return(NULL)}
  fm <- suppressWarnings(fitMeasures(f))
  get <- function(k) if(k %in% names(fm)) unname(fm[k]) else NA_real_
  st <- standardizedSolution(f)
  cv <- lavInspect(f,'cov.lv');if(is.matrix(cv)) cv <- list(cv)
  co <- lavInspect(f,'cor.lv');if(is.matrix(co)) co <- list(co)
  mineig <- min(sapply(cv,function(m)min(eigen(m,symmetric=TRUE)$values)))
  maxcor <- max(c(0,unlist(lapply(co,function(m)abs(m[upper.tri(m)])))))
  adm <- isTRUE(suppressWarnings(lavInspect(f,'post.check'))) && mineig>0 && maxcor<1 && all(abs(st$est.std[st$op=='=~'])<=1)
  rows[[key]] <<- data.frame(key=key,converged=conv,admissible=adm,n=nrow(data),chisq=get('chisq.scaled'),df=get('df.scaled'),p=get('pvalue.scaled'),CFI=get('cfi.scaled'),TLI=get('tli.scaled'),RMSEA=get('rmsea.scaled'),RMSEA_low=get('rmsea.ci.lower.scaled'),RMSEA_high=get('rmsea.ci.upper.scaled'),SRMR=get('srmr'),CFI_robust=get('cfi.robust'),TLI_robust=get('tli.robust'),RMSEA_robust=get('rmsea.robust'),min_latent_eigen=mineig,max_factor_correlation=maxcor,error='')
  write.csv(parameterEstimates(f,standardized=TRUE,ci=TRUE),file.path(out,'parameters',paste0(key,'_parameters.csv')),row.names=FALSE)
  write.csv(st,file.path(out,'parameters',paste0(key,'_standardized.csv')),row.names=FALSE)
  capture.output(summary(f,fit.measures=TRUE,standardized=TRUE),file=file.path(out,'parameters',paste0(key,'_summary.txt')))
  fits[[key]] <<- f
  f
}
# Full-sample and sensitivity CFA, with the same prespecified item allocation.
for(s in paste0('S',0:4)) for(m in names(models)) {
  dd <- d[d[[s]]==1,]
  # std.lv identification is explicit in the generated single-group syntax.
  sy <- measEq.syntax(configural.model=models[[m]],data=dd,ordered=items,parameterization='theta',ID.fac='std.lv',ID.cat='Wu.Estabrook.2016')
  runfit(paste(s,m,sep='_'),as.character(sy),dd)
}
# Proper ordinal sequence: thresholds, then loadings, then intercepts.
# Diagnostic only: original RU Q01 category 0 is empty. Both groups receive
# the identical adjacent-category merge; this does not modify primary scores.
d$Q01 <- ifelse(d$Q01<=1,1,d$Q01)
levels <- list(configural=character(),thresholds='thresholds',metric=c('thresholds','loadings'),scalar=c('thresholds','loadings','intercepts'))
for(m in names(models)) {
  prev <- NULL
  for(lv in names(levels)) {
    key <- paste('MG',m,lv,sep='_')
    sy <- measEq.syntax(configural.model=models[[m]],data=d,ordered=items,group='language',parameterization='theta',ID.fac='std.lv',ID.cat='Wu.Estabrook.2016',group.equal=levels[[lv]])
    writeLines(as.character(sy),file.path(out,'parameters',paste0(key,'_syntax.txt')))
    f <- runfit(key,as.character(sy),d,group='language')
    if(!is.null(f) && !is.null(prev)) {
      z <- tryCatch(suppressWarnings(lavTestLRT(fits[[prev]],f,method='satorra.2000')),error=function(e)e)
      if(!inherits(z,'error')) {zz <- as.data.frame(z);zz$model <- rownames(zz);write.csv(zz,file.path(out,paste0(key,'_LRT.csv')),row.names=FALSE)}
      a <- rows[[prev]];b <- rows[[key]]
      comparisons[[key]] <- data.frame(model=m,step=lv,delta_CFI=b$CFI-a$CFI,delta_RMSEA=b$RMSEA-a$RMSEA,delta_SRMR=b$SRMR-a$SRMR,both_admissible=isTRUE(a$admissible)&&isTRUE(b$admissible))
    }
    if(!is.null(f)) prev <- key
  }
}
cols <- unique(unlist(lapply(rows,names)))
rows <- lapply(rows,function(x){for(n in setdiff(cols,names(x))) x[[n]] <- NA;x[,cols]})
write.csv(do.call(rbind,rows),file.path(out,'CFA_ALL_MODELS.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,comparisons),file.path(out,'INVARIANCE_DELTAS.csv'),row.names=FALSE,na='')
writeLines(unlist(lapply(names(warnings.all),function(k)paste(k,warnings.all[[k]],sep=': '))),file.path(out,'warnings.txt'))
pk <- installed.packages(lib.loc=.libPaths(),fields='License')
write.csv(pk[,c('Package','Version','License')],file.path(out,'R_PACKAGES.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_SESSION.txt'))

