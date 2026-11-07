$ErrorActionPreference = 'Stop'
$runtime = Join-Path $PSScriptRoot 'runtime'
$rscript = Join-Path $runtime 'R\bin\Rscript.exe'
$library = Join-Path $runtime 'R-library'
if (-not (Test-Path -LiteralPath $rscript)) {
    throw 'Bundled R runtime missing. Restore runtime/R from the section-17 package. This script does not modify the system R installation.'
}
New-Item -ItemType Directory -Force -Path $library | Out-Null
$env:BDSF_Q3_R_LIBRARY = $library
$env:BDSF_Q3_R_DOWNLOADS = Join-Path $runtime 'downloads'
$env:LC_ALL = 'C'
& $rscript --vanilla -e 'p <- list.files(Sys.getenv("BDSF_Q3_R_DOWNLOADS"), pattern="[.]zip$", full.names=TRUE); if(length(p)!=5) stop("Expected five pinned package archives"); install.packages(p, lib=Sys.getenv("BDSF_Q3_R_LIBRARY"), repos=NULL, type="win.binary")'
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
