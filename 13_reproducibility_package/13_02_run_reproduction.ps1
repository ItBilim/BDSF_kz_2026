param([string]$ProjectRoot=(Split-Path $PSScriptRoot -Parent))
& (Join-Path $PSScriptRoot '13_run_pipeline.ps1') -ProjectRoot $ProjectRoot
exit $LASTEXITCODE
