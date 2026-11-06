param([string]$ProjectRoot=(Split-Path $PSScriptRoot -Parent),[switch]$BuildOnly)
$ErrorActionPreference='Stop'
try {
    $ProjectRoot=(Resolve-Path -LiteralPath $ProjectRoot).Path
    $python=Join-Path $ProjectRoot 'runtime\python\python.exe'
    if(-not(Test-Path -LiteralPath $python)){throw 'Bundled runtime missing: restore the complete reviewer package.'}
    $arguments=@('-I','-X','utf8',(Join-Path $ProjectRoot '13_reproducibility_package\13_pipeline.py'),'--project-root',$ProjectRoot)
    if($BuildOnly){$arguments+='--build-only'}
    & $python @arguments
    exit $LASTEXITCODE
} catch { Write-Error $_; exit 1 }
