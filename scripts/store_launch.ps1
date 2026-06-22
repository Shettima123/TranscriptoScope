$ErrorActionPreference = "Stop"

$appDir = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot ".."))
$runtimeDir = Join-Path $appDir "runtime"
$bundledRscript = Join-Path $runtimeDir "R\bin\Rscript.exe"
$libraryRoot = Join-Path $runtimeDir "R-library"
$logDir = Join-Path $env:LOCALAPPDATA "TranscriptoScope\logs"
$outLog = Join-Path $logDir "store-shiny-out.log"
$errLog = Join-Path $logDir "store-shiny-err.log"

if (-not (Test-Path $logDir)) {
  New-Item -Path $logDir -ItemType Directory -Force | Out-Null
}

if (-not (Test-Path $bundledRscript)) {
  throw "Bundled Rscript was not found at $bundledRscript. Reinstall TranscriptoScope."
}

$libraryCandidates = @()
if (Test-Path $libraryRoot) {
  $libraryCandidates = Get-ChildItem -Path $libraryRoot -Directory -ErrorAction SilentlyContinue |
    Sort-Object Name -Descending
}
if ($libraryCandidates.Count -lt 1) {
  throw "Bundled R package library was not found at $libraryRoot. Reinstall TranscriptoScope."
}

$env:R_LIBS_USER = $libraryCandidates[0].FullName
$env:R_LIBS_SITE = ""
$env:PATH = @(
  (Join-Path $runtimeDir "R\bin"),
  (Join-Path $runtimeDir "R\bin\x64"),
  $env:PATH
) -join ";"

Set-Location $appDir
& $bundledRscript (Join-Path $appDir "scripts\launch_app.R") 1> $outLog 2> $errLog
if ($LASTEXITCODE -ne 0) {
  throw "TranscriptoScope exited with code $LASTEXITCODE. See $errLog"
}
