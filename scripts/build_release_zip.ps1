param(
  [string] $OutputDir = ""
)

. "$PSScriptRoot\windows_common.ps1"

$appDir = Get-AppDir
$version = Get-AppVersion
if (-not $OutputDir) {
  $OutputDir = Join-Path $appDir "..\..\outputs"
}
$outputDir = [System.IO.Path]::GetFullPath($OutputDir)
$packageName = "TranscriptoScope_Windows_v$version"
$stagingRoot = Join-Path ([System.IO.Path]::GetTempPath()) "deseq2_dge_workbench_release"
$stagingDir = Join-Path $stagingRoot $packageName
$zipPath = Join-Path $outputDir "$packageName.zip"

if (Test-Path $stagingRoot) {
  Assert-PathInsideDirectory -Path $stagingRoot -ParentDirectory ([System.IO.Path]::GetTempPath())
  Remove-Item -LiteralPath $stagingRoot -Recurse -Force
}
New-Item -Path $stagingDir -ItemType Directory -Force | Out-Null
New-Item -Path $outputDir -ItemType Directory -Force | Out-Null

$excludeNames = @(
  ".git",
  ".Rproj.user",
  "__pycache__",
  "shiny-err.log",
  "shiny-out.log"
)

Get-ChildItem -Path $appDir -Force |
  Where-Object { $excludeNames -notcontains $_.Name } |
  ForEach-Object {
    Copy-Item -Path $_.FullName -Destination $stagingDir -Recurse -Force
  }

if (Test-Path $zipPath) {
  Remove-Item -LiteralPath $zipPath -Force
}

Compress-Archive -Path (Join-Path $stagingRoot $packageName) -DestinationPath $zipPath -CompressionLevel Optimal

Write-Host $zipPath
