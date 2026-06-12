. "$PSScriptRoot\windows_common.ps1"

$appDir = Get-AppDir
$installScript = Join-Path $appDir "scripts\install_packages.R"

Write-Host "Checking and installing required R packages..."
Write-Host "This can take several minutes the first time, especially for DESeq2."

Set-Location $appDir
Invoke-AppRScript -ScriptPath $installScript

Write-Host ""
Write-Host "Package setup complete."
