. "$PSScriptRoot\windows_common.ps1"

$appDir = Get-AppDir
$launchScript = Join-Path $appDir "scripts\launch_app.R"

Write-Host "Starting TranscriptoScope..."
Write-Host "App folder: $appDir"

Set-Location $appDir
Invoke-AppRScript -ScriptPath $launchScript
