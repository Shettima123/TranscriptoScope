param(
  [switch] $Force
)

. "$PSScriptRoot\windows_common.ps1"

$appName = "TranscriptoScope"
$defaultInstallDir = Join-Path $env:LOCALAPPDATA "TranscriptoScope"
$scriptInstallDir = Get-AppDir
if (Test-Path (Join-Path $scriptInstallDir "app.R")) {
  $installDir = $scriptInstallDir
} else {
  $installDir = $defaultInstallDir
}
$startMenuProgramsDir = Join-Path $env:APPDATA "Microsoft\Windows\Start Menu\Programs"
$startMenuDir = Join-Path $startMenuProgramsDir $appName
$startMenuShortcut = Join-Path $startMenuProgramsDir "$appName.lnk"
$desktopShortcut = Join-Path ([Environment]::GetFolderPath("Desktop")) "$appName.lnk"
$uninstallKey = "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\TranscriptoScope"
$localAppData = [System.IO.Path]::GetFullPath($env:LOCALAPPDATA)
$appDataPrograms = [System.IO.Path]::GetFullPath((Join-Path $env:APPDATA "Microsoft\Windows\Start Menu\Programs"))

if (-not $Force) {
  $answer = Read-Host "Uninstall $appName from $installDir? Type Y to continue"
  if ($answer -ne "Y" -and $answer -ne "y") {
    Write-Host "Uninstall cancelled."
    exit 0
  }
}

Set-Location $env:TEMP

try {
  Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
    Where-Object {
      $_.CommandLine -and
      $_.CommandLine.IndexOf($installDir, [System.StringComparison]::OrdinalIgnoreCase) -ge 0
    } |
    ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
  Start-Sleep -Milliseconds 500
} catch {
  Write-Host "Could not stop running app processes. Continuing with uninstall."
}

if (Test-Path $desktopShortcut) {
  Remove-Item -LiteralPath $desktopShortcut -Force
}

if (Test-Path $startMenuShortcut) {
  Remove-Item -LiteralPath $startMenuShortcut -Force
}

if (Test-Path $startMenuDir) {
  Assert-PathInsideDirectory -Path $startMenuDir -ParentDirectory $appDataPrograms
  Remove-Item -LiteralPath $startMenuDir -Recurse -Force
}

if (Test-Path $uninstallKey) {
  Remove-Item -LiteralPath $uninstallKey -Recurse -Force
}

if (Test-Path $installDir) {
  Assert-PathInsideDirectory -Path $installDir -ParentDirectory $localAppData
  Remove-Item -LiteralPath $installDir -Recurse -Force
}

Write-Host "$appName has been uninstalled."
