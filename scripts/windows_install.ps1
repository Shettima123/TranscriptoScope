param(
  [string] $InstallDir = (Join-Path $env:LOCALAPPDATA "TranscriptoScope"),
  [switch] $SkipPackageInstall,
  [switch] $NoDesktopShortcut
)

. "$PSScriptRoot\windows_common.ps1"

$sourceDir = Get-AppDir
$installDir = [System.IO.Path]::GetFullPath($InstallDir)
$appName = "TranscriptoScope"
$appVersion = Get-AppVersion
$startMenuProgramsDir = Join-Path $env:APPDATA "Microsoft\Windows\Start Menu\Programs"
$startMenuDir = Join-Path $startMenuProgramsDir $appName
$startMenuShortcut = Join-Path $startMenuProgramsDir "$appName.lnk"
$desktopShortcut = Join-Path ([Environment]::GetFolderPath("Desktop")) "$appName.lnk"
$launcherScript = Join-Path $installDir "Launch_TranscriptoScope.vbs"
$wscriptPath = Join-Path $env:WINDIR "System32\wscript.exe"
$uninstallBatch = Join-Path $installDir "Uninstall_TranscriptoScope.bat"
$localAppData = [System.IO.Path]::GetFullPath($env:LOCALAPPDATA)
$legacyAppName = "DESeq2 DGE Workbench"
$legacyInstallDir = Join-Path $env:LOCALAPPDATA "DESeq2DGEWorkbench"
$legacyStartMenuDir = Join-Path $startMenuProgramsDir $legacyAppName
$legacyStartMenuShortcut = Join-Path $startMenuProgramsDir "$legacyAppName.lnk"
$legacyDesktopShortcut = Join-Path ([Environment]::GetFolderPath("Desktop")) "$legacyAppName.lnk"
$legacyUninstallKey = "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\DESeq2DGEWorkbench"

Write-Host "Installing $appName $appVersion..."
Write-Host "Source: $sourceDir"
Write-Host "Install folder: $installDir"

Write-Host "Checking for old beta install names..."
try {
  Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
    Where-Object { $_.CommandLine -and $_.CommandLine -like "*DESeq2DGEWorkbench*" } |
    ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
} catch {
  Write-Host "Could not check running legacy app processes. Continuing."
}

foreach ($legacyShortcut in @($legacyDesktopShortcut, $legacyStartMenuShortcut)) {
  if (Test-Path $legacyShortcut) {
    Remove-Item -LiteralPath $legacyShortcut -Force -ErrorAction SilentlyContinue
  }
}
if (Test-Path $legacyStartMenuDir) {
  Remove-Item -LiteralPath $legacyStartMenuDir -Recurse -Force -ErrorAction SilentlyContinue
}
if (Test-Path $legacyUninstallKey) {
  Remove-Item -LiteralPath $legacyUninstallKey -Recurse -Force -ErrorAction SilentlyContinue
}
if (Test-Path $legacyInstallDir) {
  try {
    Assert-PathInsideDirectory -Path $legacyInstallDir -ParentDirectory $localAppData
    Remove-Item -LiteralPath $legacyInstallDir -Recurse -Force -ErrorAction Stop
  } catch {
    Write-Host "Old beta folder could not be removed yet: $legacyInstallDir"
  }
}

if (-not (Test-Path $installDir)) {
  New-Item -Path $installDir -ItemType Directory -Force | Out-Null
}

if ($sourceDir -ne $installDir) {
  $excludeNames = @(".git", ".Rproj.user", "__pycache__", "shiny-err.log", "shiny-out.log")
  Get-ChildItem -Path $sourceDir -Force |
    Where-Object { $excludeNames -notcontains $_.Name } |
    ForEach-Object {
      Copy-Item -Path $_.FullName -Destination $installDir -Recurse -Force
    }
}

New-AppShortcut `
  -ShortcutPath (Join-Path $startMenuDir "$appName.lnk") `
  -TargetPath $wscriptPath `
  -Arguments "`"$launcherScript`"" `
  -WorkingDirectory $installDir `
  -Description "Launch $appName" `
  -IconLocation "$wscriptPath,2"

New-AppShortcut `
  -ShortcutPath $startMenuShortcut `
  -TargetPath $wscriptPath `
  -Arguments "`"$launcherScript`"" `
  -WorkingDirectory $installDir `
  -Description "Launch $appName" `
  -IconLocation "$wscriptPath,2"

New-AppShortcut `
  -ShortcutPath (Join-Path $startMenuDir "Uninstall $appName.lnk") `
  -TargetPath $uninstallBatch `
  -WorkingDirectory $installDir `
  -Description "Uninstall $appName"

if (-not $NoDesktopShortcut) {
  New-AppShortcut `
    -ShortcutPath $desktopShortcut `
    -TargetPath $wscriptPath `
    -Arguments "`"$launcherScript`"" `
    -WorkingDirectory $installDir `
    -Description "Launch $appName" `
    -IconLocation "$wscriptPath,2"
}

$uninstallKey = "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\TranscriptoScope"
if (-not (Test-Path $uninstallKey)) {
  New-Item -Path $uninstallKey -Force | Out-Null
}
Set-ItemProperty -Path $uninstallKey -Name "DisplayName" -Value $appName
Set-ItemProperty -Path $uninstallKey -Name "DisplayVersion" -Value $appVersion
Set-ItemProperty -Path $uninstallKey -Name "Publisher" -Value "Dr. Abubakar Abdulkadir, Southern University A and M"
Set-ItemProperty -Path $uninstallKey -Name "InstallLocation" -Value $installDir
Set-ItemProperty -Path $uninstallKey -Name "UninstallString" -Value "`"$uninstallBatch`""
Set-ItemProperty -Path $uninstallKey -Name "DisplayIcon" -Value $launcherScript

if (-not $SkipPackageInstall) {
  Write-Host ""
  Write-Host "Checking R and installing required packages..."
  & (Join-Path $installDir "scripts\windows_install_packages.ps1")
  if ($LASTEXITCODE -ne 0) {
    throw "Package installation failed. You can rerun Install_Packages.bat from $installDir."
  }
}

Write-Host ""
Write-Host "$appName $appVersion is installed."
Write-Host "Use the Desktop or Start Menu shortcut to launch it."
