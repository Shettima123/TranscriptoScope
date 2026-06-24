param(
  [string] $OutputDir = (Join-Path (Join-Path $PSScriptRoot "..") "outputs\store_installer"),
  [string] $BuildRoot = (Join-Path $env:TEMP "TranscriptoScopeStoreBuild"),
  [switch] $SkipCompile,
  [switch] $SkipRuntimeCopy
)

$ErrorActionPreference = "Stop"

function Assert-PathInsideDirectory {
  param(
    [Parameter(Mandatory = $true)]
    [string] $Path,

    [Parameter(Mandatory = $true)]
    [string] $ParentDirectory
  )

  $resolvedPath = [System.IO.Path]::GetFullPath($Path)
  $resolvedParent = [System.IO.Path]::GetFullPath($ParentDirectory).TrimEnd('\')
  if (-not $resolvedPath.StartsWith($resolvedParent + "\", [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to modify path outside expected directory: $resolvedPath"
  }
}

function Find-Iscc {
  $command = Get-Command "iscc.exe" -ErrorAction SilentlyContinue
  if ($command) {
    return $command.Source
  }

  $candidates = @(
    (Join-Path $env:LOCALAPPDATA "Programs\Inno Setup 6\ISCC.exe"),
    "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
    "$env:ProgramFiles\Inno Setup 6\ISCC.exe"
  )
  foreach ($candidate in $candidates) {
    if ($candidate -and (Test-Path $candidate)) {
      return $candidate
    }
  }

  return $null
}

$appDir = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot ".."))
$appVersion = (Get-Content -Path (Join-Path $appDir "VERSION") -Raw).Trim()
if (-not $appVersion) {
  throw "VERSION file is empty."
}

$outputDirFull = [System.IO.Path]::GetFullPath($OutputDir)
$buildRoot = [System.IO.Path]::GetFullPath($BuildRoot)
$buildRootParent = Split-Path -Parent $buildRoot
if (-not $buildRootParent -or $buildRoot -eq [System.IO.Path]::GetPathRoot($buildRoot)) {
  throw "Unsafe build root: $buildRoot"
}
$stageRoot = Join-Path $buildRoot "stage"
$stageAppDir = Join-Path $stageRoot "TranscriptoScope"
$installerScript = Join-Path $buildRoot "TranscriptoScope_Store.iss"
$rootPackages = "shiny,ggplot2,DESeq2,SummarizedExperiment,S4Vectors,fgsea,WGCNA,impute,preprocessCore"

if (Test-Path $buildRoot) {
  Assert-PathInsideDirectory -Path $buildRoot -ParentDirectory $buildRootParent
  Remove-Item -LiteralPath $buildRoot -Recurse -Force
}
New-Item -Path $stageAppDir -ItemType Directory -Force | Out-Null
New-Item -Path $outputDirFull -ItemType Directory -Force | Out-Null

Write-Host "Staging TranscriptoScope $appVersion..."
Write-Host "Source: $appDir"
Write-Host "Stage:  $stageAppDir"

$excludeNames = @(
  ".git",
  ".Rproj.user",
  "__pycache__",
  "outputs",
  "runtime",
  "shiny-err.log",
  "shiny-out.log"
)

Get-ChildItem -Path $appDir -Force |
  Where-Object { $excludeNames -notcontains $_.Name } |
  ForEach-Object {
    Copy-Item -LiteralPath $_.FullName -Destination $stageAppDir -Recurse -Force
  }

$rscript = (Get-Command "Rscript.exe" -ErrorAction Stop).Source
$rHome = & $rscript -e "cat(normalizePath(R.home(), winslash='/', mustWork=TRUE))"
$rHome = $rHome.Trim()
if (-not (Test-Path $rHome)) {
  throw "Could not resolve R.home(): $rHome"
}

$runtimeDir = Join-Path $stageAppDir "runtime"
$runtimeRDir = Join-Path $runtimeDir "R"
$runtimeLibRoot = Join-Path $runtimeDir "R-library"
$runtimeLibVersion = & $rscript -e "cat(paste(R.version`$major, sub('\\..*$', '', R.version`$minor), sep='.'))"
$runtimeLib = Join-Path $runtimeLibRoot $runtimeLibVersion.Trim()

if (-not $SkipRuntimeCopy) {
  Write-Host "Copying R runtime from $rHome..."
  Copy-Item -LiteralPath $rHome -Destination $runtimeRDir -Recurse -Force
}

Write-Host "Staging R package dependency library..."
& $rscript (Join-Path $stageAppDir "scripts\stage_store_r_library.R") $runtimeLib $rootPackages
if ($LASTEXITCODE -ne 0) {
  throw "R package staging failed."
}

$packageManifest = Join-Path $runtimeLib "TRANSCRIPTOSCOPE_STORE_R_PACKAGES.csv"
if (-not (Test-Path $packageManifest)) {
  throw "Package manifest was not created: $packageManifest"
}
$packagesToCopy = Import-Csv -LiteralPath $packageManifest
foreach ($package in $packagesToCopy) {
  if (-not (Test-Path $package.SourcePackageDir)) {
    throw "Package source directory does not exist: $($package.SourcePackageDir)"
  }
  Write-Host "Copying R package $($package.Package) $($package.Version)..."
  Copy-Item -LiteralPath $package.SourcePackageDir -Destination $runtimeLib -Recurse -Force
}
Write-Host "Copied $($packagesToCopy.Count) R package(s)."

Write-Host "Running bundled runtime/package preflight..."
$bundledRscript = Join-Path $runtimeRDir "bin\Rscript.exe"
if (-not (Test-Path $bundledRscript)) {
  throw "Bundled Rscript not found after staging: $bundledRscript"
}
$env:R_LIBS_USER = $runtimeLib
$env:R_LIBS_SITE = ""
& $bundledRscript (Join-Path $stageAppDir "scripts\store_preflight.R") $stageAppDir
if ($LASTEXITCODE -ne 0) {
  throw "Store preflight failed."
}

$escapedStage = $stageAppDir.Replace("\", "\\")
$escapedOutput = $outputDirFull.Replace("\", "\\")
$escapedIcon = (Join-Path $stageAppDir "www\favicon.ico").Replace("\", "\\")
$setupGuid = "{{7B311F6D-71D4-4DA8-9BBA-4976EE502CA9}"

$iss = @"
#define MyAppName "TranscriptoScope"
#define MyAppVersion "$appVersion"
#define MyAppPublisher "Dr. Abubakar Abdulkadir | Dr. Rosby's Lab, Southern University A and M"
#define MyAppURL "https://github.com/Shettima123/TranscriptoScope"

[Setup]
AppId=$setupGuid
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={localappdata}\TranscriptoScope
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=$escapedOutput
OutputBaseFilename=TranscriptoScope_Windows_Store_v$appVersion
SetupIconFile=$escapedIcon
UninstallDisplayIcon={app}\www\favicon.ico
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Files]
Source: "$escapedStage\\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\TranscriptoScope"; Filename: "{sys}\wscript.exe"; Parameters: """{app}\Launch_TranscriptoScope_Store.vbs"""; WorkingDir: "{app}"; IconFilename: "{app}\www\favicon.ico"
Name: "{autodesktop}\TranscriptoScope"; Filename: "{sys}\wscript.exe"; Parameters: """{app}\Launch_TranscriptoScope_Store.vbs"""; WorkingDir: "{app}"; IconFilename: "{app}\www\favicon.ico"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: checkedonce

[Run]
Filename: "{app}\runtime\R\bin\Rscript.exe"; Parameters: """{app}\scripts\store_preflight.R"" ""{app}"""; WorkingDir: "{app}"; Flags: runhidden

[UninstallDelete]
Type: filesandordirs; Name: "{localappdata}\TranscriptoScope\logs"
"@

New-Item -Path $buildRoot -ItemType Directory -Force | Out-Null
Set-Content -Path $installerScript -Value $iss -Encoding UTF8
Write-Host "Inno Setup script written to $installerScript"

if ($SkipCompile) {
  Write-Host "SkipCompile was supplied; installer compilation was not run."
  exit 0
}

$iscc = Find-Iscc
if (-not $iscc) {
  throw "Inno Setup compiler ISCC.exe was not found. Install Inno Setup 6, then rerun this script."
}

Write-Host "Compiling installer with $iscc..."
& $iscc $installerScript
if ($LASTEXITCODE -ne 0) {
  throw "Inno Setup compilation failed with exit code $LASTEXITCODE."
}

$installerPath = Join-Path $outputDirFull "TranscriptoScope_Windows_Store_v$appVersion.exe"
if (-not (Test-Path $installerPath)) {
  throw "Expected installer was not created: $installerPath"
}

$hash = Get-FileHash -Algorithm SHA256 -LiteralPath $installerPath
Write-Host ""
Write-Host "Store installer created:"
Write-Host $installerPath
Write-Host "SHA256: $($hash.Hash)"
