$ErrorActionPreference = "Stop"

function Get-AppDir {
  return (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
}

function Get-AppVersion {
  $appDir = Get-AppDir
  $versionFile = Join-Path $appDir "VERSION"
  if (Test-Path $versionFile) {
    $version = (Get-Content -Path $versionFile -Raw).Trim()
    if ($version) {
      return $version
    }
  }
  return "0.2.0"
}

function Find-Rscript {
  $pathCommand = Get-Command "Rscript.exe" -ErrorAction SilentlyContinue
  if ($pathCommand) {
    return $pathCommand.Source
  }

  $searchRoots = @(
    $env:ProgramFiles,
    ${env:ProgramFiles(x86)},
    $env:LOCALAPPDATA
  ) | Where-Object { $_ -and (Test-Path $_) }

  $candidates = New-Object System.Collections.Generic.List[string]
  foreach ($root in $searchRoots) {
    $rRoot = Join-Path $root "R"
    if (Test-Path $rRoot) {
      Get-ChildItem -Path $rRoot -Directory -Filter "R-*" -ErrorAction SilentlyContinue |
        Sort-Object Name -Descending |
        ForEach-Object {
          $candidates.Add((Join-Path $_.FullName "bin\Rscript.exe"))
        }
    }
  }

  $found = $candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
  if ($found) {
    return $found
  }

  throw @"
Rscript.exe was not found.

Install R for Windows from https://cran.r-project.org/bin/windows/base/
Then run this installer again.
"@
}

function Invoke-AppRScript {
  param(
    [Parameter(Mandatory = $true)]
    [string] $ScriptPath,

    [string[]] $Arguments = @()
  )

  $rscript = Find-Rscript
  & $rscript $ScriptPath @Arguments
  if ($LASTEXITCODE -ne 0) {
    throw "Rscript failed with exit code $LASTEXITCODE."
  }
}

function New-AppShortcut {
  param(
    [Parameter(Mandatory = $true)]
    [string] $ShortcutPath,

    [Parameter(Mandatory = $true)]
    [string] $TargetPath,

    [string] $Arguments = "",
    [string] $WorkingDirectory = "",
    [string] $Description = "TranscriptoScope",
    [string] $IconLocation = "",
    [int] $WindowStyle = 1
  )

  $shortcutDir = Split-Path -Parent $ShortcutPath
  if (-not (Test-Path $shortcutDir)) {
    New-Item -Path $shortcutDir -ItemType Directory -Force | Out-Null
  }

  $shell = New-Object -ComObject WScript.Shell
  $shortcut = $shell.CreateShortcut($ShortcutPath)
  $shortcut.TargetPath = $TargetPath
  $shortcut.Arguments = $Arguments
  $shortcut.WorkingDirectory = $WorkingDirectory
  $shortcut.Description = $Description
  $shortcut.WindowStyle = $WindowStyle
  if ($IconLocation) {
    $shortcut.IconLocation = $IconLocation
  }
  $shortcut.Save()
}

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
