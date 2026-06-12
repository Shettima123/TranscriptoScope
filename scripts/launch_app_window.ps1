param(
  [int] $Port = 7353,
  [switch] $NoBrowser
)

. "$PSScriptRoot\windows_common.ps1"

function Test-AppUrl {
  param([Parameter(Mandatory = $true)][string] $Url)

  try {
    $response = Invoke-WebRequest -Uri $Url -UseBasicParsing -TimeoutSec 2
    return (
      $response.StatusCode -ge 200 -and
      $response.StatusCode -lt 500 -and
      $response.Content -match "TranscriptoScope"
    )
  } catch {
    return $false
  }
}

function Test-PortInUse {
  param([Parameter(Mandatory = $true)][int] $PortNumber)

  $connection = Get-NetTCPConnection -LocalAddress "127.0.0.1" -LocalPort $PortNumber -ErrorAction SilentlyContinue
  return [bool]$connection
}

function Get-FreeTcpPort {
  $listener = [System.Net.Sockets.TcpListener]::new([System.Net.IPAddress]::Parse("127.0.0.1"), 0)
  try {
    $listener.Start()
    return ([System.Net.IPEndPoint]$listener.LocalEndpoint).Port
  } finally {
    $listener.Stop()
  }
}

function Find-AppBrowser {
  $pathCommand = Get-Command "msedge.exe" -ErrorAction SilentlyContinue
  if ($pathCommand) {
    return $pathCommand.Source
  }

  $candidates = @(
    (Join-Path ${env:ProgramFiles(x86)} "Microsoft\Edge\Application\msedge.exe"),
    (Join-Path $env:ProgramFiles "Microsoft\Edge\Application\msedge.exe"),
    (Join-Path $env:LOCALAPPDATA "Microsoft\Edge\Application\msedge.exe"),
    (Join-Path $env:ProgramFiles "Google\Chrome\Application\chrome.exe"),
    (Join-Path ${env:ProgramFiles(x86)} "Google\Chrome\Application\chrome.exe"),
    (Join-Path $env:LOCALAPPDATA "Google\Chrome\Application\chrome.exe")
  )

  $found = $candidates | Where-Object { $_ -and (Test-Path $_) } | Select-Object -First 1
  if ($found) {
    return $found
  }

  return $null
}

function Open-AppWindow {
  param([Parameter(Mandatory = $true)][string] $Url)

  $browser = Find-AppBrowser
  if ($browser) {
    Start-Process -FilePath $browser -ArgumentList @("--app=$Url", "--new-window") | Out-Null
    return
  }

  Start-Process $Url | Out-Null
}

$appDir = Get-AppDir
$serverScript = Join-Path $appDir "scripts\run_local_dev_server.R"
$url = "http://127.0.0.1:$Port/"

if (-not (Test-AppUrl -Url $url)) {
  if (Test-PortInUse -PortNumber $Port) {
    $Port = Get-FreeTcpPort
    $url = "http://127.0.0.1:$Port/"
  }

  $rscript = Find-Rscript
  $outLog = Join-Path $appDir "shiny-app-$Port-out.log"
  $errLog = Join-Path $appDir "shiny-app-$Port-err.log"
  if (Test-Path $outLog) { Remove-Item -LiteralPath $outLog -Force -ErrorAction SilentlyContinue }
  if (Test-Path $errLog) { Remove-Item -LiteralPath $errLog -Force -ErrorAction SilentlyContinue }

  Start-Process `
    -FilePath $rscript `
    -ArgumentList @($serverScript, $Port) `
    -WorkingDirectory $appDir `
    -WindowStyle Hidden `
    -RedirectStandardOutput $outLog `
    -RedirectStandardError $errLog | Out-Null

  $ready = $false
  for ($i = 0; $i -lt 60; $i++) {
    Start-Sleep -Seconds 1
    if (Test-AppUrl -Url $url) {
      $ready = $true
      break
    }
  }

  if (-not $ready) {
    throw "TranscriptoScope did not start. Check $errLog for details."
  }
}

if (-not $NoBrowser) {
  Open-AppWindow -Url $url
}
