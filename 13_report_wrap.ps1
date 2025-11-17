[CmdletBinding()]
Param(
  [string]$Root      = $(Get-Location).Path,
  [string]$ToolsDir  = "",
  [string]$CertsDir  = "",
  [string]$InputsDir = "",
  [int]$Dps          = 900,
  [switch]$UseVenv   = $true,
  [string]$VenvDir   = ""
)

function Msg($m) { Write-Host "[13_report_wrap] $m" -ForegroundColor Cyan }
function Die($m) { throw "[13_report_wrap] FATAL: $m" }

function Resolve-PathSafe($p) {
  if (-not $p) { return "" }
  return (Resolve-Path -LiteralPath $p).Path
}

function Get-PythonExe {
  param(
    [switch]$UseVenv,
    [string]$VenvDir
  )
  if ($UseVenv -and $VenvDir) {
    $py = Join-Path $VenvDir "Scripts\python.exe"
    if (Test-Path $py) { return $py }
  }
  $pyGlobal = "python"
  return $pyGlobal
}

if (-not $ToolsDir)  { $ToolsDir  = Join-Path $Root "tools" }
if (-not $CertsDir)  { $CertsDir  = Join-Path $Root "PROOF_PACKET" }
if (-not $InputsDir) { $InputsDir = $CertsDir }

$ToolsDir  = Resolve-PathSafe $ToolsDir
$CertsDir  = Resolve-PathSafe $CertsDir
$InputsDir = Resolve-PathSafe $InputsDir

if (-not (Test-Path $ToolsDir))  { Die "ToolsDir not found: $ToolsDir" }
if (-not (Test-Path $CertsDir))  { Die "CertsDir not found: $CertsDir" }
if (-not (Test-Path $InputsDir)) { Die "InputsDir not found: $InputsDir" }

$PythonExe = Get-PythonExe -UseVenv:$UseVenv -VenvDir $VenvDir
Msg "Using Python: $PythonExe"

$Script = Join-Path $ToolsDir "better_report_wrap.py"
if (-not (Test-Path $Script)) { Die "better_report_wrap.py not found at: $Script" }

$outJson = Join-Path $CertsDir "report_wrap.json"
$outMd   = Join-Path $CertsDir "report_wrap.md"

Msg "Generating final report wrap..."
$argsList = @(
  $Script,
  "--certs-dir", $CertsDir,
  "--inputs-dir", $InputsDir,
  "--out-md", $outMd,
  "--out-json", $outJson,
  "--dps", $Dps
)

& $PythonExe @argsList
if ($LASTEXITCODE -ne 0) {
  Die "better_report_wrap.py failed with exit code $LASTEXITCODE"
}

Msg "Report wrap written:"
Msg "  JSON : $outJson"
Msg "  Markdown : $outMd"
