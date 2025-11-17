[CmdletBinding()]
Param(
  # Paths
  [string]$Root      = $(Get-Location).Path,
  [string]$ToolsDir  = "",
  [string]$CertsDir  = "",

  # Rolling-T settings (single sweep interval)
  [string]$T0        = "1e16",   # start T
  [string]$T1        = "1e20",   # end T
  [string]$DeltaTarget = "1e-12",
  [int]$MeshInitial  = 256,
  [int]$MeshMax      = 131072,

  # Precision / formatting
  [int]$Dps          = 950,
  [int]$Digits       = 220,

  # Python / venv
  [switch]$UseVenv   = $true,
  [string]$VenvDir   = ""
)

function Msg($m)  { Write-Host "[12_rolling_t] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[12_rolling_t] WARNING: $m" -ForegroundColor Yellow }
function Assert-File($p, $m) { if (-not (Test-Path $p)) { Die $m } }

if (-not $ToolsDir)  { $ToolsDir  = Join-Path $Root "tools" }
if (-not $CertsDir)  { $CertsDir  = Join-Path $Root "PROOF_PACKET" }
if (-not $VenvDir)   { $VenvDir   = Join-Path $Root "rh-env" }

$ErrorActionPreference = "Stop"

$Py = "python"
if ($UseVenv) {
  $PyCandidate = Join-Path $VenvDir "Scripts\python.exe"
  if (Test-Path $PyCandidate) {
    $Py = $PyCandidate
  } else {
    Warn "venv requested but $PyCandidate not found, using python on PATH"
  }
}

Msg "Rolling-T uniform sweep stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, CertsDir=$CertsDir, Dps=$Dps"
Msg "T in [$T0, $T1], delta_target=$DeltaTarget, mesh_initial=$MeshInitial, mesh_max=$MeshMax"

$Tool = Join-Path $ToolsDir "rolling_T_uniform_cert_v3.py"
Assert-File $Tool "rolling_T_uniform_cert_v3.py not found in $ToolsDir"

$OutFile = Join-Path $CertsDir "rolling_T_uniform_cert.json"

# Build args to match the current Python CLI:
# rolling_T_uniform_cert_v3.py --packet-dir PACKET_DIR --T0 T0 --T1 T1
#   [--delta-target ...] [--mesh-initial ...] [--mesh-max ...] [--dps ...] [--digits ...] --out OUT
$ArgsList = @(
  $Tool,
  "--packet-dir", $CertsDir,
  "--T0",         $T0,
  "--T1",         $T1,
  "--delta-target", $DeltaTarget,
  "--mesh-initial", ("{0}" -f $MeshInitial),
  "--mesh-max",     ("{0}" -f $MeshMax),
  "--dps",          ("{0}" -f $Dps),
  "--digits",       ("{0}" -f $Digits),
  "--out",          $OutFile
)

Msg "python $($ArgsList -join ' ')"
& $Py @ArgsList
$ec = $LASTEXITCODE
if ($ec -ne 0) {
  Die "rolling_T_uniform_cert_v3.py exited with code $ec"
}

if (Test-Path $OutFile) {
  Msg "rolling_T_uniform_cert.json created"
} else {
  Die "Expected output not found: $OutFile"
}
