Param(
  [string]$Root = $(Get-Location).Path,
  [string]$ToolsDir = "",
  [string]$CertsDir = "",
  [int]$Dps = 900,
  [string]$T0 = "10000000000000000",   # 1e16
  [string]$X0 = "10000000000000000",   # kept for compatibility; may be unused
  [switch]$UseVenv = $true,
  [string]$VenvDir = ""
)

function Msg($m)  { Write-Host "[10_infinite_analysis] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[10_infinite_analysis] WARNING: $m" -ForegroundColor Yellow }
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

Msg "Infinite analysis stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, CertsDir=$CertsDir, Dps=$Dps, T0=$T0, X0=$X0"

# 1) Riemann-von Mangoldt / infinitude witness -> rv_mangoldt_bounds.json
$RvMain   = Join-Path $CertsDir "rv_mangoldt_bounds.json"
$RvTheory = Join-Path $CertsDir "rv_mangoldt_bounds.theory.json"
$rv = Join-Path $ToolsDir "rv_mangoldt_bounds.py"

if (Test-Path $rv) {
  Msg "python rv_mangoldt_bounds.py --T0 $T0 --dps $Dps --out $RvMain --theory-out $RvTheory"
  & $Py $rv `
    --T0 "$T0" `
    --dps "$Dps" `
    --out $RvMain `
    --theory-out $RvTheory
} else {
  Warn "rv_mangoldt_bounds.py not found in $ToolsDir; skipping infinite-analysis witness"
}

if (Test-Path $RvMain) {
  Msg "[rv_mangoldt_bounds] ok -> $RvMain"
} else {
  if (Test-Path $rv) {
    Warn "rv_mangoldt_bounds.json not created (check rv_mangoldt_bounds.py CLI flags)"
  }
}

Msg "Infinite-analysis stage complete."
