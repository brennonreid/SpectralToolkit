Param(
  [string]$Root = $(Get-Location).Path,
  [string]$ToolsDir = "",
  [string]$InputsDir = "",
  [string]$CertsDir = "",
  [int]$Dps = 220,
  [switch]$UseVenv = $true,
  [string]$VenvDir = ""
)

function Msg($m)  { Write-Host "[04_weil_psd] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[04_weil_psd] WARNING: $m" -ForegroundColor Yellow }
function Assert-File($p, $m) { if (-not (Test-Path $p)) { Die $m } }

if (-not $ToolsDir)  { $ToolsDir  = Join-Path $Root "tools" }
if (-not $InputsDir) { $InputsDir = Join-Path $Root "packs\rh\inputs" }
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

Msg "Weil PSD stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, CertsDir=$CertsDir, Dps=$Dps"

$WindowJson = Join-Path $CertsDir "window.json"
Assert-File $WindowJson "window.json not found in PROOF_PACKET (run 01_inputs.ps1 first)"

$Bochner = Join-Path $CertsDir "weil_psd_bochner.json"

# Prefer analytic Weil kernel, fall back to numeric Bochner PSD
$weilk = Join-Path $ToolsDir "weil_kernel.py"
if (Test-Path $weilk) {
  Msg "python weil_kernel.py --window-config $WindowJson --out $Bochner --dps $Dps --method bochner"
  & $Py $weilk `
    --window-config $WindowJson `
    --out $Bochner `
    --dps "$Dps" `
    --method "bochner"
} else {
  $bo = Join-Path $ToolsDir "bochner_psd_cert.py"
  if (Test-Path $bo) {
    Msg "python bochner_psd_cert.py --window-config $WindowJson --out $Bochner --dps $Dps"
    & $Py $bo `
      --window-config $WindowJson `
      --out $Bochner `
      --dps "$Dps"
  } else {
    Die "No PSD tool found (neither weil_kernel.py nor bochner_psd_cert.py exists in $ToolsDir)"
  }
}

Assert-File $Bochner "Weil PSD must produce weil_psd_bochner.json"

Msg "Weil PSD stage complete."
