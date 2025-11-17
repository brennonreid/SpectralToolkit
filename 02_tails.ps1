Param(
  [string]$Root      = $(Get-Location).Path,
  [string]$ToolsDir  = "",
  [string]$InputsDir = "",
  [string]$CertsDir  = "",
  [int]$Dps          = 220,
  [double]$Sigma     = 6.0,
  [double]$K0        = 0.25,
  [string]$T0        = "1000000",
  [switch]$UseVenv   = $true,
  [string]$VenvDir   = ""
)

function Msg($m)  { Write-Host "[02_tails] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[02_tails] WARNING: $m" -ForegroundColor Yellow }
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
    Warn "venv requested but $PyCandidate not found, falling back to python on PATH."
  }
}

Msg "Tails stage (gamma tail + prime tail)"
Msg "Root=$Root, ToolsDir=$ToolsDir, CertsDir=$CertsDir, Dps=$Dps, T0=$T0, Sigma=$Sigma, K0=$K0"

# 1) Gamma tail envelope -> gamma_tail.json
$GammaTail = Join-Path $CertsDir "gamma_tail.json"
$te = Join-Path $ToolsDir "tail_envelope.py"
if (Test-Path $te) {
  Msg "python tail_envelope.py --sigma $Sigma --T0 $T0 --out $GammaTail --dps $Dps"
  & $Py $te `
    --sigma "$Sigma" `
    --T0 "$T0" `
    --out $GammaTail `
    --dps "$Dps"
}
Assert-File $GammaTail "tail_envelope must produce gamma_tail.json"

# 2) Prime tail envelope -> prime_tail_envelope.json
$PrimeTail = Join-Path $CertsDir "prime_tail_envelope.json"
$pte = Join-Path $ToolsDir "prime_tail_envelope.py"
if (Test-Path $pte) {
  Msg "python prime_tail_envelope.py --sigma $Sigma --k0 $K0 --T0 $T0 --out $PrimeTail --dps $Dps --tqdm"
  & $Py $pte `
    --sigma "$Sigma" `
    --k0 "$K0" `
    --T0 "$T0" `
    --out $PrimeTail `
    --dps "$Dps" `
}
Assert-File $PrimeTail "prime_tail_envelope must produce prime_tail_envelope.json"

Msg "Tails stage complete."
