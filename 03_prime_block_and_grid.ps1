Param(
  [string]$Root = $(Get-Location).Path,
  [string]$ToolsDir = "",
  [string]$InputsDir = "",
  [string]$CertsDir = "",
  [int]$Dps = 950,
  [double]$Sigma = 6.0,
  [double]$K0 = 0.25,
  [int]$Grid = 8000,
  [switch]$UseVenv = $true,
  [string]$VenvDir = ""
)

function Msg($m)  { Write-Host "[03_prime_block_and_grid] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[03_prime_block_and_grid] WARNING: $m" -ForegroundColor Yellow }
function Assert-File($p, $m) { if (-not (Test-Path $p)) { Die $m } }

if (-not $ToolsDir)  { $ToolsDir  = Join-Path $Root "tools" }
if (-not $InputsDir) { $InputsDir = Join-Path $Root "" }
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

Msg "Prime block and grid stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, InputsDir=$InputsDir, CertsDir=$CertsDir, Dps=$Dps, Grid=$Grid, Sigma=$Sigma, K0=$K0"

# Window config (for grid error)
$WindowJson = Join-Path $CertsDir "window.json"
Assert-File $WindowJson "window.json not found in PROOF_PACKET (run 01_inputs.ps1 first)"

# zeros.txt input
$ZerosTxt = Join-Path $InputsDir "zeros.txt"
Assert-File $ZerosTxt "zeros.txt not found in inputs directory ($InputsDir)"

# 1) prime_block_norm -> prime_block_norm.json
$PrimeBlock = Join-Path $CertsDir "prime_block_norm.json"
$pb = Join-Path $ToolsDir "prime_block_norm.py"
if (Test-Path $pb) {
  Msg "python prime_block_norm.py --zeros $ZerosTxt --sigma $Sigma --k0 $K0 --out $PrimeBlock --dps $Dps"
  & $Py $pb `
    --zeros "$ZerosTxt" `
    --sigma "$Sigma" `
    --k0 "$K0" `
    --out $PrimeBlock `
    --dps "$Dps"
} else {
  Die "prime_block_norm.py not found in $ToolsDir"
}
Assert-File $PrimeBlock "prime_block_norm must produce prime_block_norm.json"

# 2) op_grid_error_bound -> grid_error_bound.json
$GridErr = Join-Path $CertsDir "grid_error_bound.json"
$geb = Join-Path $ToolsDir "op_grid_error_bound.py"
if (Test-Path $geb) {
  Msg "python op_grid_error_bound.py --grid $Grid --out $GridErr --dps $Dps"
  & $Py $geb `
    --grid "$Grid" `
    --out $GridErr `
    --dps "$Dps"
} else {
  Warn "op_grid_error_bound.py not found in $ToolsDir; skipping grid error bound"
}

Msg "Prime block and grid stage complete."
