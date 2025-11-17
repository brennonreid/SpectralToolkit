Param(
  [string]$Root      = $(Get-Location).Path,
  [string]$ToolsDir  = "",
  [string]$CertsDir  = "",
  [int]$Dps          = 950,
  [switch]$UseVenv   = $true,
  [string]$VenvDir   = ""
)

function Msg($m)  { Write-Host "[05_continuum_rollup] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[05_continuum_rollup] WARNING: $m" -ForegroundColor Yellow }
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

Msg "Continuum / operator rollup stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, CertsDir=$CertsDir, Dps=$Dps"

# Core inputs that should already exist from earlier stages
$BandCert    = Join-Path $CertsDir "band_cert.json"
$GammaTail   = Join-Path $CertsDir "gamma_tail.json"
$PrimeTail   = Join-Path $CertsDir "prime_tail_envelope.json"
$PrimeBlock  = Join-Path $CertsDir "prime_block_norm.json"
$GridErr     = Join-Path $CertsDir "grid_error_bound.json"

Assert-File $BandCert   "band_cert.json not found (run 01_inputs.ps1)"
Assert-File $GammaTail  "gamma_tail.json not found (run 02_tails.ps1)"
Assert-File $PrimeTail  "prime_tail_envelope.json not found (run 02_tails.ps1)"
Assert-File $PrimeBlock "prime_block_norm.json not found (run 03_prime_block_and_grid.ps1)"
if (-not (Test-Path $GridErr)) {
  Warn "grid_error_bound.json not found; continuum rollup will proceed without a grid error term"
}

# Continuum operator rollup -> continuum_operator_cert.json
$ContinuumJson = Join-Path $CertsDir "continuum_operator_cert.json"
$cor = Join-Path $ToolsDir "continuum_operator_rollup.py"
Assert-File $cor "continuum_operator_rollup.py not found in $ToolsDir"

# Build args per normalization spec:
#   --band-cert, --prime-block, --prime-tail, --gamma-tails, [--grid-error], --dps, --out
$contArgs = @(
  "--band-cert",   $BandCert,
  "--prime-block", $PrimeBlock,
  "--prime-tail",  $PrimeTail,
  "--gamma-tails", $GammaTail,
  "--dps",         "$Dps",
  "--out",         $ContinuumJson
)
if (Test-Path $GridErr) {
  $contArgs = @(
    "--band-cert",   $BandCert,
    "--prime-block", $PrimeBlock,
    "--prime-tail",  $PrimeTail,
    "--gamma-tails", $GammaTail,
    "--grid-error",  $GridErr,
    "--dps",         "$Dps",
    "--out",         $ContinuumJson
  )
}

Msg "python continuum_operator_rollup.py $($contArgs -join ' ')"
& $Py $cor @contArgs
Assert-File $ContinuumJson "continuum_operator_rollup did not produce continuum_operator_cert.json"

Msg "Continuum rollup stage complete."
