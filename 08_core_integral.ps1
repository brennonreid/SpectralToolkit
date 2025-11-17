Param(
  [string]$Root      = $(Get-Location).Path,
  [string]$ToolsDir  = "",
  [string]$InputsDir = "",
  [string]$CertsDir  = "",
  [int]$Dps          = 950,
  [string]$T0        = "10000000000000000",  # base T0 for certificates (1e16)
  [switch]$UseVenv   = $true,
  [string]$VenvDir   = ""
)

function Msg($m)  { Write-Host "[08_core_integral] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[08_core_integral] WARNING: $m" -ForegroundColor Yellow }
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

Msg "Core integral / finite-interval and uniform certificate stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, CertsDir=$CertsDir, Dps=$Dps, T0=$T0"

# Core inputs this stage expects from earlier steps
$BandCert    = Join-Path $CertsDir "band_cert.json"
$GammaTail   = Join-Path $CertsDir "gamma_tail.json"
$PrimeTail   = Join-Path $CertsDir "prime_tail_envelope.json"
$PrimeBlock  = Join-Path $CertsDir "prime_block_norm.json"
$GridErr     = Join-Path $CertsDir "grid_error_bound.json"
$WeilPsd     = Join-Path $CertsDir "weil_psd_bochner.json"
$Continuum   = Join-Path $CertsDir "continuum_operator_cert.json"
$DensityJson = Join-Path $CertsDir "density_prover.json"
$WindowJson  = Join-Path $CertsDir "window.json"

Assert-File $BandCert    "band_cert.json not found (run 01_inputs.ps1)"
Assert-File $GammaTail   "gamma_tail.json not found (run 02_tails.ps1)"
Assert-File $PrimeTail   "prime_tail_envelope.json not found (run 02_tails.ps1)"
Assert-File $PrimeBlock  "prime_block_norm.json not found (run 03_prime_block_and_grid.ps1)"
Assert-File $GridErr     "grid_error_bound.json not found (run 03_prime_block_and_grid.ps1)"
Assert-File $WeilPsd     "weil_psd_bochner.json not found (run 04_weil_psd.ps1)"
Assert-File $Continuum   "continuum_operator_cert.json not found (run 05_continuum_rollup.ps1)"
Assert-File $WindowJson  "window.json not found in PROOF_PACKET (run 01_inputs.ps1)"

if (-not (Test-Path $DensityJson)) {
  Warn "density_prover.json missing; density info may be absent"
}

# 1) Uniform rollup certificate -> uniform_certificate.json
$UniformCert = Join-Path $CertsDir "uniform_certificate.json"
$urc = Join-Path $ToolsDir "uniform_rollup_cert.py"

if (Test-Path $urc) {
  Msg "python uniform_rollup_cert.py --T0 $T0 --certs-dir $CertsDir --out $UniformCert --dps $Dps"

  & $Py $urc `
    --T0 "$T0" `
    --certs-dir $CertsDir `
    --out $UniformCert `
    --dps "$Dps"

  if ($LASTEXITCODE -ne 0) {
    Warn "uniform_rollup_cert.py exited with code $LASTEXITCODE"
  } elseif (Test-Path $UniformCert) {
    Msg "uniform_certificate.json created -> $UniformCert"
  } else {
    Warn "uniform_certificate.json not created (expected at $UniformCert)"
  }
} else {
  Warn "uniform_rollup_cert.py not found in $ToolsDir; skipping uniform certificate"
}

# 2) Core finite-interval integral proof -> core_interval_prover.json
$CoreJson = Join-Path $CertsDir "core_interval_prover.json"
$cip = Join-Path $ToolsDir "core_interval_prover.py"
if (Test-Path $cip) {
  # CLI: core_interval_prover.py --T0 T0 --window-config WINDOW_CONFIG [--dps DPS] --out OUT
  Msg "python core_interval_prover.py --T0 $T0 --window-config $WindowJson --out $CoreJson --dps $Dps"
  & $Py $cip `
    --T0 "$T0" `
    --window-config $WindowJson `
    --out $CoreJson `
    --dps "$Dps"
} else {
  Warn "core_interval_prover.py not found in $ToolsDir; skipping core finite-interval proof"
}
if (Test-Path $CoreJson) {
  Msg "core_interval_prover.json created"
} else {
  if (Test-Path $cip) {
    Warn "core_interval_prover.json not created (check core_interval_prover.py CLI flags)"
  }
}

Msg "Core integral / finite-interval and uniform certificate stage complete."
