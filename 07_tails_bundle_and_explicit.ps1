Param(
  [string]$Root      = $(Get-Location).Path,
  [string]$ToolsDir  = "",
  [string]$InputsDir = "",
  [string]$CertsDir  = "",
  [int]$Dps          = 950,
  [double]$Sigma     = 6.0,
  [double]$K0        = 0.25,
  [string]$T0        = "10000000000000000",   # 1e16, now just context for logs
  [double]$X0        = 40.0,
  [double]$Ap        = 1.0,
  [double]$Ag        = 1.0,
  [switch]$UseVenv   = $true,
  [string]$VenvDir   = ""
)

function Msg($m)  { Write-Host "[07_tails_bundle_and_explicit] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[07_tails_bundle_and_explicit] WARNING: $m" -ForegroundColor Yellow }
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

Msg "Tails bundle and explicit-formula stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, InputsDir=$InputsDir, CertsDir=$CertsDir, Dps=$Dps, Sigma=$Sigma, K0=$K0, T0=$T0, X0=$X0, Ap=$Ap, Ag=$Ag"

# Core certs this stage ultimately depends on
$BandCert       = Join-Path $CertsDir "band_cert.json"
$WeilPsd        = Join-Path $CertsDir "weil_psd_bochner.json"
$ContinuumCert  = Join-Path $CertsDir "continuum_operator_cert.json"
$GridErr        = Join-Path $CertsDir "grid_error_bound.json"
$PrimeTailEnv   = Join-Path $CertsDir "prime_tail_envelope.json"
$GammaTail      = Join-Path $CertsDir "gamma_tail.json"

Assert-File $BandCert      "band_cert.json not found (run 01_inputs.ps1)"
Assert-File $WeilPsd       "weil_psd_bochner.json not found (run 04_weil_psd.ps1)"
Assert-File $ContinuumCert "continuum_operator_cert.json not found (run 05_continuum_rollup.ps1)"
Assert-File $GridErr       "grid_error_bound.json not found (run 03_prime_block_and_grid.ps1)"
Assert-File $PrimeTailEnv  "prime_tail_envelope.json not found (run 02_tails.ps1)"
Assert-File $GammaTail     "gamma_tail.json not found (run 02_tails.ps1)"

# 1) Operator prime tail bound -> op_prime_tail_bound.json
$PrimeOpTail = Join-Path $CertsDir "op_prime_tail_bound.json"
$opTail = Join-Path $ToolsDir "op_prime_tail_bound.py"
if (Test-Path $opTail) {
  # CLI: --x0, --A-prime, --K, [--scale-by-log], --out, [--dps]
  Msg "python op_prime_tail_bound.py --x0 $X0 --A-prime $Ap --K 3 --out $PrimeOpTail --dps $Dps"
  & $Py $opTail `
    --x0 "$X0" `
    --A-prime "$Ap" `
    --K "3" `
    --out $PrimeOpTail `
    --dps "$Dps"
} else {
  Warn "op_prime_tail_bound.py not found in $ToolsDir; skipping op_prime_tail_bound"
}
if (Test-Path $PrimeOpTail) {
  Msg "op_prime_tail_bound.json created"
} else {
  Warn "op_prime_tail_bound.json not created"
}

# 2) Analytic bounds -> analytic_bounds.json
$AnalyticBounds = Join-Path $CertsDir "analytic_bounds.json"
$ab = Join-Path $ToolsDir "analytic_bounds.py"
if (Test-Path $ab) {
  # CLI: analytic_bounds.py [--sigma SIGMA] [--A-prime A_PRIME] [--K K] --out OUT [--dps DPS]
  Msg "python analytic_bounds.py --sigma $Sigma --A-prime $Ap --K 3 --out $AnalyticBounds --dps $Dps"
  & $Py $ab `
    --sigma "$Sigma" `
    --A-prime "$Ap" `
    --K "3" `
    --out $AnalyticBounds `
    --dps "$Dps"
} else {
  Warn "analytic_bounds.py not found in $ToolsDir; skipping analytic bounds"
}
if (Test-Path $AnalyticBounds) {
  Msg "analytic_bounds.json created"
} else {
  Warn "analytic_bounds.json not created"
}

# 3) Analytic tail fit bundle -> analytic_tail_fit.json  (tails JSON)
$TailsJson = Join-Path $CertsDir "analytic_tail_fit.json"
$atf = Join-Path $ToolsDir "analytic_tail_fit.py"
if (Test-Path $atf) {
  # CLI: analytic_tail_fit.py --packet-dir PROOF_PACKET --Ap Ap --Ag Ag [--out] [--dps]
  Msg "python analytic_tail_fit.py --packet-dir $CertsDir --Ap $Ap --Ag $Ag --out $TailsJson --dps $Dps"
  & $Py $atf `
    --packet-dir $CertsDir `
    --Ap "$Ap" `
    --Ag "$Ag" `
    --out $TailsJson `
    --dps "$Dps"
} else {
  Warn "analytic_tail_fit.py not found in $ToolsDir; skipping analytic tail bundle"
}
if (Test-Path $TailsJson) {
  Msg "analytic_tail_fit.json created"
} else {
  Warn "analytic_tail_fit.json not created; explicit_formula tails input will be missing"
}

# 4) Explicit formula rollup -> explicit_formula.json
$ExplicitJson = Join-Path $CertsDir "explicit_formula.json"
$ef = Join-Path $ToolsDir "explicit_formula.py"

if (-not (Test-Path $ef)) {
  Warn "explicit_formula.py not found in $ToolsDir; skipping explicit formula"
} elseif (-not (Test-Path $TailsJson)) {
  Warn "tails JSON (analytic_tail_fit.json) not found; skipping explicit_formula.py"
} else {
  # CLI: explicit_formula.py --band-cert BAND_CERT --weil-psd WEIL_PSD --tails TAILS [--continuum-cert CONTINUUM_CERT] [--dps] --out OUT
  Msg "python explicit_formula.py --band-cert $BandCert --weil-psd $WeilPsd --tails $TailsJson --continuum-cert $ContinuumCert --dps $Dps --out $ExplicitJson"
  & $Py $ef `
    --band-cert $BandCert `
    --weil-psd $WeilPsd `
    --tails $TailsJson `
    --continuum-cert $ContinuumCert `
    --dps "$Dps" `
    --out $ExplicitJson

  if (Test-Path $ExplicitJson) {
    Msg "explicit_formula.json created"
  } else {
    Warn "explicit_formula.json not created"
  }
}

Msg "Tails bundle and explicit-formula stage complete."
