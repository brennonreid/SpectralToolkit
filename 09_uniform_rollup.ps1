Param(
  [string]$Root      = $(Get-Location).Path,
  [string]$ToolsDir  = "",
  [string]$InputsDir = "",
  [string]$CertsDir  = "",
  [int]$Dps          = 950,

  # Uniform sweep / base heights
  [string]$T0        = "10000000000000000",   # base reference T0 (1e16)
  [string]$TMin      = "10000000000000000",   # sweep start  (1e16)
  [string]$TMax      = "100000000000000000",  # sweep end    (1e17)

  # Rolling-T uniform cert controls
  [double]$DeltaTarget  = 1e-12,
  [int]$MeshInitial     = 256,
  [int]$MeshMax         = 131072,
  [int]$Digits          = 220,

  [switch]$UseVenv   = $true,
  [string]$VenvDir   = ""
)

function Msg($m)  { Write-Host "[09_uniform_rollup] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[09_uniform_rollup] WARNING: $m" -ForegroundColor Yellow }
function Assert-File($p, $m) { if (-not (Test-Path $p)) { Die $m } }

if (-not $ToolsDir)  { $ToolsDir  = Join-Path $Root "tools" }
if (-not $InputsDir) { $InputsDir = Join-Path $Root "" }
if (-not $CertsDir)  { $CertsDir  = Join-Path $Root "PROOF_PACKET" }
if (-not $VenvDir)   { $VenvDir   = Join-Path $Root "rh-env" }

$ErrorActionPreference = "Stop"

# Python resolver (same pattern as other stages)
$Py = "python"
if ($UseVenv) {
  $PyCandidate = Join-Path $VenvDir "Scripts\python.exe"
  if (Test-Path $PyCandidate) {
    $Py = $PyCandidate
  } else {
    Warn "venv requested but $PyCandidate not found, using python on PATH"
  }
}

Msg "Uniform rollup stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, CertsDir=$CertsDir, Dps=$Dps"
Msg "T0=$T0, TMin=$TMin, TMax=$TMax, DeltaTarget=$DeltaTarget"

# ---------------------------------------------------------------------------
# Core inputs this stage expects from earlier steps
# ---------------------------------------------------------------------------

# From early stages (for STP + rollup)
$BandCert    = Join-Path $CertsDir "band_cert.json"
$GammaTail   = Join-Path $CertsDir "gamma_tail.json"
$PrimeTail   = Join-Path $CertsDir "prime_tail_envelope.json"
$PrimeBlock  = Join-Path $CertsDir "prime_block_norm.json"
$WeilPsd     = Join-Path $CertsDir "weil_psd_bochner.json"
$Continuum   = Join-Path $CertsDir "continuum_operator_cert.json"

# From tails/explicit stage
$TailsJson   = Join-Path $CertsDir "analytic_tail_fit.json"
$Explicit    = Join-Path $CertsDir "explicit_formula.json"

# Basic assertions for things we *know* must already exist
Assert-File $BandCert   "band_cert.json not found (run 01_inputs.ps1)"
Assert-File $GammaTail  "gamma_tail.json not found (run 02_tails.ps1)"
Assert-File $PrimeTail  "prime_tail_envelope.json not found (run 02_tails.ps1)"
Assert-File $PrimeBlock "prime_block_norm.json not found (run 03_prime_block_and_grid.ps1)"
Assert-File $WeilPsd    "weil_psd_bochner.json not found (run 04_weil_psd.ps1)"
Assert-File $Continuum  "continuum_operator_cert.json not found (run 05_continuum_rollup.ps1)"

# Tails + explicit may or may not be present depending on how far pipeline has run
if (-not (Test-Path $TailsJson)) {
  Warn "analytic_tail_fit.json not found (run 07_tails_bundle_and_explicit.ps1); rollup_uniform.py and stp_test.py may be skipped"
}
if (-not (Test-Path $Explicit)) {
  Warn "explicit_formula.json not found (run 07_tails_bundle_and_explicit.ps1); rollup_uniform.py and stp_test.py may be skipped"
}

# ---------------------------------------------------------------------------
# 1) Rolling-T uniform certificate sweep -> rolling_uniform_cert.json
# ---------------------------------------------------------------------------

$RollingJson = Join-Path $CertsDir "rolling_uniform_cert.json"
$rt = Join-Path $ToolsDir "rolling_T_uniform_cert_v3.py"

if (Test-Path $rt) {
  Msg "python rolling_T_uniform_cert_v3.py --packet-dir $CertsDir --T0 $TMin --T1 $TMax --delta-target $DeltaTarget --mesh-initial $MeshInitial --mesh-max $MeshMax --digits $Digits --dps $Dps --out $RollingJson"
  & $Py $rt `
    --packet-dir $CertsDir `
    --T0 "$TMin" `
    --T1 "$TMax" `
    --delta-target "$DeltaTarget" `
    --mesh-initial "$MeshInitial" `
    --mesh-max "$MeshMax" `
    --digits "$Digits" `
    --dps "$Dps" `
    --out $RollingJson

  if (Test-Path $RollingJson) {
    Msg "rolling_uniform_cert.json created"
  } else {
    Warn "rolling_uniform_cert.json not created (check rolling_T_uniform_cert_v3.py CLI flags or write logic)"
  }
} else {
  Warn "rolling_T_uniform_cert_v3.py not found in $ToolsDir; skipping rolling uniform certificate"
}

# ---------------------------------------------------------------------------
# 2) Uniform rollup over gamma/prime/envelope + explicit formula
#    -> rollup_uniform.json (+ optional theory file)
# ---------------------------------------------------------------------------

$RollupJson  = Join-Path $CertsDir "rollup_uniform.json"
$TheoryJson  = Join-Path $CertsDir "rollup_uniform.theory.json"
$ru = Join-Path $ToolsDir "rollup_uniform.py"

if (Test-Path $ru) {
  if ((Test-Path $GammaTail) -and (Test-Path $PrimeTail) -and (Test-Path $Explicit)) {
    Msg "python rollup_uniform.py --T0 $T0 --gamma-envelope $GammaTail --prime-envelope $PrimeTail --explicit-formula $Explicit --out $RollupJson --theory-out $TheoryJson --dps $Dps"
    & $Py $ru `
      --T0 "$T0" `
      --gamma-envelope $GammaTail `
      --prime-envelope $PrimeTail `
      --explicit-formula $Explicit `
      --out $RollupJson `
      --theory-out $TheoryJson `
      --dps "$Dps"

    if (Test-Path $RollupJson) {
      Msg "rollup_uniform.json created"
    } else {
      Warn "rollup_uniform.json not created (check rollup_uniform.py CLI flags)"
    }
  } else {
    Warn "Missing one of gamma_tail.json / prime_tail_envelope.json / explicit_formula.json; skipping rollup_uniform.py"
  }
} else {
  Warn "rollup_uniform.py not found in $ToolsDir; skipping rollup uniform rollup"
}

# ---------------------------------------------------------------------------
# 3) STP consistency test -> diagnostic only (no JSON file)
# ---------------------------------------------------------------------------

$stp = Join-Path $ToolsDir "stp_test.py"

if (Test-Path $stp) {
  if ((Test-Path $TailsJson) -and (Test-Path $Explicit)) {
    Msg "python stp_test.py --band-cert $BandCert --weil-psd $WeilPsd --tails $TailsJson --prime-norm $PrimeBlock --explicit $Explicit --continuum-cert $Continuum --dps $Dps"
    & $Py $stp `
      --band-cert $BandCert `
      --weil-psd $WeilPsd `
      --tails $TailsJson `
      --prime-norm $PrimeBlock `
      --explicit $Explicit `
      --continuum-cert $Continuum `
      --dps "$Dps"

    if ($LASTEXITCODE -ne 0) {
      Warn "stp_test.py reported STP FAIL (see diagnostic output above)"
    } else {
      Msg "stp_test.py completed successfully"
    }
  } else {
    Warn "analytic_tail_fit.json or explicit_formula.json missing; skipping stp_test.py"
  }
} else {
  Warn "stp_test.py not found in $ToolsDir; skipping STP test"
}

Msg "Uniform rollup stage complete."
