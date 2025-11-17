Param(
  [string]$Root      = $(Get-Location).Path,
  [string]$ToolsDir  = "",
  [string]$InputsDir = "",
  [string]$CertsDir  = "",
  [int]$Dps          = 950,
  [double]$Sigma     = 6.0,
  [double]$K0        = 0.25,
  # Phase-9 knobs (can be overridden on the command line)
  [string]$T0        = "10000000000000000",   # 1e16
  [double]$X0        = 40.0,
  [double]$ACenter   = 0.30,
  [double]$BCenter   = 2.80,
  [switch]$UseVenv   = $true,
  [string]$VenvDir   = ""
)

function Msg($m)  { Write-Host "[06_lipschitz_and_density] $m" -ForegroundColor Cyan }
function Die($m)  { throw $m }
function Warn($m) { Write-Host "[06_lipschitz_and_density] WARNING: $m" -ForegroundColor Yellow }
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

Msg "Lipschitz and density stage"
Msg "Root=$Root, ToolsDir=$ToolsDir, InputsDir=$InputsDir, CertsDir=$CertsDir, Dps=$Dps, Sigma=$Sigma, K0=$K0"
Msg "T0=$T0, X0=$X0, ACenter=$ACenter, BCenter=$BCenter"

# 1) Lipschitz/q-bound -> lipschitz_q_bound.json
$LipsJson = Join-Path $CertsDir "lipschitz_q_bound.json"
$lq = Join-Path $ToolsDir "lipschitz_q_bound.py"
if (Test-Path $lq) {
  Msg "python lipschitz_q_bound.py --T0 $T0 --x0 $X0 --A-prime 1.0 --sigma-scale 1.0 --k0-scale 1.0 --K 1 --out $LipsJson --dps $Dps"
  & $Py $lq `
    --T0 "$T0" `
    --x0 "$X0" `
    --A-prime "1.0" `
    --sigma-scale "1.0" `
    --k0-scale "1.0" `
    --K "1" `
    --out $LipsJson `
    --dps "$Dps"
} else {
  Warn "lipschitz_q_bound.py not found in $ToolsDir; skipping Lipschitz bound"
}
if (Test-Path $LipsJson) {
  Msg "lipschitz_q_bound.json created"
} else {
  Warn "lipschitz_q_bound.json not created"
}

# 2) Density prover -> density_prover.json
$DensityJson = Join-Path $CertsDir "density_prover.json"
$dp = Join-Path $ToolsDir "density_prover.py"
if (Test-Path $dp) {
  Msg "python density_prover.py --a-center $ACenter --b-center $BCenter --out $DensityJson --dps $Dps"
  & $Py $dp `
    --a-center "$ACenter" `
    --b-center "$BCenter" `
    --out $DensityJson `
    --dps "$Dps"
} else {
  Warn "density_prover.py not found in $ToolsDir; skipping density certificate"
}
if (Test-Path $DensityJson) {
  Msg "density_prover.json created"
} else {
  Warn "density_prover.json not created"
}

Msg "Lipschitz and density stage complete."
